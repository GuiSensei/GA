# -*- coding: utf-8 -*-

from math import log2
from itertools import accumulate
from bisect import bisect_right
from random import uniform, random


class GAIndividual():
    def __init__(self,ranges,encoding = 'binary', eps = 0.001):#初始化
        self.ranges = ranges      #取值范围
        self.encoding = encoding  #编码方式
        self.eps = eps            #精度
        self._check_parameters()  #参数检查
        self.lengths, self.precisions = [],[] #染色体的长度、精度

        
        for (a,b), eps in zip(self.ranges,self.eps):
            length = int(log2((b-a)/eps)+1) #基因长度
            precision = (b - a)/(2**length-1) #基因的实际精度
            self.lengths.append(length) #所有基因的长度
            self.precisions.append(precision) #所有基因的实际精度

        self.gene_indices = self._get_gene_indices() #获取基因索引
        self.variants = self._init_variants() #初始化变量
        self.chromsome = self.encode() #生成染色体/个体
        self.real_number = self.decode(self.chromsome) #染色体解码
        self.each_fitness = fitness(self.real_number) #计算适应度

    
    #参数检查
    def _check_parameters(self): #类中方法名前一条下滑线表示私有方法，外部不要访问它（其实是可以访问的）
        if isinstance(self.eps,float):#确定精度是不是只有一个值，而且是float，不能是其他类型
            self.eps = [self.eps]*len(self.ranges)#如果只有一个值，将其转换成与ranges同一长度的列表
        else:                         #确定精度是否为list/tuple
            if len(self.eps) != len(self.ranges): #判断精度长度与ranges长度是否一致
                raise ValueError("变量个数与精度个数不一致！")
            for eps, (a,b) in zip(self.eps, self.ranges):  #一致的话，判断精度是否超出了ranges区间
                if eps > (b - a):
                    msg = "精度{}在范围({},{})有误！".format(eps,a,b)
                    raise ValueError(msg)

    
    #获取基因索引    
    def _get_gene_indices(self):
        end_indices = list(accumulate(self.lengths)) #每个基因结束的索引号保存在list中（因为最后一个数字不包括）
        start_indices = [0] + end_indices[:-1] #每个基因开始的索引号
        return list(zip(start_indices,end_indices))#[(gene1start,gene1end),(gene2start,gene2end),...]
    
    #初始化变量
    def _init_variants(self):
        variants = []
        for eps,(a,b) in zip(self.precisions,self.ranges):
            n_intervals = (b - a)//eps  #[a,b]的中间有n_intervals个数
            n = int(uniform(0,n_intervals + 1)) #基于均匀分布随机生成[0,n_intervals+1]中的一个数，取整
            variants.append(a+n*eps) #转换回来，添加到list中
        return variants
    
    #编码
    def encode(self):
        if self.encoding == 'decimal':  #十进制编码直接返回十进制数
            return self.variants
        
        #二进制编码
        chromsome = []
        for var, (a, _), length, eps in zip(self.variants,self.ranges,self.lengths,self.precisions):
            chromsome.extend(self.binarize(var - a, eps, length))
        return chromsome
    
    #转换二进制数
    @staticmethod #静态方法，既可以类访问，也可以实例访问;@classmethod动态方法/类方法，不通过实例化类来访问，跟类外声明函数是一样的
    def binarize(decimal, eps, length):#静态方法不用加self，动态方法第一个必须是cls（其他名称也可以）表示类。
        n = int(decimal/eps)
        bin_str = '{:0>{}b}'.format(n, length)#对n进行转进制，长度为{length}，b为二进制，>为右对齐，0表示前面用0补全，:前面为序号，可以不写
        return [int(i) for i in bin_str]
    
    #对编码进行解码
    def decode(self,chromsome):
        if self.encoding =='decimal':
            return self.variants
        
        variants = [self.decimalize(chromsome[start:end],eps,lower_bound)
                    for (start,end),(lower_bound, _), eps in 
                    zip(self.gene_indices,self.ranges,self.precisions)]
        return variants
    
    #转换十进制数
    @staticmethod
    def decimalize(binary,eps,lower_bound):
        bin_str = "".join([str(bit) for bit in binary])#将列表转换成字符串
        return lower_bound + int(bin_str,2)*eps
    

class GAPopulation():
    def __init__(self,indv_template,size=50): 
        if size % 2 != 0:
            raise ValueError('Population size must be an even number')
        self.size = size
        self.indv_template = indv_template
        self.individuals, self.indvs_decode = [],[] #编码、解码
        self.all_fitness = []
        for _ in range(self.size):
            indv = GAIndividual(ranges=self.indv_template.ranges,
                                    encoding=self.indv_template.encoding,
                                    eps=self.indv_template.eps)
            self.individuals.append(indv.chromsome)
            self.indvs_decode.append(indv.real_number)
            self.all_fitness.append(indv.each_fitness)

class Selection():
    def __init__(self,selection,individual,all_fitness):
        self.selection = selection
        self.individual = individual
        self.all_fitness = all_fitness
        if self.selection == 'rw':
            self.selected = self.roulettewheelselect()
        elif self.selection == 'er':
            pass
        elif self.selection == 'lr':
            pass
        elif self.selection =='tour':
            pass

    def roulettewheelselect(self):
        sum_fit = sum(self.all_fitness)
        wheel = list(accumulate([i/sum_fit for i in self.all_fitness]))    
        # Select a father and a mother.
        father_idx = bisect_right(wheel, random())
        father = self.individual[father_idx]
        mother_idx = (father_idx + 1) % len(wheel)
        mother = self.individual[mother_idx]
        return father, mother
    
class Crossover():
    def __init__(self, two_indv, pc=0.25,pe=0.5):
        if pc <= 0.0 or pc > 1.0:
            raise ValueError('错误的交叉概率')
        self.pc = pc
        if pe <= 0.0 or pe > 1.0:
            raise ValueError('Invalid genome exchange probability')
        self.pe = pe
        self.father = two_indv[0]
        self.mother = two_indv[1]
        self.crossed = self.uniformcrossover()

    def uniformcrossover(self):
        do_cross = True if random() <= self.pc else False
        if not do_cross:
            return self.father, self.mother
        chrom1 = self.father.copy()
        chrom2 = self.mother.copy()
        for i, (g1, g2) in enumerate(zip(chrom1, chrom2)):
            do_exchange = True if random() < self.pe else False
            if do_exchange:
                chrom1[i], chrom2[i] = g2, g1
        return chrom1, chrom2
            
class Mutation():
    def __init__(self,indv,pm=0.02):
        if pm <= 0.0 or pm > 1.0:
            raise ValueError('错误的变异概率')
        self.pm = pm
        self.indv = indv
        self.mutated = self.mutate()
        
    def mutate(self):
        do_mutation = True if random() <= self.pm else False
        individual = self.indv.copy()
        if do_mutation:
            for i, genome in enumerate(individual):
                do_flip = True if random() <= self.pm else False
                if do_flip:
                    individual[i] = genome^1
        return individual


#适应度公式
def fitness(real_number):#real_number为解码后的数字组成的list
    return real_number[0] + real_number[1]+ real_number[2]+ real_number[3]

if __name__ == "__main__":
    temp = GAIndividual(ranges = [(35,40),(20,30),(30,40),(8,15)],eps = 0.1)
    population = GAPopulation(temp) 
    selection = Selection('rw',population.individuals,population.all_fitness)
    crossover = Crossover(selection.selected)
    mutation = Mutation(crossover.crossed[1])
    print(population.individuals,population.indvs_decode,population.all_fitness)
    print(selection.selected)
    print(crossover.crossed)
    print(mutation.mutated)
