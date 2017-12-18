# -*- coding: utf-8 -*-

from math import log2
from itertools import accumulate
from random import uniform

class GAIndivadual():
    def __init__(self,ranges,encoding = 'binary', eps = 0.001):#初始化
        """
        1. 变量需要用list/tuple嵌套，如ranges = [(1,2),(3,4),(5,6)]
        2. 精度可以只有一个，如果只有一个，表示所有变量均为同一精度，默认0.001（保留三位小数）；
        3. 如果不同变量取不同精度，那么也需要用list/tuple，且其长度应与ranges长度一致。
        4. 编码方式默认二进制binary，也可十进制decimal
        """
        self.ranges = ranges      #取值范围
        self.encoding = encoding  #编码方式
        self.eps = eps            #精度
        self._check_parameters()  #参数检查
        self.lengths, self.precisions = [],[] #染色体的长度、精度
        
        
        for (a,b), eps in zip(self.ranges,self.eps):
            length = int(log2((b-a)/eps)) #基因长度（此处长度可能偏小）
            precision = (b - a)/(2**length) #基因的精度（长度小了之后，精度可能达不到）
            
            #此处省略了判断精度的冗余
            
            self.lengths.append(length) #所有基因的长度
            self.precisions.append(precision) #所有基因的精度
            
        self.gene_indices = self._get_gene_indices() #获取基因尺寸
        self.variants = self._init_variants() #初始化变量
        self.chromosome = self.encode() #生成染色体/个体
        
    
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
    
    #获取基因尺寸    
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
        chromosome = []
        for var, (a, _), length, eps in zip(self.variants,self.ranges,self.lengths,self.precisions):
            chromosome.extend(self.binarize(var - a, eps, length))
        return chromosome
    
    #转换二进制数
    @staticmethod #静态方法，既可以类访问，也可以实例访问;@classmethod动态方法/类方法，不通过实例化类来访问，跟类外声明函数是一样的
    def binarize(decimal, eps, length):#静态方法不用加self，动态方法第一个必须是cls（其他名称也可以）表示类。
        n = int(decimal/eps)
        bin_str = '{:0>{}b}'.format(n, length)#对n进行转进制，长度为{length}，b为二进制，>为右对齐，0表示前面用0补全，:前面为序号，可以不写
        return [int(i) for i in bin_str]



a = GAIndivadual(ranges = [(-3.0,12.1),(4.1,5.8)],eps = 0.0001)
b = a.chromosome
print(b)
print(len(b))
