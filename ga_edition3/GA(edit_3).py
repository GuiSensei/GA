#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
# Genetic Algorithm
参数修改：只需要修改run和fitness中的参数
变量：a,b
目标： max f(a,b)=21.5+ a*sin(4*pi*a)+b*sin(20*pi*b)
约束：   -3.0 <= a <= 12.1
         4.1 <= b <=5.8
'''

import heapq  #精英保留中使用
from math import log2,log10,sin,pi #sin和pi是目标函数中出现的，如无请删除
import random
from itertools import accumulate #轮盘赌选择中使用
from bisect import bisect_right #轮盘赌选择中使用


#基因长度、精度、索引号
def gene(ranges,eps):
    if isinstance(eps,float):#确定精度是不是只有一个值，而且是float，不能是其他类型
        eps = [eps]*len(ranges)#如果只有一个值，将其转换成与ranges同一长度的列表
    else:                         #确定精度是否为list/tuple
        if len(eps) != len(ranges): #判断精度长度与ranges长度是否一致
            raise ValueError("变量个数与精度个数不一致！")
        for eps, (a,b) in zip(eps, ranges):  #一致的话，判断精度是否超出了ranges区间
            if eps > (b - a):
                msg = "精度{}在范围({},{})有误！".format(eps,a,b)
                raise ValueError(msg)

    lengths,precisions=[],[]
    for (a, b), eps in zip(ranges, eps):
        length = int(log2((b - a)/eps)+1)
        precision = (b - a)/(2**length-1)
        lengths.append(length)
        precisions.append(precision)
    
    #基因索引号
    end_indices = list(accumulate(lengths)) #每个基因结束的索引号保存在list中（因为最后一个数字不包括）
    start_indices = [0] + end_indices[:-1] #每个基因开始的索引号
    gene_indices = list(zip(start_indices,end_indices))#[(gene1start,gene1end),(gene2start,gene2end),...]
    return lengths,precisions,gene_indices

#初始化变量，十进制数
def init_variants(precisions,ranges,eps):
    variants = []
    for eps,(a,b) in zip(precisions,ranges):
        n_intervals = (b - a)//eps  #[a,b]的中间有n_intervals个数
        n = int(random.uniform(0,n_intervals + 1)) #基于均匀分布随机生成[0,n_intervals+1]中的一个数，取整
        variants.append(a+n*eps) #转换回来，添加到list中
    return variants

#将初始变量编码
def encode(variants,ranges,lengths,precisions):
    each_chrom = []
    for var, (a, _), length, eps in zip(variants,ranges,lengths,precisions):
        each_chrom.extend(binarize(var - a, eps, length))
    return each_chrom

#十进制转换为二进制
def binarize(decimal, eps, length):
    n = int(decimal/eps)
    bin_str = '{:0>{}b}'.format(n, length)#对n进行转进制，长度为{length}，b为二进制，>为右对齐，0表示前面用0补全，:前面为序号，可以不写
    return [int(i) for i in bin_str]

#初始化种群
def populate(ranges,eps,indv_number):
    each_gene = gene(ranges,eps)
    gene_length = each_gene[0] #基因长度list
    gene_precisions = each_gene[1]#基因精度list
    population = []
    for _ in range(indv_number):
        variants = init_variants(gene_precisions,ranges,eps)
        chromsome = encode(variants,ranges,gene_length,gene_precisions)
        population.append(chromsome)        
    return population


#解码
def decode(ranges,eps,chromsome):
    each_gene = gene(ranges,eps)
    gene_indices = each_gene[2]
    precisions = each_gene[1]
    variants = [decimalize(chromsome[start:end],eps,lower_bound)
                for (start,end),(lower_bound, _), eps in 
                zip(gene_indices,ranges,precisions)]
    return variants

#将二进制数转换十进制数
def decimalize(binary,eps,lower_bound):
    bin_str = "".join([str(bit) for bit in binary])#将列表转换成字符串
    return lower_bound + int(bin_str,2)*eps

#精英保留，默认2个
def elit(popu,fit,elit_num=2):
    if elit_num%2 != 0:
        raise ValueError("保留的精英数应为偶数")
    #elite_fit = heapq.nlargest(n,fit)#最大n个适应度
    elite_index = heapq.nlargest(elit_num,range(len(fit)),fit.__getitem__)#最大适应度的索引
    elite_popu = [popu[index]for index in elite_index] #最大n个个体
    return elite_popu

#移除劣解，默认2个
def removebad(popu,fit,n=2):
    popu = list(popu)
    fit = list(fit)
    for i in range(n):
        popu.pop(fit.index(min(fit)))
        fit.pop(fit.index(min(fit)))
    return tuple(popu)

#选择（轮盘赌RouletteWheelSelection）
def rwselect(popu,fit):
    sum_fit = sum(fit)
    each_fit = [i/sum_fit for i in fit]
    wheel = list(accumulate([i for i in each_fit]))
    popu_encode_select = []
    #popu_decode_select = []
    for i in range(2):
        index = bisect_right(wheel, random.random())
        popu_encode_select.append(popu[index])
    return popu_encode_select

#均匀交叉
def uniformcrossover(two_indv, pc=0.25,pe=0.5):
    if pc <= 0.0 or pc > 1.0:
        raise ValueError('错误的交叉概率')
    if pe <= 0.0 or pe > 1.0:
        raise ValueError('无效的基因点交换概率')
    do_cross = True if random.random() <= pc else False
    if not do_cross:
        return two_indv
    father = two_indv[0]
    mother = two_indv[1]
    chrom1 = father.copy()
    chrom2 = mother.copy()
    for i, (g1, g2) in enumerate(zip(chrom1, chrom2)):
        do_exchange = True if random.random() < pe else False
        if do_exchange:
            chrom1[i], chrom2[i] = g2, g1
    return [chrom1, chrom2]

#变异
def mutation(indv,pm=0.05):
    if pm <= 0.0 or pm > 1.0:
        raise ValueError('错误的变异概率')
    do_mutation = True if random.random() <= pm else False
    individual = indv.copy()
    if do_mutation:
        for i, genome in enumerate(individual):
            do_flip = True if random.random() <= pm else False
            if do_flip:
                individual[i] = genome^1
    return individual

#选择、交叉、变异操作
def operate(population,cal_fit,pc,pm): 
    step1 = rwselect(population,cal_fit)#从population中选择出2个
    step2 = uniformcrossover(step1,pc) #对上述两个进行交叉
    step3 = [mutation(each,pm) for each in step2] #对交叉结果进行变异
    return step3


#数据存储到当前目录下的all_fit.py文件中
def datastore(iter_num,mean,variants,fitness_values):    
    with open('all_fit.py', 'w', encoding='utf-8') as f:
        f.write('all_fit = [\n')
        for ng,m, x, y in zip(iter_num,mean, variants, fitness_values):
            f.write('    ({}, {}, {}, {}),\n'.format(ng,m, x, y))
        f.write(']\n\n')

#绘图
def draw_fig():
    import matplotlib.pyplot as plt
    
    from all_fit import all_fit    
    steps, mean, fits, variants = list(zip(*all_fit))
    best_step, best_v, best_f = steps[-1], variants[-1], fits[-1]

    fig = plt.figure()   
    ax = fig.add_subplot(111)
    ax.plot(steps, fits)    
    ax.plot(steps, mean)
    plt.legend(['Max Fitness','Mean Fitness'])
    ax.set_xlabel('Generation')
    ax.set_ylabel('Fitness')
    
    # 标注最大值点
    ax.scatter([best_step], [best_f], color='r',marker = '^')
    ax.annotate(s='x: [{:.2f}, {:.2f}]\ny:{:.2f}'.format(*best_v, best_f),
                xy=(best_step, best_f),
                xytext=(best_step-10, best_f))
       
    plt.savefig('fit.png',dpi=600,transparent=True)    #保存图片
    plt.show()  #图片显示

        
#适应度函数
def fitness(indv):   
    a,b = indv[0],indv[1]
    return 21.5+ a*sin(4*pi*a)+b*sin(20*pi*b)

def run(n=1):#n为迭代次数
    indv_number = 20   #种群中个体数
    ranges=[(-3.0,12.1),(4.1,5.8)] #取值范围
    eps = 0.0001    #精度要求
    pc = 0.35       #交叉概率
    pm = 0.15       #变异概率
    point = int(-log10(eps))
    init_population = populate(ranges,eps,indv_number)
    init_decoded = [decode(ranges,eps,chromsome) for chromsome in init_population]
    cal_fit = [fitness(i) for i in init_decoded]
    elited = elit(init_population,cal_fit,elit_num=2)    
    allkids = (init_population).copy()
    population = init_population
    iter_list = []
    best_solution = []
    best_fitness = []
    all_mean_fit = []
    
    for iter_num in range(n):
        new_kids,kids = [],[]
        while True:
            #print('执行循环...')
            for _ in range(indv_number//2):
                kids += operate(population,cal_fit,pc,pm)
            kids += (elited + new_kids) #精英保留的+生成的
            new_kids = [list(t) for t in set(tuple(_) for _ in kids)]#去除重复个体
            #new_kids = [list(t) for t in (set(tuple(_) for _ in new_kids)-set(tuple(_) for _ in allkids))]#去除已出现过个体
            allkids = [list(t) for t in set(tuple(_) for _ in allkids+new_kids)] #所有个体           
            if len(new_kids) - indv_number < 0:#个数不足继续执行
                continue
            if len(new_kids) - indv_number > 0:#个数过多删除劣解
                indvs_fits=[]
                for kid in new_kids:                
                    indv = decode(ranges,eps,kid)                
                    indv_fit = fitness(indv)
                    indvs_fits.append(indv_fit)
                new_kids = removebad(new_kids,indvs_fits,n=len(new_kids) - indv_number)
            if len(new_kids) - indv_number == 0:
                children = new_kids
                break        
        population = children
        decoded_population = [decode(ranges,eps,chromsome) for chromsome in population]
        cal_fit = [fitness(i) for i in decoded_population]
        elited = elit(population,cal_fit,elit_num=2)
        optimal_fitness = max(cal_fit) #最大适应度值
        mean_fitness = sum(cal_fit)/len(cal_fit)
        optimal_index = cal_fit.index(optimal_fitness)
        #optimal_individual = population[optimal_index]
        optimal_solution = [round(i,point) for i in decoded_population[optimal_index]] #最好的解码后的染色体
        print("第{0}代:平均适应度值为{1:.{4}f};最优解:{2:.{4}f},其值分别为:{3}".format(
                iter_num,mean_fitness,optimal_fitness,optimal_solution,point))
        iter_list.append(iter_num)
        all_mean_fit.append(mean_fitness)
        best_solution.append(optimal_solution)
        best_fitness.append(optimal_fitness)        
    datastore(iter_list,all_mean_fit,best_fitness,best_solution)
    return 0

if __name__ == '__main__':
    run(100) #执行次数
    draw_fig() #图形显示，需要安装matplotlib，如未安装请注释


