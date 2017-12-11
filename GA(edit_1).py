# Genetic Algorithm
#max f(a,b)=21.5+ a*sin(4*pi*a)+b*sin(20*pi*b)
#    -3.0 <= a <= 12.1
#    4.1 <= b <=5.8


import math
import random
from itertools import accumulate
from bisect import bisect_right


def initial(a,b,*code):
    '''
    a: the individual number of population
    b: the length of a(n) chromosome/individual
    '''
    L = []
    for j in range(a):
        LL=[]
        for i in range(b):
            LL.append(random.choice(code))
        L.append(tuple(LL))
    return tuple(L)


def decode(start,stop,L,k):
    return start + sum(j*2**(len(L)-i-1) for i,j in enumerate(L))*(stop-start)/(2**k-1)

def fitness(popu):
    A = []
    B = []
    for each in popu:
        a = decode(-3.0,12.1,each[:18],18)
        b = decode(4.1,5.8,each[18:],15)
        A.append(a)
        B.append(b)
    Y = []
    for i in range(len(A)):
         y = 21.5 + A[i] * math.sin(4 * math.pi * A[i]) + B[i] * math.sin(20 * math.pi * B[i])
         Y.append(y)
    return Y


def selection(popu,fit):
    sum_fit = sum(fit)
    each_fit = [i/sum_fit for i in fit]
    wheel = list(accumulate([i for i in each_fit]))
    #print(wheel)
    popu_select = []
    for i in range(len(popu)):
        index = bisect_right(wheel, random.random())
        popu_select.append(popu[index])
    return tuple(popu_select)


def crossover(popu,pc):
    cross_popu = []
    while True:
        mother,father = list(popu[random.randrange(len(popu))]),list(popu[random.randrange(len(popu))])
        if mother == father:
            continue
        if random.random() <= pc:
            index = random.randrange(len(popu[0]))
            mother[index:],father[index:] = father[index:],mother[index:]
        cross_popu.append(tuple(mother))
        cross_popu.append(tuple(father))
        if len(cross_popu) == 20:
            break
    return tuple(cross_popu)


def mutation(popu,pm):
    mut_popu = []
    for each in popu:
        each = list(each)
        if random.random() <= pm:
            index = random.randrange(len(each))
            each[index] = int(not each[index])
        mut_popu.append(tuple(each))
    return tuple(mut_popu)


def run(n=1):
    '''
    n: the number of iterations
    '''
    max_fit = []
    max_ind = []
    ind_num = 20 #个体数
    chrom_length = 33 #染色体长度
    code = [0,1]  #编码
    pc = 0.25   #交叉概率
    pm = 0.01   #变异概率
    population = initial(ind_num,chrom_length,*code)
    cal_fit = fitness(population)
    max_fit.append(max(cal_fit))
    max_ind.append(population[cal_fit.index(max(cal_fit))])
    for i in range(n):        
        step1 = selection(population,cal_fit)
        step2 = crossover(step1,pc)
        step3 = mutation(step2,pm)
        population = step3
        cal_fit = fitness(population)
        max_fit.append(max(cal_fit))
        max_ind.append(population[cal_fit.index(max(cal_fit))])
    optimal_result = max(max_fit)
    optimal_index = max_fit.index(optimal_result)
    optimal_individual = max_ind[optimal_index]
    a = decode(-3.0,12.1,optimal_individual[:18],18)
    b = decode(4.1,5.8,optimal_individual[18:],15)
    print("最优解：",optimal_result)    
    print("其值为：a=",a,'b=',b)
    print("最优代：",optimal_index)

if __name__ == '__main__':
    run(1000)


