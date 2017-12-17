# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 19:44:18 2017

@author: gfz
"""
import copy
import numpy as np
import random
import pandas as pd
from bisect import bisect_right
from itertools import accumulate
import matplotlib.pyplot as plt

#归一化方法，ways=0按列归一化，1按行归一化，2全局归一化，默认全局归一化
def normalization(dataframe, ways = 2):
    if isinstance(dataframe,pd.DataFrame):
        columns_max = np.max(dataframe,axis = 0)
        columns_min = np.min(dataframe,axis = 0)
        rows_max = np.max(dataframe,axis = 1)
        rows_min = np.min(dataframe,axis = 1)
        if ways == 0:
            return (dataframe - columns_min)/(columns_max - columns_min)
        elif ways == 1:
            return ((dataframe.T - rows_min)/(rows_max - rows_min)).T
        elif ways == 2:
            return (dataframe - min(rows_min))/(max(rows_max) - min(rows_min))
        else:
            return "归一化参数2输入有误！！！"
    else:
        return "归一化参数1输入有误！！！"
#成本数据
cost = normalization(pd.read_excel('myDataFrame.xlsx',sheetname = 'cost'))
#碳足迹数据
cfp = normalization(pd.read_excel('myDataFrame.xlsx',sheetname = 'cfp'))
#零部件
part = pd.read_excel('myDataFrame.xlsx',sheetname = 'part', header = 0)
#可合并零部件
combine_part = pd.read_excel('myDataFrame.xlsx',sheetname = 'combine')
#可分解零部件
decompose_part = pd.read_excel('myDataFrame.xlsx',sheetname = 'decompose')
#材料类型数据
material = list(cost.columns)




#A是DataFrame，B是一维list
def chromesome(A,B):
    chromesomeList = []
    for i in range(len(list(A.columns))):
        chromesomeList.append(random.choice(((np.array(A.T)).tolist())[i]))
    for i in range(len(list(A.columns))):
        chromesomeList.append(random.choice(B))
    return chromesomeList

#种群初始化
def population(*args):#args中包含3个参数，按顺序为：种群个体数目，零件参数、材料参数
    try:
        for i in range(args[0]):
            #populationList.append(chromesome(args[1],args[2],args[3]))
            yield chromesome(args[1],args[2])
    except:
        return "初始种群生成失败！"

def fit_fun(C):
    return 1 - (sum(cost[C[i+int(len(C)/2)]][C[i]] for i in range(int(len(C)/2))) +\
           sum(cfp[C[i+int(len(C)/2)]][C[i]] for i in range(int(len(C)/2))))/(len(C))

def elite(population, fitness, num):
    popu = population.copy()
    fit = fitness.copy()
    max_pop = []
    for i in range(num):
        idx = fit.index(max(fit))
        max_pop.append(popu.pop(idx))
        fit.pop(idx)
    return max_pop

def select(population, fitness):

    fit = fitness
    sum_fit = sum(fit)
    wheel = list(accumulate([i/sum_fit for i in fit]))
    father_idx = bisect_right(wheel, random.random())
    father = population[father_idx]
    mother_idx = (father_idx + 1) % len(wheel)
    mother = population[mother_idx]
    return father, mother

def cross(father, mother, pc, pe = 0.5):
    chrom1 = father.copy()
    chrom2 = mother.copy()

    do_cross = True if random.random() <= pc else False
    if do_cross:
        for i, (g1, g2) in enumerate(zip(chrom1, chrom2)):
            do_exchange = True if random.random() < pe else False
            if do_exchange:
                chrom1[i], chrom2[i] = g2, g1
    return chrom1, chrom2

def mutate(individual, pm):
    do_mutation = True if random.random() <= pm else False
    if do_mutation:
        for idx, gene in enumerate(individual):
            do_flip = True if random.random() <= pm else False
            if do_flip:
                A = list(part['A'])
                B = list(part['B'])
                C = list(part['C'])
                D = list(part['D'])
                E = list(part['E'])
                F = list(part['F'])
                M = material.copy()
                if gene in list(part['A']):
                    A.remove(gene)
                    gene = random.choice(A)
                if gene in list(part['B']):
                    B.remove(gene)
                    gene = random.choice(B)
                if gene in list(part['C']):
                    C.remove(gene)
                    gene = random.choice(C)
                if gene in list(part['D']):
                    D.remove(gene)
                    gene = random.choice(D)
                if gene in list(part['E']):
                    E.remove(gene)
                    gene = random.choice(E)
                if gene in list(part['F']):
                    F.remove(gene)
                    gene = random.choice(F)
                if gene in material:
                    M.remove(gene)
                    gene = random.choice(M)
            individual[idx] = gene
    return individual

def run(popu,fit):
    allpopu = []
    newfit = []
    step_0 = elite(popu,fit,6)
    #print(step_0)
    allpopu = copy.deepcopy(step_0)
    #print(allpopu)
    for i in range(len(popu)):
        step_1 = select(popu, fit)
        step_2 = cross(step_1[0],step_1[1],0.25)#交叉概率一般0.25~1.0
        for each in step_2:
            if each in allpopu:
                continue
            else:
                allpopu.append(each)
    #print(allpopu)

    for each in popu:
        #print(each)
        step_3 = mutate(each,0.5)#交叉概率
        #print(step_3)
        allpopu.append(step_3)
        #print(allpopu)
    allfit = [fit_fun(i) for i in allpopu]
    #print(allpopu)
    newdata = elite(allpopu,allfit,20)#进入下一代个数
    #print(newdata)
    for i in newdata:
        newfit.append(fit_fun(i))

    return newdata,newfit


if __name__ == '__main__':
    iteration = 200 #迭代次数

    data = []
    fdata = []
    for i in population(20, part, material):
        fdata.append(fit_fun(i))
        data.append(i)

    with open('data.txt','w') as f:
        a,b = data,fdata
        mean = []
        Max = []
        for count in range(iteration+1):
            f.write(("第{0}代").format(count)+'\n')
            for index,each in enumerate(a):
                f.write(str(each)+"\'s fitness = "+str(b[index])+'\n')
            mean.append(sum(b)/len(b))
            Max.append(max(b))
            f.write("均值："+str(mean[-1])+'\n')
            f.write("最值："+str(Max[-1])+'\n')
            if count == iteration:
                break
            else:
                a,b = run(a,b)

    try:
        plt.plot(list(np.arange(len(mean))),mean,marker='.', mec='r',label="mean")
        plt.plot(list(np.arange(len(Max))),Max,marker='^', mec='b',label="max")
        plt.xlabel("n") #横坐标
        plt.ylabel("Fitness") #纵坐标
        plt.legend() #显示图例
        plt.savefig('GA', dpi = 600)#保存图片
        plt.show()#图片显示
        print("图片已保存，程序运行结束！")
    except Exception as e:
        print("图片有误！")
