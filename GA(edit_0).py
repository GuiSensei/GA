# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 18:24:52 2017

@author: Francis
"""

# Genetic Algorithm
# max f(a,b)=21.5+ a*sin(4*pi*a)+b*sin(20*pi*b)
#    -3.0 <= a <= 12.1
#    4.1 <= b <=5.8


import math
import random
from itertools import accumulate
from bisect import bisect_right


def initial():
    L = []
    for j in range(20):
        LL=[]
        for i in range(33):
            LL.append(random.choice([0,1]))
        L.append(LL)
    return L


def decode(start,stop,L,k):
    return start + sum(j*2**(len(L)-i-1) for i,j in enumerate(L))*(stop-start)/(2**k-1)

def fitness(A,B):
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
    return popu_select


def crossover(popu,pc):
    cross_popu = []
    while True:
        mother,father = popu[random.randrange(len(popu))],popu[random.randrange(len(popu))]
        if mother == father:
            continue
        if random.random() <= pc:
            index = random.randrange(len(popu[0]))
            mother[index:],father[index:] = father[index:],mother[index:]
        cross_popu.append(mother)
        cross_popu.append(father)
        if len(cross_popu) == 20:
            break
    return cross_popu


def mutation(popu,pm):
    mut_popu = []
    for each in popu:
        if random.random() <= pm:
            index = random.randrange(len(each))
            each[index] = int(not each[index])
        mut_popu.append(each)
    return mut_popu


def run():
    pc = 0.25
    pm = 0.01
    population = initial()
    Max = []
    for i in range(10000):
        A = []
        B = []
        for each in population:
            a = decode(-3.0,12.1,each[:18],18)
            b = decode(4.1,5.8,each[18:],15)
            A.append(a)
            B.append(b)
        cal_fit = fitness(A,B)
        step1 = selection(population,cal_fit)
        step2 = crossover(step1,pc)
        step3 = mutation(step2,pm)
        Max.append(max(cal_fit))
    print(max(Max))

if __name__ == '__main__':
    run()
