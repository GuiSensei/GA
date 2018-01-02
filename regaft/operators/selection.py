#!/usr/bin/env python
# -*- coding: utf-8 -*-

from components.population import GAPopulation

class GASelection():
    def __init__(self,population,fitness):
        self.population = population
        self.fitness = fitness
        

    def RouletteWheelSelection(self):
        fit = population.all_fit(fitness)#种群所有的适应度值
        min_fit = min(fit)#最小的适应度值
        fit = [(i-minfit) for i in fit]
        
        sum_fit = sum(fit)
        
    
    
a = RouletteWheelSelection(1,2)
print(a.fitness)