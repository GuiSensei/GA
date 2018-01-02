# -*- coding: utf-8 -*-

#from math import pi

from components.individual import GAIndividual
from components.population import GAPopulation
from components.fitness import GAFitness

def func():
    y = a+b+c+d
    return y


fitness = GAFitness(population,func)

print(population)

