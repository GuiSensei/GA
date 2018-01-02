# -*- coding: utf-8 -*-

from .individual import GAIndividual


    
  
        
if __name__ == "__main__":
    temp = GAIndividual(ranges = [(35,40),(20,30),(30,40),(8,15)],eps = 0.1)
    population = GAPopulation(temp).result
    print(population)
    
