
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:08:58 2015

@author: maxwell
"""

import sys
import subprocess as sp
import numpy as np

def column(matrix, i):
    return [row[i] for row in matrix]


class Simulation():
    
    def __init__(self,parameters):
        self.parameters = parameters
        self.real = self.get_real_data()
        self.model = self.call_simulation()
        self.mse = self.get_MSE()        
        
    def get_real_data(self):
        
        real = np.loadtxt('data/real_data.txt')[:6000]
        real = np.transpose(np.array([column(real,1),column(real,2)]))    
     
        filter_range = 5 
     
        for i in range(len(real)):    
            if real[i][1] == 0 :
                real[i][1] = real[i-1][1]
                
            if i > filter_range:
                avg_sum = 0   
                for j in range(filter_range):
                    avg_sum = avg_sum + real[i-j][1]
                real[i][1] = avg_sum/filter_range    
        
        return real
       
    def call_simulation(self):
       
        f = open('parms','w')  
        for item in self.parameters:
            f.write('%f\n'%item)    
        f.close()        
            
        sp.call(["make"])    
        sp.call(["./modell","250","parms"])
        
        return np.loadtxt('model_data')
    
    def get_MSE(self):
    
        if len(self.model) != len(self.real):
            print "Error! Data not the same size"
            sys.exit()
        else :
            data_len = len(self.model)
           
        mse = 0           
        for i in range(data_len):
            mse += np.power(self.real[i][1] - self.model[i][1],2)
    
        return mse/data_len



