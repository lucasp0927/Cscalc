#!/usr/bin/python
from __future__ import division
import sys
from numpy import *
from scipy import linalg
from ElectricField import ElectricField
import matplotlib.pyplot as plt
import pickle

class Pulse(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        file_in = sys.argv[1]        
        #        dictf = open(file_in,'r')
        self.parameter = pickle.load( open( file_in, "rb" ) )
         #self.parameter = eval(dictf.read())
        self.T = self.parameter['T']
        self.P = self.parameter['P']
        self.n = self.parameter['n']
        self.N = self.parameter['N']
        self.group = self.parameter['group']                
        self.cutoff = self.parameter['cutoff']        

    def ij2idx(self,i,j):
        """
        0 1 2
        3 4 5
        6 7 8
        """
        #idx = (i*(2*self.n-i+1))/2+(j-i)
        idx = self.n*i+j
        return idx
        
    def time_plot(self,rep,num):
        r_t = rep - self.cutoff
        expT = linalg.expm(self.T*r_t)
        start = 1
        state = zeros(self.N,complex)
        time_arr = []
        time = 0.0
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = zeros((3,num*2))
        for i in range(num):
            state = dot(self.P,state.T)
            time += self.cutoff
            time_arr.append(time)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][2*i] += real(state[self.ij2idx(k,k)])
            state = dot(expT,state.T)
            print state
            time += r_t
            time_arr.append(time)            
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][2*i+1] += real(state[self.ij2idx(k,k)])                    
        plt.figure(1)                    
        fig = plt.subplot(1,1,1)
        plt.title("test")
        plt.ylim(0,1)
        plt.xlabel('time')
        plt.ylabel('population')
        for i in range(3):
            fig.plot(time_arr,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.show()
        
    def correct(self,arr):
        # TODO: write in better numpy way
        for i in range(self.n):
            if (real(arr[self.ij2idx(i,i)])<0.0 or real(arr[self.ij2idx(i,i)])>1.0):
                return False
        else:
            return True

    def freq_plot(self,start,end,number):
        rept = linspace(start,end,number)
        for t in rept:
            M = dot(self.P,linalg.expm(self.T*(t-self.cutoff)))
            print linalg.eig(M)
if __name__ == '__main__':
    p = Pulse()
    p.time_plot(1e-8,5000)
    M = dot(linalg.expm(p.T*(1e-8-p.cutoff)),p.P)
    W,V = linalg.eig(M)
    for v in V.T:
        if p.correct(v):
            print v
    
    #p.freq_plot(1e-7,2e-7,2)
