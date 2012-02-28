#!/usr/bin/python
from __future__ import division
import sys
import numpy as np
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
        state = np.zeros(self.N,complex)
        time_arr = []
        time = 0.0
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,num*2))
        for i in range(num):
            state = np.dot(self.P,state.T)
            time += self.cutoff
            time_arr.append(time)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][2*i] += np.real(state[self.ij2idx(k,k)])
            state = np.dot(expT,state.T)
            #print state
            time += r_t
            time_arr.append(time)            
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][2*i+1] += np.real(state[self.ij2idx(k,k)])                    
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

    def freq_plot(self,start,end,number,pnum):
        repf = np.linspace(start,end,number)
        rept = 1.0/repf
        start = 1
        state = np.zeros(self.N,complex)
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,number))        
        for t in enumerate(rept):
            M = np.dot(linalg.expm(self.T*(t[1]-self.cutoff)),self.P)
            M = np.linalg.matrix_power(M,pnum)
            state1 = np.dot(M,state.T)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][t[0]] += np.real(state1[self.ij2idx(k,k)])
        plt.figure(1)                    
        fig = plt.subplot(1,1,1)
        plt.title("test")
        plt.ylim(-1,2)
        plt.xlabel('repetition rate(Hz)')
        plt.ylabel('population')
        for i in range(3):
            fig.plot(repf,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.show()                            

            
if __name__ == '__main__':
    p = Pulse()
    #p.time_plot(1.7e-9,10000)
    p.freq_plot(1e8,2e8,10000,10000)    
    #p.freq_plot(1e-9,2e-9,10000,20000)
