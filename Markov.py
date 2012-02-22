#!/usr/bin/python2
from __future__ import division
import sys
from ElectricField import ElectricField
import numpy as np
from scipy import linalg

class Markov(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        file_in = sys.argv[1]
        file_out = sys.argv[2]
        dictf = open(file_in,'r')
        self.parameter = eval(dictf.read())
        self.omega = self.parameter['omega']
        self.group = self.parameter['level_group']                
        self.n = self.parameter['n']# number of levels
        self.N = ((self.n+1)*self.n)/2# the number of independent terms in  markov matrix
        self.decoherence = self.parameter['decoherence_matrix']
        
        self.T = np.zeros((self.N,self.N),complex) # time independent part d rho/ dt = T rho
        self.final = np.zeros((self.N,self.N),complex) # final markov matrix
        
        self.EField = ElectricField()
        self.tsample = np.linspace(0,self.EField.cutoff,self.EField.sample)
        #self.Dfunction = [np.zeros((self.N,self.N),complex) for i in range(self.EField.sample)] # time dependent part
        self.Dfunction=[]
        self.order = 0
        dictf.close()

    def ij2idx(self,i,j):
        """
        0 1 2
          3 4
            5
        """
        assert i <= j             # i row j column
        idx = (i*(2*self.n-i+1))/2+(j-i)
        assert idx < self.N
        return idx

    def rotate_omega(self,i,j):
        if i in self.group[0]:
            if j not in self.group[0]:
                return self.omega[i]-self.omega[j]-self.EField.carrier_freq
            else:
                return self.omega[i]-self.omega[j]
        else:
            return self.omega[i]-self.omega[j]
        
    def prepareT(self):
        for i in range(self.n):
            for j in range(i,self.n):
                for k in self.decoherence[i][j]:
                    self.T[self.ij2idx(i,j)][self.ij2idx(k[0],k[1])]+=k[2]
        for i in range(self.n):
             for j in range(i+1,self.n):
                 self.T[self.ij2idx(i,j)][self.ij2idx(i,j)]+= -1.0j*self.rotate_omega(i,j)
                 
    def zeroOrder(self):
        for i in enumerate(self.tsample):
            print i[0]
            self.Dfunction.append(linalg.expm(markov.T*i[1]))

    def addOrder(self):
        self.order += 1
        for i in range(self.n):
            for j in range(i,self.n):
                # plus H rho
                # minus rho H
                pass
    
    def finalResult(self):
        time = 2*np.pi/self.EField.repetition_freq-self.EField.cutoff
        print time
        self.final = np.dot(self.Dfunction[-1],linalg.expm(markov.T*time))
    
if __name__ == '__main__':
    print 'markov test'
    markov = Markov()
    markov.prepareT()
    # for i in range(markov.n):
    #     print markov.omega[i]
    #print markov.decoherence[0][0]
    state = np.zeros((markov.N,1),complex)
    state [0][0] = 1.0
    result = linalg.expm(markov.T*1e-9,20)
    markov.zeroOrder()

    markov.finalResult()
    print markov.final
