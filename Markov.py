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
        self.EField = ElectricField()
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
             
if __name__ == '__main__':
    print 'markov test'
    markov = Markov()
    print "markov.n = ",markov.n
    print "markov.N = ",markov.N
    markov.prepareT()
    # for i in range(markov.n):
    #     print markov.omega[i]
    #print markov.decoherence[0][0]
    print markov.T
    state = np.zeros((markov.N,1),complex)
    state [0][0] = 1.0
    print state
    result = linalg.expm(markov.T*1e-9,20)

    for i in range(100):
        state =  np.dot(result,state)
        print np.real(state[markov.N-1][0])

    
