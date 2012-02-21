#!/usr/bin/python2
from __future__ import division
import sys
from ElectricField import ElectricField
import numpy as np

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
        self.n = self.parameter['n']# number of levels
        self.N = ((self.n+1)*self.n)/2# the number of independent terms in  markov matrix
        self.decoherence = self.parameter['decoherence_matrix']
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

if __name__ == '__main__':
    print 'markov test'
    markov = Markov()
    print "markov.n = ",markov.n
    print "markov.N = ",markov.N
    print markov.decoherence[0][0]
