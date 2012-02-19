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
        self.n = self.parameter['n']
        self.N = ((self.n+1)*self.n)/2# the size of markov matrix
        dictf.close()

    def ij2idx(self,i,j):
        """
        0 1 2
          3 4
            5
        """
        assert i <= j             # i row j column
        idx = (i*(2*self.n-i+1))/2+(j-i)
        return idx

if __name__ == '__main__':
    print 'markov test'
    markov = Markov()
    print markov.n
    print markov.parameter['decoherence_matrix'][0][0]
    markov.n = 4
    for i in range(4):
        for j in range(i,4):
            print i,j
            print "idx",markov.ij2idx(i,j)
