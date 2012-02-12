import sys
from ElectricField import ElectricField
import numpy as np

class Markov(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        self.parameter = {}

if __name__ == '__main__':
    print 'markov test'
    markov = Markov()
    file_in = sys.argv[1]
    file_out = sys.argv[2]
    dictf = open(file_in,'r')
    markov.parameter = eval(dictf.read()) # watch out for security issue
    dictf.close()
    print markov.parameter['level_group']
    print markov.parameter['n']
    
