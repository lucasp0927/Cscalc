#!/usr/bin/python2
from ElectricField import ElectricField
from Pulse import Pulse
from Markov import Markov
import sys
import time
import numpy as np
import argparse

def creat_matrix_file(input,output,ef):
    markov = Markov(input,output,ef)
    markov.prepareT()
    markov.prepareD()
    markov.zeroOrder()
    for i in xrange(50):
        print "-------------------------"
        print "order ",markov.order+1
        t1 = time.time()
        norm = markov.addOrder2()
        print "difference norm %e" %norm        
        t2 = time.time()
        print 'took %0.3f ms' % ((t2-t1)*1000.0)
        # markov.prepareT()
        # markov.prepareD()        
        #markov.plotGraph(title=str(i)+"th order")
        if norm == 0:
            break
        #markov.pp.close()
    markov.write()
    
def freq_data(input,ef):
    p = Pulse(input,ef)
    M = p.P - np.identity(p.N)
    p.freq_plot(1e5,500)#better be multiple of process number
    p.file_out.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("name", type=str,help="calculate")
    parser.add_argument("-m", "--matrix", action="store_true",help="calculate transition matrix file")
    parser.add_argument("-f", "--freq", action="store_true",help="plot frequency domain")    
    args = parser.parse_args()
    input_name = args.name
    base_name = str.split(args.name,'.')[0]
    factor = [1,2,3,4,5,6,7,8,9,10,100,200,300,400,500,600,700,800,900,1000]
    ef = ElectricField()
    
    for f in factor:
        ef.setfactor(f)
        print 'factor',f
        if args.matrix:
            creat_matrix_file(input_name,base_name+str(f),ef)
        if args.freq:
            freq_data(base_name+str(f)+".p",ef)
