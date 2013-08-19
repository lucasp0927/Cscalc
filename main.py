#!/usr/bin/python2
from ElectricField import ElectricField
from Pulse import Pulse
from Markov import Markov
import sys
import time
import numpy as np
import argparse

def create_matrix_file(input,output,ef):
"""
Create the .p files contain supermatrix and other parameters.
"""
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
        if norm == 0:
            break
        sys.stdout.flush()
    markov.write()
    
def freq_data(input,ef):
"""
Create frequency domain plot
"""
    p = Pulse(input,ef)
    M = p.P - np.identity(p.N)
    p.freq_plot(2000,500)#better be multiple of process number
    p.file_out.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("name", type=str,help="calculate")
    parser.add_argument("jobs", type=str)
    parser.add_argument("-m", "--matrix", action="store_true",help="calculate transition matrix file")
    parser.add_argument("-f", "--freq", action="store_true",help="plot frequency domain")    
    args = parser.parse_args()
    input_name = args.name
    jobs = args.jobs
    base_name = str.split(args.name,'.')[0]
    factor = eval(jobs)
    ef = ElectricField()
    
    for f in factor:
        ef.setfactor(f)
        print 'factor',f
        if args.matrix:
            create_matrix_file(input_name,base_name+str(f),ef)
        if args.freq:
            freq_data(base_name+str(f)+".p",ef)
