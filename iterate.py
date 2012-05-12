#!/usr/bin/python2
from ElectricField import ElectricField
from Pulse import Pulse
from Markov import Markov
import sys
import time
import numpy as np

def creat_matrix_file(input,output,ef):
    markov = Markov(input,output,ef)
    markov.prepareT()
    markov.prepareD()
    markov.zeroOrder()
    for i in xrange(50):
        print "-------------------------"
        print "order ",markov.order+1
        t1 = time.time()
        norm = markov.addOrder()
        print "difference norm %e" %norm        
        t2 = time.time()
        print 'took %0.3f ms' % ((t2-t1)*1000.0)
        markov.prepareT()
        markov.prepareD()        
        markov.plotGraph(title=str(i)+"th order")
        if norm == 0:
            break
    markov.pp.close()
    markov.write()
    
def freq_data(input,ef):
    p = Pulse(input,ef)
    M = p.P - np.identity(p.N)
    p.freq_plot(1e3,200)
    p.file_out.close()
    
if __name__ == '__main__':
    input_name = sys.argv[1]
    base_name = str.split(sys.argv[1],'.')[0]
    factor = [500,1000]
    ef = ElectricField()
    for f in factor:
        ef.setfactor(f)
        print 'power',ef.calpower()
        creat_matrix_file(input_name,base_name+str(f),ef)
        freq_data(base_name+str(f)+".p",ef)
