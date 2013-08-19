#!/usr/bin/python2
from ctypes import *
import time
import numpy as np
import copy
if __name__ == '__main__':
    libmatrixMul = CDLL("./obj/libmatrixMul.so")
    
    def cuda_zdot(A,B):
        """
        Only for square matrix!!
        """
        N = A.shape[0]
        result = (c_double*(2*N**2))()
        t1 = time.time()                                        
        libmatrixMul.zdot(A.ctypes.data_as(POINTER(c_double)),B.ctypes.data_as(POINTER(c_double)),N,result)
        t2 = time.time()
        print 'function took %0.3f ms' % ((t2-t1)*1000.0)                
        C = np.frombuffer(result, dtype='complex', count=N**2)
        C.shape = (N,N)
        del result
        return C
    

    N = 10000
    def testz():
        print '-------------------------------------------'
        A = np.random.random((N,N))+1.0j*np.random.random((N,N))   
        B = np.random.random((N,N))+1.0j*np.random.random((N,N))   
        t1 = time.time()                            
        C_h = A+B
        t2 = time.time()
        print 'numpy took %0.3f ms' % ((t2-t1)*1000.0)    
        t1 = time.time()                                
        C = cuda_zdot(A,B)
        t2 = time.time()
        print 'ctype_dot took %0.3f ms' % ((t2-t1)*1000.0)        
        print 'norm',np.linalg.norm(C-C_h)
        
    for i in range(10):
        testz()
    print 'finish'

