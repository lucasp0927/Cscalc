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
        libmatrixMul.zdot(A.ctypes.data_as(POINTER(c_double)),B.ctypes.data_as(POINTER(c_double)),N,result)
        C = np.frombuffer(result, dtype='complex', count=N**2)
        C.shape = (N,N)
        del result
        return C
    
    def cuda_cdot(A,B):
        """
        Only for square matrix!!
        """
        N = A.shape[0]
        result = (c_float*(2*N**2))()
        libmatrixMul.cdot(A.ctypes.data_as(POINTER(c_float)),B.ctypes.data_as(POINTER(c_float)),N,result)
        C = np.frombuffer(result, dtype='complex64', count=N**2)
        C.shape = (N,N)
        del result
        return C
    #libmatrixMul.runTest()
    #    tenDouble = c_double * 10
    #    data = tenDouble(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0)
    libmatrixMul.deviceVerify()
    N = 1024
    def testz():
        print '-------------------------------------------'
        A = np.random.random((N,N))+1.0j*np.random.random((N,N))   
        B = np.random.random((N,N))+1.0j*np.random.random((N,N))   
        t1 = time.time()                            
        C_h = np.dot(A,B)
        t2 = time.time()
        print 'np.dot took %0.3f ms' % ((t2-t1)*1000.0)    
        t1 = time.time()                                
        C = cuda_zdot(A,B)
        t2 = time.time()
        print 'cuda_dot took %0.3f ms' % ((t2-t1)*1000.0)        
        print 'norm',np.linalg.norm(C-C_h)
        
    def testc():
        print '-------------------------------------------'
        A = np.random.random((N,N))+1.0j*np.random.random((N,N))
        B = np.random.random((N,N))+1.0j*np.random.random((N,N))
        A = np.complex64(A)
        B = np.complex64(B)        
        t1 = time.time()                            
        C_h = np.dot(A,B)
        t2 = time.time()
        print 'np.dot took %0.3f ms' % ((t2-t1)*1000.0)    
        t1 = time.time()                                
        C = cuda_cdot(A,B)
        t2 = time.time()
        print 'cuda_dot took %0.3f ms' % ((t2-t1)*1000.0)        
        print 'norm',np.linalg.norm(C-C_h)        

    for i in range(10):
        testc()

    print 'finish'

