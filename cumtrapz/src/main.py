#!/usr/bin/python2
from ctypes import *
import time
import numpy as np
from scipy import linalg
import copy

libcumtrapz = CDLL("./obj/libcumtrapz.so")
libcumtrapz.matrix_power.restype = None
libcumtrapz.matrix_power.argtypes = [np.ctypeslib.ndpointer(c_double),
        c_int,
        c_int,
        np.ctypeslib.ndpointer(c_double)]
def ctype_matrix_power(A,N,p,B):
    A = A.view('float64')
    A = np.ascontiguousarray(A)
    B = B.view('float64')
    B = np.ascontiguousarray(B)
    libcumtrapz.matrix_power(A,N,p,B)
    A = A.view('complex')
    B = B.view('complex')

N = 100
A = np.random.random((N,N))+1j*np.random.random((N,N))
B = np.zeros((N,N),dtype='complex')

for p in range(10):
    print "---------------"
    print "p =",p
    t1 = time.time()
    ctype_matrix_power(A,N,p,B)
    t2 = time.time()
    print 'ctype took %0.3f ms' % ((t2-t1)*1000.0)
    t1 = time.time()
    B_np = np.linalg.matrix_power(A,p)
    t2 = time.time()
    print 'numpy took %0.3f ms' % ((t2-t1)*1000.0)
    # print B
    # print B_np
    print "B diff", np.linalg.norm(B-B_np) 

