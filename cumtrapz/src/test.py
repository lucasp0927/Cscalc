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
libcumtrapz.matrix_vector_power.restype = None
libcumtrapz.matrix_vector_power.argtypes = [np.ctypeslib.ndpointer(c_double),
        c_int,
        c_int,
        np.ctypeslib.ndpointer(c_double),
        np.ctypeslib.ndpointer(c_double)]
def ctype_matrix_power(A,N,p,B):
    A = A.view('float64')
    A = np.ascontiguousarray(A)
    B = B.view('float64')
    B = np.ascontiguousarray(B)
    libcumtrapz.matrix_power(A,N,p,B)
    A = A.view('complex')
    B = B.view('complex')

def ctype_matrix_vector_power(A,N,p,v_r,v):
    A = A.view('float64')
    A = np.ascontiguousarray(A)
    v_r = v_r.view('float64')
    v_r = np.ascontiguousarray(v_r)
    v = v.view('float64')
    v = np.ascontiguousarray(v)
    libcumtrapz.matrix_vector_power(A,N,p,v_r,v)
    A = A.view('complex')
    v_r = v_r.view('complex')
    v = v.view('complex')

N = 100
A = np.random.random((N,N))+1j*np.random.random((N,N))
v_r = np.zeros(N,dtype='complex')
B = np.zeros((N,N),dtype='complex')
v = np.random.random((N,))+1j*np.random.random((N,))
for p in range(4,20):
    print "---------------"
    print "p =",p
    t1 = time.time()
    ctype_matrix_power(A,N,2**p,B)
    v_r_c1 = np.dot(B,v)
    t2 = time.time()
    print 'ctype_power took %0.3f ms' % ((t2-t1)*1000.0)
    t1 = time.time()
#    ctype_matrix_vector_power(A,N,p,v_r,v)
    t2 = time.time()
    print 'ctype_power2 took %0.3f ms' % ((t2-t1)*1000.0)
    t1 = time.time()
    B_np = np.linalg.matrix_power(A,2**p)
    v_r_np = np.dot(B_np,v)
    t2 = time.time()
    print 'numpy took %0.3f ms' % ((t2-t1)*1000.0)
    # print B
    # print B_np
    print "v_r diff", np.linalg.norm(v_r_np-v_r_c1) 

