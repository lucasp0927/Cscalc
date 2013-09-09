#!/usr/bin/python2
from __future__ import division
import sys,gc
from ElectricField import ElectricField
import numpy as np
from scipy import linalg,integrate
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import time
import pickle
from ctypes import *
from copy import copy
import mkl
mkl.set_num_threads(12)

#from scipy.interpolate import interp1d,UnivariateSpline
HBAR = 1.05457148e-34
#HBAR = 1.0

class SuperMatrix(object):
    """
    """
    def __init__(self,file_in,file_out,ef):
        """
        """
        self.libcumtrapz = CDLL("./cumtrapz/src/obj/libcumtrapz.so")
        self.initlibcumtrapz()
        self.file_out = file_out
        #self.pp = PdfPages(self.file_out+".pdf")
        dictf = open(file_in,'r')
        self.parameter = eval(dictf.read())
        self.omega = self.parameter['omega']
        self.gamma = self.parameter['gamma']
        self.group = self.parameter['level_group']
        self.dipole = self.parameter['dipole'][0]
        self.n = self.parameter['n']# number of levels
        self.N = self.n**2 # the number of independent terms in  density matrix
        self.decoherence = self.parameter['decoherence_matrix']
        self.T = [] # time independent part d rho/ dt = T rho
        self.D = [] # time dependent part d rho/ dt = T rho
        self.final = np.zeros((self.N,self.N),complex) # final markov matrix
        self.EField = ef
        self.smpnum = self.EField.sample
        self.cutoff = self.EField.cutoff
        self.tsample = np.linspace(0,self.EField.cutoff,self.smpnum)
        env_vec = np.vectorize(self.EField.envelope)
        self.envelope = env_vec(self.tsample)
        self.dt = self.tsample[1]-self.tsample[0]
        self.Dfunction = np.empty((self.smpnum,self.N,self.N),complex)
        self.DfunctionTemp =[] #np.empty((self.N,self.N,self.smpnum),complex)
        self.order = 0
        self.parameter = {}
        dictf.close()

    def ij2idx(self,i,j):
        """
        0 1 2
        3 4 5
        6 7 8
        """
        idx = self.n*i+j
        return idx

    def rotate_omega(self,i,j):
        if j<i:
            print "error"
        if i in self.group[0]:
            if j in self.group[0]:
                return self.omega[i]-self.omega[j]
            else:
                return self.omega[i]-self.omega[j]-self.EField.carrier_freq
        else:
            return self.omega[i]-self.omega[j]

    def prepareT(self):
        self.T = np.zeros((self.N,self.N),complex) # time independent part d rho/ dt = T rho
        for i in xrange(self.n):
            for j in xrange(i,self.n):
                for k in self.decoherence[i][j]:
                    self.T[self.ij2idx(i,j)][self.ij2idx(k[0],k[1])]+=k[2]
                    if i != j:
                        self.T[self.ij2idx(j,i)][self.ij2idx(k[1],k[0])]+=k[2]
        for i in xrange(self.n):
             for j in xrange(i+1,self.n):
                 if self.rotate_omega(i,j) > 1e11:
                     print i,j
                     print "error"
                 self.T[self.ij2idx(i,j)][self.ij2idx(i,j)]+= -1.0j*self.rotate_omega(i,j)
                 if i!= j:
                     self.T[self.ij2idx(j,i)][self.ij2idx(j,i)]+= 1.0j*self.rotate_omega(i,j)

    def prepareD(self):
        self.D = np.zeros((self.N,self.N),complex) # time dependent part d rho/ dt = T rho
        for i in xrange(self.n):
            for j in xrange(self.n):
                for k in xrange(self.n):
                    self.D[self.ij2idx(i,j)][self.ij2idx(k,j)] += -1.0j*(self.dipole[i][k] )/ HBAR
                    self.D[self.ij2idx(i,j)][self.ij2idx(i,k)] -= -1.0j*(self.dipole[k][j] )/ HBAR

    def zeroOrder(self):
        print "zero order"
        for i in enumerate(self.tsample):
            # sys.stdout.write('%s\r' % i[0])
            print i[0]
            sys.stdout.flush()
            self.Dfunction[i[0],:,:]=linalg.expm(self.T*i[1],15)

    # def addOrder(self):
    #     last = self.Dfunction[...,-1]
    #     self.DfunctionTemp = np.empty((self.smpnum,self.N,self.N),complex)
    #     for i in xrange(self.smpnum):
    #         self.DfunctionTemp[i,...] = np.dot(self.T+self.D*self.EField.envelope(self.tsample[i]),self.Dfunction[i,...])
    #     # self.Dfunction = []
    #     # self.T = []
    #     # self.D = []
    #     del self.Dfunction
    #     del self.T
    #     del self.D
    #     gc.collect()                    # clean up memory
    #     #should rewrite cumtrapz with C to use less memory
    #     #self.Dfunction = np.zeros((self.N,self.N,self.smpnum),complex)
    #     self.Dfunction = integrate.cumtrapz(self.DfunctionTemp,self.tsample)
    #     del self.DfunctionTemp
    #     gc.collect()
    #     self.Dfunction = np.concatenate((np.zeros((self.N,self.N,1),complex),self.Dfunction),axis=-1)
    #     for i in xrange(self.N):
    #         self.Dfunction[i,i,:] += np.ones(self.smpnum,complex)
    #     self.order += 1
    #     now = self.Dfunction[...,-1]
    #     return linalg.norm(now-last)
    # #print "difference norm %f" %linalg.norm(now-last)

    def addOrder2(self):
        last = copy(self.Dfunction[-1,...])
        #self.DfunctionTemp = np.empty((self.smpnum,self.N,self.N),complex)
        #self.DfunctionTemp = np.empty((self.N,self.N,self.smpnum),complex)

        # for i in xrange(self.smpnum):
        #      tmp = np.dot(self.T+self.D*self.envelope[i],self.Dfunction[i,...])
        #      self.Dfunction[i,...] = copy(tmp)

        #del self.DfunctionTemp
        # del self.T
        # del self.D
        #gc.collect()                    # clean up memory
        self.ctype_addorder()
        # self.ctype_cumtrapz()

        # self.Dfunction = integrate.cumtrapz(self.DfunctionTemp,self.tsample)
        # self.Dfunction = np.concatenate((np.zeros((self.N,self.N,1),complex),self.Dfunction),axis=-1)
        # for i in xrange(self.N):
        #     self.Dfunction[i,i,:] += np.ones(self.smpnum,complex)
        self.order += 1
        now = self.Dfunction[-1,...]
        return linalg.norm(now-last)

    def write(self):
        data={}
        data['T'] = self.T
        data['P'] = self.Dfunction[-1,:,:]
        data['cutoff'] = self.cutoff
        data['n'] = self.n
        data['N'] = self.N
        data['group'] = self.group
        data['power'] = self.EField.calpower()
        data['sigma'] = self.EField.sigma
        data['maxima'] = self.EField.maxima
        data['factor'] = self.EField.factor
        data['gamma'] = self.gamma
        pickle.dump( data, open( self.file_out+".p", "wb" ) )

    def initlibcumtrapz(self):
        self.libcumtrapz.cumtrapz.restype = None
        self.libcumtrapz.cumtrapz.argtypes = [np.ctypeslib.ndpointer(c_double),
                c_double,
                c_int,
                c_int]

        self.libcumtrapz.addorder.restype = None
        self.libcumtrapz.addorder.argtypes = [np.ctypeslib.ndpointer(c_double),
                np.ctypeslib.ndpointer(c_double),
                np.ctypeslib.ndpointer(c_double),
                np.ctypeslib.ndpointer(c_double),
                c_int,
                c_int,
                c_double]

    def ctype_cumtrapz(self):
        #result = (c_double*(2*N**2))()
        self.Dfunction = self.Dfunction.view('float64')
        self.Dfunction = np.ascontiguousarray(self.Dfunction)
        self.libcumtrapz.cumtrapz(self.Dfunction,self.dt,self.N,self.smpnum)
        self.Dfunction = self.Dfunction.view('complex')

    def ctype_addorder(self):
        #result = (c_double*(2*N**2))()

        self.Dfunction = self.Dfunction.view('float64')
        self.Dfunction = np.ascontiguousarray(self.Dfunction)
        self.T = self.T.view('float64')
        self.T = np.ascontiguousarray(self.T)
        self.D = self.D.view('float64')
        self.D = np.ascontiguousarray(self.D)
        self.envelope = np.ascontiguousarray(self.envelope,dtype = 'float64')

        self.libcumtrapz.addorder(self.T,self.D,self.envelope,self.Dfunction,self.N,self.smpnum,self.dt)

        self.Dfunction = self.Dfunction.view('complex')
        self.T = self.T.view('complex')
        self.D = self.D.view('complex')


if __name__ == '__main__':
    ef = ElectricField()
    markov = SuperMatrix(sys.argv[1],sys.argv[2],ef)
    # markov.prepareT()
    # markov.prepareD()
    # markov.zeroOrder()
    markov.ctype_test()
    # for i in xrange(50):
    #     print "-------------------------"
    #     print "order ",markov.order+1
    #     t1 = time.time()
    #     norm = markov.addOrder2()
    #     print "difference norm %e" %norm
    #     t2 = time.time()
    #     print 'took %0.3f ms' % ((t2-t1)*1000.0)
    #     markov.prepareT()
    #     markov.prepareD()
    #     #markov.plotGraph(title=str(i)+"th order")
    #     if norm == 0:
    #         break
    # markov.pp.close()
    # markov.write()
    # print "See the output PDF file to check if purturbation converge."
