#!/usr/bin/python2
from __future__ import division
import sys,gc
from ElectricField import ElectricField
import numpy as np
from scipy import linalg,integrate
#from scipy.interpolate import interp1d,UnivariateSpline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import pickle

HBAR =  1.05457148e-34
class Markov(object):
    """
    """
    def __init__(self, ):
        """
        """
        file_in = sys.argv[1]
        self.file_out = sys.argv[2]
        self.pp = PdfPages(self.file_out+".pdf")        
        dictf = open(file_in,'r')
        self.parameter = eval(dictf.read())
        self.omega = self.parameter['omega']
        self.group = self.parameter['level_group']
        self.dipole = self.parameter['dipole'][0]
        self.n = self.parameter['n']# number of levels
        self.N = self.n**2 # the number of independent terms in  density matrix
        self.decoherence = self.parameter['decoherence_matrix']
        self.T = [] # time independent part d rho/ dt = T rho
        self.D = [] # time independent part d rho/ dt = T rho
        self.final = np.zeros((self.N,self.N),complex) # final markov matrix
        self.EField = ElectricField()
        self.smpnum = self.EField.sample
        self.cutoff = self.EField.cutoff
        self.tsample = np.linspace(0,self.EField.cutoff,self.smpnum)
        self.Dfunction = np.empty((self.N,self.N,self.smpnum),complex)
        self.DfunctionTemp =[] #np.empty((self.N,self.N,self.smpnum),complex)
        self.order = 0
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
        for i in range(self.n):
            for j in range(i,self.n):
                for k in self.decoherence[i][j]:
                    self.T[self.ij2idx(i,j)][self.ij2idx(k[0],k[1])]+=k[2]
                    if i != j:
                        self.T[self.ij2idx(j,i)][self.ij2idx(k[1],k[0])]+=k[2]
        for i in range(self.n):
             for j in range(i+1,self.n):
                 if self.rotate_omega(i,j) > 1e11:
                     print i,j
                     print "error"
                 self.T[self.ij2idx(i,j)][self.ij2idx(i,j)]+= -1.0j*self.rotate_omega(i,j)
                 if i!= j:
                     self.T[self.ij2idx(j,i)][self.ij2idx(j,i)]+= 1.0j*self.rotate_omega(i,j)

    def prepareD(self):
        self.D = np.zeros((self.N,self.N),complex) # time independent part d rho/ dt = T rho        
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    self.D[self.ij2idx(i,j)][self.ij2idx(k,j)] += -1.0j*(self.dipole[i][k] )/ HBAR
                    self.D[self.ij2idx(i,j)][self.ij2idx(i,k)] -= -1.0j*(self.dipole[k][j] )/ HBAR

    def zeroOrder(self):
        for i in enumerate(self.tsample):
            print i[0]
            self.Dfunction[:,:,i[0]]=linalg.expm(self.T*i[1])

    def addOrder(self):
        self.DfunctionTemp = np.empty((self.N,self.N,self.smpnum),complex)        
        for i in range(self.smpnum):
            self.DfunctionTemp[...,i] = np.dot(self.T+self.D*self.EField.envelope(self.tsample[i]),self.Dfunction[...,i])
        self.Dfunction = []
        self.T = []
        self.D = []
        gc.collect()        
        self.Dfunction = integrate.cumtrapz(self.DfunctionTemp,self.tsample)
        self.DfunctionTemp = []
        gc.collect()        
        self.Dfunction = np.concatenate((np.zeros((self.N,self.N,1),complex),self.Dfunction),axis=-1)
        for i in range(self.N):
            self.Dfunction[i,i,:] += np.ones(self.smpnum,complex)
        self.order += 1
        
    def plotGraph(self,title=""):
        start = 1
        state = np.zeros(self.N,complex)
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,self.smpnum))
        for i in range(self.smpnum):
            state1 = np.dot(self.Dfunction[:,:,i],state.T)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][i] += np.real(state1[self.ij2idx(k,k)])
        plt.figure(1)                    
        fig = plt.subplot(1,1,1)
        plt.title(title)
        plt.ylim(0,1)
        plt.xlabel('time')
        plt.ylabel('population')
        for i in range(3):
            fig.plot(self.tsample,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.savefig(self.pp,format='pdf')
        plt.clf()
        #show()
        
    def write(self):
        #        txtf = open(self.file_out+'.txt','w')
        data={}
        data['T'] = self.T
        data['P'] = self.Dfunction[:,:,-1]
        data['cutoff'] = self.cutoff
        data['n'] = self.n
        data['N'] = self.N
        data['group'] = self.group
        #txtf.write(str(data))
        #txtf.close()
        pickle.dump( data, open( self.file_out+".p", "wb" ) )
        
    def tidy(self):
        # remove imag part in diagnopart
        pass
    
if __name__ == '__main__':
    markov = Markov()
    markov.prepareT()
    markov.prepareD()
    markov.zeroOrder()
    for i in range(20):
        print "order ",i
        t1 = time.time()                        
        markov.addOrder()
        t2 = time.time()
        print 'took %0.3f ms' % ((t2-t1)*1000.0)
        markov.prepareT()
        markov.prepareD()        
        markov.plotGraph(title=str(i)+"th order")
    markov.pp.close()
    markov.write()
    
