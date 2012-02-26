#!/usr/bin/python2
from __future__ import division
# threading
import sys
import threading
from itertools import izip, count

from ElectricField import ElectricField
import numpy as np
from scipy import linalg,integrate
from scipy.interpolate import interp1d,UnivariateSpline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time

HBAR =  1.05457148e-34
class Markov(object):
    """
    """
    def __init__(self, ):
        """
        """
        file_in = sys.argv[1]
        file_out = sys.argv[2]
        self.pp = PdfPages(file_out)        
        dictf = open(file_in,'r')
        self.parameter = eval(dictf.read())
        self.omega = self.parameter['omega']
        self.group = self.parameter['level_group']
        self.dipole = self.parameter['dipole'][0]
        self.n = self.parameter['n']# number of levels
        self.N = self.n**2 # the number of independent terms in  density matrix
        self.decoherence = self.parameter['decoherence_matrix']
        self.T = np.zeros((self.N,self.N),complex) # time independent part d rho/ dt = T rho
        self.D = np.zeros((self.N,self.N),complex) # time independent part d rho/ dt = T rho
        self.final = np.zeros((self.N,self.N),complex) # final markov matrix
        self.EField = ElectricField()
        self.smpnum = self.EField.sample
        self.cutoff = self.EField.cutoff
        self.tsample = np.linspace(0,self.EField.cutoff,self.smpnum)
        self.Dfunction = np.empty((self.N,self.N,self.smpnum),complex)
        self.DfunctionTemp = np.empty((self.N,self.N,self.smpnum),complex)
        self.order = 0
        dictf.close()

    def ij2idx(self,i,j):
        """
        0 1 2
        3 4 5
        6 7 8
        """
        #idx = (i*(2*self.n-i+1))/2+(j-i)
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
        for i in range(self.n):
            for j in range(i,self.n):
                for k in self.decoherence[i][j]:
                    self.T[self.ij2idx(i,j)][self.ij2idx(k[0],k[1])]+=k[2]
                    self.T[self.ij2idx(j,i)][self.ij2idx(k[1],k[0])]+=k[2]
        for i in range(self.n):
             for j in range(i+1,self.n):
                 if self.rotate_omega(i,j) > 1e11:
                     print i,j
                     print "error"
                 self.T[self.ij2idx(i,j)][self.ij2idx(i,j)]+= -1.0j*self.rotate_omega(i,j)
                 self.T[self.ij2idx(j,i)][self.ij2idx(j,i)]+= 1.0j*self.rotate_omega(i,j)

    def prepareD(self):
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    self.D[self.ij2idx(i,j)][self.ij2idx(k,j)] += -1.0j*(self.dipole[i][k] )/ HBAR
                    self.D[self.ij2idx(i,j)][self.ij2idx(i,k)] -= -1.0j*(self.dipole[k][j] )/ HBAR


    def zeroOrder(self):
        for i in enumerate(self.tsample):
            print i[0]
            self.Dfunction[:,:,i[0]]=linalg.expm(markov.T*i[1])

    def addOrder(self):
        for i in range(self.N):
            print i
            t1 = time.time()                
            for j in range(self.N):
                self.calcDFunction(i,j)
            t2 = time.time()
            print 'took %0.3f ms' % ((t2-t1)*1000.0)
        # def calcRow(i):
        #     for j in range(self.N):
        #           self.calcDFunction(i,j)
        # foreach(calcRow,range(self.N))
        self.order += 1            
        self.Dfunction = self.DfunctionTemp.copy()

    def calcSlope(self,I,J,i):
        ans =  np.sum((self.T[I][:]+self.EField.envelope(self.tsample[i])*self.D[I][:])*self.Dfunction[:,J,i])
        #ans1 = np.sum((np.outer(self.T[I][:],np.ones(self.smpnum,complex))+np.outer(self.D[I][:],self.EField.envelope(self.tsample)))*self.Dfunction[:,J,:],axis=0)
        return ans

    def slp(self,x,t,interpolater,xs,ys):
        # xs = interpolater.x
        # ys = interpolater.y
        if t < xs[0]:
            return ys[0]+(t-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif t > xs[-1]:
            return ys[-1]+(t-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolater(t)
        
    def calcDFunction(self,I,J):
        if I==J:
            init = 1
        else:
            init = 0
        # slopfunc = np.vectorize(self.calcSlope)
        # slope = slopfunc(I,J,np.arange(self.smpnum))

        slope = np.sum((np.outer(self.T[I][:],np.ones(self.smpnum,complex))+np.outer(self.D[I][:],self.EField.envelope(self.tsample)))*self.Dfunction[:,J,:],axis=0) # 75% faster

        # interpolate slope
        slope_r = np.real(slope)
        slope_i = np.imag(slope)
        # slope_r_int = interp1d(self.tsample,slope_r,kind='cubic') 
        # slope_i_int = interp1d(self.tsample,slope_i,kind='cubic') # spline representation
        slope_r_int = UnivariateSpline(self.tsample,slope_r) # much faster than interp1d
        slope_i_int = UnivariateSpline(self.tsample,slope_i)
        result_r = integrate.odeint(self.slp,init,self.tsample,args=(slope_r_int,self.tsample,slope_r))
        result_i = integrate.odeint(self.slp,init,self.tsample,args=(slope_i_int,self.tsample,slope_i))
        self.DfunctionTemp[I,J,:] = np.transpose(result_r + 1.0j*result_i)
        # plt.plot(self.tsample,result_r,self.tsample,result_i)
        # plt.show()
        
    def finalResult(self):
        time = 2*np.pi/self.EField.repetition_freq-self.EField.cutoff
        self.final = np.dot(self.Dfunction[:,:,-1],linalg.expm(self.T*time)) # check this

    def plotGraph(self,title=""):
        state = np.zeros(self.N,complex)
        for i in self.group[1]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[0])
        data = np.zeros((3,self.smpnum))
        for i in range(self.smpnum):
            state1 = np.dot(self.Dfunction[:,:,i],state.T)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][i] += state1[self.ij2idx(k,k)]
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

def foreach(f,l,threads=8,return_=False):
    """
    Apply f to each element of l, in parallel
    """

    if threads>1:
        iteratorlock = threading.Lock()
        exceptions = []
        if return_:
            n = 0
            d = {}
            i = izip(count(),l.__iter__())
        else:
            i = l.__iter__()


        def runall():
            while True:
                iteratorlock.acquire()
                try:
                    try:
                        if exceptions:
                            return
                        v = i.next()
                    finally:
                        iteratorlock.release()
                except StopIteration:
                    return
                try:
                    if return_:
                        n,x = v
                        d[n] = f(x)
                    else:
                        f(v)
                except:
                    e = sys.exc_info()
                    iteratorlock.acquire()
                    try:
                        exceptions.append(e)
                    finally:
                        iteratorlock.release()
        
        threadlist = [threading.Thread(target=runall) for j in xrange(threads)]
        for t in threadlist:
            t.start()
        for t in threadlist:
            t.join()
        if exceptions:
            a, b, c = exceptions[0]
            raise a, b, c
        if return_:
            r = d.items()
            r.sort()
            return [v for (n,v) in r]
    else:
        if return_:
            return [f(v) for v in l]
        else:
            for v in l:
                f(v)
            return
        
if __name__ == '__main__':
    markov = Markov()
    markov.prepareT()
    markov.prepareD()
    markov.zeroOrder()
    # #markov.calcDFunction(0,0)
    for i in range(6):
        markov.addOrder()
        markov.plotGraph(title=str(i)+"th order")
    markov.pp.close()


