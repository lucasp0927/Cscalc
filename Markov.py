#!/usr/bin/python2
from __future__ import division
import sys
from ElectricField import ElectricField
import numpy as np
from scipy import linalg,integrate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

HBAR =  1.05457148e-34
class Markov(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        file_in = sys.argv[1]
        file_out = sys.argv[2]
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
        #self.Dfunction = [np.zeros((self.N,self.N),complex) for i in range(self.EField.sample)] # time dependent part
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
        if i in self.group[0]:
            if j not in self.group[0]:
                return self.omega[i]-self.omega[j]-self.EField.carrier_freq
            else:
                return self.omega[i]-self.omega[j]
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
        for i in range(int(self.N)):
            for j in range(i,int(self.N)):
                print i,j
                self.calcDFunction(i,j)

        self.order += 1            
        self.Dfunction = self.DfunctionTemp
        
    def calcSlope(self,I,J,i):
        #ans = 0.0 + 0.0j
        # for k in range(self.N):
        #     ans += (self.T[I][k]+self.EField.envelope(self.tsample[i])*self.D[I][k])*self.Dfunction[k,J,i]
        ans =  np.sum((self.T[I][:]+self.EField.envelope(self.tsample[i])*self.D[I][:])*self.Dfunction[:,J,i])
        return ans

    def slp(self,x,t,interpolater):
        xs = interpolater.x
        ys = interpolater.y
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
        slopfunc = np.vectorize(self.calcSlope)
        slope = slopfunc(I,J,np.arange(self.smpnum))
        # interpolate slope
        slope_r = np.real(slope)
        slope_i = np.imag(slope)
        slope_r_int = interp1d(self.tsample,slope_r,kind='cubic')
        slope_i_int = interp1d(self.tsample,slope_i,kind='cubic') # spline representation
        result_r = integrate.odeint(self.slp,init,self.tsample,args=(slope_r_int,))
        result_i = integrate.odeint(self.slp,init,self.tsample,args=(slope_i_int,))
        self.DfunctionTemp[I,J,:] = np.transpose(result_r + 1.0j*result_i)

    def finalResult(self):
        time = 2*np.pi/self.EField.repetition_freq-self.EField.cutoff
        self.final = np.dot(self.Dfunction[:,:,-1],linalg.expm(self.T*time)) # check this

    def plotGraph(self):
        state = np.zeros(self.N,complex)
        for i in self.group[0]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[0])
        data = np.zeros((3,self.smpnum))
        for i in range(self.smpnum):
            state = np.dot(self.Dfunction[:,:,i],state)
            for j in range(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][i] += state[self.ij2idx(k,k)]

        plt.plot(self.tsample,data[0],self.tsample,data[1],self.tsample,data[2])
        plt.show()
    
if __name__ == '__main__':
    markov = Markov()
    markov.prepareT()
    markov.prepareD()    
    markov.zeroOrder()
    markov.addOrder()
    markov.plotGraph()
    #    markov.addOrder()


