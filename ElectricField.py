#!/usr/bin/python2
import numpy as np
#import matplotlib.pyplot as plt
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, integrate

class ElectricField(object):
    """
    """
    
    def __init__(self,):
        """
        """
        #self.carrier_freq = (1000000000000000)*2*np.pi # carrier frequency in rad
        self.carrier_freq = 2105628723506709.8 #(335.116048807e12+5.170855370625e9)*2*np.pi # carrier frequency in rad
        #self.carrier_freq = 2105631933334755.8 # carrier frequency in rad
        self.repetition_freq =  57759008871.57628/10
        #self.repetition_freq = 91.016156e6*2*np.pi # repetition frequency in rad
        #self.repetition_freq = 91.9257518e6*2*np.pi # repetition frequency in rad        
        self.factor = 1
        self.setfactor(self.factor)        
        self.cutoff = self.sigma*10.0#where electric field start consider to be zero
        self.sample = 100
        
    def setfactor(self,f):
        self.factor = float(f)
        self.maxima = 1e5*np.sqrt(self.factor)
        self.sigma = 20e-15

    def envelope(self,t):
        """
        the envelope function must start from t = 0
        """
        return self.maxima*np.exp(-(t-self.sigma*5)**2/(2*self.sigma**2))

    def plotSpectrum(self,y,Ts):
        """
        Plots a Single-Sided Amplitude Spectrum of y(t)
        """
        Fs = 1/Ts
        n = len(y) # length of the signal
        k = np.arange(n)
        T = n/Fs
        frq = k/T # two sides frequency range
        frq = frq[range(n/2)] # one side frequency range
    
        Y = fft(y)/n # fft computing and normalization
        Y = Y[range(n/2)]
        
        plot(frq,abs(Y),'pr') # plotting the spectrum
        xlabel('Freq (Hz)')
        ylabel('|Y(freq)|')

    def calpower(self,):
        return '{:e}'.format(8.85e-12/2*integrate.romberg(lambda x:(self.envelope(x)*np.sin(self.carrier_freq*x))**2,0,self.cutoff,divmax=20)/(2*np.pi/self.repetition_freq))+"w/m^2"
    
    def check(self,):
        print "rabi frequency area: "
        print 4e5*integrate.quad(self.envelope,0,self.cutoff)[0] #4e5 is dipole moment / hbar
        print "average power:"
        print self.calpower()
        print str(8.85e-12/4*integrate.romberg(lambda x:(self.envelope(x))**2,0,self.cutoff)/(2*np.pi/self.repetition_freq))+"w/m^2"        
        Ts = self.cutoff/500.0
        x = np.arange(0, self.cutoff, Ts)
        env_vec = np.vectorize(self.envelope)
        y = env_vec(x)
        print "please check if rotating wave approximation is valid."
        # plot envelope and its spectrum
        subplot(2,1,1)
        plot(x,y)
        xlabel('Time')
        ylabel('Amplitude')
        subplot(2,1,2)
        self.plotSpectrum(y,Ts)
        show()

if __name__ == '__main__':
    ef = ElectricField()
    ef.check()

