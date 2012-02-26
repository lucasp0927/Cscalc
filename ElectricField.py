#!/usr/bin/python2
import numpy as np
#import matplotlib.pyplot as plt
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, integrate

class ElectricField(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        self.carrier_freq = (335.116048807e12-4.021776399375e9)*2*np.pi # carrier frequency in rad
        self.repetition_freq = 100e6*2*np.pi # repetition frequency in rad
        self.cutoff = 2e-13#where electric field start consider to be zero
        self.sample = 10
        
    def envelope(self,t):
        """
        the envelope function must start from t = 0
        """
        sigma = 20e-15
        maxima = 1e8
        return maxima*np.exp(-(t-sigma*5)**2/(2*sigma**2))

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

    def check(self,):
        print "rabi frequency area: "
        print 4e5*integrate.quad(ef.envelope,0,ef.cutoff)[0] #4e5 is dipole moment / hbar
        print "average power:"
        print str(8.85e-12/2*integrate.quad(lambda x:ef.envelope(x)**2,0,ef.cutoff)[0]/(2*np.pi/ef.repetition_freq))+"w/m^2"
        Ts = ef.cutoff/500.0
        x = np.arange(0, ef.cutoff, Ts)
        env_vec = np.vectorize(ef.envelope)
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

