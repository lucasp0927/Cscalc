import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

class ElectricField(object):
    """
    """
    
    def __init__(self, ):
        """
        """
        self.carrier_freq = (335.116048807e12-4.021776399375e9)*2*np.pi # carrier frequency in rad
        self.repetition_freq = 100e6*2*np.pi # repetition frequency in rad
        self.cutoff = 2e-13#where electric field start consider to be zero
    def envelope(self,t):
        """
        the envelope function must start from t = 0
        """
        sigma = 20e-15
        maxima = 1e6
        return maxima*np.exp(-(t-sigma*5)**2/(2*sigma**2))

    def check(self,):
        print "rabi frequency area: "
        print 4e5*integrate.quad(ef.envelope,0,ef.cutoff)[0] #4e5 is dipole moment / hbar
        print "average power:"
        print str(8.85e-12/2*integrate.quad(lambda x:ef.envelope(x)**2,0,ef.cutoff)[0]/(2*np.pi/ef.repetition_freq))+"w/m^2"
        x = np.arange(0, ef.cutoff, ef.cutoff/500.0)
        env_vec = np.vectorize(ef.envelope)
        y = env_vec(x)
        plt.plot(x, y)
        plt.show()
    
if __name__ == '__main__':
    ef = ElectricField()
    ef.check()









