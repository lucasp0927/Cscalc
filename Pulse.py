#!/usr/bin/python2
from __future__ import division
import sys
import numpy as np
from scipy import linalg
from ElectricField import ElectricField
import matplotlib.pyplot as plt
import pickle

class Pulse(object):
    """
    """

    def __init__(self, file_in):
        """
        """
        self.file_out=open(str.split(file_in,'.')[0]+"_freq.dat","w")
        self.parameter = pickle.load( open( file_in, "rb" ) )
        self.T = self.parameter['T']
        self.P = self.parameter['P']
        self.n = self.parameter['n']
        self.N = self.parameter['N']
        self.group = self.parameter['group']
        self.cutoff = self.parameter['cutoff']
        self.lastrow = np.zeros(self.N,complex)
        self.con = np.zeros(self.N,complex)
        self.con[-1] = 1.0
        self.ef = ElectricField()
        p2 = [self.ij2idx(x,x) for x in range(self.n)]
        for i in p2:
            self.lastrow[i] = 1.0
        self.ii2idxv = np.vectorize(self.ii2idx)
        self.dump_header()

    def ij2idx(self,i,j):
        """
        0 1 2
        3 4 5
        6 7 8
        """
        #idx = (i*(2*self.n-i+1))/2+(j-i)
        idx = self.n*i+j
        return idx

    def ii2idx(self,i):
        return (self.n+1)*i

    def dump_header(self,):
        self.file_out.write("#carrier freq: "+str(self.ef.carrier_freq)+" rad\n")
        self.file_out.write("#centeral repetitoin freq: "+str(self.ef.repetition_freq)+" rad\n")
        self.file_out.write("#sigma: "+str(self.ef.sigma)+"\n")
        self.file_out.write("#maxima: "+str(self.ef.maxima)+"\n")
        self.file_out.write("#average power: "+str(self.ef.calpower())+"\n")
        self.file_out.write("#factor: "+str(self.ef.factor)+"\n")        
        self.file_out.write("#{:-<80}\n".format(''))
        self.file_out.write("#{:<20} {:<20} {:<20} {:<20}\n".format("rep_freq(Hz)","population(0)","population(1)","population(2)"))

    def time_plot(self,rep,num):
        print "plot time domain, total",num,"points."
        r_t = rep - self.cutoff
        expT = linalg.expm(self.T*r_t)
        start = 1
        state = np.zeros(self.N,complex)
        time_arr = []
        time = 0.0
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,num+num))
        for i in xrange(num):
            sys.stdout.write('%s\r' % i)
            sys.stdout.flush()
            state = np.dot(self.P,state.T)
            time += self.cutoff
            time_arr.append(time)
            for j in xrange(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][i+i] += np.real(state[self.ij2idx(k,k)]) # i+i is slightly faster than i*2
            state = np.dot(expT,state.T)
            #print state
            time += r_t
            time_arr.append(time)
            for j in xrange(3):           # make this more elegent
                for k in self.group[j]:
                    data[j][i+i+1] += np.real(state[self.ij2idx(k,k)])
        plt.figure(1)
        fig = plt.subplot(1,1,1)
        plt.title("test")
        plt.ylim(0,1)
        plt.xlabel('time')
        plt.ylabel('population')
        for i in xrange(3):
            fig.plot(time_arr,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.show()

    # def correct(self,arr):
    #     # TODO: write in better numpy way
    #     for i in xrange(self.n):
    #         if (real(arr[self.ij2idx(i,i)])<0.0 or real(arr[self.ij2idx(i,i)])>1.0):
    #             return False
    #     else:
    #         return True

    def freq_plot(self,freq_range,number):#,pnum):
        print "plot frequency domain, total",number,"points."
        rf = self.ef.repetition_freq/(2*np.pi)
        repf = np.linspace(rf-freq_range,rf+freq_range,number)
        rept = 1.0/repf
        start = 1
        state = np.zeros(self.N,complex)
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,number))
        for t in enumerate(rept):
            sys.stdout.write('%s\r' % t[0])
            sys.stdout.flush()
            #print t[0]
            M = np.dot(linalg.expm(self.T*(t[1]-self.cutoff)),self.P)

            M = np.linalg.matrix_power(M,200000)
            state1 = np.dot(M,state.T)

            # M = M - np.identity(p.N)
            # M[-1,...] = self.lastrow
            # state1 = linalg.solve(M,self.con)
            for g in enumerate(self.group):
                data[g[0],t[0]] = np.sum(np.real(state1[self.ii2idxv(g[1][:])]))

        for rf in enumerate(repf):
            self.file_out.write('{0:<20} {1[0]:<20} {1[1]:<20} {1[2]:<20}\n'.format(rf[1],data[:,rf[0]]))
        # plt.figure(1)
        # fig = plt.subplot(1,1,1)
        # plt.title("population vs repetition rate")
        # #plt.ylim(-0.1,1.1)
        # plt.xlabel('repetition rate(Hz)')
        # plt.ylabel('population')
        # for i in xrange(0,3):
        #     fig.plot(repf,data[i],label=str(i))
            # if i == 0:
            #     mean = np.min(data[0]) + (np.max(data[0])-np.min(data[0]))/2.0
            #     for f in xrange(len(repf)-1):
            #         if data[0][f] == np.min(data[0]):
            #             print "min is at:",repf[f],"Hz"
            #         if data[0][f+1]<mean and data[0][f]>mean:
            #             low = repf[f]
            #         if data[0][f+1]>mean and data[0][f]<mean:
            #             high = repf[f]
            #     print(high - low)
        # handles, labels = fig.get_legend_handles_labels()
        # fig.legend(handles[::-1], labels[::-1])
        # plt.show()

if __name__ == '__main__':
    p = Pulse(sys.argv[1])
    M = p.P - np.identity(p.N)
    #p.time_plot(1.67e-8,100)
    p.freq_plot(2e3,100)
    #p.freq_plot(1e-9,2e-9,10000,20000)
    p.file_out.close()
