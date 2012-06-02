#!/usr/bin/python2
from __future__ import division
import sys
import numpy as np
from scipy import linalg
from ElectricField import ElectricField
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pickle
from multiprocessing import Process, Queue
from math import ceil
import mkl
mkl.set_num_threads(12)

class Pulse(object):
    """
    """

    def __init__(self, file_in, ef):
        """
        """
        self.filename =str.split(file_in,'.')[0]+"_freq"
        self.file_out = ''
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
        self.ef = ef
        self.process = 1
        p2 = [self.ij2idx(x,x) for x in range(self.n)]
        for i in p2:
            self.lastrow[i] = 1.0
        self.ii2idxv = np.vectorize(self.ii2idx)


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
        self.file_out.write("#sigma: "+str(self.parameter['sigma'])+"\n")
        self.file_out.write("#maxima: "+str(self.parameter['maxima'])+"\n")
        self.file_out.write("#average power: "+str(self.parameter['power'])+"\n")
        self.file_out.write("#factor: "+str(self.parameter['factor'])+"\n")
        self.file_out.write("#{:-<80}\n".format(''))
        self.file_out.write("#{:<20} {:<20} {:<20} {:<20}\n".format("rep_freq(Hz)","population(0)","population(1)","population(2)"))

    def time_plot(self,num,step):
        self.file_out=open(self.filename+".dat","w")
        print "plot time domain, total",num,"points."
        rep = 1.0/(self.ef.repetition_freq/(2*np.pi))
        r_t = rep - self.cutoff
        M = np.dot(linalg.expm(self.T*r_t),self.P)                
        start = 1
        state = np.zeros(self.N,complex)
        time_arr = []
        time = 0.0
        stepM = np.linalg.matrix_power(M, step)        
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        data = np.zeros((3,num/step))
        for i in xrange(0,num,step):
            for j in xrange(3):      # make this more elegent
                for k in self.group[j]:
                    data[j][i/step] += np.real(state[self.ij2idx(k,k)]) # i+i is slightly faster than i*2
            # sys.stdout.write('%s\r' % i)
            print i
            sys.stdout.flush()
            state = np.dot(stepM,state.T)
            time += rep*step
            time_arr.append(time)

        for t in enumerate(time_arr):
            self.file_out.write('{0:<20} {1[0]:<20} {1[1]:<20} {1[2]:<20}\n'.format(t[1],data[:,t[0]])) # output to log file
        self.file_out.close()

        plt.figure(1)
        fig = plt.subplot(1,1,1)
        plt.title("test")
        plt.ylim(-0.1,1.1)
        plt.xlabel('time')
        plt.ylabel('population')
        for i in xrange(3):
            fig.plot(time_arr,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.savefig(self.filename+"_time")
        plt.clf()
        # plt.show()

    # def correct(self,arr):
    #     # TODO: write in better numpy way
    #     for i in xrange(self.n):
    #         if (real(arr[self.ij2idx(i,i)])<0.0 or real(arr[self.ij2idx(i,i)])>1.0):
    #             return False
    #     else:
    #         return True

    def freq_plot(self,freq_range,number):
        self.file_out=open(self.filename+".dat","w")
        self.dump_header()        
        print "plot frequency domain, total",number,"points."
        rf = self.ef.repetition_freq/(2*np.pi)
        repf = np.linspace(rf-freq_range,rf+freq_range,number)
        rept = 1.0/repf
        data = np.zeros((3,number))

        def chunks(l, n):
            n = int(ceil(float(len(l))/float(n)))
            l = list(enumerate(l))
            return [l[i:i+n] for i in range(0, len(l), n)]
        q = Queue()
        process = []
        for i in range(self.process):
            process.append(Process(target=self.plot_worker, args=(q,chunks(rept,self.process)[i])))
        for p in process:
            p.start()
        n = 0
        while n < len(rept)*len(self.group): # better way?
            try:
                d = q.get()
                data[d[0],d[1]] = d[2]
                n += 1
                #sys.stdout.write('%s\r' % int(n/3))
                if n%30 == 0:
                    print int(n/3)
                sys.stdout.flush()
            except:
                pass
        for p in process:
            p.join()

        for rf in enumerate(repf):
            self.file_out.write('{0:<20} {1[0]:<20} {1[1]:<20} {1[2]:<20}\n'.format(rf[1],data[:,rf[0]])) # output to log file

        plt.figure(1)
        fig = plt.subplot(1,1,1)
        plt.title("population vs repetition rate")
        plt.xlabel('repetition rate(Hz)')
        plt.ylabel('population')
        for i in xrange(0,1):           # plot only highest level
            fig.plot(repf,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.savefig(self.filename)

        for i in xrange(1,3):           # plot only highest level
            fig.plot(repf,data[i],label=str(i))
        handles, labels = fig.get_legend_handles_labels()
        fig.legend(handles[::-1], labels[::-1])
        plt.savefig(self.filename+"_all")
        plt.clf()
        
        self.file_out.close()
        
    # def matrix_vector_power(self,M,v,n):
    #     ###
    #     #find closest log 2
    #     ###
    #     new_n = np.ceil(np.log2(n))
    #     part = np.floor(np.log2(self.N/np.log(2)))
    #     partM = np.linalg.matrix_power(M,int(2**(new_n-part)))
    #     for i in range(int(2**part)):
    #         v = np.dot(partM,v.T)
    #     return v

    def plot_worker(self,q,job):
        state = np.zeros(self.N,complex)
        start = 1
        for i in self.group[start]:
            state[self.ij2idx(i,i)] = 1.0/len(self.group[start])
        for t in job:
            M = np.dot(linalg.expm(self.T*(t[1]-self.cutoff)),self.P)
            #state1 = self.matrix_vector_power(M,state.T,2**26)
            M = np.linalg.matrix_power(M,100000000)
            state1 = np.dot(M,state.T)            
            
            # M = M - np.identity(self.N)
            # M[-1,...] = self.lastrow
            # state1 = linalg.solve(M,self.con)
            for g in enumerate(self.group):
                q.put([g[0],t[0],np.sum(np.real(state1[self.ii2idxv(g[1][:])]))])

       
if __name__ == '__main__':
    ef = ElectricField()
    p = Pulse(sys.argv[1],ef)
    M = p.P - np.identity(p.N)
    p.time_plot(100000000,100000)
    #p.freq_plot(1e6,100)
    #p.freq_plot(1e-9,2e-9,10000,20000)

