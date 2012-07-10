#!/usr/bin/python2
from __future__ import division
import sys
from parameter_common import Parameter
import numpy as np
# D2 line Energy level diagram
# refer to page 25 of D line data IMPORTANT the graph underneath is not accurate
#     _F=5
# _<B _F=4 C
# |   _F=3 D  32 levels
# |   _F=2 E
# |A
# |
# |  _F=4 F  9 levels
# _< _
#    _F=3 G  7 levels
# total 48 levels
# in Hz
twopi = 2*np.pi
#following freq are in rad
A=351.72571850e12*twopi
B=12.79851e6*twopi
C=263.8906e6*twopi
D=188.4885e6*twopi
E=399.7128e6*twopi
F=4.021776399375e9*twopi
G=5.170855370625e9*twopi

if __name__ == '__main__':
    filename =  sys.argv[1]
    arg_dict = {
        'l1f':(5,4,3,2),
        'l0f':(4,3),
        'B':0.0,
        'd1':0,
        'gamma':0.0,
        'egpair':(((1,5),(0,4)),((1,5),(0,3)),((1,4),(0,4)),((1,4),(0,3)),((1,3),(0,4)),((1,3),(0,3)),((1,2),(0,4)),((1,2),(0,3))),
        'omega_list':(G+A+C,G+A+B,G+A-D,G+A-E,G+F,0),
        'parameter':{'nup': [((1,3,0),(0,4,-1)),((1,3,0),(0,3,1))],
                      'e_amp': [(1,(1,)), (5,(-1,1))],
                      'sweep_profile':[0,-4E5,4E5,500],
                      },
        'filename':filename
        }
    para = Parameter(**arg_dict)
    para.write()
