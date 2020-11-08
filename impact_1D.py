#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:02:35 2019

@author: dingding
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt
import scipy.io as scio
import random
import os.path
from scipy import signal
from scipy import integrate
import scipy.stats as ss


# Choose Parameter
delta_t = 0.1     # single vibration duration
t = 3    # exposure time
T = 3600    # total observation time
N = int(t/delta_t)
n =  int(T/t)   # test counts
d = 5    # vibration change range
a = 4.0    # intra pixel 2a*2
sigma = 1.0/4    # vibration distribution std


# Gaussian Distribution Function
def normal_distribution(x, mean, sigma):
    y = np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)
    return y

def square(x):
    for i in range(0,len(x)):        
        if x[i]>=-0.5 and x[i]<=0.5:
            y = 1
        else:
            y = 0
    return y
#    return 1

def trapezoid_distribution(x, mean, sigma):
    a = mean-3.0*sigma
    b = mean+3.0*sigma
    c = mean-1.0*sigma
    d = mean+1.0*sigma
    h = 2.0/((b-a)+(d-c))
    f = np.where((x>=c)&(x<=d), h, np.where((x>=a)&(x<c), h/(c-a)*(x-a), np.where((x>d)&(x<=b), h/(d-b)*(x-b), 0)))
    return f

# Point Spread Function (PSF) as Gaussian Distribution
y0 = np.linspace(-a, a, num = int(90*a))
#x0 = normal_distribution(y0, 0, 1)/(45)
#x0 = square(y0)
x0 = trapezoid_distribution(y0, 0, 1)/(45)
#print x0.shape
#print np.sum(x0)


## Normalize Matrix
#def normalize(data):
#    m = np.mean(data)
#    mx = np.max(data)
#    mn = np.min(data)
#    return [(float(u) - m) / (mx - mn) for u in data]


# Intra Pixel Matrix
names = locals()
fn = scio.loadmat('450nm_shift.mat')
orimat = fn['pixel22']
#fig1 = plt.figure()
#plt.imshow(orimat)
#plt.savefig('orimat.png')
normmat = [u/(np.max(orimat)/2.0) for u in orimat]
np.random.seed(0)
for i in range(0,40):
    for j in range(0,2):
        names['randmat'+str(i)+'_'+str(j)] =  (2*np.random.random(size=(45,45))-1)*0.01
        names['randmat'+str(i)+'_'+str(j)] = names['randmat'+str(i)+'_'+str(j)] + normmat

intramat0 =  np.hstack((randmat00,randmat01))   
intramat = intramat0   
for i in range(1,40):
    names['intramat'+str(i)] = names['randmat'+str(i)+'_'+'0']
    for j in range(1,2):
        names['intramat'+str(i)] = np.hstack((names['intramat'+str(i)],names['randmat'+str(i)+'_'+str(j)]))
    intramat = np.vstack((intramat,names['intramat'+str(i)]))
#fig0 = plt.figure()
#plt.imshow(intramat)
#plt.savefig('intramat_10*2.png')
#intramat = np.ones((450,90))


# Repeatedly Test
f_tot = []

for i in range(0,n):
    if i%100 == 0:
        print i
    # Vibration array
    #np.random.seed(i)
    #gau_array = np.random.normal(0, 0.5, 75000)
    #del_yi_array = [u for u in gau_array if u>=-d and u<=d]
    #print len(del_yi_array)
    #del_yi = np.random.choice(del_yi_array, int(N))
        
    x = np.arange(-d, d, 1.0/45)
    xU, xL = x+1.0/90, x - 1.0/90 
    prob = ss.norm.cdf(xU, scale=sigma) - ss.norm.cdf(xL, scale=sigma)
    prob = prob / prob.sum()    #normalize the probabilities so their sum is 1
    np.random.seed(i) 
    del_yi = np.random.choice(x, size = N, p = prob)
    
    for j in range(0,N):
        names['y'+str(j+1)] = np.linspace(-a+del_yi[j], a+del_yi[j], num = int(90*a))
#        names['x'+str(j+1)] = normal_distribution(names['y'+str(j)], 0, 1)
#        names['x'+str(j+1)] = square(names['y'+str(j+1)])
        names['x'+str(j+1)] = trapezoid_distribution(names['y'+str(j)], 0, 1)
    
    # Dot Intra Pixel Matrix
    dotmat0 = intramat[int(900-45*a):int(900+45*a)]
    for j in range(0,N):
        names['dotmat'+str(j+1)] = intramat[int(900-45*a+45*del_yi[j]):int(900+45*a+45*del_yi[j])]
    
    
    # Flux / Photons Counts
    f0 = np.dot(x0,dotmat0)
    f = []
    for j in range(0,N):
        f.append(np.sum(np.dot(names['x'+str(j+1)],names['dotmat'+str(j+1)])))
    f_tot.append(np.sum(f)*delta_t)
        
    
# Test Error
f_tot = [u/np.mean(f_tot) for u in f_tot]    
m = np.mean(f_tot)
s = np.std(f_tot)
pr = s/m
print m,s,pr
