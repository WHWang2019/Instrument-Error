#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:02:35 2019

@author: dingding
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy
import scipy.optimize as opt
import scipy.io as scio
import random
import os.path
from scipy import signal
from scipy import integrate
import scipy.stats as ss
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


# Choose Parameter
delta_t = 0.1     # single vibration duration
t = 3    # exposure time
T = 60    # total observation time
N = int(t/delta_t)
n =  int(T/t)   # test counts
a1 = 4.0    # pixel number of x axis, intra pixel 2a1*2a2
a2 = 4.0    # pixel number of y axis, intra pixel 2a1*2a2
d = 10    # vibration change range
si = 1.0    # vibration distribution sigma
rnum = 4    #reference star number


# 1D Gaussian Distribution Function
def normal_distribution_1d(x, mean, sigma):
    f = np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)
    return f

# 2D Gaussian Distribution Function
def normal_distribution_2d(x, y, mean1, mean2, sigma1, sigma2):
    f = np.exp(-1.0/2*(((x-mean1)/sigma1)**2+((y-mean2)/sigma2)**2))/(2*np.pi*sigma1*sigma2)
    return f

# 1D Trapezoid Distribution Function
def trapezoid_distribution_1d(x, mean, sigma):
    a = mean-3.0*sigma
    b = mean+3.0*sigma
    c = mean-1.0*sigma
    d = mean+1.0*sigma
    h = 2.0/((b-a)+(d-c))
    f = np.where((x>=c)&(x<=d), h, np.where((x>=a)&(x<c), h/(c-a)*(x-a), np.where((x>d)&(x<=b), h/(d-b)*(x-b), 0)))
    return f

# 2D Trapezoid Distribution Function
def trapezoid_distribution_2d(x, y, mean1, mean2, sigma1, sigma2):    
    f = trapezoid_distribution_1d(x, mean1, sigma1)
    g = trapezoid_distribution_1d(y, mean2, sigma2)
    return f*g


# Point Spread Function (PSF) as Gaussian Distribution
x0 = np.linspace(-a1, a1, num = int(45*2*a1))
y0 = np.linspace(-a2, a2, num = int(45*2*a2))
x0, y0 = np.meshgrid(x0, y0)
#z0 = normal_distribution_2d(x0, y0, 0, 0, 1, 1)
z0 = trapezoid_distribution_2d(x0, y0, 0, 0, 1, 1)
print z0.shape
print np.sum(z0/(45*45))

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_surface(x0, y0, z0, rstride=1, cstride=1, cmap='rainbow')
#plt.show()
#plt.savefig('2d_x0y0z0.png')


## Normalize Matrix
#def normalize(data):
#    m = np.mean(data)
#    mx = np.max(data)
#    mn = np.min(data)
#    return [(float(u) - m) / (mx - mn) for u in data]



## Choose Reference Star
#np.random.seed(rnum)
#rstarx =  (2*np.random.random(size=rnum)-1)*0.25*45.0
#rstarx = [int(rstarx[i])/45.0 for i in range(len(rstarx))]
#rstary =  (2*np.random.random(size=rnum)-1)*0.25*45.0
#rstary = [int(rstary[i])/45.0 for i in range(len(rstary))]
#print rstarx,rstary        


# Intra Pixel Matrix
names = locals()
fn = scio.loadmat('450nm_shift.mat')
orimat = fn['pixel22']
#fig1 = plt.figure()
#plt.imshow(orimat)
#plt.savefig('orimat.png')
normmat = [np.sqrt(u) for u in orimat]
normmat = [u/(np.max(normmat)*2.0/3) for u in normmat]
np.random.seed(0)
for i in range(0,40):
    for j in range(0,40):
        names['randmat'+str(i)+str(j)] =  (2*np.random.random(size=(45,45))-1)*0.01
        names['randmat'+str(i)+str(j)] = names['randmat'+str(i)+str(j)] + normmat

intramat0 = randmat00
for j in range(1,40):
    intramat0 =  np.hstack((intramat0,names['randmat'+'0'+str(j)]))   
intramat = intramat0   
for i in range(1,40):
    names['intramat'+str(i)] = names['randmat'+str(i)+'0']
    for j in range(1,40):
        names['intramat'+str(i)] = np.hstack((names['intramat'+str(i)],names['randmat'+str(i)+str(j)]))
    intramat = np.vstack((intramat,names['intramat'+str(i)]))
#fig0 = plt.figure()
#plt.imshow(intramat)
#plt.savefig('intramat_20*20.png')
#intramat = np.ones((1800,1800))


# Repeatedly Test
fl_tot = []
for i in range(0,rnum):
    names['fl_rstar_tot_'+str(i)] = []   

for k in range(0,n):
    if k%100 == 0:
        print k
    
    # Vibration array
    #np.random.seed(i)
    #gau_array = np.random.normal(0, 1, 300)
    #del_xi_array = [u for u in gau_array if u>=-d and u<=d]
    #del_xi_array = [u for u in gau_array if u>=-d and u<=d]
    #del_yi = np.random.choice(del_xi_array, int(N))
    #del_yi = np.random.choice(del_yi_array, int(N))

    u = np.arange(-d, d, 1.0/45)
    uU, uL = u+1.0/90, u-1.0/90 
    prob = ss.norm.cdf(uU, scale=si) - ss.norm.cdf(uL, scale=si)
    prob = prob / prob.sum()  
    np.random.seed(k)       
    del_xi = np.random.choice(u, size = N, p = prob)
    del_yi = np.random.choice(u, size = N, p = prob)
    
    for h in range(0,N):
        names['x'+str(h+1)] = np.linspace(-a1+del_xi[h], a1+del_xi[h], num = int(45*2*a1))
        names['y'+str(h+1)] = np.linspace(-a2+del_yi[h], a2+del_yi[h], num = int(45*2*a2))
        names['x'+str(h+1)], names['y'+str(h+1)] = np.meshgrid(names['x'+str(h+1)], names['y'+str(h+1)])
#        names['z'+str(h+1)] = normal_distribution_2d(names['x'+str(h+1)], names['y'+str(h+1)], 0, 0, 1, 1)
        names['z'+str(h+1)] = trapezoid_distribution_2d(names['x'+str(h+1)], names['y'+str(h+1)], 0, 0, 1, 1)
#        names['x'+str(j+1)] = square(names['y'+str(j+1)])
#        print names['z'+str(h+1)].shape
#        print np.sum(names['z'+str(h+1)])
        
#        for i in range(0,rnum):
#            names['x_rstar_'+str(i)+'_'+str(h+1)] = np.linspace(rstarx[i]-a1+del_xi[h], rstarx[i]+a1+del_xi[h], num = int(45*2*a1))
#            names['y_rstar_'+str(i)+'_'+str(h+1)] = np.linspace(rstary[i]-a2+del_yi[h], rstary[i]+a2+del_yi[h], num = int(45*2*a2))
#            names['x_rstar_'+str(i)+'_'+str(h+1)], names['y_rstar_'+str(i)+'_'+str(h+1)] = np.meshgrid(names['x_rstar_'+str(i)+'_'+str(h+1)], names['y_rstar_'+str(i)+'_'+str(h+1)])
#            names['z_rstar_'+str(i)+'_'+str(h+1)] = normal_distribution_2d(names['x_rstar_'+str(i)+'_'+str(h+1)], names['y_rstar_'+str(i)+'_'+str(h+1)], 0, 0, 1, 1)
    
    
    # Dot Intra Pixel Matrix
    dotmat0 = intramat[int(900-45*a2):int(900+45*a2),int(900-45*a1):int(900+45*a1)]
    for h in range(0,N):
        names['dotmat'+str(h+1)] = intramat[int(900-45*a2+45*del_yi[h]):int(900+45*a2+45*del_yi[h]),int(900-45*a1+45*del_xi[h]):int(900+45*a1+45*del_xi[h])]
#        for i in range(0,rnum):
#            names['dotmat_'+str(i)+'_'+str(h+1)] = intramat[int(900+45*rstary[i]-45*a2+45*del_yi[h]):int(900+45*rstary[i]+45*a2+45*del_yi[h]),int(900+45*rstarx[i]-45*a1+45*del_xi[h]):int(900+45*rstarx[i]+45*a1+45*del_xi[h])]

    
    # Flux / Photons Counts
    fl0 = np.dot(z0,dotmat0)
    fl = []
#    for i in range(0,rnum):
#        names['fl_rstar_'+str(i)+'_rp'+str(k)] = []
    for h in range(0,N):
        fl.append(np.sum(np.dot(names['z'+str(h+1)],names['dotmat'+str(h+1)])/(45*45)))
#        for i in range(0,rnum):
#            names['fl_rstar_'+str(i)+'_rp'+str(k)].append(np.sum(np.dot(names['z_rstar_'+str(i)+'_'+str(h+1)],names['dotmat_'+str(i)+'_'+str(h+1)])/(45*45)))
#    for i in range(0,rnum):
#        names['fl_rstar_tot_'+str(i)].append(np.sum(names['fl_rstar_'+str(i)+'_rp'+str(k)])*delta_t)
        
    fl_tot.append(np.sum(fl)*delta_t)
        
#fl_tstar_tot = [u/np.mean(fl_tot) for u in fl_tot]
#for i in range(0,rnum):
#    names['fl_rstar_tot_'+str(i)] = [u/np.mean(names['fl_rstar_tot_'+str(i)]) for u in names['fl_rstar_tot_'+str(i)]]

#fl_rstar_mean = []
#for k in range(0,n):
#    fl_rstar_sum = 0
#    for i in range(0,rnum):
#        fl_rstar_sum = fl_rstar_sum + names['fl_rstar_tot_'+str(i)][k]
#    fl_rstar_mean.append(fl_rstar_sum/rnum)
#    
#fl_del = [fl_rstar_mean[i] - fl_tstar_tot[i] for i in range(len(fl_tstar_tot))]
    
    
    names['x_'+str(k+1)], names['y_'+str(k+1)], names['z_'+str(k+1)] = x1, y1, z1
    for h in range(1,N):
        names['x_'+str(k+1)] = np.hstack((names['x_'+str(k+1)],names['x'+str(h+1)]))
        names['y_'+str(k+1)] = np.hstack((names['y_'+str(k+1)],names['y'+str(h+1)]))
        names['z_'+str(k+1)] = np.hstack((names['z_'+str(k+1)],names['z'+str(h+1)]))
        

# save data
#fi = open('2d_data.txt','w')
#np.savetxt(fi, fl_tot, fmt='%.6f',newline='\n')
    

# Test Error
#m = np.mean(fl_tstar_tot)
#s = np.std(fl_tstar_tot)
#print m,s
#m_com = np.mean(fl_rstar_mean)
#s_com = np.std(fl_rstar_mean)
#print m_com,s_com
#m_del = np.mean(fl_del)
#s_del = np.std(fl_del)
#print m_del,s_del
#
#plt.figure()
#plt.plot(fl_rstar_mean,'b')
#plt.plot(fl_tstar_tot,'g')
#plt.savefig('flux_2star.png')
#plt.figure()
#plt.plot(fl_del)
#plt.savefig('flux_del.png')


fig = plt.figure(figsize=(4,4))
plt.xlim(-4,4)
plt.ylim(-4,4)
plt.contourf(x_20, y_20, z_20, cmap='Greys')
currentAxis=plt.gca()
rect1=patches.Rectangle((-3, -3),6,6,linewidth=1,edgecolor='b',facecolor='none')
rect2=patches.Rectangle((-2, -2),4,4,linewidth=1,edgecolor='g',facecolor='none')
currentAxis.add_patch(rect1)
currentAxis.add_patch(rect2)
#plt.savefig('contourf_4_4.png')