#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2015.03.23
###############################################

from matplotlib import pylab
import scipy as sp
import scipy.optimize
import numpy as np
import math 

def expFunc1(x, alpha, beta, gamma):
    return alpha * sp.exp(-(x-beta)/gamma)
def expFunc2(x, alpha1, beta, gamma1, alpha2, gamma2):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  alpha2 * sp.exp(-(x-beta)/gamma2)
def expFunc3(x, alpha1, beta, gamma1, a, b):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  a*x+b
def cost_fn1(param):
    return sum((expFunc1(x, *param)-y)**2 for x,y in zip(X,Y))
def cost_fn2(param):
    return sum((expFunc2(x, *param)-y)**2 for x,y in zip(X,Y))
def cost_fn3(param):
    return sum((expFunc3(x, *param)-y)**2 for x,y in zip(X,Y))
def residue(param, y, x):
    return y-expFunc3(x,*param)
def gamma(x, k, theta, alpha, beta):
    gamma = alpha*pow(x-beta, k-1) * sp.exp(-(x-beta)/theta)/(pow(theta,k) * math.gamma(k))
    return gamma
def cost_gamma(g_para):
    return sum((gamma(x,*g_para)-y)**2 for x, y in zip(hx[1:40],hy[1:]))

def gaussian(x, sigma, mu,alpha):
    return alpha/(sigma*np.sqrt(np.pi*2))*sp.exp(-pow((x-mu),2)/(2*pow(sigma,2)))

def cost_gaussian(gaus_para):
    return sum((gaussian(x, *gaus_para)-y)**2 for x, y in zip(GX, GY))

readfilelist = open("listoffiles.txt",'r')
filenames = readfilelist.readlines()
readfilelist.close()
filename = [' ' for i in range(6)]
dummy    = [' ' for i in range(6)]

for i in range(6):
    filename[i], dummy[i] = filenames[i].split(' ')
    filename[i] = filename[i].rstrip('.atf')
    print filename[i]

readdatfile = open("average.dat",'r')
data = readdatfile.readlines()
readdatfile.close()

X = sp.zeros(81)
Y = sp.zeros(81)
for i in range(0,81):
    X[i], Y[i] = data[i+71].split('\t')
    X[i] = float(X[i])
    Y[i] = float(Y[i])
"""    
print 'X'
print X 
print 'Y'
print Y
"""
param_fin = [0.0001333, 12.91, 0.5157,5.639, 4.726]
#param_fin = [0.25, 11.9, 0.7, -1.06, 19.7]

param_fin = scipy.optimize.fmin(cost_fn2, param_fin,xtol=1e-8)
print param_fin

"""
print 'leastsq---------------------------------------------------------------------------'
param_leastsq = scipy.optimize.leastsq(residue, param_ini, args=(Y,X),full_output=True)
print param_leastsq
"""
graph = pylab.figure()
pylab.plot(X,Y,'g', label='data', linewidth=1.0)
pylab.bar(X,Y,width =  0.1)
fX = X
fY1 = [expFunc2(x, *param_fin) for x in fX]
pylab.plot(fX, fY1, 'r--', label = 'F2', linewidth = 2.0)

pylab.xlabel('Time[s]')
pylab.ylabel('Frequency[Hz]')
pylab.title('Parameter optimization')
pylab.legend()
pylab.grid(True)
pylab.savefig("fitting.png")

residual_Y = [None for i in range(81*6)]
cnt = 0
cnt2 = 0
for i in range(6):
    readbinfile= open(filename[i]+'bin.dat','r')
    bindata = readbinfile.readlines()
    steps, delay = bindata[0].split(' ')
    steps = int(steps)
    delay = int(delay)

    Xe = sp.zeros(steps+10+delay)
    Ye = sp.zeros(steps+10+delay)
    print steps
    print filename[i]
    for j in range(steps+10+delay):
        Xe[j],Ye[j] = bindata[j+1].split('\t')
        Xe[j] = float(Xe[j])
        Ye[j] = float(Ye[j])

    for k in range(81):
        if (k+70)<steps+10+delay:
            if(Ye[k+70] != 0):
                residual_Y[cnt] = (Ye[k+70]-fY1[k])
                #if (residual_Y[cnt]>-8.0)&(residual_Y[cnt]<-4.0):
                #print repr(cnt) + ' '+repr(Ye[k+70])+ '  '+repr(fY1[k]) + ' '+repr(residual_Y[cnt])
                #cnt2 +=1
                # if (residual_Y[cnt]==-1.0):
                #  print repr(cnt) + ' '+repr(Ye[k+70])+ '  '+repr(fY1[k]) + ' '+repr(residual_Y[cnt])
                #    cnt2 +=1
                cnt +=1
            


hist = pylab.figure()
text = 'a '+ repr(param_fin[0])+ ' b '+ repr(param_fin[1])+ ' c '+ repr(param_fin[2])+'\n'+'d ' + repr(param_fin[3])+ ' e '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*exp(-(x-b)/e)'

#pylab.text(-1.0,13,text,fontsize=10)
histreturn = pylab.hist(residual_Y[0:cnt], bins=40,range=(-100,100))
pylab.xlim = (-100,100)
GX = histreturn[1]
GY = histreturn[0]
#hy = histreturn[0]
#hx = histreturn[1]
#print hx
#print hy
pylab.savefig("fitting_residualhistogram.png")
"""
g_para_ini = [13.0, 0.16, 68.7, -1.9]
g_para_fin = scipy.optimize.fmin(cost_gamma, g_para_ini)
stats = pylab.figure()
pylab.bar(hx[1:40],hy[1:],width = 0.25)
hx_fit = sp.arange(-1.0,8.0,0.1)
hy_fit = [gamma(x,*g_para_fin) for x in hx_fit]
print g_para_fin
pylab.plot(hx_fit, hy_fit,'r--',linewidth=1.0)
pylab.savefig("fitting_gamma.png")
"""
gaus_para_ini = [2.0, 0.7, 10]
gaus_para_fin = scipy.optimize.fmin(cost_gaussian, gaus_para_ini)
gausgraph = pylab.figure()
pylab.bar(GX[0:40],GY,width = 0.25)
gx_fit = sp.arange(-100.0,100.0,1.0)
gy_fit = [gaussian(x,*gaus_para_fin) for x in gx_fit]
print gaus_para_fin
pylab.plot(gx_fit, gy_fit,'r--',linewidth=1.0)
pylab.savefig("fitting_gaus.png")


pylab.show()
