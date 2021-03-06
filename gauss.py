#!/usr/bin/python
# coding: UTF-8

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy as sp
import os.path

def gaussian(x, mu , sigma):
    return 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-np.power(x-mu, 2.)/(2*np.power(sigma,2.)))
"""
def gaussian2(x, mu , sigma):
    return np.exp(-np.power(x-mu, 2.)/(2*np.power(sigma,2.)))
"""
def kernel(x, mus, sigma):
    y = 0
    for mu in mus:
        print "mu",mu
        y+=gaussian(x,mu,sigma)
    return y

print np.pi
readfilelist = open("ListOfData.txt","r")
Mfreq = open("./analyzed_data/maximum_frequency.dat","w")
lines = readfilelist.readlines()
for line in lines:
    filename, dummy = line.split(' ')
    filename = os.path.basename(filename)
    filename2 = filename.rstrip('.atf')
    f = open("./analyzed_data/"+filename2+"spiketiming.dat","r")
    print filename2
    spiketiming = f.readlines()
    f.close()
    ST = sp.zeros(len(spiketiming))
    for i in range(len(ST)):
        ST[i] = float(spiketiming[i])

    print spiketiming
    print ST
    #kde = stats.gaussian_kde(ST)
    fig = plt.figure()
    X = np.linspace(0,17,num=1700)
    Y = kernel(X,ST,0.200)
    maximum = 0
    for y1 in Y:
        if(y1>maximum):
            maximum = y1
    plt.plot(ST,np.zeros(ST.shape),'b+', ms=20)
    plt.plot(X, kernel(X,ST, 0.200),'r-', label="sigma=0.200")
    #plt.plot(x, kde(x), 'r-', label="Scotts's Rule")
    plt.text(0,1,"maximum = %f"%maximum)
    plt.legend(loc=1)
    plt.savefig("./graph/gk/"+filename2+"_gaussian-kernel.png")
    Mfreq.write(filename2+"\t%f\n"%maximum)
