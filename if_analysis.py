from matplotlib import pylab
import scipy as sp
import numpy as np

filenames = open("listoffilenames.txt",'r')
fname = filenames.readlines()
filenames.close()

for i in range(6):
    fname[i] = fname[i].rstrip('.atf\n')



def if_histogram(filename):
    fp = open(filename+'ifreq.dat','r')
    data = fp.readlines()
    data[0] = int(data[0].rstrip('\n'))
    print data[0]
    time = [0 for i in range(data[0])]
    ifreq = [0 for i in range(data[0])]
    
    for i in range(1,data[0]):
        #print i
        time[i], ifreq[i] = data[i].split('\t')
        time[i] = float(time[i])
        ifreq[i] = float(ifreq[i].rstrip('\n'))
    fig = pylab.figure()
    pylab.hist(ifreq, bins= 250)
    pylab.xlim(0,600)
    pylab.savefig(filename+'ifreq_hist.png')

def if_histogram_whole():
    fp = open('ifreq_whole6.csv','r')
    data = fp.readlines()
    print len(data)
    time = [0 for i in range(len(data))]
    ifreq = [0 for i in range(len(data))]
    
    for i in range(len(data)):
        print i
        time[i], ifreq[i] = data[i].split(',')
        time[i] = float(time[i])
        ifreq[i] = float(ifreq[i].rstrip('\n'))
    fig = pylab.figure()
    pylab.plot(time, ifreq,'b')
    pylab.ylim(0,300)
    pylab.savefig('aver_if6.png')
"""
for i in range(6):
    if_histogram(fname[i])
"""

if_histogram_whole()
