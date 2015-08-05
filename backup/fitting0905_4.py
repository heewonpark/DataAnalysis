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
    return sum((gamma(x,*g_para)-y)**2 for x, y in zip(hx[0:20],hy))

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
fX = X
fY1 = [expFunc2(x, *param_fin) for x in fX]
pylab.plot(fX, fY1, 'b--', label = 'fitting', linewidth = 1.0)
"""
fY2 = expFunc3(X, *param_leastsq[0])
pylab.plot(fX, fY2, 'r--', label = 'leastsq', linewidth = 1.0)
"""
pylab.xlabel('X axis')
pylab.ylabel('Y axis')
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
            residual_Y[cnt] = (Ye[k+70]-fY1[k])/fY1[k]
            #if (residual_Y[cnt]>-8.0)&(residual_Y[cnt]<-4.0):
            #print repr(cnt) + ' '+repr(Ye[k+70])+ '  '+repr(fY1[k]) + ' '+repr(residual_Y[cnt])
            #cnt2 +=1
            if (residual_Y[cnt]==-1.0):
                print repr(cnt) + ' '+repr(Ye[k+70])+ '  '+repr(fY1[k]) + ' '+repr(residual_Y[cnt])
                cnt2 +=1
            cnt +=1
            


hist = pylab.figure()
text = 'a '+ repr(param_fin[0])+ ' b '+ repr(param_fin[1])+ ' c '+ repr(param_fin[2])+'\n'+'d ' + repr(param_fin[3])+ ' e '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*exp(-(x-b)/e)'

pylab.text(-1.0,13,text,fontsize=10)
histreturn = pylab.hist(residual_Y[0:cnt+1], bins=40,range=(-2,8))
hy = histreturn[0]
hx = histreturn[1]
print hx
print hy
pylab.savefig("fitting_residualhistogram.png")
"""
g_para_ini = [5.05125, 1.059, 77.887, -5.0002]
g_para_fin = scipy.optimize.fmin(cost_gamma, g_para_ini)
stats = pylab.figure()
pylab.bar(hx[0:20],hy,width = 1.0)
hx_fit = sp.arange(-5.0,15.0,0.1)
hy_fit = [gamma(x,*g_para_fin) for x in hx_fit]
print g_para_fin
pylab.plot(hx_fit, hy_fit,'r--',linewidth=1.0)
pylab.savefig("fitting_gamma.png")
"""
pylab.show()
