#! /usr/bin/python
# coding: UTF-8

### HEEWON PARK
### 2016.01.23

### OPTIMIZATION PROGRAM FOR MICHAELIS-MENTEN EQUATION ###


import matplotlib.pyplot as plt
import scipy.optimize as optimize
import scipy as sp
import numpy as np
import os.path

def Michaelis_Menten_1(c,K,Tau_max):
    n = 1.0
    r = Tau_max/(1+np.power(K/c,n))
    return r

#def Michaelis_Menten_2(c,parameter):
def Michaelis_Menten_2(c,K,Tau_max,n):
    #K = parameter[0]
    #Tau_max = parameter[1]
    #n = parameter[2]
    r = Tau_max/(1+np.power(K/c,n))
    return r

def fit_func(parameter,x,y):
    residual = y-Michaelis_Menten_2(x,parameter)
    return residual

def cost_mm2(param):
    return sum((Michaelis_Menten_2(x, *param)-y)**2 for x,y in zip(X,Y))


save_dir = "./bombykol_200ms/mm_optimize/"
save_param_dir = "./bombykol_200ms/mm_optimize/parameters"

def optimization(parameters):
    data = np.loadtxt("./bombykol_200ms/analyzed_data/psth/psth-dose-response.dat",delimiter=",")
    global X,Y
    X = data[0]
    Y = data[1]
    print X, Y
    
    print "\n"
    print "::::: Curve_fit :::::"
    init_curve_fit = parameters[0]
    print "Initial Value : ", init_curve_fit
    popt,pcov = optimize.curve_fit(Michaelis_Menten_2,X,Y,p0=init_curve_fit)
    print"[CURVE_FIT REUSLT] k = %f, Tau_max = %f, n = %f"%(popt[0],popt[1],popt[2])
    print pcov
    np.savetxt(save_param_dir+"/mm-curve_fit.txt",popt,fmt='%.6f',delimiter=',')        

    print "\n"
    print "::::: fmin :::::"
    init_fmin = parameters[1]
    print "Initial Value : ", init_fmin
    xopt = optimize.fmin(cost_mm2,init_fmin)
    print"[FMIN REUSLT] k = %f, Tau_max = %f, n = %f"%(xopt[0],xopt[1],xopt[2])
    np.savetxt(save_param_dir+"/mm-fmin.txt",popt,fmt='%.6f',delimiter=',')        

    fig = plt.figure(figsize=(10,8),dpi=400)
    plt.rcParams['font.size']=15
    plt.plot(X,Y,'go',label="PSTH Maximum Frequency")
    
    n = 1000
    Xmm = np.logspace(1,5,num=100)
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*popt),'r-',label=r"$M(c)$,Curve_fit")
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*xopt),'b-',label=r"$M(c)$,fmin")
    plt.xscale('log')
    plt.xlabel('Dose[ng]',fontsize=30)
    plt.ylabel('Maximum frequency[Hz]',fontsize=30)
    plt.legend(loc='upper left')
    #plt.savefig('Michaelis-menten_fitted_curve2.png')
    plt.savefig(save_dir+'Michaelis-menten_fitted.png',transparent=True)
    plt.show()

init_curve_fit = [1,10,0.5]
init_fmin = [1,10,0.5]
parameters = [init_curve_fit, init_fmin]

optimization(parameters)
