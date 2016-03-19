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

import scipy.stats as stats

def Michaelis_Menten_1(c,K,Tau_max):
    n = 1.0
    r = Tau_max/(1+np.power(K/c,n))
    return r

#def Michaelis_Menten_2(c,parameter):
"""
def Michaelis_Menten_2(c,K,Tau_max,n):
    #K = parameter[0]
    #Tau_max = parameter[1]
    #n = parameter[2]
    Tau_max = 60
    n = 1
    r = Tau_max/(1+np.power(K/c,n))
    return r
"""

def Michaelis_Menten_2(c,K):
    #K = parameter[0]
    #Tau_max = parameter[1]
    #n = parameter[2]
    Tau_max = 60
    n = 1
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
    data_avg = np.loadtxt("./bombykol_200ms/analyzed_data/psth/psth-dose-response.dat",delimiter=",")
    data = np.loadtxt("./bombykol_200ms/analyzed_data/psth-dose-response-all.csv",delimiter=",")
    global X,Y
    X = data[0]
    Y = data[1]
    print X, Y
    
    print "\n"
    print "::::: Curve_fit :::::"
    init_curve_fit = parameters[0]
    print "Initial Value : ", init_curve_fit
    popt,pcov = optimize.curve_fit(Michaelis_Menten_2,X,Y,p0=init_curve_fit)
    print"[CURVE_FIT REUSLT]", popt
    print pcov
    np.savetxt(save_param_dir+"/mm-curve_fit.txt",popt,fmt='%.6f',delimiter=',')        

    print "\n"
    print "::::: fmin :::::"
    init_fmin = parameters[1]
    print "Initial Value : ", init_fmin
    xopt = optimize.fmin(cost_mm2,init_fmin)
    print"[FMIN REUSLT]",xopt#     k = %f, Tau_max = %f, n = %f"%(xopt[0],xopt[1],xopt[2])
    np.savetxt(save_param_dir+"/mm-fmin.txt",popt,fmt='%.6f',delimiter=',')        

    fig = plt.figure(figsize=(10,8),dpi=400)
    plt.rcParams['font.size']=15

    Dose_tau = [2000,5000,10000]
    tau_rise = [-23,-71.6,-50.4]
    tau_fall = [80,125,175]
    #plt.plot(Dose_tau, tau_rise,label="tau_rise")
    #plt.plot(Dose_tau, tau_fall,label="tau_fall")

    plt.plot(X,Y,'g+',label="PSTH Maximum Frequency")
    plt.plot(data_avg[0],data_avg[1],'go',label="PSTH Average Maximum Frequency")
    n = 1000
    Xmm = np.logspace(1,5,num=100)
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*popt),'r-',label=r"$M(c)$,Curve_fit")
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*xopt),'b-',label=r"$M(c)$,fmin")
    plt.xscale('log')
    plt.xlabel('Dose[ng]',fontsize=30)
    plt.ylabel('Maximum frequency[Hz]',fontsize=30)
    #plt.legend(loc='upper left')
    #plt.savefig('Michaelis-menten_fitted_curve2.png')
    plt.savefig(save_dir+'Michaelis-menten_fitted.png',transparent=True)
    #plt.show()

    #### Graph for Thesis small
    fig = plt.figure(figsize=(4,3),dpi=250)
    fig.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 10 #フォントサイズを設定
    #plt.plot(X,Y,'g+',label="PSTH Maximum Frequency")
    plt.plot(data_avg[0],data_avg[1],'go',label="PSTH Maximum Frequency")
    n = 1000
    Xmm = np.logspace(0.1,5,num=100)
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*popt),'r-',label=r"Fitted Curve")
    #plt.plot(Xmm,Michaelis_Menten_2(Xmm,*xopt),'b-',label=r"$M(c)$,fmin")
    plt.xscale('log')
    plt.xlabel('Dose[ng]',fontsize=10)
    plt.ylabel('Maximum frequency[Hz]',fontsize=10)
    plt.legend(loc=2,frameon=False,fontsize=10)
    #plt.legend(loc='upper left')
    #plt.savefig('Michaelis-menten_fitted_curve2.png')
    plt.savefig(save_dir+'Michaelis-menten_fitted_forThesis.png',transparent=True)
    #plt.show()

    psth_data = [[10,0,0,0,0,0],
                 [10,0,0,0,0,0],
                 [20,0,0,0,0,0],
                 [0,80,0,10,10,0],
                 [40,90,60,0,70,60],
                 [30,70,10,50,60]]
    print psth_data
    psth_data = np.array(psth_data)
    psth_sem = [None for i in range(6)]
    psth_sem[0] = stats.sem(psth_data[0])
    psth_sem[1] = stats.sem(psth_data[1])
    psth_sem[2] = stats.sem(psth_data[2])
    psth_sem[3] = stats.sem(psth_data[3])
    psth_sem[4] = stats.sem(psth_data[4])
    psth_sem[5] = stats.sem(psth_data[5])
    print psth_sem
    
    #### Graph for Thesis big
    fig = plt.figure(figsize=(10,8),dpi=400)
    fig.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    #plt.plot(X,Y,'g+',label="PSTH Maximum Frequency")
    plt.errorbar(data_avg[0],data_avg[1],yerr=psth_sem,fmt='ko',label="PSTH Maximum Frequency")
    n = 1000
    Xmm = np.logspace(0.1,5,num=100)
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*popt),'r-',label=r"Fitted Curve",linewidth=2.0)
    #plt.plot(Xmm,Michaelis_Menten_2(Xmm,*xopt),'b-',label=r"$M(c)$,fmin")
    plt.xscale('log')
    plt.xlabel('Dose[ng]')
    plt.ylabel('Maximum frequency[Hz]')
    plt.ylim(0,70)
    plt.legend(loc=2,frameon=False)
    #plt.legend(loc='upper left')
    #plt.savefig('Michaelis-menten_fitted_curve2.png')
    plt.savefig(save_dir+'Michaelis-menten_fitted_forThesis_big.png',transparent=True)
    #plt.show()


    #### Graph for Presentation
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    plt.rcParams['axes.linewidth'] = 2.0 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    fig = plt.figure(figsize=(6,4.5),dpi=250)
    fig.subplots_adjust(bottom=0.2,left=0.15)
    ax = fig.add_subplot(111)
    plt.errorbar(data_avg[0],data_avg[1],yerr=psth_sem,fmt='ko',label="PSTH Maximum Frequency")
    n = 1000
    Xmm = np.logspace(0.1,5,num=100)
    plt.plot(Xmm,Michaelis_Menten_2(Xmm,*popt),'r-',label=r"Fitted Curve",linewidth=2.0)
    #plt.plot(Xmm,Michaelis_Menten_2(Xmm,*xopt),'b-',label=r"$M(c)$,fmin")
    plt.xscale('log')
    plt.xlabel('Dose[ng]')
    plt.ylabel('Maximum frequency[Hz]')
    plt.ylim(0,70)
    plt.legend(loc=2,frameon=False,fontsize=15)
    #plt.legend(loc='upper left')
    #plt.savefig('Michaelis-menten_fitted_curve2.png')
    plt.savefig(save_dir+'Michaelis-menten_fitted_forPresentation.png',transparent=True)




init_curve_fit = [1,10,2]
init_fmin = [1,10,2]
parameters = [init_curve_fit, init_fmin]
#optimization(parameters)

init_curve_fit = [3000]
init_fmin = [3000]
parameters = [init_curve_fit, init_fmin]
optimization(parameters)
