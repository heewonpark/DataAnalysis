#! /usr/bin/python
# coding: UTF-8

"""BEGIN COMMENT
このプログラムで実装する項目:
PSTHの最大周波数を基準にしてMichaelis-Menten式でfittingをする関数->
10000ngで刺激した時のPSTHをFittingする関数．
周波数上昇部，周波数減衰部をそれぞれ違う関数でFittingする

以上．
END COMMENT"""

import numpy as np
import matplotlib.pyplot as plt
import os.path
import scipy.optimize as optimize
import scipy as sp

from SpikeAnalyzer import fitting, analyzer
#*****************************************************
# FUNCTIONS FOR FITTING

def exp_func1(x, alpha, beta, gamma):
    return alpha * sp.exp(-(x-beta)/gamma)
def exp_func2(x, alpha1, beta, gamma1, alpha2, gamma2):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  alpha2 * sp.exp(-(x-beta)/gamma2)
def exp_func3(x, alpha1, beta, gamma1, a, b):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  a*x+b

def cost_fn1(param):
    return sum((exp_func1(x, *param)-y)**2 for x,y in zip(X,Y))
def cost_fn2(param):
    return sum((exp_func2(x, *param)-y)**2 for x,y in zip(X,Y))
def cost_fn3(param):
    return sum((exp_func3(x, *param)-y)**2 for x,y in zip(X,Y))
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

def Michaelis_Menten_1(c,K,Tau_max):
    n = 1.0
    r = Tau_max/(1+np.power(K/c,n))
    return r

def Michaelis_Menten_2(c,parameter):
    K = parameter[0]
    Tau_max = parameter[1]
    n = parameter[2]
    r = Tau_max/(1+np.power(K/c,n))
    return r

def fit_func(parameter,x,y):
    residual = y-Michaelis_Menten_2(x,parameter)
    return residual

#def exp_single_pos(x, alpha, beta, gamma,c):
#    return alpha * sp.exp((x-beta)/gamma) + c

def exp_single(x, alpha, beta, gamma,c):
    c = 1.5
    #beta = 250
    return alpha * sp.exp(-(x-beta)/gamma)+c

def exp_single2(x, gamma):
    c = 1.5
    beta = BETA
    alpha = ALPHA
    return alpha * sp.exp(-(x-beta)/gamma)+c

def exp_double2(x, gamma1, gamma2, p):
    c = 1.5
    #beta = 250
    return ALPHA * p * sp.exp(-(x-beta)/gamma1)  +  ALPHA * (1-p) * sp.exp(-(x-beta)/gamma2) + c

#def exp_single_rise(x, alpha, beta, gamma):
#    return alpha*(1-sp.exp(-(x-beta)/gamma))
#    return alpha*(1-sp.exp(-(x-beta)/gamma))+3.5

def cost_es(param):
    return sum((exp_single(x, *param)-y)**2 for x,y in zip(X,Y))

def cost_es2(param):
    return sum((exp_single2(x, *param)-y)**2 for x,y in zip(X,Y))

def cost_ed2(param):
    return sum((exp_double2(x, *param)-y)**2 for x,y in zip(X,Y))

#def cost_esr(param):
#    return sum((exp_single_rise(x, *param)-y)**2 for x,y in zip(X,Y))

#*****************************************************
# MAIN FUNCTIONS

save_dir = "./bombykol_200ms/fitting/"

def fitting_PSTH(dose, parameters, bp_func,bp_func_cost, ap_func, ap_func_cost):
    global X,Y
    fit = fitting()
    fit.readPSTH("./bombykol_200ms/analyzed_data/psth/Dose%d_psth.dat"%dose)
    fit.readPSTH_ALL("./bombykol_200ms/analyzed_data/psth/Dose%d_psth_all.dat"%dose)
    #print fit.psth_f[:50]
    MAX_INDEX = fit.psth_f[:50].argmax()
    print "MAX INDEX : ",MAX_INDEX," MAX VALUE : ", fit.psth_f[MAX_INDEX]

    #### Optimization for rising phase
    X_bp = fit.psth_t[:MAX_INDEX+1] # Before Peak
    Y_bp = fit.psth_f[:MAX_INDEX+1] # Before Peak
    #X_bp = X_bp-25
    x_bp_all = []
    y_bp_all = []
    for i in range(len(fit.psth_t_all)):
        x_bp_all.append(fit.psth_t_all[i][:MAX_INDEX+1])
        y_bp_all.append(fit.psth_f_all[i][:MAX_INDEX+1])
    x_bp_all=np.array(x_bp_all)
    y_bp_all=np.array(y_bp_all)
    x_bp_all=np.array(x_bp_all)
    y_bp_all=np.array(y_bp_all)
    (a,b)=x_bp_all.shape
    x_bp_all=x_bp_all.reshape(a*b)
    #x_bp_all=x_bp_all - 25
    (a,b)=y_bp_all.shape
    y_bp_all=y_bp_all.reshape(a*b)

    #print X_bp,x_bp_all
    #print y_bp_all
    print "Parameters : ", parameters
    print "**********************************************************"
    print "*** Before Peak, Curve_fit ***\n"
    
    global BETA
    BETA = 250

    init_bp_curve_fit = parameters[0] ## initialize parameters for curve_fit of before peak
    popt,pcov = optimize.curve_fit(bp_func, x_bp_all, y_bp_all, p0=init_bp_curve_fit)
    print "Initial Value : ",init_bp_curve_fit
    print "Optimal Value : ", popt
    print pcov
    np.savetxt(save_dir+"parameters/%dng_exp_single_Rise_Curve_fit.txt"%dose,popt,fmt='%.6f',delimiter=',')        

    #X = X_bp
    #Y = Y_bp
    X = x_bp_all
    Y = y_bp_all

    print "\n\n*** Before Peak, fmin ***"
    init_bp_fmin = parameters[1]
    param_fin = optimize.fmin(bp_func_cost, init_bp_fmin, maxiter=2000,maxfun=2000)
    np.savetxt(save_dir+"/parameters/%dng_exp_single_Rise_fmin.txt"%dose,param_fin,fmt='%.6f',delimiter=',')        
    print "Initial Value : ",init_bp_fmin
    print "Optimal Value : ", param_fin
    print "Initial Cost : ",bp_func_cost(init_bp_fmin)," Optimized Cost : ",bp_func_cost(param_fin)

    #### Draw Graph of fitting result
    x_bp = np.arange(0,250.0,0.01)
    y_bp_cf = bp_func(x_bp, *popt)
    y_bp_fmin = bp_func(x_bp, *param_fin)
    fig = plt.figure()
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'o')
    plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    plt.text(250,bp_func(250,*popt)+5,bp_func(250,*popt))
    plt.text(250,bp_func(250,*param_fin)+5,bp_func(250,*param_fin))
    plt.xlim(0,300)
    plt.plot(x_bp, y_bp_cf,label="cf, exp_single")
    plt.legend()
    plt.savefig(save_dir+"fitting_before_Stim_200ms_%dng.png"%dose)

    
    print "\n\n**********************************************************"
    #### Optimization for falling phase
    
    X_ap = fit.psth_t[MAX_INDEX:50] # After Peak
    Y_ap = fit.psth_f[MAX_INDEX:50] # After Peak
    #X_ap = X_ap+25
    x_ap_all = []
    y_ap_all = []
    for i in range(len(fit.psth_t_all)):
        x_ap_all.append(fit.psth_t_all[i][MAX_INDEX:50])
        y_ap_all.append(fit.psth_f_all[i][MAX_INDEX:50])
        
    x_ap_all=np.array(x_ap_all)
    y_ap_all=np.array(y_ap_all)
    x_ap_all=np.array(x_ap_all)
    y_ap_all=np.array(y_ap_all)
    (a,b)=x_ap_all.shape
    x_ap_all=x_ap_all.reshape(a*b)
    #x_ap_all=x_ap_all + 25
    (a,b)=y_ap_all.shape
    y_ap_all=y_ap_all.reshape(a*b)
    #print X_ap, x_ap_all
    print "*** After Peak, Curve_fit ***\n"
    init_ap_curve_fit = parameters[2]
    popt,pcov = optimize.curve_fit(ap_func, x_ap_all, y_ap_all,p0=init_ap_curve_fit)
    print "Initial Value : ", init_ap_curve_fit
    print "Optimal Value : ", popt
    print pcov
    np.savetxt(save_dir+"parameters/%dng_exp_single_Fall_Curve_fit.txt"%dose,popt,fmt='%.6f',delimiter=',')        
    #X = X_ap
    #Y = Y_ap
    X = x_ap_all
    Y = y_ap_all
    print "\n\n*** After Peak, fmin ***\n"
    #parameter = [70,190,170,0]
    init_ap_fmin = parameters[3]
    param_fin = optimize.fmin(ap_func_cost, init_ap_fmin)
    
    print "Initial Value : ", init_ap_fmin
    print "Optimal Value : ", param_fin
    print "Initial Cost : ", ap_func_cost(init_ap_fmin)," Optimized Cost : ",ap_func_cost(param_fin)
    np.savetxt(save_dir+"parameters/%dng_exp_single_Fall_fmin.txt"%dose,param_fin,fmt='%.6f',delimiter=',')        
    print "\n\n*** Fitting Finish ***"
    print "**********************************************************"


    #### Draw Graphs of fitting results
    x_ap = np.arange(fit.psth_t[MAX_INDEX],5000,1)

    y_ap_fmin = ap_func(x_ap, *param_fin)
    y_ap_cf = ap_func(x_ap, *popt)
    fig = plt.figure()
    plt.plot(x_ap_all,y_ap_all,'o')
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,label="cf, exp_single")
    plt.legend()
    plt.xlim(0,5000)
    plt.savefig(save_dir+"fitting_after_Stim_200ms_%dng.png"%dose)

    fig = plt.figure()
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'o')
    plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    #plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.plot(x_bp, y_bp_cf,label="cf, exp_single")

    plt.plot(x_ap_all,y_ap_all,'o')
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,label="cf, exp_single")
    plt.legend()
    plt.xlim(0,1000)

    plt.savefig(save_dir+"fitting_Stim_200ms_%dng_10000ms.png"%dose)

    #### Graph for Thesis small
    fig = plt.figure(figsize=(4,3),dpi=250)
    fig.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 10 #フォントサイズを設定
    #plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    plt.bar(fit.psth_t-50,fit.psth_f,100,label="PSTH")
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_bp_all,y_bp_all,'o')
    #plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    #plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.plot([0,200],[58,58],'k',linewidth=4.0)
    plt.plot(x_bp, y_bp_cf,'r',label="Fitted Curve",linewidth=2.0)
    #plt.plot(x_ap_all,y_ap_all,'o')
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,'r',linewidth=2.0)
    plt.legend(frameon=False)
    plt.xlim(-100,5000)
    plt.xlabel("Time[ms]")
    plt.ylabel("Frequency[Hz]")
    
    plt.savefig(save_dir+"fitting_Stim_200ms_%dng_5s.png"%dose)

    #### Graph for Presentation
    fig = plt.figure(figsize=(6,4.5),dpi=250)
    fig.subplots_adjust(bottom=0.2,left =0.15)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    plt.rcParams['axes.linewidth'] = 2.0 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    plt.bar(fit.psth_t-50,fit.psth_f,100,label="PSTH")
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_bp_all,y_bp_all,'o')
    #plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    #plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.plot([0,200],[58,58],'k',linewidth=2.0)
    plt.plot(x_bp, y_bp_cf,'r',label="Fitted Curve",linewidth=2.0)
    #plt.plot(x_ap_all,y_ap_all,'o')
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,'r',linewidth=2.0)
    plt.legend(frameon=False,fontsize=15)
    plt.xlim(-100,5000)
    plt.xlabel("Time[ms]")
    plt.ylabel("Frequency[Hz]")
    
    plt.savefig(save_dir+"fitting_Stim_200ms_%dng_5s_fp_bar.png"%dose)


    #### Graph for Thesis big
    fig = plt.figure(figsize=(10,8),dpi=400)
    fig.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    #plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    #plt.bar(fit.psth_t-50,fit.psth_f,100,label="PSTH")
    #plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'k.')
    #plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    #plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.plot([0,200],[-3,-3],'k',linewidth=4.0)
    plt.plot(x_bp, y_bp_cf,'r',label="Fitted Curve",linewidth=2.0)
    plt.plot(x_ap_all,y_ap_all,'k.')
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,'r',linewidth=2.0)
    plt.legend(frameon=False)
    plt.xlim(-100,5000)
    plt.ylim(-5,95)
    plt.yticks(np.arange(0,100,10))
    plt.xticks(np.arange(0,6000,1000))
    plt.xlabel("Time[ms]")
    plt.ylabel("Frequency[Hz]")
    
    plt.savefig(save_dir+"fitting_Stim_200ms_%dng_5s_ft.png"%dose)

    #### Graph for Presentation
    fig = plt.figure(figsize=(6,4.5),dpi=250)
    fig.subplots_adjust(bottom=0.2,left=0.15)
    ax = fig.add_subplot(111)
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    plt.rcParams['axes.linewidth'] = 2.0 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    #plt.bar(fit.psth_t-50,fit.psth_f,100,label="PSTH")
    #plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'k.')
    #plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    #plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.plot([0,200],[-3,-3],'k',linewidth=4.0)
    plt.plot(x_bp, y_bp_cf,'r',label="Fitted Curve",linewidth=2.0)
    plt.plot(x_ap_all,y_ap_all,'k.')
    #plt.plot(fit.psth_t,fit.psth_f)
    #plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,'r',linewidth=2.0)
    plt.legend(frameon=False)
    plt.xlim(-100,5000)
    plt.ylim(-5,95)
    plt.yticks(np.arange(0,100,10))
    plt.xticks(np.arange(0,6000,1000))
    plt.xlabel("Time[ms]")
    plt.ylabel("Frequency[Hz]")
    
    plt.savefig(save_dir+"fitting_Stim_200ms_%dng_5s_fp.png"%dose)

    #plt.show()

#fitting_PSTH_10000ng()
#fitting_PSTH_5000ng()
#fitting_PSTH_2000ng()
ALPHA = 42.5
#init_bp_curve_fit = [40,-50,40,-50]
init_bp_curve_fit = [-50]
init_bp_fmin = [-50]
#init_bp_fmin      = [40,250,-150,0]
#init_ap_curve_fit = [250,50,0.7]
#init_ap_fmin      = [250,50,0.7]
init_ap_curve_fit = [175]
init_ap_fmin      = [175]
param_10000ng = [init_bp_curve_fit, init_bp_fmin, init_ap_curve_fit, init_ap_fmin]
####fitting_PSTH(10000, param_10000ng, exp_single2, cost_es2, exp_single, cost_es)
####fitting_PSTH(10000, param_10000ng, exp_single2, cost_es2, exp_double2, cost_ed2)
fitting_PSTH(10000, param_10000ng, exp_single2, cost_es2, exp_single2, cost_es2)

ALPHA = 51.833
#init_bp_curve_fit = [52,-70]
#init_bp_fmin = [40,-80]
init_bp_curve_fit = [-70]
init_bp_fmin = [-70]
#init_bp_fmin      = [40,250,-150,0]
#init_ap_curve_fit = [70,190,170,0]
#init_ap_fmin      = [70,190,170,0]
init_ap_curve_fit = [150]
init_ap_fmin      = [150]
param_5000ng = [init_bp_curve_fit, init_bp_fmin, init_ap_curve_fit, init_ap_fmin]
####fitting_PSTH(5000, param_5000ng, exp_single2, cost_es2, exp_single, cost_es)
fitting_PSTH(5000, param_5000ng, exp_single2, cost_es2, exp_single2, cost_es2)

ALPHA = 15.167
init_bp_curve_fit = [-70]
init_bp_fmin = [-70]
#init_bp_fmin      = [40,250,-150,0]
init_ap_curve_fit = [150]
init_ap_fmin      = [150]
param_2000ng = [init_bp_curve_fit, init_bp_fmin, init_ap_curve_fit, init_ap_fmin]
####fitting_PSTH(2000, param_2000ng, exp_single, cost_es, exp_single, cost_es)
fitting_PSTH(2000, param_2000ng, exp_single2, cost_es2, exp_single2, cost_es2)
