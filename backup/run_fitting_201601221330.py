#! /usr/bin/python
# coding: UTF-8

"""BEGIN COMMENT
このプログラムで実装する項目:
PSTHの最大周波数を基準にしてMichaelis-Menten式でfittingをする関数
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
    #c = 1.43
    return alpha * sp.exp(-(x-beta)/gamma)+c

def exp_single_rise(x, alpha, beta, gamma):
    return alpha*(1-sp.exp(-(x-beta)/gamma))
#    return alpha*(1-sp.exp(-(x-beta)/gamma))+3.5
def cost_es(param):
    return sum((exp_single(x, *param)-y)**2 for x,y in zip(X,Y))

def cost_esr(param):
    return sum((exp_single_rise(x, *param)-y)**2 for x,y in zip(X,Y))

#*****************************************************
# MAIN FUNCTIONS

save_dir = "./bombykol_200ms/fitting/"
def fitting_PSTH_10000ng():
    fit = fitting()
    fit.readPSTH("./bombykol_200ms/analyzed_data/psth/Dose10000_psth.dat")
    fit.readPSTH_ALL("./bombykol_200ms/analyzed_data/psth/Dose10000_psth_all.dat")
    #print fit.psth_f[:50]
    MAX_INDEX = fit.psth_f[:50].argmax()
    print MAX_INDEX, fit.psth_f[MAX_INDEX]
    X_bp = fit.psth_t[:MAX_INDEX+1] # Before Peak
    Y_bp = fit.psth_f[:MAX_INDEX+1] # Before Peak
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
    (a,b)=y_bp_all.shape
    y_bp_all=y_bp_all.reshape(a*b)

    #print x_bp_all
    #print y_bp_all
    
    init_param = [40,250,-150,0]
    #popt,pcov = optimize.curve_fit(exp_single_rise, X_bp, Y_bp, p0=init_param)
    #popt,pcov = optimize.curve_fit(exp_single_rise, x_bp_all, y_bp_all, p0=init_param)
    popt,pcov = optimize.curve_fit(exp_single, x_bp_all, y_bp_all, p0=init_param)
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/10mg_exp_single_Rise_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')



    global X,Y
    #X = X_bp
    #Y = Y_bp
    X = x_bp_all
    Y = y_bp_all
    print X,Y
    parameter = [40,250,-150,0]
    #param_fin = optimize.fmin(cost_esr, parameter)
    param_fin = optimize.fmin(cost_es, parameter)
    print parameter, param_fin
    np.savetxt(save_dir+"parameters/10mg_exp_single_Rise_fmin.txt",param_fin,fmt='%.6f',delimiter=',')

    #x_fbp = np.arange(0,fit.psth_t[MAX_INDEX],0.01)
    x_bp = np.arange(0,300.0,0.01)
    #y_bp_cf = exp_single_rise(x_bp, *popt)
    y_bp_cf = exp_single(x_bp, *popt)
    #y_fbp = exp_single_pos(x_fbp, 1.0, 1.0, 1.0)
    #y_bp_fmin = exp_single_rise(x_bp, *param_fin)
    y_bp_fmin = exp_single(x_bp, *param_fin)
    para = [40,250,-150,0]
    y_bp_custom = exp_single(x_bp, *para)
    #y_bp_cf = exp_func1(x_fbp, *popt)
    #print x_fbp, y_fbp
    fig = plt.figure()
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'o')
    plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.xlim(0,300)
    plt.plot(x_bp, y_bp_cf,label="cf, exp_single")
    plt.legend()
    plt.savefig(save_dir+"fitting_before_Stim_200ms_10000ng.png")

    print "*************************************"
    """
    popt,pcov = optimize.curve_fit(exp_func2, X_bp, Y_bp)
    print popt
    print pcov

    popt,pcov = optimize.curve_fit(exp_func3, X_bp, Y_bp)
    print popt
    print pcov
    """
    X_ap = fit.psth_t[MAX_INDEX:50] # After Peak
    Y_ap = fit.psth_f[MAX_INDEX:50] # After Peak

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
    (a,b)=y_ap_all.shape
    y_ap_all=y_ap_all.reshape(a*b)

    #popt,pcov = optimize.curve_fit(exp_single, X_ap, Y_ap,p0=[70,190,170])
    popt,pcov = optimize.curve_fit(exp_single, x_ap_all, y_ap_all,p0=[70,190,170,0])
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/10mg_exp_single_Fall_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')    
    #global X,Y
    #X = X_ap
    #Y = Y_ap
    X = x_ap_all
    Y = y_ap_all
    parameter = [70,190,170,0]
    param_fin = optimize.fmin(cost_es, parameter)
    print parameter,param_fin
    np.savetxt(save_dir+"parameters/10mg_exp_single_Fall_fmin.txt",param_fin,fmt='%.6f',delimiter=',')    
    x_ap = np.arange(fit.psth_t[MAX_INDEX],5000,1)

    y_ap_fmin = exp_single(x_ap, *param_fin)
    y_ap_cf = exp_single(x_ap, *popt)

    fig = plt.figure()
    plt.plot(x_ap_all,y_ap_all,'o')
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,label="cf, exp_single")
    plt.legend()
    plt.xlim(0,5000)
    plt.savefig(save_dir+"fitting_after_Stim_200ms_10000ng.png")

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

    plt.savefig(save_dir+"fitting_Stim_200ms_10000ng_10000ms.png")
    plt.show()

def fitting_PSTH_5000ng():
    fit = fitting()
    fit.readPSTH("./bombykol_200ms/analyzed_data/psth/Dose5000_psth.dat")
    fit.readPSTH_ALL("./bombykol_200ms/analyzed_data/psth/Dose5000_psth_all.dat")
    #print fit.psth_f[:50]
    MAX_INDEX = fit.psth_f[:50].argmax()
    print MAX_INDEX, fit.psth_f[MAX_INDEX]
    X_bp = fit.psth_t[:MAX_INDEX+1] # Before Peak
    Y_bp = fit.psth_f[:MAX_INDEX+1] # Before Peak
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
    (a,b)=y_bp_all.shape
    y_bp_all=y_bp_all.reshape(a*b)

    #print x_bp_all
    #print y_bp_all
    
    init_param = [40,250,-150,0]
    #popt,pcov = optimize.curve_fit(exp_single_rise, X_bp, Y_bp, p0=init_param)
    #popt,pcov = optimize.curve_fit(exp_single_rise, x_bp_all, y_bp_all, p0=init_param)
    popt,pcov = optimize.curve_fit(exp_single, x_bp_all, y_bp_all, p0=init_param)
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/5000ng_exp_single_Rise_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')        


    global X,Y
    #X = X_bp
    #Y = Y_bp
    X = x_bp_all
    Y = y_bp_all
    print X,Y
    parameter = [40,250,-150,0]
    #param_fin = optimize.fmin(cost_esr, parameter)
    param_fin = optimize.fmin(cost_es, parameter)
    np.savetxt(save_dir+"/parameters/5000ng_exp_single_Rise_fmin.txt",param_fin,fmt='%.6f',delimiter=',')        
    print parameter, param_fin
    #x_fbp = np.arange(0,fit.psth_t[MAX_INDEX],0.01)
    x_bp = np.arange(0,300.0,0.01)
    #y_bp_cf = exp_single_rise(x_bp, *popt)
    y_bp_cf = exp_single(x_bp, *popt)
    #y_fbp = exp_single_pos(x_fbp, 1.0, 1.0, 1.0)
    #y_bp_fmin = exp_single_rise(x_bp, *param_fin)
    y_bp_fmin = exp_single(x_bp, *param_fin)
    para = [40,250,-150,0]
    y_bp_custom = exp_single(x_bp, *para)
    #y_bp_cf = exp_func1(x_fbp, *popt)
    #print x_fbp, y_fbp
    fig = plt.figure()
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'o')
    plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.xlim(0,300)
    plt.plot(x_bp, y_bp_cf,label="cf, exp_single")
    plt.legend()
    plt.savefig(save_dir+"fitting_before_Stim_200ms_5000ng.png")

    print "*************************************"
    """
    popt,pcov = optimize.curve_fit(exp_func2, X_bp, Y_bp)
    print popt
    print pcov

    popt,pcov = optimize.curve_fit(exp_func3, X_bp, Y_bp)
    print popt
    print pcov
    """
    X_ap = fit.psth_t[MAX_INDEX:50] # After Peak
    Y_ap = fit.psth_f[MAX_INDEX:50] # After Peak

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
    (a,b)=y_ap_all.shape
    y_ap_all=y_ap_all.reshape(a*b)

    #popt,pcov = optimize.curve_fit(exp_single, X_ap, Y_ap,p0=[70,190,170])
    popt,pcov = optimize.curve_fit(exp_single, x_ap_all, y_ap_all,p0=[70,190,170,0])
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/5000ng_exp_single_Fall_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')        
    #global X,Y
    #X = X_ap
    #Y = Y_ap
    X = x_ap_all
    Y = y_ap_all
    parameter = [70,190,170,0]
    param_fin = optimize.fmin(cost_es, parameter)
    print parameter,param_fin
    np.savetxt(save_dir+"parameters/5000ng_exp_single_Fall_fmin.txt",param_fin,fmt='%.6f',delimiter=',')        
    x_ap = np.arange(fit.psth_t[MAX_INDEX],5000,1)

    y_ap_fmin = exp_single(x_ap, *param_fin)
    y_ap_cf = exp_single(x_ap, *popt)

    fig = plt.figure()
    plt.plot(x_ap_all,y_ap_all,'o')
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,label="cf, exp_single")
    plt.legend()
    plt.xlim(0,5000)
    plt.savefig(save_dir+"fitting_after_Stim_200ms_5000ng.png")


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

    plt.savefig(save_dir+"fitting_Stim_200ms_5000ng_10000ms.png")
    plt.show()

def fitting_PSTH_2000ng():
    fit = fitting()
    fit.readPSTH("./bombykol_200ms/analyzed_data/psth/Dose2000_psth.dat")
    fit.readPSTH_ALL("./bombykol_200ms/analyzed_data/psth/Dose2000_psth_all.dat")
    #print fit.psth_f[:50]
    MAX_INDEX = fit.psth_f[:50].argmax()
    print MAX_INDEX, fit.psth_f[MAX_INDEX]
    X_bp = fit.psth_t[:MAX_INDEX+1] # Before Peak
    Y_bp = fit.psth_f[:MAX_INDEX+1] # Before Peak
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
    (a,b)=y_bp_all.shape
    y_bp_all=y_bp_all.reshape(a*b)

    #print x_bp_all
    #print y_bp_all
    
    init_param = [8,200,-50,0]
    #popt,pcov = optimize.curve_fit(exp_single_rise, X_bp, Y_bp, p0=init_param)
    #popt,pcov = optimize.curve_fit(exp_single_rise, x_bp_all, y_bp_all, p0=init_param)
    popt,pcov = optimize.curve_fit(exp_single, x_bp_all, y_bp_all, p0=init_param)
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/2000ng_exp_single_Rise_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')        
    global X,Y
    #X = X_bp
    #Y = Y_bp
    X = x_bp_all
    Y = y_bp_all
    print X,Y
    parameter = [40,250,-150,0]
    #param_fin = optimize.fmin(cost_esr, parameter)
    param_fin = optimize.fmin(cost_es, parameter)
    np.savetxt(save_dir+"/parameters/2000ng_exp_single_Rise_fmin.txt",param_fin,fmt='%.6f',delimiter=',')        
    print parameter, param_fin
    #x_fbp = np.arange(0,fit.psth_t[MAX_INDEX],0.01)
    x_bp = np.arange(0,300.0,0.01)
    #y_bp_cf = exp_single_rise(x_bp, *popt)
    y_bp_cf = exp_single(x_bp, *popt)
    #y_fbp = exp_single_pos(x_fbp, 1.0, 1.0, 1.0)
    #y_bp_fmin = exp_single_rise(x_bp, *param_fin)
    y_bp_fmin = exp_single(x_bp, *param_fin)
    para = [40,250,-150,0]
    y_bp_custom = exp_single(x_bp, *para)
    #y_bp_cf = exp_func1(x_fbp, *popt)
    #print x_fbp, y_fbp
    fig = plt.figure()
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_bp_all,y_bp_all,'o')
    plt.plot(x_bp, y_bp_fmin,'-',label="fmin, exp_single")
    plt.plot(x_bp, y_bp_custom,'-',label="custom,exp_single")
    plt.xlim(0,300)
    plt.plot(x_bp, y_bp_cf,label="cf, exp_single")
    plt.legend()
    plt.savefig(save_dir+"fitting_before_Stim_200ms_2000ng.png")

    print "*************************************"
    """
    popt,pcov = optimize.curve_fit(exp_func2, X_bp, Y_bp)
    print popt
    print pcov

    popt,pcov = optimize.curve_fit(exp_func3, X_bp, Y_bp)
    print popt
    print pcov
    """
    X_ap = fit.psth_t[MAX_INDEX:50] # After Peak
    Y_ap = fit.psth_f[MAX_INDEX:50] # After Peak

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
    (a,b)=y_ap_all.shape
    y_ap_all=y_ap_all.reshape(a*b)

    #popt,pcov = optimize.curve_fit(exp_single, X_ap, Y_ap,p0=[70,190,170])
    popt,pcov = optimize.curve_fit(exp_single, x_ap_all, y_ap_all,p0=[70,190,170,0])
    print popt
    print pcov
    np.savetxt(save_dir+"parameters/2000ng_exp_single_Fall_Curve_fit.txt",popt,fmt='%.6f',delimiter=',')        
    #global X,Y
    #X = X_ap
    #Y = Y_ap
    X = x_ap_all
    Y = y_ap_all
    parameter = [70,190,170,0]
    param_fin = optimize.fmin(cost_es, parameter)
    print parameter,param_fin
    np.savetxt(save_dir+"parameters/2000ng_exp_single_Fall_fmin.txt",param_fin,fmt='%.6f',delimiter=',')        
    x_ap = np.arange(fit.psth_t[MAX_INDEX],5000,1)

    y_ap_fmin = exp_single(x_ap, *param_fin)
    y_ap_cf = exp_single(x_ap, *popt)

    fig = plt.figure()
    plt.plot(x_ap_all,y_ap_all,'o')
    plt.plot(fit.psth_t,fit.psth_f)
    plt.plot(x_ap, y_ap_fmin,label="fmin, exp_single")
    plt.plot(x_ap, y_ap_cf,label="cf, exp_single")
    plt.legend()
    plt.xlim(0,5000)
    plt.savefig(save_dir+"fitting_after_Stim_200ms_2000ng.png")

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

    plt.savefig(save_dir+"fitting_Stim_200ms_2000ng_10000ms.png")
    plt.show()


#fitting_PSTH_10000ng()
#fitting_PSTH_5000ng()
fitting_PSTH_2000ng()
    
