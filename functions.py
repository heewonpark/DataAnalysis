#!#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# 2015.08.04
###############################################

class functions:
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
    def cost_gaussian(gaus_para,X,Y):
        return sum((gaussian(x, *gaus_para)-y)**2 for x, y in zip(X, Y))
    def distribution_exp1(x,alpha,beta,gamma):
        return alpha * sp.exp((x-beta)/gamma)
    def cost_de1(de1_para):
        return sum((distribution_exp1(x,*de1_para)-y)**2 for x,y in zip(DX1, DY1))
    def distribution_exp2(x,alpha1,beta,gamma1,alpha2,gamma2,a,b):
        return alpha1 *sp.exp(-(x-beta)/gamma1)+alpha2 * sp.exp(-(x-beta)/gamma2)+a*x+b
    def cost_de2(de1_para):
        return sum((distribution_exp2(x,*de1_para)-y)**2 for x,y in zip(DX2, DY2))
    def distribution_exp3(x,alpha,beta,gamma,a,b):
        return alpha * sp.exp(-(x-beta)/gamma)+a*x+b
    def cost_de3(de3_para):
        return sum((distribution_exp3(x,*de3_para)-y)**2 for x,y in zip(DX2, DY2))
    def distribution_exp4(x,alpha,beta,gamma,c):
        return alpha * sp.exp(-(x-beta)/gamma)+c
    def cost_de4(de4_para):
        return sum((distribution_exp4(x,*de4_para)-y)**2 for x,y in zip(DX2, DY2))
