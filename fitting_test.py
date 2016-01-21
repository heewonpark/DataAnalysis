#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2015.03.23
###############################################

from matplotlib import pylab
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize
import numpy as np
import math 
import csv

#from functions import functions as fc
def exp_func1(x, alpha, beta, gamma):
    return alpha * sp.exp(-(x-beta)/gamma)
def exp_func2(x, alpha1, beta, gamma1, alpha2, gamma2):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  alpha2 * sp.exp(-(x-beta)/gamma2)
"""
def exp_func3(x, alpha1, beta, gamma1, a, b):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  a*x+b
"""

def exp_func3(x, alpha1, beta, gamma1, a, b):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  + b*x +a

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


def cost_fn3_www(param):
    sum_ = 0
    for i in range(len(X)):
        sum_+=(exp_func3(X[i], *param)-Y[i])**2
    return sum_

def cost_fn2_www(param):
    sum_ = 0
    for i in range(len(X)):
        sum_+=(exp_func2(X[i], *param)-Y[i])**2
    return sum_

def cost_fn1_www(param):
    sum_ = 0
    for i in range(len(X)):
        sum_+=(exp_func1(X[i], *param)-Y[i])**2
    return sum_

#---------------------------------------------
# READ FILE LIST

# lof List of Files
data = csv.reader(open("List_Of_Files.csv",'r'))
filename = []
for D in data:
    filename.append(D[0].rstrip('.atf'))
print filename

#---------------------------------------------
# READ average.dat FILE

ReadFile = open("average.dat",'r')
data = ReadFile.readlines()
ReadFile.close()

X = sp.zeros(81)
Y = sp.zeros(81)
for i in range(0,81):
    X[i], Y[i] = data[i+71].split('\t')
    #X[i] = float(X[i])
    X[i] = float(X[i]-6.1)
    Y[i] = float(Y[i])

print 'X'
print X 
print 'Y'
print Y


#parameter1 = [5,0,10]
#parameter2 = [1,1,1,1,1]
#parameter3 = [1,1,1,1,1]
"""
#---------------------------------------------
# LOAD PARAMETERS FOR FUNCTION
#parameter = np.loadtxt("exponential_function_parameters.txt",float)
parameter1 = np.loadtxt("exponential_function_parameters_result_func1.txt",float)
#param_fin1 = scipy.optimize.fmin(cost_fn1, parameter1,xtol=1e-8,ftol=0.0001,retall=True)
param_fin1 = scipy.optimize.fmin(cost_fn1, parameter1)
print param_fin1
print "result: ", cost_fn1(param_fin1)
np.savetxt("exponential_function_parameters_result_func1.txt",param_fin1)

parameter2 = np.loadtxt("exponential_function_parameters_result_func2.txt",float)
parameter2 = [0.5,1,1,10,1]
param_fin2 = scipy.optimize.fmin(cost_fn2, parameter2)
print param_fin2
print "result: ", cost_fn2(param_fin2)
np.savetxt("exponential_function_parameters_result_func2.txt",param_fin2)

parameter3 = np.loadtxt("exponential_function_parameters_result_func3.txt",float)
param_fin3 = scipy.optimize.fmin(cost_fn3, parameter3)
print param_fin3
print "result: ", cost_fn3(param_fin3)
np.savetxt("exponential_function_parameters_result_func3.txt",param_fin3)
"""
param_fin1 = [10,20,10]
#param_fin2 = [0.00005,7.0,0.501,5.199,4.17]
param_fin2 = [26,0.3,4.0,37,0.5]
param_fin3 = [1,1,1,1,1]
print "result: ", cost_fn1(param_fin1)
print "result: ", cost_fn2(param_fin2)
print "result: ", cost_fn3(param_fin3)


print "--------Single----------"
parameter1 = np.loadtxt("exponential_function_parameters_result_func1.txt",float)
param_fin1,pcov1 = scipy.optimize.curve_fit(exp_func1,X,Y,p0=parameter1)
print param_fin1
print "result: ", cost_fn1(param_fin1)
print "result_www: ", cost_fn1_www(param_fin1)
perr1 = np.sqrt(np.diag(pcov1))
#print perr1
np.savetxt("exponential_function_parameters_result_func1.txt",param_fin1)

print "--------Double----------"
parameter2 = np.loadtxt("exponential_function_parameters_result_func2.txt",float)
#param_fin2,pcov2 = scipy.optimize.curve_fit(exp_func2,X,Y,p0=parameter2)
param_fin2,pcov2 = scipy.optimize.curve_fit(exp_func2,X,Y)
print param_fin2
print "result: ", cost_fn2(param_fin2)
print "result_www: ", cost_fn2_www(param_fin2)
#perr2 = np.sqrt(np.diag(pcov2))
#print perr2
np.savetxt("exponential_function_parameters_result_func2.txt",param_fin2)

print "--------Single+linear----------"
parameter3 = np.loadtxt("exponential_function_parameters_result_func3.txt",float)
#param_fin3,pcov3 = scipy.optimize.curve_fit(exp_func3,X,Y,p0=parameter3)
param_fin3,pcov3 = scipy.optimize.curve_fit(exp_func3,X,Y)
print param_fin3
print "result: ", cost_fn3(param_fin3)
print "result_www: ", cost_fn3_www(param_fin3)
#perr3 = np.sqrt(np.diag(pcov3))
#print perr3
np.savetxt("exponential_function_parameters_result_func3.txt",param_fin3)


#graph = plt.figure(figsize=(10,8),dpi=400)
graph = plt.figure()
#plt.rcParams['font.size']=20
#graph.patch.set_alpha(0.0)
#pylab.plot(X,Y,'g', label='data', linewidth=1.0)
plt.bar(X,Y,width=0.1,color="#3b3b3b",label='PSTH')
fX = X
fY1 = [exp_func1(x, *param_fin1) for x in fX]
fY2 = [exp_func2(x, *param_fin2) for x in fX]
fY3 = [exp_func3(x, *param_fin3) for x in fX]

plt.plot(fX, fY1, 'r-', label = r"$F1$", linewidth = 3.0)
plt.plot(fX, fY2, 'b-', label = r"$F2$", linewidth = 3.0)
plt.plot(fX, fY3, 'g-', label = r"$F3$", linewidth = 3.0)
plt.xlabel('Time[s]')
plt.ylabel('Frequency[Hz]')
#pylab.title('Parameter optimization')
plt.xlim(0,8)
plt.ylim(0,120)
plt.legend(fontsize=20)
plt.grid(True)
plt.savefig("fitting_test.png")
plt.show()
