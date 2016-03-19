import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import math
#data = np.loadtxt("MaximumFreq_ByGaussianKernel.txt",float)
data = np.loadtxt("./saved_parameter/MaximumFreq_ByPSTH.txt",float)
x = data[:,0]
y = data[:,1]
print x
print y

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

fig = plt.figure(figsize=(10,8),dpi=400)
plt.rcParams['font.size']=15
plt.plot(x,y,'go',label="Maximum Frequency")
x_tmp = x
y_tmp = y
x = [float(xn) for xn in x]
y = [float(yn) for yn in y]
x = np.array(x)
y = np.array(y)

popt,pcov = curve_fit(Michaelis_Menten_1,x,y)
#Parameter = popt
print"[CURVE_FIT REUSLT] k = %f, Tau_max = %f"%(popt[0],popt[1])
#np.savetxt("Michaelis-Menten_Parameter.txt",popt)
np.savetxt("Michaelis-Menten_ParameterPSTH.txt",popt)
Parameter=np.append(popt,1.0)
print Parameter
result = leastsq(fit_func,Parameter,args=(x,y))
print"[LEASTSQ RESULT] k = %f, Tau_max = %f, n = %f\n"%(result[0][0],result[0][1],result[0][2])
#np.savetxt("Michaelis-Menten_Parameter_n.txt",result[0])
np.savetxt("Michaelis-Menten_Parameter_nPSTH.txt",result[0])
x = [10 for _ in range(100)]
x_exponent = [float(i/14.5) for i in range(100)]
x = np.array(x)
x_exponent = np.array(x_exponent)
#print x, x_exponent
x_ = np.power(x,x_exponent)
#print x_
#fig1 =plt.figure()
print "30 ",Michaelis_Menten_2(30,result[0])
print "100 ",Michaelis_Menten_2(100,result[0])
print "300 ",Michaelis_Menten_2(300,result[0])
print "1000 ",Michaelis_Menten_2(1000,result[0])
print "3000 ",Michaelis_Menten_2(3000,result[0])
print "10000 ",Michaelis_Menten_2(10000,result[0])

#plt.plot(x_,Michaelis_Menten_2(x_,result[0]),label="Michaelis-Menten Fitted Curve(n)")
plt.plot(x_,Michaelis_Menten_2(x_,result[0]),'r-',label=r"$M(c)$")
#plt.plot(x_,Michaelis_Menten_1(x_,*popt),label="Michaelis-Menten Fitted Curve")
#plt.legend(loc='upper left')
plt.xscale('log')
plt.xlabel('Dose[ng]',fontsize=30)
plt.ylabel('Maximum frequency[Hz]',fontsize=30)
plt.legend(loc='upper left')
#plt.savefig('Michaelis-menten_fitted_curve2.png')
plt.savefig('Michaelis-menten_fitted_curve2PSTH.png',transparent=True)

#plt.show()

#### Graph for Presentation
#plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 2.0
#plt.rcParams['xtics.major.size'] = 10
#plt.rcParams['xtics.major.width'] = 1.5

fig = plt.figure(figsize=(5,4),dpi=250)
fig.subplots_adjust(bottom=0.2,left =0.20)
ax = fig.add_subplot(111)
fig.patch.set_alpha(0.0)

plt.plot(x_tmp,y_tmp,'ko')
plt.plot(x_,Michaelis_Menten_2(x_,result[0]),'r-',label=r"$M(c)$")
plt.xscale('log')
plt.xlabel('Dose[ng]')
plt.ylabel('Maximum frequency[Hz]')
print np.logspace(1,7,num=3)
plt.xticks(np.logspace(0,6,num=3))
plt.yticks(np.arange(0,210,50))
plt.ylim(0,200)
plt.legend(loc='upper left',frameon=False,fontsize=15)
#plt.savefig('Michaelis-menten_fitted_curve2.png')
plt.savefig('Michaelis-menten_fitted_curve2PSTH_fp.png',transparent=True)
