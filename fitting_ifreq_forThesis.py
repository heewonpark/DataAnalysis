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
    return sum((gamma(x,*g_para)-y)**2 for x, y in zip(hx[1:40],hy[1:]))
def gaussian(x, sigma, mu,alpha):
    return alpha/(sigma*np.sqrt(np.pi*2))*sp.exp(-pow((x-mu),2)/(2*pow(sigma,2)))
def cost_gaussian(gaus_para):
    return sum((gaussian(x, *gaus_para)-y)**2 for x, y in zip(GX, GY))
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


readfilelist = open("listoffiles.txt",'r')
filenames = readfilelist.readlines()
readfilelist.close()
filename = [' ' for i in range(6)]
dummy    = [' ' for i in range(6)]

for i in range(6):
    filename[i], dummy[i] = filenames[i].split(' ')
    filename[i] = filename[i].rstrip('.atf')
    print filename[i]

readdatfile = open("./analyzed_data/average.dat",'r')
data = readdatfile.readlines()
readdatfile.close()

X = sp.zeros(81)
Y = sp.zeros(81)
for i in range(0,81):
    X[i], Y[i] = data[i+71].split('\t')
    X[i] = float(X[i])
    Y[i] = float(Y[i])

param_fin = [0.0001333, 12.91, 0.5157,5.639, 4.726]
param_fin = scipy.optimize.fmin(cost_fn2, param_fin,xtol=1e-8)
print param_fin

graph = pylab.figure()
pylab.plot(X,Y,'g', label='data', linewidth=1.0)
fX = X
fY1 = [expFunc2(x, *param_fin) for x in fX]
pylab.plot(fX, fY1, 'b--', label = 'fitting', linewidth = 1.0)
graph.set_ylim = (0,300)
pylab.xlabel('X axis')
pylab.ylabel('Y axis')
pylab.title('Parameter optimization')
pylab.legend()
pylab.grid(True)
pylab.savefig("fitting.png")

residual_Y = [None for i in range(1000)]
residual_kikaku = [None for i in range(1000)]
residual_kikaku_time = [None for i in range(1000)]
residual_time = [None for i in range(1000)]
cnt = 0
for i in range(6):
    readifreqfile= open('./analyzed_data/'+filename[i]+'ifreq.dat','r')
    ifreqdata = readifreqfile.readlines()
    steps = int(ifreqdata[0].rstrip('\n'))
    delay = 0.0
    Xe = sp.zeros(steps)
    Ye = sp.zeros(steps)
    print steps
    print filename[i]
    for j in range(0,steps):
        Xe[j],Ye[j] = ifreqdata[j+1].split('\t')
        Xe[j] = float(Xe[j])
        Ye[j] = float(Ye[j])

    delay = Xe[0]-6.0
    for k in range(steps):
        if Ye[k]<300:
            funcresult = expFunc2(Xe[k]-delay, *param_fin)
            residual_Y[cnt] = Ye[k]-funcresult
            residual_time[cnt] = 1/Ye[k]-1/funcresult
            residual_kikaku[cnt] =  (Ye[k]-funcresult)/funcresult
            residual_kikaku_time[cnt] =  (1/Ye[k]-1/funcresult)*funcresult
            if (residual_kikaku[cnt]<-0.5):
                print cnt, residual_kikaku[cnt], Ye[k], funcresult
            #print Xe[k], Ye[k], funcresult
            cnt +=1
      
    fig3 = pylab.figure()
    Xee=  [Xe[k]-delay for k in range(steps)]
    Yee = [expFunc2(Xe[k],*param_fin) for k in range(steps)]
    pylab.plot(Xee,Yee)
    pylab.plot(Xee,Ye)
    pylab.ylim(0,300)
    pylab.savefig(filename[i]+'zansa.png')

    #### Graph for Presentation
    #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
    plt.rcParams['font.size'] = 20 #フォントサイズを設定
    plt.rcParams['axes.linewidth'] = 2.0 #軸の太さを設定。目盛りは変わらない
    #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
    #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    
    fig = plt.figure(figsize=(5,4),dpi=250)
    fig.subplots_adjust(bottom=0.2,left =0.20)
    ax = fig.add_subplot(111)
    fig.patch.set_alpha(0.0)
    plt.plot(Xee,Yee,'r',label='Fitted Curve')
    plt.plot(Xee,Ye,'g',label='ISF')
    plt.ylim(0,300)
    plt.ylabel('ISF[Hz]')
    plt.xlabel('Time[s]')
    plt.legend(frameon=False,fontsize=15)
    plt.savefig(filename[i]+'zansa_fp.png')
    


text = 'a '+ repr(param_fin[0])+ ' b '+ repr(param_fin[1])+ ' c '+ repr(param_fin[2])+'\n'+'d ' + repr(param_fin[3])+ ' e '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*exp(-(x-b)/e)'

hist = pylab.figure()
#pylab.text(-1.0,13,text,fontsize=10)
histreturn = pylab.hist(residual_Y[0:cnt], bins=100, range=(-100,250))
pylab.savefig("fitting_ifreqresidualhistogram.png")

BINS = 125
hist_kikaku = pylab.figure()
histreturn = pylab.hist(residual_kikaku[0:cnt], bins=BINS, range=(-1.0,24.0))
DX1 = histreturn[1][0:5]
DY1 = histreturn[0][0:5]
pylab.savefig("fitting_ifreqresidual_kikaku.png")

hist_kikaku_time = pylab.figure()
histreturn_time = pylab.hist(residual_kikaku_time[0:cnt], bins=BINS)
pylab.savefig("fitting_ifreqresidual_kikaku_time.png")


print DX1
print DY1


fig = plt.figure(figsize=(10,8),dpi=400) 
plt.rcParams['font.size']=15
de1_para_ini = [3.4, -2.05, 0.55]
de1_para_fin = scipy.optimize.fmin(cost_de1, de1_para_ini)
de1graph = pylab.figure(figsize=(10,7),dpi=400)
plt.bar(DX1,DY1,width = 0.2,color="#3b3b3b",alpha=0.7,label=r'Residual histogram')
dx1_fit = sp.arange(-1.0, -0.14, 0.05)
dy1_fit = [distribution_exp1(x,*de1_para_fin) for x in dx1_fit]
dx1_fit2 = sp.append(-1.0,dx1_fit)
dy1_fit2 = sp.append(0,dy1_fit)
print de1_para_fin
plt.plot(dx1_fit2, dy1_fit2,'g-',linewidth=2.0)
#plt.fill(dx1_fit2, dy1_fit2,'g',alpha=0.8)

DX2 = histreturn[1][5:BINS]
DY2 = histreturn[0][5:]
plt.bar(DX2,DY2,width = 0.2,color="#3b3b3b",alpha=0.7)
de3_para_fin =[2.5,3.4,0.96,-0.09,1.99]
de3_para_fin = scipy.optimize.fmin(cost_de3, de3_para_fin)
dx3_fit = sp.arange(-0.14, 24.0, 0.05)
dy3_fit = [distribution_exp3(x,*de3_para_fin) for x in dx3_fit]
print de3_para_fin
plt.plot(dx3_fit, dy3_fit,'g-',linewidth=2.0,label=r'$R(x)$')
#plt.fill(dx3_fit, dy3_fit,'g',alpha=0.8)
plt.fill_between(dx1_fit2,0,dy1_fit2,color='g',alpha=0.7)
plt.fill_between(dx3_fit,0,dy3_fit,color='g',alpha=0.7)
plt.xlim(-2,15)
plt.ylim(0,120)
plt.legend()
plt.xlabel('Normalized residual',fontsize=30)
plt.ylabel('Counts',fontsize=30)
plt.savefig("fitting_de.png",transparent=True)


#### Graph for Presentation
#plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['font.size'] = 20 #フォントサイズを設定
plt.rcParams['axes.linewidth'] = 2.0 #軸の太さを設定。目盛りは変わらない
#plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                                         
#plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     

#graph = plt.figure(figsize=(6,4.5),dpi=250)

#fig = plt.figure(figsize=(10,8),dpi=400) 
#plt.rcParams['font.size']=15
de1_para_ini = [3.4, -2.05, 0.55]
de1_para_fin = scipy.optimize.fmin(cost_de1, de1_para_ini)
#de1graph = pylab.figure(figsize=(10,7),dpi=400)
fig = plt.figure(figsize=(5,4),dpi=250)
fig.subplots_adjust(bottom=0.2,left =0.20)
ax = fig.add_subplot(111)
fig.patch.set_alpha(0.0)

plt.bar(DX1,DY1,width = 0.2,color="#3b3b3b",alpha=0.7,label=r'Residual')
dx1_fit = sp.arange(-1.0, -0.14, 0.05)
dy1_fit = [distribution_exp1(x,*de1_para_fin) for x in dx1_fit]
dx1_fit2 = sp.append(-1.0,dx1_fit)
dy1_fit2 = sp.append(0,dy1_fit)
print de1_para_fin
plt.plot(dx1_fit2, dy1_fit2,'g-',linewidth=2.0)
#plt.fill(dx1_fit2, dy1_fit2,'g',alpha=0.8)

DX2 = histreturn[1][5:BINS]
DY2 = histreturn[0][5:]
plt.bar(DX2,DY2,width = 0.2,color="#3b3b3b",alpha=0.7)
de3_para_fin =[2.5,3.4,0.96,-0.09,1.99]
de3_para_fin = scipy.optimize.fmin(cost_de3, de3_para_fin)
dx3_fit = sp.arange(-0.14, 24.0, 0.05)
dy3_fit = [distribution_exp3(x,*de3_para_fin) for x in dx3_fit]
print de3_para_fin
plt.plot(dx3_fit, dy3_fit,'g-',linewidth=2.0,label=r'$R(x)$')
#plt.fill(dx3_fit, dy3_fit,'g',alpha=0.8)
plt.fill_between(dx1_fit2,0,dy1_fit2,color='g',alpha=0.7)
plt.fill_between(dx3_fit,0,dy3_fit,color='g',alpha=0.7)
plt.xlim(-2,15)
plt.ylim(0,120)
plt.legend(frameon=False,fontsize=15)
plt.xlabel('Normalized residual')
plt.ylabel('Counts')
plt.savefig("fitting_de_fp.png",transparent=True)


xarr = sp.arange(-0.5,24.0,0.1)
for i in xarr:
    if distribution_exp3(i,*de3_para_fin) < 0:
        print i

save_para = open("saveParameter.dat",'w')
for i in range(len(param_fin)):
    save_para.writelines(str(param_fin[i])+'\t')
save_para.write('\n')
for i in range(len(de1_para_fin)):
    save_para.writelines(str(de1_para_fin[i])+'\t')
save_para.write('\n')
for i in range(len(de3_para_fin)):
    save_para.writelines(str(de3_para_fin[i])+'\t')
save_para.write('\n')
save_para.close()
#pylab.show()
