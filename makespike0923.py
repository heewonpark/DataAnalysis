from matplotlib import pylab
import scipy as sp
import scipy.optimize
import numpy as np
import math


read_num  = open("save_filenumber.dat",'r')
num_dat   = read_num.readlines()
read_num.close()
file_num  = num_dat[0].rsplit('\n')
file_num  = int(file_num[0])
print file_num

read_para = open("saveParameter.dat",'r')
para_dat = read_para.readlines()
read_para.close()

C = 0.0

de_para = [0 for _ in range(5)]
de_para = para_dat[0].split('\t')
de_para.remove('\n')
for i in range(5):
    de_para[i] = float(de_para[i])
    print de_para
    
cf1_para = [0 for _ in range(3)]
cf1_para = para_dat[1].split('\t')
cf1_para.remove('\n')
for i in range(3):
    cf1_para[i] = float(cf1_para[i])
    print cf1_para
            
cf2_para = [0 for _ in range(5)]
cf2_para = para_dat[2].split('\t')
cf2_para.remove('\n')
for i in range(5):
    cf2_para[i] = float(cf2_para[i])
    print cf2_para
    
def double_exp(x, alpha1, beta, gamma1, alpha2, gamma2):
    return alpha1 * sp.exp(-(x-beta)/gamma1)  +  alpha2 * sp.exp(-(x-beta)/gamma2)
    
def cum_Func1(x, alpha, beta, gamma):
    return alpha * gamma * sp.exp((x-beta)/gamma) - alpha * gamma * sp.exp((-1-beta)/gamma)

def distribution_exp1(x,alpha,beta,gamma):
    return alpha * sp.exp((x-beta)/gamma)

def distribution_exp3(x, alpha, beta, gamma, a, b):
    return alpha * sp.exp(-(x-beta)/gamma)+a*x+b

def distribution_exp3_(x, de3_para):
    return distribution_exp3(x,*de3_para)


find_eq_para = []
find_eq_para.append(cf1_para)
find_eq_para.append(cf2_para[0:5])
print find_eq_para
def find_de1eqde3(x,find_eq_para):
    return distribution_exp1(x,*find_eq_para[0]) - distribution_exp3(x,*find_eq_para[1])
cf1eqcf2 = scipy.optimize.fsolve(find_de1eqde3 ,-1,args=find_eq_para)
print cf1eqcf2

C = cum_Func1(cf1eqcf2, *cf1_para)
cf2_para.append(C)
print 'C '
print C
def cum_Func2(x, alpha, beta, gamma, a, b, c):
    return -alpha * gamma * sp.exp(-(x-beta)/gamma)+alpha*gamma*sp.exp(-(cf1eqcf2-beta)/gamma) + a *(x**2)/2 + b*x + c - a*(cf1eqcf2**2)/2 - b*cf1eqcf2

print 'cf2'
print cf2_para[0:5]
x_max = scipy.optimize.fsolve(distribution_exp3_, 10, args=cf2_para[0:5])
print x_max
print cf2_para[5]
print cum_Func1(cf1eqcf2, *cf1_para),cum_Func2(cf1eqcf2,*cf2_para)
cf_para = [None for _ in range(3)]
cf_para[0]= cf1_para
cf_para[1]= cf2_para
def cum_Func(x, cf_para):
    if x<=cf1eqcf2:
        return cum_Func1(x, *cf_para[0])-cf_para[2]
    elif (x<=x_max)&(x>cf1eqcf2):
        return cum_Func2(x, *cf_para[1])-cf_para[2]
    elif x>x_max:
        return -10

cf_para[2] = 0.0
print cf_para
cumulative_area = cum_Func(x_max, cf_para)
print cumulative_area
spike_t = [None for _ in range(500)]
spike_t_now = 0.0
spike_t[0] = 6.0
spiketimef = open("spiketiming"+str(file_num)+".dat",'w')

fig =pylab.figure()
xarr = sp.arange(-1.0, 50, 0.1)
y1 = [cum_Func(i, cf_para) for i in xarr]
#y2 = [cum_Func2(i,*cf_para[1]) for i in xarr]
x1 = sp.arange(-1.0, 0, 0.05)
y3 = [distribution_exp1(i,*cf1_para) for i in x1]
y4 = [distribution_exp3(i,*cf2_para[0:5]) for i in xarr]
pylab.plot(xarr, y1)
#pylab.plot(xarr, y2)
pylab.plot(x1, y3)
pylab.plot(xarr,y4)
pylab.grid()

pylab.savefig("makespikecheck.png")

cnt = 0
for i in range(500):
    if spike_t[i]>14.5:
        break

    if i<500-1:
        random = cumulative_area * np.random.random()
        cf_para[2] = random
       
        err = scipy.optimize.fsolve(cum_Func,-1,args=cf_para)
        print i, random, err
        
        if err<-1:
            print "XXXXXXXXXXXXXXXXXXXXXXXXXX"
        spike_t[i+1] = spike_t[i] + 1/(double_exp(spike_t[i],*de_para)*(1+err))
        print spike_t[i+1], double_exp(spike_t[i],*de_para),cf_para[2], err
        spiketimef.writelines(str(float(spike_t[i+1]))+'\n')
        cnt +=1

spiketimef.writelines(str(cnt)+'\n')
spiketimef.close()

write_numfile = open("save_filenumber.dat",'w')
write_numfile.write(str(file_num+1)+'\n')
write_numfile.close()
#pylab.show()
