from matplotlib import pylab
import scipy as sp
import scipy.optimize

readdatfile = open("average.dat",'r')
data = readdatfile.readlines()
readdatfile.close()

X = sp.zeros(89)
Y = sp.zeros(89)
for i in range(0,89):
    X[i], Y[i] = data[i+71].split('\t')
    X[i] = float(X[i])
    Y[i] = float(Y[i])
 
print 'X'
print X 
print 'Y'
print Y

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


param_fin = [0.0000795, 12.8708, 0.4938, 5.1405, 3.819]
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

residual_Y1 = [(Y[i]-fY1[i]) for i in range(89)]
print residual_Y1
hist = pylab.figure()
#text = 'alpha '+ repr(param_fin[0])+ ' beta '+ repr(param_fin[1])+ ' gamma '+ repr(param_fin[2])+'\n'+'a ' + repr(param_fin[3])+ ' b '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*x+e'
text = 'a '+ repr(param_fin[0])+ ' b '+ repr(param_fin[1])+ ' c '+ repr(param_fin[2])+'\n'+'d ' + repr(param_fin[3])+ ' e '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*exp(-(x-b)/e)'

print text
pylab.text(-1.0,13,text,fontsize=10)
pylab.hist(residual_Y1, bins=15)
pylab.savefig("residual_histogram.png")
pylab.show()
