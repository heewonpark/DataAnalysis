from matplotlib import pylab
import scipy as sp
import scipy.optimize

readdatfile = open("average.dat",'r')
data = readdatfile.readlines()
readdatfile.close()

X = sp.zeros(89)
Y = sp.zeros(89)
for i in range(71,160):
    X[i-71], Y[i-71] = data[i].split('\t')
    X[i-71] = float(X[i-71])
    Y[i-71] = float(Y[i-71])

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


param_ini = [0.25, 11.9, 0.7, -1.06, 19.7]
param_fin = scipy.optimize.fmin(cost_fn3, param_ini)
print param_fin

pylab.plot(X,Y,'g', label='data', linewidth=1.0)

fX = X
fY = [expFunc3(x, *param_fin) for x in fX]
pylab.plot(fX, fY, 'b--', label = 'fitting', linewidth = 1.0)

pylab.xlabel('X axis')
pylab.ylabel('Y axis')
pylab.title('Parameter optimization')
pylab.legend()
pylab.grid(True)

pylab.savefig("fitting.png")
pylab.show()
