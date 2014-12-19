from matplotlib import pylab
import scipy as sp
import scipy.optimize

X = sp.arange(3,19)
Y = sp.zeros((len(X),))

for i in range(3000):
    Y[sum(sp.random.randint(1,7) for i in range(3))-3] +=1

def gauss(x, amp, mean, sigma):
    return amp*sp.exp(-(mean - x)**2/sigma)

def cost_fn(param):
    return sum((gauss(x, *param)-y)**2 for x,y in zip(X,Y))

param_ini = [100.0, 5.0, 3.0]
param_fin = scipy.optimize.fmin(cost_fn, param_ini, xtol=1e-8)

def cost_fn_w(param):
    return sum(((gauss(x,*param)-y)/y)**2 for x,y in zip(X,Y))

param_fin_w = scipy.optimize.fmin(cost_fn_w, param_ini,xtol=1e-8)

mean_min, mean_max = 7.0, 9.0

def limit_param(param):
    amp, mean, sigma = param
    mean = mean_min + (mean_max - mean_min)/(1+sp.exp(mean))
    return (amp, mean, sigma)

def cost_fn_l(param):
    param = limit_param(param)
    return sum(((gauss(x, *param)-y)/y)**2 for x,y in zip(X,Y))

param_fin_l = scipy.optimize.fmin(cost_fn_l, param_ini,xtol=1e-8)
param_fin_l = limit_param(param_fin_l)

pylab.plot(X,Y,'go', label='data', linewidth=1.0)

fX = sp.arange(2.0, 19.0, 0.1)
fY = [gauss(x, *param_fin) for x in fX]
pylab.plot(fX, fY, 'b--', label = 'fitting', linewidth = 1.0)
fY = [gauss(x, *param_fin_w) for x in fX]
pylab.plot(fX, fY, 'r'  , label = 'weighted', linewidth = 1.0)
fY = [gauss(x, *param_fin_l) for x in fX]
pylab.plot(fX, fY, 'b-.', label = 'limited', linewidth = 1.0)

pylab.xlabel('X axis')
pylab.ylabel('Y axis')
pylab.title('Parameter optimization')
pylab.legend()
pylab.grid(True)

pylab.semilogy()
pylab.axis([2.0,19.0,1.0,500])
pylab.savefig("fittingSample.png")
pylab.show()
