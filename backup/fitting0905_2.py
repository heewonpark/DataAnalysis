from matplotlib import pylab
import scipy as sp
import scipy.optimize


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


readfilelist = open("listoffiles.txt",'r')
filenames = readfilelist.readlines()
readfilelist.close()
filename = [' ' for i in range(6)]
dummy    = [' ' for i in range(6)]

for i in range(6):
    filename[i], dummy[i] = filenames[i].split(' ')
    filename[i] = filename[i].rstrip('.atf')
    print filename[i]

readdatfile = open("average.dat",'r')
data = readdatfile.readlines()
readdatfile.close()

X = sp.zeros(81)
Y = sp.zeros(81)
for i in range(0,81):
    X[i], Y[i] = data[i+71].split('\t')
    X[i] = float(X[i])
    Y[i] = float(Y[i])
"""    
print 'X'
print X 
print 'Y'
print Y
"""
param_fin = [0.0001333, 12.91, 0.5157,5.639, 4.726]
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

residual_Y1 = [0 for i in range(81)]

for i in range(6):
    readbinfile= open(filename[i]+'bin.dat','r')
    bindata = readbinfile.readlines()
    steps, delay = bindata[0].split(' ')
    steps = int(steps)
    delay = int(delay)

    Xe = sp.zeros(steps+10+delay)
    Ye = sp.zeros(steps+10+delay)
    print steps
    print filename[i]
    for j in range(steps+10+delay):
        Xe[j],Ye[j] = bindata[j+1].split('\t')
        Xe[j] = float(Xe[j])
        Ye[j] = float(Ye[j])
        print str(Xe[j]) +'\t'+str(Ye[j])+'\n'
    

    for k in range(81):
        if (k+71)<steps+10+delay:
            residual_Y1[k] += (Ye[k+71]-fY1[k])
  
hist = pylab.figure()
#text = 'alpha '+ repr(param_fin[0])+ ' beta '+ repr(param_fin[1])+ ' gamma '+ repr(param_fin[2])+'\n'+'a ' + repr(param_fin[3])+ ' b '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*x+e'
text = 'a '+ repr(param_fin[0])+ ' b '+ repr(param_fin[1])+ ' c '+ repr(param_fin[2])+'\n'+'d ' + repr(param_fin[3])+ ' e '+repr(param_fin[4])+'\nf(x)=a*exp(-(x-b)/c)+d*exp(-(x-b)/e)'

print text
pylab.text(-1.0,13,text,fontsize=10)
pylab.hist(residual_Y1, bins=50)
pylab.savefig("fitting_residualhistogramKikaku.png")
pylab.show()
