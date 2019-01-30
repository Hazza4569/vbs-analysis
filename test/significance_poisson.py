import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit

sig,diff = np.loadtxt('/home/user108/y4p/graph_logs/30-01/significance_dsc_sigs.dat',skiprows=0,unpack=True,delimiter=',')
n=len(sig)

def func(x,a,b,c,d):
    return a*(b**(c*x)) + d

x = np.linspace(1,11,11)
x2 = np.linspace(1,11,100)
y = diff[1:n]

mpl.plot(x,y,'k.')

popt, pcov = curve_fit(func,x,y)
print popt

mpl.plot(x,func(x,*popt), 'r-')

mpl.show()

