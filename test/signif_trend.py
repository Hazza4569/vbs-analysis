import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
import glob

#data_files = glob.glob('/home/user108/y4p/graph_logs/24-02_optsigs/opt_sig_*.dat');
data_files = glob.glob('/home/user108/y4p/graph_logs/15-03_optsigs/opt_sig_*.dat');
data_files = glob.glob('/home/user108/y4p/graph_logs/27-03/opt_sig_*.dat');
lum = []
sigf = []
sigf_std = []
mjj = []
njj = []
ns = []
nb = []

for iFile in data_files:
    a,b,c,d,e,f,g = np.loadtxt(iFile,skiprows=1,unpack=True)
    lum.append(a)
    sigf.append(b)
    sigf_std.append(c)
    mjj.append(d)
    njj.append(e)
    ns.append(f)
    nb.append(g)

print lum
print sigf
print sigf_std
print mjj
print njj

xmin = 2e1
xmax = 5000
ymin = 0
ymax = 12

hel = {'fontname':''}

lum_label_size = 9
lumlinestyle = {'color':'black','linestyle':'dotted','linewidth':1.5}
siglinestyle = {'color':'grey','linestyle':'dashed','alpha':0.5,'linewidth':1.5}
datastyle = {'marker':'o','color':'black','linestyle':'','markersize':6}
errstyle = {'marker':'o','color':'black','linestyle':'','markersize':6,'capsize':2}
legstyle = {'loc':2,'bbox_to_anchor':(.05,.95),'framealpha':0}

mpl.figure()
mpl.loglog(1,10,label='Data',**datastyle);
#mpl.semilogx(1,10,label='Data',**datastyle);
mpl.semilogx([xmin,xmax],[5,5],**siglinestyle)
mpl.semilogx([149,149],[ymin,ymax],**lumlinestyle)
mpl.plot([300,300],[ymin,ymax],**lumlinestyle)
mpl.plot([3000,3000],[ymin,ymax],**lumlinestyle)
mpl.errorbar(lum,sigf,yerr=sigf_std,**errstyle)

mpl.ylim(ymin,ymax)
mpl.xlim(xmin,xmax)
mpl.xlabel('Integrated luminosity [$\mathsf{fb^{-1}}$]',horizontalalignment='right',x=1.0)
mpl.ylabel('Significance [$\mathsf{\sigma}$]',horizontalalignment='right',y=1.0)

mpl.gca().tick_params(direction='in',axis='both',which='both',bottom=True,top=True,left=True,right=True)

mpl.text(130,7.6,'Current (Run 2) ATLAS dataset',rotation = 90, fontsize=lum_label_size,**hel)
mpl.text(265,6.25,'Projected Run 3 dataset',rotation = 90, fontsize=lum_label_size,**hel)
mpl.text(2650,4.5,'Projected HL-LHC dataset',rotation = 90, fontsize=lum_label_size,**hel)

def func(x,a,b):
    return a*np.sqrt(x)+b

popt, pcov = curve_fit(func,lum,sigf)
print popt

print ((5-popt[1])/popt[0])**2

x = np.linspace(xmin,xmax,10000)
mpl.plot(x,func(x,*popt), 'r-',label="$\sqrt{N}$ fit")

mpl.legend(**legstyle)

mpl.show()
