import numpy as np
import matplotlib.pyplot as mpl

f_l = 0.338
f_p = 0.27
f_m = 0.392

x = np.linspace(-1,1,1000)

A = 100.

p_l = 3*(1-x**2)/4.
p_p = 3*(1-2*A*x+x**2)/8.
p_m = (3.*(1+2*A*x+x**2))/8.

#mpl.rcParams['axes.linewidth'] = 1.5
#mpl.rcParams['xtick.major.size'] = 8
#mpl.rcParams['xtick.minor.size'] = 4
#mpl.rcParams['ytick.major.size'] = 6
#mpl.rcParams['ytick.minor.size'] = 3
#mpl.rcParams['ytick.minor.size'] = 3

mpl.figure()
mpl.plot(x,f_l*p_l,color='blue')
#mpl.plot(x,(f_p*p_p+f_m*p_m)/(f_p+f_m),color='red')
mpl.plot(x,f_p*p_p,color='limegreen')
mpl.plot(x,f_m*p_m,color='red')
mpl.plot(x,f_l*p_l+f_p*p_p+f_m*p_m,color='black')

mpl.xlim(-1,1)
mpl.ylim(0,1)

ax = mpl.gca()
ax.tick_params(direction='in',which='both',bottom=True,top=True,left=True,right=True,)
ax.set_xticks([-1,-0.5,0,0.5,1])
#ax.legend(['$V_LV_L\\to V_LV_L$','$V_TV_T\\to V_TV_T$','Combined'],framealpha=0)
ax.legend(['Longitudinal','Right Polarised', 'Left Polarised','Combined'],framealpha=0)
ax.minorticks_on()
mpl.xlabel('$\\mathsf{\\cos\\,\\theta^*}$',horizontalalignment='right',x=1.0)
mpl.ylabel('Probability',horizontalalignment='right',y=1.0)

mpl.text(-0.8,0.9,'A = ' + str(A))

mpl.show()
