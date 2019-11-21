import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
     

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.markeredgewidth'] = 1
mpl.rcParams['font.size'] = 16 # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 24
mpl.rcParams['figure.figsize'] = [6.,5]

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 


x_min = -np.pi/2
x_max = +np.pi/2

file_BANDS = open('bands0.txt','r')
MAT_BANDS = (np.loadtxt(file_BANDS))
file_BANDS.close()


N_BAND = np.size(MAT_BANDS[:,0])
k=np.linspace(x_min, x_max, N_BAND)     
    

fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)
#fig1.suptitle(r"SSH bands", y=1.00)

ax11 = fig1.add_subplot(gs1[0,0])
ax11.set_xlabel('$\mathrm{k}$ ($\mathrm{a_0^{-1}}$)')
ax11.set_ylabel('$\mathrm{Energy}$ $\mathrm{(eV)}$')
ax11.set_xticks([x_min, 0.0, x_max])
ax11.set_xticklabels([r'$-\mathrm{\pi/2}$', '0' , '$+\mathrm{\pi/2}$'])
ax11.plot(k,MAT_BANDS[:,0], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,1], 'k--', linewidth=2.0)
#ax11.plot(k,[0]*N_BAND, 'k--', linewidth=1.0)
ax11.set_xlim(x_min, x_max)
ax11.set_xticks(np.linspace(x_min,x_max,3))

plt.subplots_adjust(top=0.95, bottom=0.15, left=0.15, right=0.95, hspace=0.15, wspace=0.15)

plt.show()
