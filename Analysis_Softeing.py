import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.signal import argrelextrema

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

tc = 0.658 # 1/eV = 0.658 fs


LIST_AA = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.092, 0.094, 0.096, 0.098, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16])#, 0.17, 0.18, 0.19, 0.20])
LIST_EE = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.092, 0.094, 0.096, 0.098, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16])*188./3.84*0.1# 0.17, 0.18, 0.19, 0.20])*188./3.84*0.1
print(LIST_EE)
LIST_AMPS = np.zeros(np.size(LIST_AA))


for i in range(np.size(LIST_AA)):
    print(LIST_AA[i])
    file_ORDER = open(str(LIST_AA[i])+'/ORDER_t.txt','r')
    ORDER_ION = np.loadtxt(file_ORDER)
    file_ORDER.close()

    t_min = 0
    t_max = 3040.0*tc
    NN = np.size(ORDER_ION[:,0])

    xs = np.linspace(0,t_max, NN)  
    data = ORDER_ION[:,0]

    u_init = ORDER_ION[0,0]
    u_min = np.amin(ORDER_ION[:,0])
    
    print(u_min)
    print(u_min/u_init)     
    LIST_AMPS[i] = u_min/u_init
    
fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)

ax00 = fig1.add_subplot(gs1[0,0])
ax00.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax00.set_yticks([-0.8, -0.4, 0.0, 0.4, 0.8])
#ax00.set_xticklabels(["0", "1"])
ax00.set_xlim(LIST_EE[0], LIST_EE[np.size(LIST_AA)-1]*1.04)
#ax00.set_ylim(LIST_AMPS[0],LIST_AMPS[7]*1.05)
ax00.set_xlabel(r'$\mathrm{E_{max}}$ $(\mathrm{MV/cm})$')
ax00.set_ylabel(r'$\mathrm{u_{min}/u_0}$')
ax00.get_yaxis().set_label_coords(-0.11,0.5)

ax00_2 = ax00.twiny()
ax00_2.plot(LIST_AA , LIST_AMPS, 'ko', label=r"$\mathrm{\Delta_K^{Floquet}}$", markersize=5)
ax00_2.plot(LIST_AA , LIST_AMPS, 'k-', label=r"$\mathrm{\Delta_K^{Floquet}}$", linewidth=2.0)
ax00_2.set_xticks([0.06, 0.08, 0.10, 0.12, 0.14, 0.16])
#ax00_2.set_xticklabels(["0", "0.02", "0.04", "0.06", "0.08", "0.1"])
ax00_2.set_xlabel(r'$\mathrm{A_{max}}$ $\mathrm{(a_0^{-1})}$')
#ax00_2.set_xlim(LIST_AA[0], LIST_AA [7]*1.04)r'$\mathrm{E_{max}}$ $(\mathrm{MV/cm})$'

plt.subplots_adjust(top=0.85, bottom=0.15, left=0.15, right=0.95, hspace=0.15, wspace=0.15)


PERIOD = np.array([3495-1646, 3521-1663, 3537-1660, 3552-1660, 3856-1856, 4037-1965, 4320-2138, 4937-2470, 2*(4775-2433), 5536-1903, 4572-2377, 3836-1973, 2591-953, 2409-884, 2289-834, 2210-795])#, 2104-760, 2033-732, 1995-711, 1973-694])
OMEGA = np.array([0.00224, 0.00223, 0.00220, 0.00219, 0.00207, 0.00200, 0.00190, 0.00168, 0.00088, 0.00114, 0.00188, 0.00222, 0.00252, 0.00271, 0.00284, 0.00292])*1000#,, 0.00308, 0.00318, 0.00322, 0.00323])*1000
print(PERIOD)

fig1 = plt.figure(2)
gs1 = gridspec.GridSpec(1, 1)

ax00 = fig1.add_subplot(gs1[0,0])
ax00.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
#ax00.set_xticklabels(["0", "1"])
ax00.set_xlim(LIST_EE[0], LIST_EE[np.size(LIST_AA)-1]*1.04)
#ax00.set_ylim(LIST_AMPS[0],LIST_AMPS[7]*1.05)
ax00.set_xlabel(r'$\mathrm{E_{max}}$ $(\mathrm{MV/cm})$')
ax00.set_ylabel(r'$\mathrm{\Omega}$ $(\mathrm{meV})$')
ax00.get_yaxis().set_label_coords(-0.13,0.5)

ax00_2 = ax00.twiny()
ax00_2.plot(LIST_AA , OMEGA, 'ko', label=r"$\mathrm{\Delta_K^{Floquet}}$", markersize=5)
ax00_2.plot(LIST_AA , OMEGA, 'k-', label=r"$\mathrm{\Delta_K^{Floquet}}$", linewidth=2.0)
#ax00_2.plot(LIST_AA , [3.47]*np.size(LIST_AA), 'k--', label=r"$\mathrm{\Delta_K^{Floquet}}$", linewidth=1.5)

ax00_2.set_xticks([0.06, 0.08, 0.10, 0.12, 0.14, 0.16])
#ax00_2.set_xticklabels(["0", "0.02", "0.04", "0.06", "0.08", "0.1"])
ax00_2.set_xlabel(r'$\mathrm{A_{max}}$ $\mathrm{(a_0^{-1})}$')
#ax00_2.set_xlim(LIST_AA[0], LIST_AA [7]*1.04)r'$\mathrm{E_{max}}$ $(\mathrm{MV/cm})$'

plt.subplots_adjust(top=0.85, bottom=0.15, left=0.15, right=0.92, hspace=0.15, wspace=0.15)

plt.show()

