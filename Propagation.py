import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
     
tc = 0.658 # 1/eV = 0.658 fs

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.markeredgewidth'] = 1
mpl.rcParams['font.size'] = 20# <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['xtick.major.size'] = 12
mpl.rcParams['ytick.major.size'] = 12
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 24
mpl.rcParams['figure.figsize'] = [10.,10.]

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 


x_min = 0.
x_max = 10000.0*tc

#LIST_AA = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.092, 0.094, 0.096, 0.098, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20])
LIST_AA = np.array([0.08])
for i in range(np.size(LIST_AA)):
    file_ORDER = open(str(LIST_AA[i])+'/ORDER_t.txt','r')
    ORDER = np.loadtxt(file_ORDER)
    file_ORDER.close()
    
    file_DRIVING = open(str(LIST_AA[i])+'/DRIVING_t.txt','r')
    DRIVING = np.loadtxt(file_DRIVING)
    file_DRIVING.close()
    
    file_ENERGY = open(str(LIST_AA[i])+'/ENERGY_t.txt','r')
    ENERGY = np.loadtxt(file_ENERGY)
    file_ENERGY.close()
    
    # =============================================================================
    # file_ENERGY = open('0.09/ENERGY_t_EXP.dat','r')
    # ENERGY_EXP = np.loadtxt(file_ENERGY)
    # file_ENERGY.close()
    # =============================================================================
    
    Nt = np.size(ORDER[:,0])
    t=np.linspace(-1500*tc, 8500*tc, Nt)
    
    xmin = -1200*tc 
    xmax = 3000
    
    fig1 = plt.figure(i)
    gs1 = gridspec.GridSpec(5, 1)
    
    ax00 = fig1.add_subplot(gs1[0,0])
    ax00.set_ylabel('$\mathrm{A_{peierls}}$ ($\mathrm{a_0^{-1}}$)')
    #ax00.plot(t, DRIVING**2/(np.amax(DRIVING**2)), color = RED, linewidth=2.0)
    ax00.plot(t, DRIVING, color = RED, linewidth=2.0)
    ax00.set_xlim(xmin,xmax)
    ax00.set_xticklabels([])
    ax00.get_yaxis().set_label_coords(-0.11,0.5)
    
    ax10 = fig1.add_subplot(gs1[1,0])
    ax10.set_ylabel('$\mathrm{u}$ ($\mathrm{\AA}$)')
    ax10.plot(t,ORDER[:,0], 'k', linewidth=2.0)
    ax10.set_xlim(xmin,xmax)
    ax10.set_xticklabels([])
    ax10.get_yaxis().set_label_coords(-0.11,0.5)
    plt.legend()
    
    ax20 = fig1.add_subplot(gs1[2,0])
    ax20.set_ylabel('$\mathrm{p}$ ($\mathrm{\AA^{-1}}$)')
    ax20.plot(t,ORDER[:,1], 'k', linewidth=2.0)
    ax20.set_xlim(xmin,xmax)
    ax20.set_xticklabels([])
    ax20.get_yaxis().set_label_coords(-0.11,0.5)
    
    ax30 = fig1.add_subplot(gs1[3,0])
    ax30.set_ylabel('$\mathrm{E_{el}}$ ($\mathrm{eV}$)')
    ax30.plot(t,ENERGY[:,0]-ENERGY[0,0], 'k', linewidth=2.0)
    #ax30.plot(t,ENERGY_EXP-ENERGY_EXP[0], 'k--', linewidth=2.0,label="EXP")
    ax30.set_xlim(xmin,xmax)
    ax30.set_xticklabels([])
    ax30.get_yaxis().set_label_coords(-0.11,0.5)
    
    ax40 = fig1.add_subplot(gs1[4,0])
    ax40.set_xlabel('$\mathrm{time}$ $\mathrm{(fs)}$')
    ax40.set_ylabel('$\mathrm{E_{nuc}}$ ($\mathrm{eV}$)')
    ax40.plot(t,ENERGY[:,1]-ENERGY[0,1], 'k', linewidth=2.0)
    ax40.set_xlim(xmin,xmax)
    ax40.get_yaxis().set_label_coords(-0.11,0.5)
    
    
    plt.subplots_adjust(top=0.95, bottom=0.10, left=0.15, right=0.95, hspace=0.15, wspace=0.15)
    plt.show()
    
    print(0.09*190./3.84*0.1)