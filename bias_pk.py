import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

Testing = True

zz = 0.987

gamma = 0.545

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h

############################################   
seldir = '/home2/vgonzalez/out/desi_samUNIT/'
############################################

if Testing: sims = [sims[0]]

logkmin = -4 ; logkmax = 3

# Read theoretical linear P(k,zz)
thfil = seldir+'linPk_z'+str(zz).replace('.','_')+'.dat'
k_th,plin_th,pnl_th = np.loadtxt(thfil,unpack=True)

# Read the galaxy P(k,zz)
for sim in sims:
    inpath = seldir+sim
    pathplot = inpath+'/plots/'
    
    # Plot for r-space
    figr = plt.figure(constrained_layout=True)
    gsr = gridspec.GridSpec(3,1) ; gsr.update(wspace=0., hspace=0.)
    axr = plt.subplot(gsr[2,0])  # Ratio plot
    axr.set_xlabel("$k$ [$h$/Mpc]")
    axr.set_ylabel("$\sqrt{P_{\rm g}/P_{\rm DM}}$")
    axr.set_autoscale_on(False) ;  axr.minorticks_on()
    xmin = 10**logkmin ; xmax = 10**logkmax
    axr.set_xlim(xmin,xmax) #; axr.set_ylim(0.8,1.2)
    axr.set_xscale('log')
        
    axp = plt.subplot(gsr[0:2,0],sharex=axr) # Pk plot
    plt.setp(axp.get_xticklabels(), visible=False)
    axp.set_autoscale_on(False) ;  axp.minorticks_on()
    #axp.set_ylim(0.7,11.)
    axp.set_yscale('log')
    axp.set_ylabel('P($k$) [Mpc/$h$)$^3$]')
    
    # r-space
    files = glob.glob(inpath+'/ascii_files/Pk/Pk*dat')

    for ff in files:
        kg,pkg = np.loadtxt(ff,unpack=True)
        axr #here make plot

