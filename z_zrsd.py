import sys,os
import numpy as np
import glob 
import matplotlib ; matplotlib.use('Agg')                                                             
from matplotlib import pyplot as plt                                                                  
import mpl_style                                                                                      
plt.style.use(mpl_style.style1)

Testing = False

zz = 0.987

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h

#############################
outdir = '/home2/vgonzalez/out/desi_samUNIT/'
plotdir = outdir+'plots/zzrsd/'
#############################

if Testing: sims = [sims[0]]

# Loop over the simulations
for ii,sim in enumerate(sims):
    files = glob.glob(outdir+sim+'/ascii_files/*.dat')

    # Loop over files
    for ff in files:
        z,zrsd = np.loadtxt(ff,usecols=(2,3),unpack=True)

        # Prepare plot
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('z')
        ax1.set_ylabel('z_rsd')
        ax1.plot(z, zrsd, 'bo')
        
        # Save plot
        plt.tight_layout()                                                                                
        plotf = plotdir+'zzrsd_'+ff.split('/')[-1].split('.dat')[0]+'.pdf'
        print(plotf)
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)
        
