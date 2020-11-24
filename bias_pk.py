import sys
import glob
from nbodykit.lab import FFTPower
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

Testing = True

zz = 0.987

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h
h0 = 0.6774

ngrid = 128

############################################   
seldir = '/home2/vgonzalez/out/desi_samUNIT/'
############################################

if Testing: sims = [sims[0]]

for iis,sim in enumerate(sims):
    knyquist = np.pi*ngrid/lboxes[iis]
    print("Nyquist frequency: {}".format(knyquist))
    
    # Path to files
    inpath = seldir+sim+'/ascii_files/'
    files = glob.glob(inpath+'*dat')
    if Testing: files = [files[0]]

    for ff in files:
        xgal,ygal,zgal,zgalz = np.loadtxt(ff,unpack=True)

