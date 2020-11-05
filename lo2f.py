import sys
import h5py
import numpy as np
from iotools import check_file

Testing = True

zz = 0.987

#sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
sims = ['UNITSIM1']

unitdir = '/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/'

min20p = 20.*1.2*10.**9 # Mpc/h

#############################
outdir = '/home2/vgonzalez/out/desi_samUNIT/'
plotdir = outdir+'plots/'
#############################

for sim in sims:
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    if (not check_file(ff)):  continue
    
    f = h5py.File(ff,'r') 
    
    ifirst=0 ; ilast=f['Mstar'].shape[0]
    if Testing:
        ilast = 50000

    fo21 = f['logFOII_3727'][ifirst:ilast]
    fo2_att1 = f['logFOII_3727_att'][ifirst:ilast]
    mhalo1 = f['Mhalo'][ifirst:ilast]
    
    f.close()

    # Remove haloes with too small mass
    ind = np.where(mhalo1 > min20p)
    fo2 = fo21[ind]
    fo2_att = fo2_att1[ind]
    fo21 = [] ; fo2_att1 = []
    #print(type(mstar))

