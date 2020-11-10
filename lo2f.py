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

# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
edges = np.array(np.arange(lmin,lmax,dl))
lhist = edges[1:]-0.5*dl

lf, lf_att = [np.zeros(shape=len(lhist)) for i in range(2)]

# Loop over the simulations
for sim in sims:
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    if (not check_file(ff)):  continue

    f = h5py.File(ff,'r') 
    
    ifirst=0 ; ilast=f['Mstar'].shape[0]
    if Testing:
        ilast = 50000

    # Luminosity of the [OII] doublet logL[erg/s]
    lum1 = 10**f['logLOII_3727'][ifirst:ilast] + 10**f['logLOII_3729'][ifirst:ilast]
    lum_att1 = 10**f['logLOII_3727_att'][ifirst:ilast] + 10**f['logLOII_3729_att'][ifirst:ilast]
    mhalo1 = f['Mhalo'][ifirst:ilast]
    f.close()

    # Remove haloes with too small mass 
    ind = np.where(mhalo1 > min20p)
    lum = lum1[ind]
    lum_att = lum_att1[ind]
    lum1 = [] ; lum_att1 = []

    # Intrinsic LF
    ind = np.where(lum>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(lum[ind])
        H, bins_edges = np.histogram(ll,bins=edges)
        lf = lf + H

    # Attenuated LF
    ind = np.where(lum_att>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(lum_att[ind])
        H, bins_edges = np.histogram(ll,bins=edges)
        lf_att = lf_att + H

    print(lf_att)
