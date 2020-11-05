import sys
import h5py
from iotools import check_file

Testing = True

zz = 0.987

#sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
sims = ['UNITSIM1']

unitdir = '/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/'

#############################
#line = 'OII3727' ; lline = '[OII]'
#outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'selections/lo2_cum_'
#plotfile = outdir+line+'.pdf'
#############################

for sim in sims:
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    if (not check_file(ff)):  continue
    
    f = h5py.File(ff,'r')
    
    ifirst=0 ; ilast=f['Mstar'].shape[0]
    if Testing:
        ilast = 500

    fo2 = f['logFOII_3727'][ifirst:ilast]
    fo2_att = f['logFOII_3727_att'][ifirst:ilast]
    mstar = f['Mstar'][ifirst:ilast]
    
    f.close()

