# conda activate nbodykit-env
import sys
import glob
from nbodykit.source.catalog import CSVCatalog
from nbodykit.lab import * #to_mesh, FFTPower, power
import numpy as np

Testing = False

zz = 0.987

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h
h0 = 0.6774

ngrid = 128
kmin = 0.01
dk = 0.005
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
        # Reading the file in nbodykit format
        names =['x', 'y', 'z', 'z_rsd']
        f = CSVCatalog(ff, names)
        f['Position'] = f['x'][:, None] * [1, 0, 0] + \
                        f['y'][:, None] * [0, 1, 0] + \
                        f['z'][:, None] * [0, 0, 1]
        f.attrs['BoxSize'] = lboxes[iis]

        # File to mesh
        mesh = f.to_mesh(compensated=True, Nmesh=ngrid, BoxSize=lboxes[iis], position='Position')
        # 1D power
        r = FFTPower(mesh, mode='1d', dk=dk, kmin=kmin)
        Pk = r.power

        # Shotnoise
        shotnoise = lboxes[iis]**3/f['x'].size
        print("Shotnoise(nbodykit) / volume/ngal = {} ({} %)".format(
            Pk.attrs['shotnoise']/shotnoise,100.*(1.-Pk.attrs['shotnoise']/shotnoise)))
        
        # Save results
        outm = ff.replace('/ascii_files/','/ascii_files/Pk/Pk_').replace(
            '.dat','_kmin'+str(kmin).replace('.','_')+\
            '_dk'+str(dk).replace('.','_')+'.dat')
        with open(outm,'w') as outf:
            outf.write('# k, Pk-shotnoise \n')

        tofile = np.array([Pk['k'],Pk['power'].real - Pk.attrs['shotnoise']])
        with open(outm,'a') as outf:
            np.savetxt(outf,tofile.T)
        print('Outfile: {}'.format(outm))
