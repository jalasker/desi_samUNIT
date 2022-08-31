import os.path, sys
import numpy as np
from scipy.interpolate import interp1d

Testing = True

zz = 0.9873

#sims = ['z0.9873']
sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']

#############################
#inpath = '/global/cscratch1/sd/jlasker/UNIT_SAM_output/'
inpath = '/global/cscratch1/sd/jlasker/UNIT_SAM_output/'
#inpath = '/home2/vgonzalez/out/desi_samUNIT/'
#############################

props = ['mass','sfr','lo2']
labelp = ['log10(M/Msun/h)','log10(SFR/Msun/Gyr)', 'log10(L[OII]/h^-2 erg/s)']

#massnd = [0.75e-4, 4.4e-4, 13.34e-4]
#sfnd = [25e-4,20e-4,5.5e-4, 13.34e-4]

massnd = [13.34e-4, 14.92e-4, 1.0E-2, 10**(-2.5)]
sfnd = [13.34e-4, 14.92e-4, 1.0E-2, 10**(-2.5)]

if Testing: sims=[sims[0]] #; props=[props[0]]

redshift = str(zz).replace('.','_')

for sim in sims:
    for ip,prop in enumerate(props):
        # Read the cumulative abundance for the stllar mass
        infile = inpath+sim+'/cumu_'+prop+'_z'+redshift+'.dat' #; print(infile)
        if (not os.path.isfile(infile)):
            print('STOP: {} not found'.format(infile)) ; sys.exit()
        data = np.loadtxt(infile, unpack=True) 
        vals = data[0,:] #; print(vals)
        lng  = data[1,:] #; print(lng)

        # Write output header
        outfile = inpath+sim+'/ngal_'+prop+'_cuts_z'+redshift+'.dat'
        ff = open(outfile,'w') ; print('Outfile: {}'.format(outfile))
        ff.write('# log(nd/(Mpc/h)^-3), Cut in '+labelp[ip]+' \n' )

        if (prop == 'mass'):
            nds = np.log10(np.array(massnd))
        else:
            nds = np.log10(np.array(sfnd))

        # Find cuts
        for nd in nds:
            if (max(lng) >= nd):
                f = interp1d(lng,vals)
                cut = f(nd)
            else:
                cut = -999.

            ff.write(' '.join(str(jj) for jj in [nd,cut]))
            ff.write(' \n')
        ff.close()
