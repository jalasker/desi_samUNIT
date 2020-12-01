#https://arxiv.org/pdf/1306.6721.pdf
# Check CAMB paper

import sys
import glob
import numpy as np
import camb as cb
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

calculate = True
plot = True

zz = [0.,0.987]

# Cosmo UNITSIMs
h0 = 0.6774
ns = 0.9667
Obh2 = 0.02234
Omh2 = 0.3089*h0**2
sigma8 = 0.8147

############################################   
outdir = '/home2/vgonzalez/out/desi_samUNIT/'
############################################

# Set Cosmology for CAMB (default As=2e-09, used for the normalisation)
params = cb.set_params(ns=ns, H0=h0*100, ombh2=Obh2, omch2=(Omh2-Obh2), WantTransfer=True)

#log-spaced array between 10^-4 to 10^3
logkmin = -4 ; logkmax=3
dk = (logkmax-logkmin)/100000
kout = 10**np.arange(logkmin,logkmax,dk)

# Get the power spectrum at z=0 to obtain sigma8
pars = params.copy()
pars.set_matter_power(redshifts=[0.], kmax=10**logkmax)
pars.NonLinear = cb.model.NonLinear_none #Same sigma8 if NonLinear_both 
results = cb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=10**logkmin, maxkh=10**logkmax, npoints = 500)
sigma8_lin = np.array(results.get_sigma8_0())

# Linear matter power spectrum interpolator
PKlin = cb.get_matter_power_interpolator(params,kmax=10**logkmax,
                                         nonlinear=False,hubble_units=True)

# Interpolator including non-linear correction from halo model
PKNL = cb.get_matter_power_interpolator(params,kmax=10**logkmax,
                                        nonlinear=True,hubble_units=True)

if plot:
    fig, ax = plt.subplots()
    xtit = '$k$ [$h$/Mpc]'
    ytit = 'P($k$) [Mpc/$h$)$^3$]'
    ax.set(xlabel=xtit,ylabel=ytit,xlim=(10**logkmin,10**logkmax),
           xscale='log',yscale='log')
    cm = plt.get_cmap('tab10') # Colour map to draw colours from
    cols = []
    
for ii,z1 in enumerate(zz):
    outfil = outdir+'linPk_z'+str(z1).replace('.','_')+'.dat'

    if calculate:        
        Poutlin = (sigma8**2/sigma8_lin**2)*PKlin.P(z1,kout)
        PoutNL = (sigma8**2/sigma8_lin**2)*PKNL.P(z1,kout)

        # Write output
        with open(outfil,'w') as ff:
            ff.write('# k [h/Mpc], P(k,z)lin, P(k,z) with NL corrections [(Mpc/h)^3] \n')
            
        tofile = np.array([kout,Poutlin,PoutNL])
        with open(outfil,'a') as ff:
            np.savetxt(ff,tofile.T)
        print('Outfile: {}'.format(outfil))

    if plot:
        col = cm(1.*ii/(len(zz)+1))
        cols.append(col) ; cols.append(col)

        k,plin,pnl = np.loadtxt(outfil,unpack=True)
        ax.plot(k,plin,color=col,linestyle='-',label='z='+str(z1))
        ax.plot(k,pnl,color=col,linestyle='--',label='NL corrections')

if plot:
    # Legend
    leg = ax.legend(loc=0,fontsize='small')
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])

    # Save plot
    plotf = outdir+'plots/linPk_z.pdf'
    print('Plot: {}'.format(plotf))
    fig.savefig(plotf)
        
