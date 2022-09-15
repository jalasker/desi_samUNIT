import sys
import glob
import numpy as np
from scipy import interpolate
from stats import chi2
import Cosmology as cosmo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import mpl_style
from desi_samUNIT import supercomputerSupport as sc

plt.style.use(mpl_style.style1)

Testing = False

zz = 0.987

ngrid = 128 # The one used in pkg.py
sims = ['all_z0.9873']

#sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h
h0 = 0.6774
omegam0 = 0.3089
omegab0 = 0.02234/h0/h0

# Limits to obtain the bias
kminb = 0.01
kmaxb = 0.1

############################################   
#seldir = '/home2/vgonzalez/out/desi_samUNIT/'
seldir = sc.scratchdir()+ '/UNIT_SAM_output/'#'/global/cscratch1/sd/jlasker/UNIT_SAM_output/'

############################################

if Testing: sims = [sims[0]]

# Read theoretical linear P(k,zz)
thfil = seldir+'linPk_z'+str(zz).replace('.','_')+'.dat'
k_th,plin_th,pnl_th = np.loadtxt(thfil,unpack=True)
fpth = interpolate.interp1d(k_th,plin_th)

# Array for trying out bias
abias = np.linspace(0.1,20.,10000)

# Obtain the derivative of the linear growth rate at the target redshift
gamma = 0.545
cosmo.set_cosmology(omega0=omegam0,omegab=omegab0,h0=h0)
omegamz = omegam0*(1+zz)**3/(cosmo.E(zz)**2)
fg = np.power(omegamz,gamma)

# Read the galaxy P(k,zz)
for iis,sim in enumerate(sims):
    inpath = seldir+sim
    knyquist = np.pi*ngrid/lboxes[iis]
    volumen = lboxes[iis]**3
    
    # Plot for r-space
    figr = plt.figure() 
    gsr = gridspec.GridSpec(3,1) ; gsr.update(wspace=0., hspace=0.)
    cm = plt.get_cmap('tab10') # Colour map to draw colours from  
    cols = ['k','k']
    
    axr = plt.subplot(gsr[2,0])  # Ratio plot
    axr.set_xlabel("$k$ [$h$/Mpc]")
    axr.set_ylabel("$\\sqrt{P^s_{\\rm i}/P_{\\rm DM}}$")
    axr.set_autoscale_on(False) ;  axr.minorticks_on()
    axr.set_xlim(0.01,knyquist) ; axr.set_ylim(0.8,2)
    axr.set_xscale('log')
    axr.plot(k_th,plin_th/plin_th,color=cols[0])
    axr.plot(k_th,pnl_th/plin_th,color=cols[1],linestyle='--')
    
    axp = plt.subplot(gsr[0:2,0],sharex=axr) # Pk plot
    plt.setp(axp.get_xticklabels(), visible=False)
    axp.set_autoscale_on(False) ;  axp.minorticks_on()
    axp.set_yscale('log') ; axp.set_ylim(11,100000.)
    axp.set_ylabel('P$^s$($k$) [Mpc/$h$)$^3$]')
    axp.axvline(x=knyquist,color=cols[1],linestyle=':')
    axp.plot(k_th,plin_th,color=cols[0],label='DM, z='+str(zz))
    axp.plot(k_th,pnl_th,color=cols[1],linestyle='--',label='DM, NL')
    
    # z-space
    files = sorted(glob.glob(inpath+'/ascii_files/Pkz/Pkz*dat'))
    if Testing: files = [files[0]]

    bfile = inpath+'/ascii_files/Pkz/bias_kmin'+str(kminb).replace('.','_')+\
            '_kmax'+str(kmaxb).replace('.','_')+\
            '_z'+str(zz).replace('.','_')+'.dat'
    with open(bfile,'w') as outf:
        outf.write('# bias selection \n')

    for ii,ff in enumerate(files):
        col = cm(1.*ii/len(files))
        cols.append(col)

        # Get number density
        root = ff.split('/Pkz/Pkz_')[1].split('_kmin')[0]
        galff = ff.split('Pkz/')[0] + root +'.dat'
        xgal = np.loadtxt(galff,usecols=(0),unpack=True)
        nd = len(xgal)/volumen ; xgal = []
        
        # Read the power spectrum and calculate the corresponding error
        kg,pkg = np.loadtxt(ff,unpack=True)
        dk = kg[1] - kg[0] # Sim. bin
        errorPk = np.sqrt(((2*np.pi)**2/(kg**2*dk*volumen))*(pkg + 1/nd)**2)
        pth = fpth(kg)

        # Bias calculation z-space
        chis = np.zeros((len(abias))) ; chis.fill(999.) ; bias = -999.
        ind = np.where((kg>=kminb) & (kg<kmaxb))
        if (np.shape(ind)[1]>1):
            for ib,bb in enumerate(abias):
                obs = pkg[ind]
                model = pth[ind]*(bb*bb + fg*bb*2/3 + fg*fg/5)
                error = errorPk[ind]
                chis[ib] = chi2(obs,model,error)
                #chis[ib] = (obs-model)**2/error

            ib = np.where(chis == np.nanmin(chis))
            bias = abias[ib][0]
            with open(bfile,'a') as outf:
                outf.write(str(bias)+'  '+ff.split('/Pkz/Pkz_')[1]+' \n')
        print('nd={:.3f} for file {}, bias={}'.format(np.log10(nd),root,bias))
            
        # Plot bias z-space
        axr.plot(kg,np.sqrt(pkg/pth),color=col)
        axr.axhline(np.sqrt(bias**2 + fg*bias*2/3 + fg*fg/5),color=col,linestyle=':')

        # Plot bias P(k)
        label = ff.split('Pkz_')[1].split('_z')[0]
        axp.fill_between(kg,pkg-errorPk,pkg+errorPk,color=col,alpha=0.4)
        axp.plot(kg,pkg,color=col,label=label)
        axp.axhline(1/nd,color=col,linestyle=':')
        
    # Legend
    leg = axp.legend(loc=0,fontsize='small')
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])

    # Save figure
    plotfile = inpath+'/plots/Pkz_zspace_z'+str(zz).replace('.','_')+'.pdf'
    figr.savefig(plotfile,constrained_layout=True)
    print('Output: ',plotfile)
    print('   r-space bias in ',bfile)
