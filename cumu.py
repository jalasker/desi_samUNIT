import sys,os
import h5py
import numpy as np
from iotools import check_file
from stats import n_gt_x
import matplotlib ; matplotlib.use('Agg')                                                             
from matplotlib import pyplot as plt                                                                  
import mpl_style                                                                                      
plt.style.use(mpl_style.style1)

Testing = False

zz = 0.987

sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h

unitdir = '/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/'

min20p = 20.*1.2*10.**9 # Msun/h
h0 = 0.6774

nd_lrg = 4.4e-4
nd_elg1 = 25e-4
nd_elg2 = 20e-4
nd_elg3 = 5.5e-4

#############################
outdir = '/home2/vgonzalez/out/desi_samUNIT/'
plotdir = outdir+'plots/'
obsdir = '/home2/vgonzalez/lfs/'
#############################

if Testing: sims = [sims[0]]

# Set up ploting grid
fig = plt.figure(figsize=(20.,7))
gs = matplotlib.gridspec.GridSpec(1,3)

axm = plt.subplot(gs[0])

axs = plt.subplot(gs[1], sharey=axm)
plt.setp(axs.get_yticklabels(), visible=False)

axl = plt.subplot(gs[2], sharey=axm)
plt.setp(axl.get_yticklabels(), visible=False)

ytit="${\\rm log}_{10}(n_{\\rm gal}(>X)/Mpc^{-3}h^3)$"

cm = plt.get_cmap('tab10') # Colour map to draw colours from
nsims = len(sims)
ocol = 'grey'                                                                                     

# Initialize the Mass cum arrays and plot
mmin = 8.5 ; mmax = 16. ; dm = 0.05
medges = np.array(np.arange(mmin,mmax,dm))
mhist = medges[1:]-0.5*dm
xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot}h^{-1})$"
axm.set_xlim(9.,11.6) ;  axm.set_ylim(-5.,-0.7)
axm.set_xlabel(xtit)  ;  axm.set_ylabel(ytit)

# Initialize the SFR cum plot
smin = 7. ; smax = 14. ; ds = 0.02
sedges = np.array(np.arange(smin,smax,ds))
shist = sedges[1:]-0.5*ds
xtit = "${\\rm log}_{10}(SFR/{\\rm M}_{\odot}h^{-1}{\\rm Gyr}^{-1})$"
axs.set_xlim(8.,11.1) 
axs.set_xlabel(xtit)

# Initialize the LO2 cum plot
lmin = 37.5 ; lmax = 44. ; dl = 0.1
ledges = np.array(np.arange(lmin,lmax,dl))
lhist = ledges[1:]-0.5*dl
xtit = "${\\rm log}_{10}(L\\rm{[OII]}/h^{-2}{\\rm erg\,s}^{-1})$"
axl.set_xlim(38.5,42.5) 
axl.set_xlabel(xtit)

# Loop over the simulations
cols =[]
for ii,sim in enumerate(sims):
    ntot = 0
    volume = lboxes[ii]**3
    col = cm(ii) ; cols.append(col)

    # Initialize arrays
    mcum = np.full((len(mhist)),0.)
    scum = np.full((len(shist)),0.)
    lcum = np.full((len(lhist)),0.)

    # File to read
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    if (not check_file(ff)):  continue

    f = h5py.File(ff,'r') 
    
    ifirst=0 ; ilast=f['Mstar'].shape[0]
    if Testing:
        ilast = 50000

    # Luminosity of the [OII] doublet logL[erg/s]
    lum_att1 = 10**f['logLOII_3727_att'][ifirst:ilast] + 10**f['logLOII_3729_att'][ifirst:ilast]
    mass1 = f['Mstar'][ifirst:ilast] # Msun/h
    sfr1 = f['SFR'][ifirst:ilast]*10**9 # Msun/h/Gyr
    mhalo1 = f['Mhalo'][ifirst:ilast] # Msun/h
    f.close()

    ntot += len(mhalo1)
    
    # Remove haloes with too small mass 
    ind = np.where(mhalo1 > min20p)
    mass = mass1[ind]
    sfr = sfr1[ind]
    lum_att = lum_att1[ind]
    mass1 = [] ; sfr1 = [] ; lum_att1 = []

    # Mass
    ind = np.where(mass>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(mass[ind])
        H = n_gt_x(medges[:-1],ll)
        mcum += H

    # SFR
    ind = np.where(sfr>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(sfr[ind])
        H = n_gt_x(sedges[:-1],ll)
        scum += H

    # L[OII]
    ind = np.where(lum_att>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(lum_att[ind]) + 2*np.log10(h0) # erg/s/h**2
        H = n_gt_x(ledges[:-1],ll)
        lcum += H

    # Normalize the cumulative function
    mcum = mcum/volume
    scum = scum/volume
    lcum = lcum/volume
    print('{}: log10(num.den./(Mpc/h)^-3) = {}'.format(sim,np.log10(ntot/volume)))

    # Plot cumulative mass function
    ind = np.where(mcum > 0)
    x = mhist[ind]
    y = np.log10(mcum[ind])
    axm.plot(x,y,color=col,label=sim)

    axm.axhline(y=np.log10(nd_lrg),linestyle='--',color=ocol)
    
    # Plot cumulative SFR function
    ind = np.where(scum > 0)
    x = shist[ind]
    y = np.log10(scum[ind])
    axs.plot(x,y,color=col)

    axs.axhline(y=np.log10(nd_elg1),linestyle='--',color=ocol)
    axs.axhline(y=np.log10(nd_elg2),linestyle='--',color=ocol)
    axs.axhline(y=np.log10(nd_elg3),linestyle='--',color=ocol)

    # Plot cumulative L[OII] function
    ind = np.where(lcum > 0)
    x = lhist[ind]
    y = np.log10(lcum[ind])
    axl.plot(x,y,color=col)

    axl.axhline(y=np.log10(nd_elg1),linestyle='--',color=ocol)
    axl.axhline(y=np.log10(nd_elg2),linestyle='--',color=ocol)
    axl.axhline(y=np.log10(nd_elg3),linestyle='--',color=ocol)


    # Write output for each property
    props = ['mass','sfr','lo2']
    labelp = ['log10(M/Msun/h)','log10(SFR/Msun/Gyr)', 'log10(L[OII]/h^-2 erg/s)']
    for j, prop in enumerate(props):
        dirsim = outdir+sim
        if not os.path.exists(dirsim):
            os.makedirs(dirsim)
        outfile = dirsim+'/cumu_'+prop+'_z'+str(zz).replace('.','_')+'.dat'

        with open(outfile,'w') as outf:
            outf.write('# '+labelp[j]+', log10(n>X/(Mpc/h)^-3) [Sim:'+sim+'] \n')

            if (prop == 'mass'):
                ncum = mcum ; pbins = medges[:-1]
            elif (prop == 'sfr'):
                ncum = scum ; pbins = sedges[:-1]
            elif (prop == 'lo2'):
                ncum = lcum ; pbins = ledges[:-1]

            tofile1 = np.copy(np.transpose(ncum))
            tofile1[tofile1<1e-8] = -999.
            
            ind = np.where(tofile1>0.)
            tofile1[ind] = np.log10(tofile1[ind])

            tofile = np.column_stack((pbins,tofile1))
            np.savetxt(outf, tofile, fmt='%.5f')
        print('Output: ',outfile)
        
# Legend    
leg = axm.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)                         
leg.draw_frame(False)                                                                             
for ii,text in enumerate(leg.get_texts()):                                                        
    text.set_color(cols[ii])                                                                     
for item in leg.legendHandles:                                                                
    item.set_visible(False)                                                                   
                                                                                                  
# Save plot
plt.tight_layout()                                                                                
plotf = plotdir+'cumu_UNIT_'+str(zz)+'.pdf'  
print('Plot: {}'.format(plotf))
fig.savefig(plotf)

