import sys
import h5py
import numpy as np
#from iotools import check_file
import os.path
import read_jc_obs as jc
import matplotlib ; matplotlib.use('Agg')  
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

zz = 0.987
sims = ['all_z0.9873']
#sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h


#unitdir = '/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/'
unitdir = '/global/project/projectdirs/desi/mocks/UNIT/SAM_madrid/'

min20p = 20.*1.2*10.**9 # Msun/h
h0 = 0.6774

#############################
#outdir = '/home2/vgonzalez/out/desi_samUNIT/'
outdir = '/global/cscratch1/sd/jlasker/UNIT_SAM_output/'
plotdir = outdir+'plots/'
obsdir = '/global/project/projectdirs/desi/mocks/'

#############################

if Testing: sims = [sims[0]] 


# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
edges = np.array(np.arange(lmin,lmax,dl))
lhist = edges[1:]-0.5*dl

lf, lf_att = [np.zeros(shape=len(lhist)) for i in range(2)]

# Initialize plot
fig, ax = plt.subplots()
xtit = "${\\rm log}_{10}(L\\rm{[OII]}/h^{-2}erg\, s^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dex}^{-1})$"

xmin = 40.2 ; xmax = 43.7
ymin = -5.9 ; ymax = -1.

ax.set(xlabel=xtit, ylabel=ytit,xlim=(xmin,xmax),ylim=(ymin,ymax))

nsims = len(sims)
ocol = 'grey'                                                                                     
cm = plt.get_cmap('tab10') # Colour map to draw colours from

# Loop over the simulations
cols =[]
for ii,sim in enumerate(sims):
    col = cm(1.*ii/(nsims+1)) ; cols.append(col)

    volume = lboxes[ii]**3
    
    # File to read
    #ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    ff = unitdir+sim+'/'+'UNITSIM1_model_z'+str(zz)+'_ELGs.h5'

    #if (not check_file(ff)):  continue
    if (not os.path.exists(ff)):
        print(ff)
        continue

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
        ll = np.log10(lum[ind]) + 2*np.log10(h0) # erg/s/h**2
        H, bins_edges = np.histogram(ll,bins=edges)
        lf = lf + H

    # Attenuated LF
    ind = np.where(lum_att>0.)
    if (np.shape(ind)[1] > 0.):
        ll = np.log10(lum_att[ind]) + 2*np.log10(h0) # erg/s/h**2
        H, bins_edges = np.histogram(ll,bins=edges)
        lf_att = lf_att + H

    # Normalize the LF
    lf = lf/dl/volume
    lf_att = lf_att/dl/volume


    # Plot all observations
    ox, oy, el, eh = jc.read_jc_lf(obsdir+'OIIlf_may16_comparat/',zz,
                                   infile='O2_3728-data-summary-Planck15.txt')
    ind = np.where(oy>-5) 
    oxr = ox[ind] ; oyr = oy[ind]
    arrinds = oxr.argsort()
    oxr = oxr[arrinds]
    oyr = oyr[arrinds]

    if(isinstance(ox, (np.ndarray))):
        ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',ecolor=ocol,color=ocol,mec=ocol)
            
    # Plot intrinsic model LF
    ind = np.where(lf > 0)
    x = lhist[ind]
    y = np.log10(lf[ind])
    ax.plot(x,y,color=col,linestyle='-',label=sim)

    # Plot attenuated model LF
    ind = np.where(lf_att > 0)
    x = lhist[ind]
    y = np.log10(lf_att[ind])
    ax.plot(x,y,color=col,linestyle='--')


# Legend    
leg = ax.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)                         
leg.draw_frame(False)                                                                             
for ii,text in enumerate(leg.get_texts()):                                                        
    text.set_color(cols[ii])                                                                     
for item in leg.legendHandles:                                                                
    item.set_visible(False)                                                                   
                                                                                                  
# Save plot
plt.tight_layout()                                                                                
plotf = plotdir+'lf02_UNIT_'+str(zz)+'.pdf'  
print('Plot: {}'.format(plotf))
fig.savefig(plotf)

