sat_dis/                                                                                            0000775 0001750 0001750 00000000000 14315340437 012402  5                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                sat_dis/nsat.py                                                                                     0000664 0001750 0001750 00000011774 14315340437 013733  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                import glob
import h5py
import numpy as np
import bahamas as b
from iotools import is_sorted,create_dir

Testing = False

z_list = [0.75, 1.]

# ARI
sims = ['AGN_TUNED_nu0_L400N1024_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
envs = ['arilega']*len(sims)
outdir = '/hpcdata4/arivgonz/BAHAMAS/'

# COSMA
#sims = ['L400N1024/WMAP9/Sims/BAHAMAS']
#envs = ['cosmalega']
#outdir = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'

############################################
ftype = '.hdf5'

if (len(sims)==1):
    sim_label = b.get_simlabels(sims)[0]
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir,sim_label=sim_label)
else:
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir)

dirplots = dirplots+'satpdf/' ; create_dir(dirplots)

if Testing: sims=[sims[-1]] ; z_list=[z_list[-1]]
############################################
print('\n')

# File to save percentages
perff = dirplots+'nsat_percent.txt'
with open(perff, 'w') as outf:
    header = "# sim file ntot per_cen per_sat per_nocen ncenmin ncenmax nsatmin nsatmax"
    outf.write(header + "\n")

    # Loop over redshifts
    for ii, zz in enumerate(z_list):
        zmins,zmaxs = b.get_zminmaxs([zz])
    
        for iis,sim in enumerate(sims):
            # Get the closest snapshot to the given redshift
            env = envs[iis]
            snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env,dirz=dirz)
            redshift = str(z_snap).replace('.','_')
    
            # Get the selections files
            files = glob.glob(outdir+sim+'/sel_*_z'+redshift+ftype)
            if (len(files)<1): continue
            files.sort()
            if Testing: files = [files[0]]

            for ff in files:
                f = h5py.File(ff, 'r') #; print(ff)
                hid = f['data/hid'][:]
                sat = f['data/sat'][:]
                if (not is_sorted(hid)):
                    ind = np.argsort(hid)
                    hid = hid[ind]
                    sat = sat[ind]
                f.close()
                    
                i=0 # Initialize arrays with first values
                uhid = np.asarray([hid[i]]) ; lastid = uhid[0]
                if (sat[i]==0):
                    nsat = np.asarray([0])
                    ncen = np.asarray([1])
                else:
                    nsat = np.asarray([1])
                    ncen = np.asarray([0])
                
                iiu = 0; #Index for unique haloes
                for i in range(1, len(hid), 1):
                    if (hid[i] == lastid):
                        if (sat[i]==0):
                            ncen[iiu] += 1 
                        else:
                            nsat[iiu] += 1
                    else:
                        iiu += 1
                        uhid = np.append(uhid,hid[i])
                        if (sat[0]==0):
                            nsat = np.append(nsat,0)
                            ncen = np.append(ncen,1)
                        else:
                            nsat = np.append(nsat,0)
                            ncen = np.append(ncen,1)
                        lastid = uhid[-1]
    
                if (len(np.unique(hid)) != len(uhid)):
                    print('STOP: mismatch of unique haloes.'); exit()
                if (max(ncen)>1):
                    print('STOP: halo with more than 1 central? Check halo masses and radius'); exit()
    
                # Write percentages to file
                header = "# file ntot per_cen per_sat per_nocen ncenmin ncenmax nsatmin nsatmax"
                dataf = sim+" "+ff.split('sel_')[-1].split(ftype)[0]

                ntot = sum(ncen)+sum(nsat)
                dataf += "  " + str(ntot)

                dataf += "  {:.1f}%".format(sum(ncen)*100./ntot)
                dataf += "  {:.1f}%".format(sum(nsat)*100./ntot)
                
                if (min(ncen)<1):  # per_nocen
                    dataf += "  {:.1f}%".format(len(ncen[ncen<1])*100./len(ncen))
                else:
                    dataf += "  0%"

                dataf += "  " + str(min(ncen))
                dataf += "  " + str(max(ncen))
                dataf += "  " + str(min(nsat))
                dataf += "  " + str(max(nsat))

                outf.write(dataf + "\n")
    
                # Write output into a new file
                outfile = ff.replace('sel_','nsat_')
                hf = h5py.File(outfile, 'w')
    
                # Data
                hfdat = hf.create_group('data')
                
                hfdat.create_dataset('hid',data=uhid)
                hfdat['hid'].dims[0].label = 'Halo ID: FOF Group Number subhalo belongs to'
    
                hfdat.create_dataset('nsat',data=nsat)
                hfdat['nsat'].dims[0].label = 'Number of selected satellite galaxies in the halo.'
    
                hfdat.create_dataset('ncen',data=ncen)
                hfdat['ncen'].dims[0].label = 'Number of selected central galaxies in the halo.'
                
                hf.close()
                print('  Output: ',outfile)

    outf.close
print('File with percentages: {}'.format(perff)) 
    sat_dis/pdf.py                                                                                      0000664 0001750 0001750 00000006637 14315340437 013541  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                import glob
import h5py
import numpy as np
import bahamas as b
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

z_list = [0.75, 1.]

# ARI
sims = ['AGN_TUNED_nu0_L400N1024_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
envs = ['arilega']*len(sims)
outdir = '/hpcdata4/arivgonz/BAHAMAS/'

# COSMA
#sims = ['L400N1024/WMAP9/Sims/BAHAMAS']
#envs = ['cosmalega']
#outdir = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'

############################################
ftype = '.hdf5'

if (len(sims)==1):
    sim_label = b.get_simlabels(sims)[0]
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir,sim_label=sim_label)
else:
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir)

if Testing: sims=[sims[-1]] ; z_list=[z_list[-1]]
############################################

# Bins for histogram
nmin = 0 ; nmax = 50.; dn = 1
nbins = np.arange(nmin,nmax,dn)
nhist = nbins + dn*0.5

cm = plt.get_cmap('tab20c') # Colour map to draw colours from  

# Loop over redshift
for ii, zz in enumerate(z_list):
    zmins,zmaxs = b.get_zminmaxs([zz])

    # Loop over the files
    for iis,sim in enumerate(sims):        
        # Get the closest snapshot to the given redshift
        env = envs[iis]
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env,dirz=dirz)
        redshift = str(z_snap).replace('.','_')

        # Figure
        fig = plt.figure()
        ytit = "Normalised counts (Area = 1)"
        xtit = "Number of Satellites"
        xmin = 0. ; xmax = 5.
        ymin = 0. ; ymax = 0.1

        ax = fig.add_subplot(111)
        ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.title.set_text(sim+", z="+str(z_snap))
        
        # Get the files with the number of satellites
        files = glob.glob(outdir+sim+'/nsat_*_z'+redshift+ftype)
        if (len(files)<1): continue
        #if Testing: files = [files[0]]
        files.sort()

        cols = []; colsfr = 0; colm = 4
        for ii,ff in enumerate(files):
            # Label
            leg = ff.split('nsat_')[1].split('_z')[0]
            if ('sfr' in leg):
                col =cm(colsfr) ; cols.append(col)
                colsfr += 1
            elif ('mass' in leg):
                col =cm(colm) ; cols.append(col)
                colm += 1
            
            # Open and read the file
            hf = h5py.File(ff, 'r') #; print(ff)
            hid = hf['data/hid'][:]
            nsat = hf['data/nsat'][:]
            ncen = hf['data/ncen'][:]            

            # Histogram
            pnsat, bin_edges = np.histogram(nsat, bins=np.append(nbins,nmax))
            intpn = np.sum(pnsat)*dn
            if(intpn<1e-7): continue
            pn = pnsat/intpn

            # Step plot
            tmp = np.insert(pn,0,pn[0])
            yy = np.asarray(tmp,dtype=float)
            tmp = np.insert(nhist,0,nhist[0]-1)
            xx = np.asarray(tmp,dtype=float)

            ax.step(xx,yy,label=leg,color=col)

        # Legend
        leg = ax.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)
        leg.draw_frame(False)
        for ii,text in enumerate(leg.get_texts()):
            text.set_color(cols[ii])
            
        # Save figure        
        plotf = dirplots+'satpdf/satpdf_'+sim.split('/')[-1]+'_z'+redshift+'.pdf'
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)

                                                                                                 sat_dis/vrsat_dis.py                                                                                0000664 0001750 0001750 00000011360 14315340437 014753  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                import glob
import h5py
import numpy as np
import bahamas as b
from iotools import is_sorted,create_dir
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

z_list = [0.75, 1.]

# ARI
sims = ['AGN_TUNED_nu0_L400N1024_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
envs = ['arilega']*len(sims)
outdir = '/hpcdata4/arivgonz/BAHAMAS/'

# COSMA
#sims = ['L400N1024/WMAP9/Sims/BAHAMAS']
#envs = ['cosmalega']
#outdir = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'

############################################
ftype = '.hdf5'

if (len(sims)==1):
    sim_label = b.get_simlabels(sims)[0]
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir,sim_label=sim_label)
else:
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir)

dirplots = dirplots+'satpdf/' ; create_dir(dirplots)

if Testing: sims=[sims[-1]] ; z_list=[z_list[0]]
############################################
print('\n')

cm = plt.get_cmap('tab20c')  # Colour map to draw colours from

# Loop over redshifts
for ii, zz in enumerate(z_list):
    zmins,zmaxs = b.get_zminmaxs([zz])

    for iis,sim in enumerate(sims):
        # Bins in distance (Mpc/h)
        vmin = -3000.; vmax = -vmin; dv = 100.
        if ('HIRES' in sim): dv = 200.
        vbins = np.arange(vmin,vmax,dv)
        vhist = vbins +dv*0.5

        # Get the closest snapshot to the given redshift
        env = envs[iis]
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env,dirz=dirz)
        redshift = str(z_snap).replace('.','_')

        # Get the selections files
        files = glob.glob(outdir+sim+'/sel_*_z'+redshift+ftype)
        if (len(files)<1): continue
        files.sort()
        if Testing: files = [files[0]]
        
        # Prep. figure
        fig = plt.figure()
        ytit = r"P$_{\rm v}$"
        xtit = r"$v_{\rm r}$(km/s)"
        xmin = -1500. ; xmax = 1499.
        ymin = 1e-12 ; ymax = 2.5e-11
        if ('HIRES' in sim): ymin = 1e-11 ; ymax = 2e-9
        
        ax = fig.add_subplot(111)
        ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        #ax.set_yscale('log')
        ax.title.set_text(sim+", z="+str(z_snap))
        cols = []; colsfr = 0; colm = 4

        f2plot = 0
        for ff in files:
            leg = ff.split('sel_')[-1].split(ftype)[0]
            if ('nd5_0' in leg): continue
            
            # Read file
            f = h5py.File(ff, 'r') #; print(ff)
            header = f['header'] 
            volume = (f['header'].attrs['boxsize'])**3
            
            dpos = f['data/dPos'][:]
            dvel = f['data/dVel'][:]
            sat = f['data/sat'][:]
            f.close()

            # Get only satellite galaxies
            ind1 = np.where(sat > 0)
            if (np.shape(ind1)[1] < 2): continue
            ind = np.squeeze(ind1)
            ntot = ind.shape[0]
            
            if (ntot>1):
                dx = dpos[ind,0]
                dy = dpos[ind,1]
                dz = dpos[ind,2]

                dvx = dvel[ind,0]
                dvy = dvel[ind,1]
                dvz = dvel[ind,2]

                # Distance to the central galaxy
                rsat = np.sqrt(dx*dx + dy*dy + dz*dz)

                # Define the radial velocity and check limiting values 
                vrsat = (dx*dvx + dy*dvy + dz*dvz)/rsat
                if ((max(vrsat) > vmax) or (min(vrsat < vmin))):
                    print('{}, {}: vrmin={:.2f}, vrmax={:.2f}'.format(sim,leg,min(vrsat),max(vrsat)))

                # Histogram
                H, bin_edges = np.histogram(vrsat, bins=np.append(vbins,vbins[-1]+dv))
                nsat = H/volume

                # Label
                if ('sfr' in leg):
                    col = cm(colsfr); cols.append(col); colsfr +=1
                elif('mass' in leg):
                    col = cm(colm); cols.append(col); colm +=1
                else: continue

                # Area to normalise the plot                
                area = np.sum(nsat)*dv
                
                # Plot
                ind = np.where(nsat>0)
                if(np.shape(ind)[1]>1):
                    f2plot += 1
                    y = nsat[ind]/volume/area
                    print("Pv_min={}, Pv_max={}".format(min(y),max(y)))
                    ax.plot(vhist[ind],y,label=leg,color=col)

        # Legend
        if f2plot < 1: continue
        leg = ax.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)
        leg.draw_frame(False)
        for ii, text in enumerate(leg.get_texts()):
            text.set_color(cols[ii])

        # Save figure
        plotf = dirplots+'vrsat_dis_'+sim.split('/')[-1]+'_z'+redshift+'.pdf'
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)
                                                                                                                                                                                                                                                                                sat_dis/README.md                                                                                   0000664 0001750 0001750 00000001022 14315340437 013654  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                Exploration of the distribution of satellites in Bahamas, selected using cuts in different properties.

**nsat.py** Generate a file with the number of sat and centrals for each selection.

**pdf.py** PDF of satellite galaxies using the information generated with nsat.py.

**rsat_dis.py** Radial distribution of satellite galaxies as a function of log r(kpc)

**rnorm_dis.py** Radial distribution of satellite galaxies as a function of both log r/rvir with rvir=R200Crit.

**vsat_dis.py** Satellite's radial velocity distribution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              sat_dis/rsat_dis.py                                                                                 0000664 0001750 0001750 00000010233 14315340437 014563  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                import glob
import h5py
import numpy as np
import bahamas as b
from iotools import is_sorted,create_dir
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

z_list = [0.75, 1.]

# ARI
sims = ['AGN_TUNED_nu0_L400N1024_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
envs = ['arilega']*len(sims)
outdir = '/hpcdata4/arivgonz/BAHAMAS/'

# COSMA
#sims = ['L400N1024/WMAP9/Sims/BAHAMAS']
#envs = ['cosmalega']
#outdir = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'

############################################
ftype = '.hdf5'

if (len(sims)==1):
    sim_label = b.get_simlabels(sims)[0]
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir,sim_label=sim_label)
else:
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir)

dirplots = dirplots+'satpdf/' ; create_dir(dirplots)

if Testing: sims=[sims[-1]] ; z_list=[z_list[0]]
############################################
print('\n')

cm = plt.get_cmap('tab20c')  # Colour map to draw colours from

# Loop over redshifts
for ii, zz in enumerate(z_list):
    zmins,zmaxs = b.get_zminmaxs([zz])

    for iis,sim in enumerate(sims):
        # Bins in distance (Mpc/h)
        dr = 0.1
        if ('HIRES' in sim): dr = 0.2 
        rmin = 0. ; rmax = 5.
        rbins = np.arange(rmin,rmax,dr)
        rhist = rbins +dr*0.5

        # Get the closest snapshot to the given redshift
        env = envs[iis]
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env,dirz=dirz)
        redshift = str(z_snap).replace('.','_')

        # Get the selections files
        files = glob.glob(outdir+sim+'/sel_*_z'+redshift+ftype)
        if (len(files)<1): continue
        files.sort()
        if Testing: files = [files[0]]
        
        # Prep. figure
        fig = plt.figure()
        ytit = r"log$_{10}(N_{\rm sat,dex}/N_{\rm sat,tot})$"
        xtit = r"$r$($h^{-1}$Mpc)"
        xmin = 0; xmax = 5
        ymin = -5. ; ymax = 0.

        ax = fig.add_subplot(111)
        ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.title.set_text(sim+", z="+str(z_snap))
        cols = []; colsfr = 0; colm = 4

        f2plot = 0
        for ff in files:
            leg = ff.split('sel_')[-1].split(ftype)[0]
            if ('nd5_0' in leg): continue
            
            # Read file
            f = h5py.File(ff, 'r') #; print(ff)
            dpos = f['data/dPos'][:]
            sat = f['data/sat'][:]
            f.close()

            # Get only satellite galaxies
            ind1 = np.where(sat > 0)
            if (np.shape(ind1)[1] < 2): continue
            ind = np.squeeze(ind1)
            ntot = ind.shape[0]
            
            if (ntot>1):
                dx = dpos[ind,0]
                dy = dpos[ind,1]
                dz = dpos[ind,2]

                # Get r in kpc/h for satellite galaxies
                rsat = np.sqrt(dx*dx + dy*dy + dz*dz)

                # Show if impose limits are not adequate
                if ((max(rsat) > rmax) or (min(rsat < rmin))):
                    print('{}, {}: rmin={:.2f}, rmax={:.2f}'.format(sim,leg,min(rsat),max(rsat)))

                # Histogram
                H, bin_edges = np.histogram(rsat, bins=np.append(rbins,rbins[-1]+dr))
                nsat = H/ntot

                # Label
                if ('sfr' in leg):
                    col = cm(colsfr); cols.append(col); colsfr +=1
                elif('mass' in leg):
                    col = cm(colm); cols.append(col); colm +=1
                else: continue
                
                # Plot
                ind = np.where(nsat>0)
                if(np.shape(ind)[1]>1):
                    f2plot += 1
                    ax.plot(rhist[ind],np.log10(nsat[ind]),label=leg,color=col)

        # Legend
        if f2plot < 1: continue
        leg = ax.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)
        leg.draw_frame(False)
        for ii, text in enumerate(leg.get_texts()):
            text.set_color(cols[ii])

        # Save figure
        plotf = dirplots+'rsat_dis_'+sim.split('/')[-1]+'_z'+redshift+'.pdf'
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)
                                                                                                                                                                                                                                                                                                                                                                     sat_dis/rnorm_dis.py                                                                                0000664 0001750 0001750 00000010340 14315340437 014746  0                                                                                                    ustar   violeta                         violeta                                                                                                                                                                                                                import glob
import h5py
import numpy as np
import bahamas as b
from iotools import is_sorted,create_dir
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

z_list = [0.75, 1.]

# ARI
sims = ['AGN_TUNED_nu0_L400N1024_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
envs = ['arilega']*len(sims)
outdir = '/hpcdata4/arivgonz/BAHAMAS/'

# COSMA
#sims = ['L400N1024/WMAP9/Sims/BAHAMAS']
#envs = ['cosmalega']
#outdir = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'

############################################
ftype = '.hdf5'

if (len(sims)==1):
    sim_label = b.get_simlabels(sims)[0]
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir,sim_label=sim_label)
else:
    outdir, dirz, dirplots = b.get_outdirs(envs[0],outdir=outdir)

dirplots = dirplots+'satpdf/' ; create_dir(dirplots)

if Testing: sims=[sims[-1]] ; z_list=[z_list[0]]
############################################
print('\n')

cm = plt.get_cmap('tab20c')  # Colour map to draw colours from

# Loop over redshifts
for ii, zz in enumerate(z_list):
    zmins,zmaxs = b.get_zminmaxs([zz])

    for iis,sim in enumerate(sims):
        # Bins in distance (Mpc/h)
        dr = 0.1
        if ('HIRES' in sim): dr = 0.2
        rmin = 0. ; rmax = 6.
        rbins = np.arange(rmin,rmax,dr)
        rhist = rbins +dr*0.5
        
        # Get the closest snapshot to the given redshift
        env = envs[iis]
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env,dirz=dirz)
        redshift = str(z_snap).replace('.','_')

        # Get the selections files
        files = glob.glob(outdir+sim+'/sel_*_z'+redshift+ftype)
        if (len(files)<1): continue
        files.sort()
        if Testing: files = [files[0]]

        # Prep. figure
        fig = plt.figure()
        ytit = r"log$_{10}(N_{\rm sat,dex}/N_{\rm sat,tot})$"
        xtit = r"$R/R_{200c}$"
        xmin = 0; xmax = 5
        ymin = -4. ; ymax = -0.5

        ax = fig.add_subplot(111)
        ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.title.set_text(sim+", z="+str(z_snap))
        cols = []; colsfr = 0; colm = 4

        f2plot = 0
        for ff in files:
            leg = ff.split('sel_')[-1].split(ftype)[0]
            if ('nd5_0' in leg): continue
            
            # Read file
            f = h5py.File(ff, 'r') #; print(ff)
            dpos = f['data/dPos'][:]
            sat = f['data/sat'][:]
            r200c = f['data/R200c'][:]
            f.close()

            # Get only satellite galaxies
            ind1 = np.where(sat > 0)
            if (np.shape(ind1)[1] < 2): continue
            ind = np.squeeze(ind1)
            ntot = ind.shape[0]
            
            if (ntot>1):
                dx = dpos[ind,0]
                dy = dpos[ind,1]
                dz = dpos[ind,2]

                # Get r in kpc/h for satellite galaxies
                rnsat = np.sqrt(dx*dx + dy*dy + dz*dz)/r200c[ind]

                # Show if impose limits are not adequate
                if ((max(rnsat) > rmax) or (min(rnsat < rmin))):
                    print('{}, {}: rmin={:.2f}, rmax={:.2f}'.format(sim,leg,min(rnsat),max(rnsat)))

                # Histogram
                H, bin_edges = np.histogram(rnsat, bins=np.append(rbins,rbins[-1]+dr))
                nsat = H/ntot
                
                # Label
                if ('sfr' in leg):
                    col = cm(colsfr); cols.append(col); colsfr +=1
                elif('mass' in leg):
                    col = cm(colm); cols.append(col); colm +=1
                else: continue
                
                # Plot
                ind = np.where(nsat>0)
                if(np.shape(ind)[1]>1):
                    f2plot += 1
                    ax.plot(rhist[ind],np.log10(nsat[ind]),label=leg,color=col)

        # Legend
        if f2plot < 1: continue
        leg = ax.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)
        leg.draw_frame(False)
        for ii, text in enumerate(leg.get_texts()):
            text.set_color(cols[ii])

        # Save figure
        plotf = dirplots+'rnormdis_'+sim.split('/')[-1]+'_z'+redshift+'.pdf'
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                