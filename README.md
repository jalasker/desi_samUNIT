To-do:
[] Plot the LF for the chosen snapshot, compared w data
[] Generate cumulative plot w L[OII], SFR and M* (to make LRG type also)
[] Obtain the cuts for the target densities
[] Measure bias from the catalogues (maybe from Pk)

Maybe:
[] HODs of the samples and other properties

In this folder there are codes to generate a [DESI-like catalogues](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/make_galaxy) from the SAGE+UNIT catalogue. Programs and subfolders here:

**lo2f.py** To plot the model [OII] luminosity function at z~1, compared with observations.

**cum_lo2.py, cum_sf.py, cum_m.py** To plot the cumulative function on the indicated property and find the cuts on it that will give the target DESI number density.

## Specifications for DESI samples

Specifications following the mocks tabulated [here](https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2).

###UNITSIM
From DESI (https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2):

*ELGs* mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}(Mpc/h)^-3.

*LRGs* mocks at z=0.9873, density=4.4e-4 (Mpc/h)^-3 and z=0.7400, density=5.5e-4 (Mpc/h)^-3.

In taurus, emission lines for the HighRes UNITSIM in [data6](/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/):

* Flux [erg/s/cm^2]
* Luminosity [erg/s]

## Output

### Format
The output should be a text file with (following [DESI wiki](https://desi.lbl.gov/trac/wiki/CosmoSimsWG/DESI_mocks)):

{x, y, z, z_RSD} in units of Mpc/h

z_RSD is including the position correction due to redshift distortion in the z-axis. 

z_RSD = z + vz*(1+redshift)/H,

where H=100*sqrt(Omega_m*(1+redshift)3 + Omega_Lambda)

### Location at NERSC

Once the catalogues are tested, they can go to NERSC:
'''
/global/project/projectdirs/desi/mocks/UNIT/SAGE_[ELG, LRG, QSO]/z*
'''
