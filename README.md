To-do:
[] Check flux cut+from Karim et al. vs. my paper
[] Plot the LF for the chosen snapshot, compared w data
[] Generate cumulative plot w L[OII] and maybe SFR
[] Obtain the cuts for the target densities


Generate a DESI-like catalogues:
https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/make_galaxy

## Specifications for DESI samples

Specifications following the mocks tabulated here:
https://desi.lbl.gov/trac/wiki/Clustering/MockChallenge/post-recon-BAO/stage2

###UNITSIM
*ELG* mocks at redshift = 0.9873 (snap 97), densities={25e-4,20e-4,5.5e-4}

## Output

### Format
The output should be a text file with (following https://desi.lbl.gov/trac/wiki/CosmoSimsWG/DESI_mocks):
{x, y, z, z_RSD} in units of Mpc/h
z_RSD is including the position correction due to redshift distortion in the z-axis. 
z_RSD = z + vz*(1+redshift)/H,
where H=100*sqrt(Omega_m*(1+redshift)3 + Omega_Lambda)

### Location at NERSC

Once the catalogues are tested, they can go to:
/global/project/projectdirs/desi/mocks/UNIT/SAGE_[ELG, LRG, QSO]/z*
