import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Testing = True

zz = 0.987

gamma = 0.545

sims = ['UNITSIM1']#,'UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
lboxes = [1000.]*len(sims) # Mpc/h

############################################   
seldir = '/home2/vgonzalez/out/desi_samUNIT/'
############################################

if Testing: sims = [sims[0]]


for sim in sims:
    inpath = seldir+sim

    # r-space
    files = glob.glob(inpath+'/ascii_files/Pk/Pk*dat')
    print(files)


    # Plot
    pathplot = inpath+'/plots/'
    print(pathplot)
