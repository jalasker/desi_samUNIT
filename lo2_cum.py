
zz = 0.987

sims = ['UNITSIM1','UNITSIM1_InvPhase','UNITSIM2','UNITSIM2_InvPhase']
#sims = ['UNITSIM1']

unitdir = '/data6/users/aknebe/Projects/UNITSIM/ELGs/'

#############################
#line = 'OII3727' ; lline = '[OII]'
#outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'selections/lo2_cum_'
#plotfile = outdir+line+'.pdf'
#############################

for sim in sims:
    ff = unitdir+sim+'/'+sim+'_model_z'+str(zz)+'_ELGs.h5'
    print(ff)

