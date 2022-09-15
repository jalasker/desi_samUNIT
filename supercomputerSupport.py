import sys
import os


def whichComputer():
    try:
        host = os.environ['NERSC_HOST']
    except:
        host = 'taurus'
    return host

def unitdir():
    host = whichComputer()
    if ('cori' in host.lower() ) or ('perlmutter' in host.lower()):
        #ud = '/global/project/projectdirs/desi/mocks/UNIT/SAM_madrid/ELGs/'
        ud = '/global/project/projectdirs/desi/mocks/UNIT/SAM_madrid/'
    elif host.lower() == 'taurus':
        #ud = '/data6/users/aknebe/Projects/UNITSIM/ELGs_DESI/'
        ud = '/data6/users/aknebe/Projects/UNITSIM/'
    return ud

def scratchdir():
    host = whichComputer()
    if 'cori' in host.lower():
        sd = os.getenv('CSCRATCH')
    elif 'perlmutter' in host.lower():
        sd = os.getenv('PSCRATCH')
    elif host.lower() == 'taurus':
        sd = '/I/made/up/a/scratch/dir/since/I/dont/know/about/taurus/'
    return sd