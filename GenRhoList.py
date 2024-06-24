# This file generates a list of different baryon densities and
# saves them in a format the slurm script expects for running
# many different hydro runs in batches

import scipy as sp
from decimal import *

###RHO IN FM-3###

minRho = .4
maxRho = 1.0
length = 61

step = round((maxRho-minRho)/length,2)
dec = abs(Decimal(str(step)).as_tuple().exponent)

rhoList = sp.linspace(minRho,maxRho,length)
listName = 'RhoList_{0}_{1}_{2}.dat'.format(minRho,maxRho,step)

with open(listName,"w") as datFile:
    for rho in rhoList:
        rho = round(rho,dec)
        datFile.write('{0}\n'.format(rho))
