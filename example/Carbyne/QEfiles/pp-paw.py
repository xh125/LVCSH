#
# Post-processing script QE --> EPW
# 14/07/2015 - Samuel Ponce
#

from builtins import input
import numpy as np
import os

# Enter the number of irr. q-points
prefix = input('Enter the prefix used for PH calculations (e.g. diam)\n')

# Enter the number of irr. q-points
nqpt = input('Enter the number of irreducible q-points\n')

# Enter the number of irr. q-points
outdir = input('Enter the outdir used for PH calculations (e.g. "./out")\n')

# Enter the number of irr. q-points
fildvscf = input('Enter the fildvscf used for PH calculations (e.g. "./out")\n')


try:
  nqpt = int(nqpt)
except ValueError:
  raise Exception('The value you enter is not an integer!')

os.system('mkdir save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp '+outdir+'/_ph0/'+prefix+'.'+fildvscf+ '1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp '+outdir+'/_ph0/'+prefix+'.'+fildvscf+ '_paw1 save/'+prefix+'.dvscf_paw_q'+label)
    os.system('cp -r '+outdir+'/_ph0/'+prefix+'.phsave save/')
  else:
    os.system('cp '+outdir+'/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.'+fildvscf+'1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp '+outdir+'/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.'+fildvscf+'_paw1 save/'+prefix+'.dvscf_paw_q'+label)
    os.system('rm '+outdir+'/_ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
