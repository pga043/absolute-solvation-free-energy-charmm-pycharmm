import os
import sys
import numpy as np
import pickle

########LOAD pyCHARMM LIBRARIES########
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.scalar as scalar
from pycharmm.lib import charmm as libcharmm

################################
################################
def set_blocks_vacuum(res,solu,lelec2=1.0,lvdw2=1.0):
    """Set-up block definitions and values:
       env <= CHARMM define for the environment selection
       solu <= CHARMM define for the solute segment to be perturbed
       lenv <= lambda scaling for the environment atoms 
       lelec, lvdw <= lambda scaling for the elec and vdW interactions
       lelec2, lvdw2 <= lambda scaling for the intra-solute elec and vdw interactions"""
    pycharmm.lingo.charmm_script('''
define solu select resname {} end
block
   clear
end
block 1
   call 1 select {} end
end
block
   coef 1 1 1.0 elec {:10.5f} vdw {:10.5f}
end '''.format(resname,
               solu,
                lelec2,lvdw2)
)
    return

################################

def set_lambda(nelec=10,nvdw=15,nshift=0):
    l = np.zeros((nvdw+nshift,2), dtype=float)
    l[:,1]=1.0
    for i in range(nelec): l[i,0] = 1-np.sin(i*np.pi/(2*(nelec-1)))
    for i in range(nvdw): l[i+nshift,1] = 1-np.sin((i)*np.pi/(2*(nvdw-1)))
    return l


##########################SETUP GLOBAL PARAMETERS##########################
T = 298         # K
k_B = 0.00198   # kcal/mol/K
ctofnb = 10
cutnb = ctofnb + 2.0
cutim = cutnb
ctonnb = ctofnb - 2.0

## Parveen Gartan, 23 Dec 2024
# Scripting starts here
solres = str(sys.argv[2])
i = int(sys.argv[1])

## make some directories and sub-directories
if not os.path.isdir(f'win{i}'): os.system(f'mkdir win{i}')
if not os.path.isdir(f'win{i}/post'): os.system(f'mkdir win{i}/post')

## read topology and paramter files
read.rtf('prep/toppar/parm14sb_all.rtf')
read.prm('prep/toppar/parm14sb_all.prm', flex=True)
read.rtf(f'prep/{solres}_charmm.rtf', append=True)
read.prm(f'prep/{solres}_charmm.prm', flex=True, append=True)

## read psf and coor
read.psf_card(f'win{i}/minimized.psf')
read.coor_card(f'win{i}/minimized.crd')

resname = psf.get_res()[0].lower() # Assumes first residue is solute residue
segid = psf.get_segid()[0].upper()

## setup box and crystal information
boxsize = 20
boxhalf = 0.0

crystal.define_cubic(boxsize)
pycharmm.lingo.charmm_script('''
open unit 10 read form name prep/cubic.xtl
crystal read card unit 10
close unit 10''')

image.setup_segment(boxhalf,boxhalf, boxhalf, 'MOL')

select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id=segid))

lambda_values = set_lambda(nelec=8,nvdw=8)                     # w/ openmm or scalecharge

## nonbonded options
my_nbonds = pycharmm.NonBondedScript(
    cutnb=12.0, ctonnb=10.0, ctofnb=11.0,
    eps=1.0,
    cdie=True,
    atom=True, vatom=True,
    switch=True, vfswitch=False,vswitch=True,
    noEwald=True) #, inbfrq=10)

#########################ANALYSIS########################
## loop over all fix lambda values using the dcd from current trajectory
for j in range(lambda_values.shape[0]):

    select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id=segid))
   
    my_nbonds.run()
    #nbonds.set_inbfrq(0)
    #nbonds.set_imgfrq(0)

    ## setup block information
    set_blocks_vacuum(res='MOL',solu='solu',lelec2=lambda_values[j,0],lvdw2=lambda_values[j,1])

    # open file for current i-loop
    print('Analyzing vacuum for lambda_elec={} and lambda_vdW={}'
            .format(lambda_values[j,0],lambda_values[j,1]))
    dcd_file = pycharmm.CharmmFile(file_name=f'win{i}/dcd/{solres}_flat.dcd',
                                           file_unit=1,
                                           formatted=False,read_only=True)

    #shake.on(bonh=True, fast=True, tol=1e-7)    
    pycharmm.lingo.charmm_script('traj query unit 1')
    nfiles = pycharmm.lingo.get_energy_value('NFILE')
    pycharmm.lingo.charmm_script('traj firstu 1 nunit 1')
    this_frame = []
    for frame in range(nfiles):
        pycharmm.lingo.charmm_script('traj read')
        energy.show()
        this_frame.append(energy.get_total())

    dcd_file.close()
    settings.set_verbosity(5)    

    print(this_frame)
    out = open(f'win{i}/post/win{j}.pkl','wb')
    pickle.dump(this_frame, out)
    out.close()
    
    #try: os.remove(f'win{i}/post/win{j}.dat')
    #except: FileNotFoundError
    #f1 = open(f'win{i}/post/win{j}.dat', 'a')
    #for fr in range(len(this_frame)):
    #    f1.write(f'{this_frame[fr]}\n')
  
    print(f'win{i}_win{j} done')
 
pycharmm.lingo.charmm_script('stop')
