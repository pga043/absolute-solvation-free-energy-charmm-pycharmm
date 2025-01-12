import os
import sys
import numpy as np

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

##======================== functions ===============##
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

#####################################
##======================================================##

##########################SETUP GLOBAL PARAMETERS##########################
T = 298         # K
k_B = 0.00198   # kcal/mol/K
ctofnb = 10
cutnb = ctofnb + 2.0
cutim = cutnb
ctonnb = ctofnb - 2.0

nequil = 20000  # number of equilibration steps
nprod = 4000000 # 2000000 production dynamics
nsavc = 20000   # 10000 save frequency

## Parveen Gartan, 23 Dec 2024
# Scripting starts here
solres = str(sys.argv[2])

## read topology and paramter files
read.rtf('prep/toppar/parm14sb_all.rtf')
read.prm('prep/toppar/parm14sb_all.prm', flex=True)
read.rtf(f'prep/{solres}_charmm.rtf', append=True)
read.prm(f'prep/{solres}_charmm.prm', flex=True, append=True)

## read psf and coor
read.psf_card('prep/{}.psf'.format(solres))
read.coor_card('prep/{}.crd'.format(solres))

coor.orient(by_rms=False,by_mass=False,by_noro=False)
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

## nonbonded options
my_nbonds = pycharmm.NonBondedScript(
    cutnb=12.0, ctonnb=10.0, ctofnb=11.0,
    eps=1.0,
    cdie=True,
    atom=True, vatom=True,
    switch=True, vfswitch=False,vswitch=True,
    noEwald=True)

my_nbonds.run()

# Loop over sets of lambda values for fixed lambda sampling
lambda_values = set_lambda(nelec=8,nvdw=8)                     # w/ openmm or scalecharge
print(lambda_values)

## win 1
i = int(sys.argv[1])

## make some directories and sub-directories
if not os.path.isdir(f'win{i}'): os.system(f'mkdir win{i}')
if not os.path.isdir(f'win{i}/res'): os.system(f'mkdir win{i}/res')
if not os.path.isdir(f'win{i}/dcd'): os.system(f'mkdir win{i}/dcd')

#print('minimizing')

#minimize.run_sd(nstep=200, tolenr=1e-3, tolgrd=1e-3)
#energy.show()

write.psf_card(f'win{i}/minimized.psf')
write.coor_card(f'win{i}/minimized.crd')

## setup block information
set_blocks_vacuum(res='MOL',solu='solu',lelec2=lambda_values[i,0],lvdw2=lambda_values[i,1])

#########################DYNAMICS########################
# Setup and run dynamics at this set of lambda values
shake.on(bonh=True, fast=True, tol=1e-7)
imgfrq = 0

dyn.set_fbetas(np.full((psf.get_natom()),1.0,dtype=float))

res_file = pycharmm.CharmmFile(file_name=f'win{i}/res/{solres}.res', file_unit=2,
                               formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name=f'win{i}/dcd/{solres}.lam', file_unit=3,
                               formatted=False,read_only=False)

my_dyn = pycharmm.DynamicsScript(leap=True, lang=False, start=True,
                                 nstep=nequil, timest=0.002,
                                 firstt=298.0, finalt=298.0, tbath=298.0,
                                 tstruc=298.0,
                                 teminc=0.0, twindh=0.0, twindl=0.0,
                                 iunwri=res_file.file_unit,
                                 iunlam=lam_file.file_unit,
                                 inbfrq=-1, imgfrq=imgfrq,
                                 iasors=0, iasvel=1, ichecw=0, iscale=0,
                                 iscvel=0,
                                 echeck=-1.0, nsavc=0, nsavv=0, nsavl=0, ntrfrq=0,
                                 isvfrq=nsavc,
                                 iprfrq=2*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                 ilbfrq=0,ihbfrq=0,)

my_dyn.run()

res_file.close()
lam_file.close()

res_file = pycharmm.CharmmFile(file_name=f'win{i}/res/{solres}_flat.res', file_unit=2,
                               formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name=f'win{i}/dcd/{solres}_flat.lam', file_unit=3,
                               formatted=False,read_only=False)
dcd_file = pycharmm.CharmmFile(file_name=f'win{i}/dcd/{solres}_flat.dcd', file_unit=1,
                               formatted=False,read_only=False)

prod_dyn = pycharmm.DynamicsScript(leap=True, lang=False, restart=False, nstep=nprod, timest=0.002,
                                 firstt=298.0, finalt=298.0, tbath=298.0,
                                tstruc=298.0,
                                teminc=0.0, twindh=0.0, twindl=0.0,
                                iunwri=res_file.file_unit,
                                iunrea=res_file.file_unit,
                                iunlam=lam_file.file_unit,
                                inbfrq=-1, imgfrq=imgfrq,
                                iuncrd=dcd_file.file_unit,
                                iasors=0, iasvel=1, ichecw=0, iscale=0, iscvel=0,
                                echeck=-1.0, nsavc=nsavc, nsavv=0, nsavl=0,
                                ntrfrq=0, isvfrq=nsavc,
                                iprfrq=5*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                ilbfrq=0, ihbfrq=0,)
                                

prod_dyn.run()

res_file.close()
lam_file.close()
dcd_file.close()

write.coor_pdb(f'win{i}/{solres}_prod.pdb')

pycharmm.lingo.charmm_script('stop')
