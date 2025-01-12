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
def set_blocks_msld(env,solu,lmbda=1.0,scalecharge=False,scalefactor=1.0):
    """Set-up msld block definitions and values:
       env <= CHARMM define for the environment selection
       solu <= CHARMM define for the solute segment to be perturbed
       lmbda <= lambda scaling for set-up of MSLD
       scalecharge <= logical specifying charge scaling in effect
       scalefactor <= scale factor to acheive separation of charge and vdW scaling 
                      MSLD"""
    
    pycharmm.lingo.charmm_script('''
define solu select resname {} end
define env select resname tip3 end
block
   clear
end
block 3
   call 1 select {} end
   call 2 select {} end
   call 3 select segid dum end

   qldm theta
   lang temp 298.0

   ldinitialize 1 1.0 0 1 0 5
   ldinitialize 2 {:8.4f} 0 1 0 5
   ldinitialize 3 {:8.4f} 0 1 0 5

   rmlambda bond angle dihe impr

    exclude 2 3
    ! elecrostatics and soft core
    pmel ex
    soft on

    msld 0 1 1 ffix
    msmatrix
    ldbi  0       ! no biasing potential

end '''.format(resname,env,solu,
               lmbda,1.0-lmbda))
    return

################################
def set_lambda(nelec=10,nvdw=15,nshift=0):
    l = np.zeros((nvdw+nshift,2), dtype=float)
    l[:,1]=1.0
    for i in range(nelec): l[i,0] = 1-np.sin(i*np.pi/(2*(nelec-1)))
    for i in range(nvdw): l[i+nshift,1] = 1-np.sin((i)*np.pi/(2*(nvdw-1)))
    return l

################################
# Ensure that FFT grid is product of small primes 2, 3, 5
def is_factor(n):
    if (n % 2 != 0): return False  # favors even number
    while n:
        flag = False
        for x in (2,3,5):
            if n % x == 0:
               n = n / x
               flag = True
               break

        if flag: continue
        break

    if n == 1: return True
    return False

def checkfft(n, margin = 5):
    n = int(n) + margin
    while 1:
        if is_factor(n): break
        else: n += 1
    return n

##########################SETUP GLOBAL PARAMETERS##########################
T = 298         # K
k_B = 0.00198   # kcal/mol/K
ctofnb = 10
cutnb = ctofnb + 2.0
cutim = cutnb
ctonnb = ctofnb - 2.0

## Parveen Gartan, 23 Dec 2024
# Scripting starts here
solres = 'lig0'
i = int(sys.argv[1])

## make some directories and sub-directories
if not os.path.isdir(f'win{i}'): os.system(f'mkdir win{i}')
if not os.path.isdir(f'win{i}/post'): os.system(f'mkdir win{i}/post')

## read topology and paramter files
read.rtf('prep/toppar/top_opls_aam.inp')
read.prm('prep/toppar/par_opls_aam.inp', flex=True)
read.stream('prep/toppar/toppar_dum_noble_gases.str')
read.rtf(f'prep/{solres}_charmm.rtf', append=True)
read.prm(f'prep/{solres}_charmm.prm', flex=True, append=True)

## read psf and coor
read.psf_card(f'win{i}/minimized.psf')
read.coor_card(f'win{i}/minimized.crd')

resname = psf.get_res()[0].lower() # Assumes first residue is solute residue
segid = psf.get_segid()[0].upper()

## setup box and crystal information
stats = coor.stat()
xsize = stats['xmax'] - stats['xmin']
ysize = stats['ymax'] - stats['ymin']
zsize = stats['zmax'] - stats['zmin']

boxsize = 26.5286831 #max(xsize, ysize, zsize)
boxhalf = boxsize / 2.0

crystal.define_cubic(boxhalf*2)
crystal.build(boxhalf)

boxhalf = 0.0 # for blade image centering
image.setup_segment(boxhalf,boxhalf, boxhalf, 'MOL')
image.setup_segment(boxhalf,boxhalf, boxhalf, 'DUM')
image.setup_residue(boxhalf, boxhalf, boxhalf, 'TIP3')

select.store_selection('ENV',pycharmm.SelectAtoms(res_name='TIP3'))
select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id=segid))

lambda_values = set_lambda(nelec=15,nvdw=15)                     # w/ openmm or scalecharge

## read nonbonds with PME
fft = checkfft(n=np.ceil(boxsize),margin=0)
nb_wPME = pycharmm.NonBondedScript(cutnb=cutnb, cutim=cutim,
                                   ctonnb=ctonnb, ctofnb=ctofnb,
                                   eps=1.0,
                                   cdie=True,
                                   atom=True, vatom=True,
                                   switch=True, vfswitch=True, vswitch=False,
                                   inbfrq=-1, imgfrq=-1,
                                   lrc_ms=False,
                                   ewald=True,pmewald=True,kappa=0.32,
                                   fftx=fft,ffty=fft,fftz=fft,order=4)

#########################ANALYSIS########################
## loop over all fix lambda values using the dcd from current trajectory
for j in range(lambda_values.shape[0]):

    select.store_selection('ENV',pycharmm.SelectAtoms(res_name='TIP3'))
    select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id=segid))
   
    nb_wPME.run() 
    ## setup block information
    set_blocks_msld(env='ENV',solu='SOLU',
                            lmbda=lambda_values[j,1])

    # open file for current i-loop
    print('Analyzing solvent for lambda_elec={} and lambda_vdW={}'
            .format(lambda_values[j,0],lambda_values[j,1]))
    dcd_file = pycharmm.CharmmFile(file_name=f'win{i}/dcd/{solres}_flat.dcd',
                                           file_unit=1,
                                           formatted=False,read_only=True)
    
    pycharmm.lingo.charmm_script('traj query unit 1')
    nfiles = pycharmm.lingo.get_energy_value('NFILE')
    pycharmm.lingo.charmm_script('traj firstu 1 nunit 1')
    pycharmm.lingo.charmm_script('domdec gpu only dlb on ndir 1 1 1')
    this_frame = []
    for frame in range(nfiles):
        pycharmm.lingo.charmm_script('traj read')
        energy.show()
        this_frame.append(energy.get_total())

    dcd_file.close()
    settings.set_verbosity(5)
    #pycharmm.charmm_script('blade off')

    #print(this_frame)
    out = open(f'win{i}/post/win{j}.pkl','wb')
    pickle.dump(this_frame, out)
    out.close()
    #try: os.remove(f'win{i}/post/win{j}.dat')
    #except: FileNotFoundError
    #f1 = open(f'win{i}/post/win{j}.dat', 'a')
    #for fr in range(len(this_frame)):
    #    f1.write(f'{this_frame[fr]}\n')
    
pycharmm.lingo.charmm_script('stop')
