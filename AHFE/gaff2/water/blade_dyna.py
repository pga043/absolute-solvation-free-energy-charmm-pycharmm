import os
import sys, subprocess
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
## Function to add dummy atom reference to psf as separate segment
def add_dummy(a=None):
    """build a dummy atom into the system and place it at the center
    of the solute.
    input: a <- vector representing geometric center of the solute"""
    
    pycharmm.lingo.charmm_script('''
* generate a dummy atom
*

set dummy = 1
!
! Add a dummy atom if asked 
set dummylabel

if @dummy .gt. 0 then 
	wrnlev 4
	
	READ SEQU DUM 1
	GENERate DUM setup 

	! 
	! Get the COM of the solute
	! to place the dummy atom initially
	coor stat sele ( .not. segid DUM )  end
	calc xdum = ?XAVE + 5
	calc ydum = ?YAVE + 5
	calc zdum = ?ZAVE + 5

	! 
	! Position anchor dummy at (xdum, ydum, zdum)
	scalar X set @xdum sele type DUM end
	scalar Y set @ydum sele type DUM end
	scalar Z set @zdum sele type DUM end

	scalar MASS set 12 sele type DUM end
	
endif
    '''.format(a.x,a.y,a.z))
    return

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

#####################################
def LRC_solute(ctofnb=None):
    """This function computes the long-range dispersion correction for the annilation of the
    ligand in solute.
    input: ctofnb <- cutoff for nonbonded interactions
    scalar.get_effect() => solute eps value
    scalar.get_radius() => solute rmin/2 value
    Dispersion coefficient C_iO => -eps_iO*rmin_iO^6
    eps_iO = sqrt(eps_ii*eps_OO) => vdW epsilon for atom i in solute and water oxygen O
    rmin_iO = rmin_ii/2 epsrmin_OO/2 => vdW rmin for atom i in solute and water oxygen O
    TIP3P Oxygen"""
    wateps = -.1521  # TIP3P emin  (kcal/mol)
    watrad = 1.7682  # TIP3P rmin/2 (A)
    watdens = 0.03222 # density of water => 1g/cc in molecules/A^3
    F_LR = -np.sum((8/3)*np.pi*watdens*np.sqrt(np.asarray(scalar.get_effect(),
                                                          dtype=float)*wateps)*\
                   np.power((np.asarray(scalar.get_radius(),
                                        dtype=float)+watrad),6)/ctofnb**3)
    print('Long-range dispersion correction to free energy: {:6.2f} kcal/mol'.format(F_LR))
    return F_LR

def get_gpu():
    gpu_info = subprocess.check_output("nvidia-smi -L", shell=True)
    #print(gpu_info.decode("ascii"))
    gpu_name = gpu_info.decode("ascii").split(':')[1].split('(')[0]
    #print(gpu_name)
    return gpu_name
##======================================================##

##########################SETUP GLOBAL PARAMETERS##########################
T = 298         # K
k_B = 0.00198   # kcal/mol/K
ctofnb = 10
cutnb = ctofnb + 2.0
cutim = cutnb
ctonnb = ctofnb - 2.0

useblade = 'prmc pref 1 iprs 100 prdv 100'
leap = 'True'

nequil = 250000  # number of equilibration steps
nprod = 2000000 # 2000000 production dynamics
nsavc = 10000   # 10000 save frequency

## Parveen Gartan, 23 Dec 2024
# Scripting starts here
solres = str(sys.argv[2])

print(f'running on GPU: {get_gpu()}')
## read topology and paramter files
read.rtf('prep/toppar/parm14sb_all.rtf')
read.prm('prep/toppar/parm14sb_all.prm', flex=True)
read.stream('prep/toppar/toppar_dum_noble_gases.str')
read.rtf(f'prep/{solres}_charmm.rtf', append=True)
read.prm(f'prep/{solres}_charmm.prm', flex=True, append=True)

## read psf and coor
read.psf_card('prep/{}-solvated.psf'.format(solres))
read.coor_card('prep/{}-solvated.crd'.format(solres))

coor.orient(by_rms=False,by_mass=False,by_noro=False)
resname = psf.get_res()[0].lower() # Assumes first residue is solute residue
segid = psf.get_segid()[0].upper()

## add dummy atom
settings.set_bomb_level(-2)
add_dummy((coor.get_positions()).mean())
settings.set_bomb_level(0)

## setup box and crystal information
stats = coor.stat()
xsize = stats['xmax'] - stats['xmin']
ysize = stats['ymax'] - stats['ymin']
zsize = stats['zmax'] - stats['zmin']

get_boxsize = f"grep 'DIMENSION' prep/{solres}-solvated.psf | awk '{{print $NF}}'"
boxsize = float(subprocess.getoutput(get_boxsize))
#boxsize = 26.5262598 #max(xsize, ysize, zsize)
boxhalf = boxsize / 2.0

crystal.define_cubic(boxhalf*2)
crystal.build(boxhalf)

boxhalf = 0.0 # for blade image centering
image.setup_segment(boxhalf,boxhalf, boxhalf, 'MOL')
image.setup_segment(boxhalf,boxhalf, boxhalf, 'DUM')
image.setup_residue(boxhalf, boxhalf, boxhalf, 'TIP3')

select.store_selection('ENV',pycharmm.SelectAtoms(res_name='TIP3'))
select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id=segid))

# Loop over sets of lambda values for fixed lambda sampling
lambda_values = set_lambda(nelec=15,nvdw=15)                     # w/ openmm or scalecharge

## win 1
i = int(sys.argv[1])

## make some directories and sub-directories
if not os.path.isdir(f'win{i}'): os.system(f'mkdir win{i}')
if not os.path.isdir(f'win{i}/res'): os.system(f'mkdir win{i}/res')
if not os.path.isdir(f'win{i}/dcd'): os.system(f'mkdir win{i}/dcd')

## long range correction
F_LRC = 0.0
#F_LRC = LRC_solute(ctofnb=ctofnb)  # Calculate long-range correction

## read nonbonds with PME
fft = checkfft(n=np.ceil(boxsize),margin=0)
nb_wPME = pycharmm.NonBondedScript(cutnb=cutnb, cutim=cutim,
                                   ctonnb=ctonnb, ctofnb=ctofnb,
                                   eps=1.0,
                                   cdie=True,
                                   atom=True, vatom=True,
                                   switch=True, vfswitch=False, vswitch=True,
                                   inbfrq=-1, imgfrq=-1,
                                   lrc_ms=True,
                                   ewald=True,pmewald=True,kappa=0.32,
                                   fftx=fft,ffty=fft,fftz=fft,order=4)
nb_wPME.run()


print('minimizing')

cons_fix.setup(selection=pycharmm.SelectAtoms(seg_id='DUM'))
cons_fix.setup(selection=pycharmm.SelectAtoms(seg_id='MOL'))
minimize.run_sd(nstep=200, tolenr=1e-3, tolgrd=1e-3)
cons_fix.turn_off()

minimize.run_sd(nstep=200, tolenr=1e-3, tolgrd=1e-3)
energy.show()

write.psf_card(f'win{i}/minimized.psf')
write.coor_card(f'win{i}/minimized.crd')

## setup block information
set_blocks_msld(env='ENV',solu='SOLU',
                            lmbda=lambda_values[i,1])

#########################DYNAMICS########################
# Setup and run dynamics at this set of lambda values
shake.on(bonh=True, fast=True, tol=1e-7)
pycharmm.lingo.charmm_script('energy blade')
imgfrq = -1

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
                                 ilbfrq=0,ihbfrq=0,
                                 blade=useblade)
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
                                ilbfrq=0, ihbfrq=0,
                                blade=useblade)

prod_dyn.run()

res_file.close()
lam_file.close()
dcd_file.close()

write.coor_pdb(f'win{i}/{solres}_prod.pdb')
pycharmm.lingo.charmm_script('blade off')

pycharmm.lingo.charmm_script('stop')
