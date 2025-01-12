#export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_charmm/
import os
import numpy as np
import pycharmm
import pycharmm.read as read
import pycharmm.settings as settings
import pycharmm.select as select
import pycharmm.scalar as scalar

# taken from:     pyCHARMM-Workshop/4AbsoluteSolvation/absolute_solvation.ipynb

solres = 'lig7'

read.rtf('prep/toppar/top_opls_aam.inp')
read.prm('prep/toppar/par_opls_aam.inp', flex=True)
pycharmm.lingo.charmm_script('stream prep/toppar/toppar_dum_noble_gases.str')

read.rtf('prep/{}_charmm.rtf'.format(solres),append=True)
read.prm('prep/{}_charmm.prm'.format(solres), flex=True, append=True)
settings.set_bomb_level(0)

read.psf_card(f'prep/{solres}-solvated.psf')
read.coor_card(f'prep/{solres}-solvated.crd')

settings.set_bomb_level(-2)
pycharmm.lingo.charmm_script(f'''
delete atom sele .not. segid MOL end
''')
settings.set_bomb_level(0)
#select.store_selection('SOLU',pycharmm.SelectAtoms(seg_id='MOL'))

#print(scalar.get_effect())
#print(scalar.get_radius())

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


lrc = LRC_solute(10)

f1 = open('../vacuum/1postprocess/results.txt', 'r')
for lines in f1.readlines():
    line = lines.split()
    if line[0] == 'MBAR':
        vac_dG = float(line[-4])
        vac_ddG = float(line[-2])

f2 = open('1postprocess/results.txt', 'r')
for lines in f2.readlines():
    line = lines.split()
    if line[0] == 'MBAR': 
        wat_dG = float(line[-4])
        wat_ddG = float(line[-2])

dG = round(vac_dG - wat_dG + lrc, 2)
ddG = round(np.sqrt(vac_ddG**2 + wat_ddG**2), 2)

print(f'Hydration free energy (kcal/mol) = {dG} +/- {ddG}')

