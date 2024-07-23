import os, sys, glob, subprocess, re
import math
import numpy as np

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

import Mol2writer

p = AllChem.ETKDGv2()
p.verbose = False

ref_mol   = sys.argv[1] 
target_mol = sys.argv[2]
out_mol = sys.argv[3]

##***************************************
## open a reference sturcture upon which others will be aligned
## preserving input coordinates (must be 3D structure data file)
ref1 = Chem.SDMolSupplier(ref_mol, removeHs=False)
for ref in ref1:
    #ref1_smi = Chem.MolToSmiles(ref)
    #ref1_mol = Chem.MolFromSmiles(ref1_smi)
    ref1_molh = ref #Chem.AddHs(ref1_mol)

#AllChem.EmbedMultipleConfs(ref1_molh)

## get MMFF parameters for the reference molecule
mmff_ref_param = AllChem.MMFFGetMoleculeProperties(ref1_molh)
ref_mol2 = ref1_molh

##--------------------------------------------------

sdf = Chem.SDMolSupplier(target_mol, removeHs=False)

tname = sdf[0].GetProp('_Name')
#for m in sdf:
#    print(m.GetProp('_Name'))
#quit()

molecules = []
smi = [Chem.MolToSmiles(x) for x in sdf]
mol = [Chem.MolFromSmiles(x) for x in smi]
molecules.append(Chem.AddHs(mol[0]))

confs = 150 ## number of conformers to generate
print(f'Generating {confs} conformers for each molecule')
for mol in molecules[0:]:
    AllChem.EmbedMultipleConfs(mol, confs, p)

## get MMFF parameters for each molecule
mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
#mmff_ref_param = mmff_params[0]
mmff_prob_params = mmff_params[0:]
#ref_mol2 = molecules[0]
prob_mols_2 = molecules[0:]


## perform alignment and save the best scored conformer
print('Performing alignment and saving the conformer with best score only.')


try:
   os.remove(str(out_mol)+'.sdf')
   os.remove(str(out_mol)+'.mol2')
except FileNotFoundError:
    print('File does not exists.')


with Chem.SDWriter(str(out_mol)+'.sdf') as writer: 
    #writer.write(ref_mol2) ## save the reference if you want to look at the alignment
    pyO3A_score = []
    rmsd = []
    for idx, mol in enumerate(prob_mols_2):
        tempscore = []
        temprmsd = []
        for cid in range(confs):
            pyO3A = rdMolAlign.GetO3A(mol, ref_mol2, mmff_prob_params[idx], mmff_ref_param, cid, 0)
            temprmsd.append(pyO3A.Align())
            tempscore.append(pyO3A.Score())
        best = np.argmax(tempscore)
        rmsd.append(temprmsd[best])
        ## save sdf file
        ## write names, to be used with network in openFE
        mol.SetProp('_Name', str(tname))
        mol.SetProp('ID', str(int(best)))
        mol.SetProp('Open3D alignment rmsd wrt reference: ', str(round(rmsd[idx], 4)))
        #Chem.MolToV3KMolFile(mol, str(target_mol) + '_wrt_' + str(ref_mol) + '.mol',  confId=int(best))
        #Chem.MolToPDBFile(mol, str(target_mol) + '_wrt_' + str(ref_mol) + '.pdb',  confId=int(best))
        #print(Mol2writer.MolToCommonMol2Block(mol))
        writer.write(mol, confId=int(best))
        Mol2writer.MolToMol2File(mol, str(out_mol) + '.mol2', confId=int(best))
        pyO3A_score.append(tempscore[best])


quit()
