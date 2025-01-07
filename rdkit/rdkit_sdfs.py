import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import glob, os
import pandas as pd
import Mol2writer

confs = 50

x=00
df = pd.read_csv('../raw_data/ct3c00691_si_002.csv', sep=',', header=1)
for i in range(10):
    name = df.values[i][2]
    print(name)
    smi = df.values[i][1] 
       
    mol = Chem.MolFromSmiles(smi)
    molh = Chem.AddHs(mol)
    confIds = Chem.AllChem.EmbedMultipleConfs(molh, confs)
    score = []
    for confId in range(confs):
        Chem.AllChem.UFFOptimizeMolecule(molh, confId=confId)
        ff = AllChem.UFFGetMoleculeForceField(molh, confId=confId)
        score.append(ff.CalcEnergy())
        #print(ff.CalcEnergy())
    minE = min(score, key=abs)
    #print(minE)
    with Chem.SDWriter(f'lig{x}.sdf') as writer:
        writer.write(molh, confId=int(score.index(minE)))
    print('-------------')
    Mol2writer.MolToMol2File(molh, f'lig{x}_tmp.mol2', confId=int(score.index(minE)))
    os.system(f'python rename_atoms.py lig{x}_tmp.mol2 lig{x}.mol2')
    os.system(f'rm lig{x}_tmp.mol2')
    x +=1


##=======================================================##
##=======================================================##

quit()

x = 8
confs = 100
for mol in suppl:
    molh = Chem.AddHs(mol)
    confIds = Chem.AllChem.EmbedMultipleConfs(molh, confs)
    score = []
    for confId in range(confs):
        Chem.AllChem.UFFOptimizeMolecule(molh, confId=confId)
        ff = AllChem.UFFGetMoleculeForceField(molh, confId=confId)
        score.append(ff.CalcEnergy())
        print(ff.CalcEnergy())
    minE = min(score, key=abs)
    print(minE)
    with Chem.SDWriter(f'sprios_dec_{x}.sdf') as writer:
         writer.write(molh, confId=int(score.index(minE)))
    print('-------------')
    with open(f'sprios_dec_{x}.pdb', 'w') as f:
         pdb = Chem.MolToPDBBlock(molh, confId=int(score.index(minE)), flavor=4)
         pdb = Chem.AtomPDBResidueInfo()
         pdb.SetResidueName('MOL')
         f.write(pdb)
    x += 1


#pdbblock = Chem.MolToPDBBlock(mh)
#open("/tmp/noh.pdb",'w').write(Chem.MolToPDBBlock(mh, flavor=4))
    
