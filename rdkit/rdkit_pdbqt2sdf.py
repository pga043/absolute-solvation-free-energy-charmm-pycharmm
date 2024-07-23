import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from networkx.algorithms import isomorphism
from rdkit.Geometry import Point3D
import sys

def mol_to_nx(mol):
      
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx())
    return G

def pdbqtToMol(fname):
    lines = [i[:54] for i in open(fname) if 'HETATM' in i]
    return Chem.MolFromPDBBlock('\n'.join(lines))


def makeMol(sdfFile, pdbqtFile):
    """
    Given an isomeric SMILES code, and it's corresponding pdbqt filename,
    infer the molecular graph from the pdbqt, then match to the actual graph
    using networkx. Use the matching to set the atom coordinates of rdkit molecule,
    resulting in an rdkit Mol object with both correct topology and correct coordinates.
    """
    sdf = Chem.SDMolSupplier(sdfFile, removeHs=False)
    smi = Chem.MolToSmiles(sdf[0])
    mol = Chem.MolFromSmiles(smi)

    pdbqt = pdbqtToMol(pdbqtFile)
    #get positions of the atoms from the pdbqt:
    positions = pdbqt.GetConformer(0).GetPositions()
    n_atoms = mol.GetNumAtoms()
    
    #get atom-to-atom matching.
    g_matching = isomorphism.GraphMatcher(
        mol_to_nx(mol), 
        mol_to_nx(pdbqt)
    )
    #if we don't have isomorphism, the inferred bonds from the pdbqt are wrong:
    assert(g_matching.is_isomorphic())
    
    #get positions of the atoms from the pdbqt:
    positions = pdbqt.GetConformer(0).GetPositions()
        
    #give the mol a blank conformer
    mol.AddConformer(Chem.Conformer(n_atoms))
    conf = mol.GetConformer(0)
    for i in range(n_atoms):
        x,y,z = positions[g_matching.mapping[i]]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    return mol


def to_sdf(original_sdf, docked_pdbqt, docked_sdf):
    mol = makeMol(original_sdf, docked_pdbqt)
    molh = Chem.AddHs(mol, addCoords=True)
    #print(Chem.MolToMolBlock(molh))
    with Chem.SDWriter(docked_sdf) as writer:
    #     Chem.AllChem.EmbedMultipleConfs(molh)
    #     Chem.AllChem.UFFOptimizeMolecule(molh, confId=-1)
         with open(docked_pdbqt) as f:
              score = f.readline()
         molh.SetProp('docking_score', score)
         writer.write(molh)


to_sdf(sys.argv[1], sys.argv[2], sys.argv[3]) 

#print(Chem.MolToMolBlock(mol))


quit()
