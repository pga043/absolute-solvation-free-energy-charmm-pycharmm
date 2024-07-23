import MDAnalysis as mda
from MDAnalysis.analysis import distances
import networkx as nx
import os

u1 = mda.Universe('mol22.psf', 'mol22.crd')
u2 = mda.Universe('sr12.psf', 'sr12.crd')

noe_out = 'noe_restraints.str'

u1_atoms = u1.select_atoms('not name H*')
u2_atoms = u2.select_atoms('not name H*')

#print(u2_atoms)
#print(len(u2_atoms))

dist_arr = distances.distance_array(u1_atoms, u2_atoms)

#print(dist_arr)

#len(dist_arr)
#print(len(dist_arr[0]))
#print(len(dist_arr[1]))

#min(dist_arr[0])

G = nx.Graph()

for i in range(len(dist_arr)):
    mindist = min(dist_arr[i])
    for j in range(len(dist_arr[i])):
        if dist_arr[i][j] == mindist:
           #print(i)
           #print(u1_atoms[i], u2_atoms[j], round(mindist,1))
           G.add_node(u1_atoms[i].name, segid=u1_atoms[i].segid)
           G.add_node(u2_atoms[j].name, segid=u2_atoms[j].segid)
           G.add_edge(u1_atoms[i].name, u2_atoms[j].name, distance=round(mindist, 1))


#for ea in G.edges(data=True):
#    print(ea[0], G.nodes[ea[0]]['segid'], ea[1], G.nodes[ea[1]]['segid'])
#    print(ea[1], G.nodes[ea[1]]['segid'])

mol1_skip = list(u1.select_atoms('(name S* or name F* or name Cl*) and not name H*').names)
mol2_skip = list(u2.select_atoms('(name S* or name F* or name Cl*) and not name H*').names)

atoms_skip = mol1_skip + mol2_skip

print(atoms_skip)

try:
   os.remove(noe_out)
except FileNotFoundError:
    print('File does not exists.')


H = nx.Graph() ## save nodes here and avoid duplication
cutoff = 1.0
kl = 40 # force constant
dl = 0.4 

with open(noe_out, 'a') as f:
     f.write(f'NOE \n')
     for ea in G.edges(data=True):
         #print(ea[0], ea[1], G.edges[ea[0], ea[1]]['distance'])
        try:
           check1 = any(ea[0] == atom for atom in atoms_skip)
           check2 = any(ea[1] == atom for atom in atoms_skip)
           if (check1 == False) and (check2 == False):
              if G.edges[ea[0], ea[1]]['distance'] < cutoff:
                 print(ea[0], ea[1], G.edges[ea[0], ea[1]]['distance'])
                 node1 = H.has_node(ea[0])
                 node2 = H.has_node(ea[1])
                 if (node1 == True) or (node2 == True):
                    print(f'One of the node already present.')
                 else:
                    try:
                       rmin = round(G.edges[ea[0], ea[1]]['distance'] - dl, 1)
                    except rmin < 0:
                          rmin = 0.0
                    rmax = round(G.edges[ea[0], ea[1]]['distance'] + dl, 1)
                    segid0 = G.nodes[ea[0]]['segid']
                    segid1 = G.nodes[ea[1]]['segid']
                    f.write(f'   assign sele atom {segid0} 1 {ea[0]} end sele atom {segid1} 1 {ea[1]} end - \n')
                    f.write(f'   kmin {kl} rmin {rmin} kmax {kl} rmax {rmax} fmax 2.0 rswitch 2.0 sexp 1.0\n')
                    #print(u, v, G.edges[u, v]['distance'])
                    H.add_node(ea[0])
                    H.add_node(ea[1])
        except:
              None

     f.write('\n')
     f.write(f'   print \n   print anal \n')
     f.write(f'END \n')
     

   

quit()

