import numpy as np
import mendeleev
import sys
from copy import deepcopy


'''
From Luis
'''

def getOPLSTypes(zMatFile,padding=4):
    """
    Get a {newAtType:oplsType} dictionary from
    LPG Z-MAT file.
    if padding == None, no padding will be done.
    if padding == int, then atom types will be padded
    with 0 until atom type have `padding` number of characters.
    """
    # Read Z-Matrix
    with open(zMatFile,'r') as f:
        lines = f.readlines()
    
    # Strip file of empty spaces
    lines = [line for line in lines if line.strip()]
    
    # Find atom type section start and get atoms section 
    ind = [i for i,line in enumerate(lines) if 'Final Non-Bonded' in line][0] + 1
    atoms = lines[ind:]
    
    # Get new atom types
    elements = [mendeleev.element(int(at.split()[1])).symbol for at in atoms]
    atomNames = [at.split()[0] for at in atoms]
    newAtTypes = [el+atName for el,atName in zip(elements,atomNames)]
    
    # Get OPLS atom types
    oplsTypes = []#[at.split()[2] for at in atoms]
    for at in atoms:
        try:
           if not at.split()[2][1] == '':
              if at.split()[2][1] == '!':
                 new = at.split()[2][0] + ''.join('X')
                 oplsTypes.append(new)
              elif at.split()[2][1] == '=':
                 new = at.split()[2][0] + ''.join('X')
                 oplsTypes.append(new)
              else:
                 new = at.split()[2]
                 oplsTypes.append(new)
        except IndexError:
              new = at.split()[2]
              oplsTypes.append(new)
    print(oplsTypes)
    if padding:
        oplsTypes = [oplsAt+''.join(['0' for _ in range(padding-len(oplsAt))]) for oplsAt in oplsTypes]
        #print(oplsTypes)
    return dict(zip(newAtTypes,oplsTypes))


def renameRTF(rtfFile, atMapping,outFile=None):
    """
    Return a new CHARMM rtf file called `outFile`
    with renamed atom types based on dictionary 
    atMapping and original `rtfFile`.

    If `outFile` is None, then it will return the
    file `f'{outFile}_renamed.rtf'`.
    """
    with open(rtfFile,'r') as f:
        lines = f.readlines()

    massStatements = [line for line in lines if line.startswith('MASS')]

    # Replace mass statements to show OPLS atom types
    for newAt, oplsAt in atMapping.items():
        for i,line in enumerate(massStatements):
            if newAt in line:
                massStatements[i] = line.replace(newAt,oplsAt)

    # Get rid of duplicates and re-enumerate
    newMass = []
    newIdx = 1
    for i,line in enumerate(massStatements):
        if sum([line.split()[2] in l for l in newMass]) == 0:
            newMass.append(f'MASS {newIdx} '+' '.join(line.split()[2:])+'\n') 
            newIdx += 1

    # Define atom statements
    atomStatements = [line for line in lines if line.startswith('ATOM')]
    for newAt, oplsAt in atMapping.items():
        for i,line in enumerate(atomStatements):
            if newAt in line:
                atomStatements[i] = line.replace(newAt,oplsAt)

    # Write new file by sections
    toWrite = [lines[0]]
    toWrite.extend(newMass)
    toWrite.extend([line for line in lines if line.startswith('AUTO')]) 
    toWrite.extend([line for line in lines if line.startswith('RESI')])
    toWrite.extend(atomStatements)
    toWrite.extend([line for line in lines if line.startswith('BOND')])
    toWrite.extend([line for line in lines if line.startswith('IMPR')])
    toWrite.extend([line for line in lines if line.startswith('PATCH')])
    toWrite.extend([line for line in lines if line.startswith('END')])

    with open(outFile,'w') as f:
        f.write(''.join(toWrite)) 

def renamePRM(prmFile,atMapping,duplicates=True,outFile=None):
    with open(prmFile,'r') as f:
        lines = f.readlines()
  
    # inds = [i for i,line in enumerate(lines) if any([line.startswith(x) for x in ['BOND','ANGLE','DIHEDRAL','IMPROPER','NONBONDED']])]
    # #BOND
    # # print(lines[inds[0]:inds[1]])
    # # #ANGLE
    # # print(lines[inds[1]:inds[2]])
    # # #DIHEDRAL
    # # print(lines[inds[2]:inds[3]])
    # # #IMPROPER
    # # print(lines[inds[3]:inds[4]])
    # # #NONBONDED 
    # # print(lines[inds[4]:])
    for i,line in enumerate(lines):
        for newAt, oplsAt in atMapping.items():
            if newAt in lines[i]:
                lines[i] = lines[i].replace(newAt,oplsAt)

    newlines = deepcopy(lines)
    # Get rid of duplicates
    if not duplicates:
        newlines = []
        for line in lines:
            if line not in newlines:
                newlines.append(line)

    with open(outFile,'w') as f:
        f.write(''.join(newlines))
        f.write('END')

if __name__ == '__main__':
    atMapping = getOPLSTypes(str(sys.argv[1]) +'.z')
    renameRTF(str(sys.argv[1]) + '.charmm.rtf',atMapping,outFile=str(sys.argv[1]) +'.tmp.rtf')
    renamePRM(str(sys.argv[1]) +'.charmm.prm',atMapping,outFile=str(sys.argv[1]) +'.tmp.prm',duplicates=False)


