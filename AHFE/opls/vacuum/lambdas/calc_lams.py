import numpy as np

def set_lambda(nelec=10,nvdw=15,nshift=0):
    l = np.zeros((nvdw+nshift,2), dtype=float)
    l[:,1]=1.0
    for i in range(nelec): l[i,0] = 1-np.sin(i*np.pi/(2*(nelec-1)))
    for i in range(nvdw): l[i+nshift,1] = 1-np.sin((i)*np.pi/(2*(nvdw-1)))
    return l

lams = set_lambda(nelec=8,nvdw=8)
for i in range(len(lams)):
    with open(f'lam{i+1}.inp', 'a') as f:
        #print(lams[i][0])
        f.write(f'set ele = {lams[i][0]} \n')
        f.write(f'set vdw = {lams[i][0]}')
