import numpy as np



def set_lambda(nelec=10,nvdw=15,nshift=0):
    l = np.zeros((nvdw+nshift,2), dtype=float)
    l[:,1]=1.0
    for i in range(nelec): l[i,0] = 1-np.sin(i*np.pi/(2*(nelec-1)))
    for i in range(nvdw): l[i+nshift,1] = 1-np.sin((i)*np.pi/(2*(nvdw-1)))
    return l


#=================================================
#if scalecharge or openmm:
    #lambda_values = set_lambda(nelec=10,nvdw=15)                     # w/ openmm or scalecharge
    # In this example script we use a smaller number of windows (8 total)
#    lambda_values = set_lambda(nelec=5,nvdw=8)                     # w/ openmm or scalecharge
#else:
#    if not vacuumonly: lambda_values = set_lambda(nelec=15,nvdw=15)  # w/ blade and not scalecharge
#    else: lambda_values = set_lambda(nelec=8,nvdw=8)

#=================================================
lams = set_lambda(nelec=15,nvdw=15)
for i in range(len(lams)):
    with open(f'lam{i+1}.inp', 'a') as f:
        #print(lams[i][0])
        #f.write(f'set lenv = 1.0 \n')
        #f.write(f'set lelec = {lams[i][0]} \n')
        #f.write(f'set lvdw = {lams[i][1]} \n')
        #f.write(f'set lelec2 = {lams[i][0]} \n')
        #f.write(f'set lvdw2 = {lams[i][1]} \n')
        f.write(f'set lam1 = {lams[i][1]} \n')
        f.write(f'set lam2 = {1-lams[i][1]} \n')
