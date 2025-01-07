import pandas as pd

df = pd.read_csv('ct3c00691_si_002.csv', sep=',', header=0)

print('name, expt_dG, gaff_dG, boresch_dG, boresch_dG_lrc, lrc')
for i in range(10):
    name = df.values[i][2]
    expt = df.values[i][3]
    gaff = df.values[i][5]
    sb   = df.values[i][10]
    sb_lrc = df.values[i][12]
    lrc  = df.values[i][13]

    print(f'{name}, {expt}, {gaff}, {sb}, {sb_lrc}, {lrc}')
     

