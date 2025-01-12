import pickle
import os

wins = 15

for win in range(wins):
    for i in range(wins):
        f = open(f'win{win}/post/win{i}.pkl', 'rb')
        eners = pickle.load(f)
        try: os.remove(f'win{win}/post/win{i}.dat')
        except: FileNotFoundError
        f1 = open(f'win{win}/post/win{i}.dat', 'a')
        for ener in range(len(eners)):
            f1.write(f'{eners[ener]}\n')
   
