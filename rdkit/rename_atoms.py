#! /usr/bin/env python

#
# Take a mol2 file and rename all the atom names to be unique!
#

import sys,os

if len(sys.argv) != 3:
    print("\n This script needs two command line argument to work:")
    print("    1st = name of the mol2 file to process")
    print("    2nd = name of renamed mol2 file to be written")
    quit()

filename = sys.argv[1]

fp=open(filename,'r')
gp=open(sys.argv[2],'w')
c_count = 1 
h_count = 1
o_count = 1
f_count = 1
n_count = 1
s_count = 1
cl_count = 1
line=fp.readline()
while line:
    if line[0:13] == '@<TRIPOS>ATOM':
        gp.write(line)
        line=fp.readline()
        while line[0:13] != '@<TRIPOS>BOND':
            tmp=line.split()
            if tmp[1] == 'C':
                tmp[1] = tmp[1]+str(c_count)
                c_count += 1
            if tmp[1] == 'O':
                tmp[1] = tmp[1]+str(o_count)
                o_count += 1
            if tmp[1] == 'N':
                tmp[1] = tmp[1]+str(n_count)
                n_count += 1
            if tmp[1] == 'S':
                tmp[1] = tmp[1]+str(s_count)
                s_count += 1
            if tmp[1] == 'F':
                tmp[1] = tmp[1]+str(f_count)
                f_count += 1
            if tmp[1] == 'H':
                tmp[1] = tmp[1]+str(h_count)
                h_count += 1
            if tmp[1] == 'Cl':
                tmp[1] = tmp[1]+str(cl_count)
                cl_count += 1
            #print(len(tmp),tmp)
            #if int(tmp[0]) < 10:
                ##tmp[1] = tmp[1]+'00'+tmp[0]
            #    tmp[1] = tmp[1]+'0'+tmp[0]
            #elif int(tmp[0]) > 9 and int(tmp[0]) < 100:
                ##tmp[1] = tmp[1]+'0'+tmp[0]
            #    tmp[1] = tmp[1]+tmp[0]
            #elif int(tmp[0]) > 99 and int(tmp[0]) < 1000:
            #    tmp[1] = tmp[1]+tmp[0]
            #else:
            #    raise ValueError("More than 1000 atoms not handled")
            gp.write("%7s %-6s  %10s%10s%10s %-6s%3s  %-4s     %9s\n" % (tmp[0],tmp[1],
                     tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8]))
            line=fp.readline()

        #gp.write(line)
    else:
        gp.write(line)
        line=fp.readline()

fp.close()
gp.close()

# overwrite the original file
#os.system('mv tmp.mol2 '+filename)

# finished
