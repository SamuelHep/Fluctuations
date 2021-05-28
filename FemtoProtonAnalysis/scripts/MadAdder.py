import glob
import math
import os
import sys

def adder(newfile,list_to_add):
    
    command = 'hadd ' + newfile + ' '
    n = len(list_to_add)
    for i in range(n):        
    
        command = command + ' ' + list_to_add[i]
        
    os.system( command )

#inputs args
prefix = sys.argv[1]
outfile = sys.argv[2]
min_add = 10

#get root files in the directory
rootfiles = list(glob.glob(prefix + '*'))
nfiles = len(rootfiles)

n = int(math.log(nfiles)/math.log(min_add))
if ( pow(min_add,n) < nfiles ):
    n = n + 1

#make names 
prefix_names = [] 

for i in range(n):

    print ( i )

    nfiles = len(rootfiles)
    ndiv = nfiles/min_add
    
    name = '_temp_adder_v' + str(i)
    prefix_names.append(name)
    newrootfiles = []
    newfile = ''

    if ( nfiles <= min_add ):
        adder(outfile,rootfiles)
        break;

    else:
        for ii in range(ndiv):
            
            newfile = prefix_names[i] + "_" + str(ii) + ".root"
            print ( newfile)
            newrootfiles.append(newfile)
            add_list = rootfiles[min_add*ii:(min_add*ii+min_add)]
            if (nfiles % min_add != 0 and ii == 0):
                add_list = add_list + rootfiles[-(nfiles % min_add):]
            
            adder(newfile,add_list)
    
        rootfiles = newrootfiles 

#clean up
remove_command = 'rm _temp_adder_v*'
os.system( remove_command )
