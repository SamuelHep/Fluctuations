#Get run#s from a list from get_.pl
#
#

import re

print "Get Run numbers from list"

st_phys = 'st_physics_'
st_raw = '_raw'
adc = '_adc_'
laser = '_laser'

outfile = open("filelist/good_3GeV.txt","w")
list_of_numbers = set([])

with open("filelist/3GeV_newProd_Fluct_GoodList.list") as run_file:
    for line in run_file:

        if re.search(adc,line):
            m1 = re.search(adc,line)
            m2 = re.search(st_raw,line)
        else:
            if not re.search(laser,line): 
                m1 = re.search(st_phys,line)
                m2 = re.search(st_raw,line)
                
                number_string = line[m1.end():m2.start()]
                print number_string
                number = [int(s) for s in re.findall(r'\b\d+\b', number_string )]
                list_of_numbers.add(number[0])

        
#print list_of_numbers

for val in list_of_numbers:

    outfile.write(str(val))
    outfile.write("\n")


outfile.close()
    
