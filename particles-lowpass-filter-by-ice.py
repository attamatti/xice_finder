#!/usr/bin/python

#filter particle stack by ice score from find xtal-ice

vers = '0.1.2'

import sys
import os
import subprocess

rih_path = 'relion_image_handler'

#-----------------------------------------------------------------------------------#
class Arg(object):
    _registry = []
    def __init__(self, flag, value, req):
        self._registry.append(self)
        self.flag = flag
        self.value = value
        self.req = req

def make_arg(flag, value, req):
    global errmsg
    errmsg = 'USAGE: particles-filter-by-ice.py --parts <particle starfile> --thresh <threshold>\nIf the log file from xice_finder is called somthing different than xice_find.log specify it with --ice <path to logfile>'
    Argument = Arg(flag, value, req)
    if Argument.req == True:
        if Argument.flag not in sys.argv:
            print(errmsg)
            sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
    if Argument.value == True:
        try:
            test = sys.argv[sys.argv.index(Argument.flag)+1]
        except ValueError:
            if Argument.req == True:
                print(errmsg)
                sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
            elif Argument.req == False:
                return False
        except IndexError:
                print(errmsg)
                sys.exit("ERROR: argument '{0}' requires a value".format(Argument.flag))
        else:
            if Argument.value == True:
                Argument.value = sys.argv[sys.argv.index(Argument.flag)+1]
        
    if Argument.value == False:
        if Argument.flag in sys.argv:
            Argument.value = True
        else:
            Argument.value = False
    return Argument.value
#-----------------------------------------------------------------------------------#
###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile(f):
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    for i in alldata:
        if '#' in i:
            labelsdic[i.split('#')[0]] = int(i.split('#')[1])-1
        if len(i.split()) > 3:
            data.append(i.split())
        if len(i.split()) < 3:
            header.append(i.strip("\n"))
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#


    
#read icefinderlog put in dict
xif_log= make_arg('--ice',True,False)
if xif_log == False:
    xif_log = 'xice_find.log'
#check for proper log file:
if os.path.isfile(xif_log) == False:
    xif_log = raw_input("ERROR: Can't find xice_find.log.  Specify its location: ")
logdata = open(xif_log,'r').readlines()
ice_scores = {}
for i in logdata:
    line = i.split()
    ice_scores[line[0]] = float(line[1])

thresh = float(make_arg('--thresh',True,True))

#read particles
starfile = make_arg('--parts',True,True)
(labels,header,data) = read_starfile(starfile)
goodparts = []
count = 0
print('Filtering particles with score > {0} to 3.7 Angstrom'.format(thresh))
for i in data:
    source  = i[labels['_rlnImageName ']].split('@')[-1].split('/')[-1]
    if float(ice_scores[source]) <= thresh:
        goodparts.append(i)
    else:
    	newline = []
        imgname = i[labels['_rlnImageName ']].split('@')[-1]
	number = i[labels['_rlnImageName ']].split('@')[0]
        newpartimg = '{0}@{1}_filt.mrcs'.format(number,imgname.split('.')[0])
	if os.path.isfile('{0}_filt.mrcs'.format(imgname.split('.')[0])) != True:
            subprocess.call('{0} --i {1} --o filt --lowpass 3.7'.format(rih_path,imgname),shell=True)
        i[labels['_rlnImageName ']] = newpartimg
        goodparts.append(i)
	count +=1

output = open('icefiltered_{0}.star'.format(str(thresh).replace('.','p')),'w')
for i in header:
    output.write('{0}\n'.format(i))
for i in goodparts:
    output.write('{0}\n'.format('\t'.join(i)))
print('{0} Iniial Particles'.format(len(data)))
print('{0} Filtered Particles'.format(count))
print('Wrote output to: icefiltered.star')
    
