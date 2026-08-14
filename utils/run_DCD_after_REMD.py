#!/usr/bin/env python

import subprocess
import os.path
import toml
import copy
import math

path_to_toml = '../input_REMD.toml'
out_dir = './'
relative_dir = './'
path_to_dcd = '../label001/'
path_to_sis = '../../unisis/unisis'

remd_toml = toml.load(path_to_toml)

try:
    nrep_temp = remd_toml['Replica']['nrep_temp']
except KeyError:
    nrep_temp = 1
try:
    nrep_force = remd_toml['Replica']['nrep_force']
except KeyError:
    nrep_force = 1
try:
    nrep_ion = remd_toml['Replica']['nrep_ion']
except KeyError:
    nrep_ion = 1

if nrep_temp < 1:
    nrep_temp = 1
if nrep_force < 1:
    nrep_force = 1
if nrep_ion < 1:
    nrep_ion = 1

print('nrep_temp = ', nrep_temp)
print('nrep_force = ', nrep_force)
print('nrep_ion = ', nrep_ion)
nrep = nrep_temp * nrep_force * nrep_ion
print('nrep = ', nrep)

irep = 0
for itemp in range(1, nrep_temp+1):
    for iforce in range(1, nrep_force+1):
        for iion in range(1, nrep_ion+1):
            irep += 1

            out_toml = copy.deepcopy(remd_toml)

            out_toml['Job']['type'] = 'DCD'
            out_toml['Files']['In']['dcd'] = path_to_dcd + remd_toml['Files']['Out']['prefix'] + f'_{irep:04d}.dcd'

            out_toml['Files']['In']['ff'] = relative_dir + remd_toml['Files']['In']['ff']
            out_toml['Files']['In']['fasta'] = relative_dir + remd_toml['Files']['In']['fasta']
            if 'pdb_ini' in remd_toml['Files']['In']:
                out_toml['Files']['In']['pdb_ini'] = relative_dir + remd_toml['Files']['In']['pdb_ini']
            if 'xyz_ini' in remd_toml['Files']['In']:
                out_toml['Files']['In']['xzy_ini'] = relative_dir + remd_toml['Files']['In']['xyz_ini']

            out_toml['Files']['Out']['prefix'] = out_dir + remd_toml['Files']['Out']['prefix'] + f'_{irep:04d}'
            out_toml['Files']['Out']['types'] = ['bp']

            if nrep_temp > 1:
                out_toml['Condition']['tempK'] = remd_toml['Replica']['Temperature'][f'{itemp}']
            if nrep_ion > 1:
                out_toml['Electrostatic']['Ionic_strength'] = remd_toml['Replica']['Ionic_strength'][f'{iion}']
            if nrep_force > 1:
                vector = remd_toml['Tweezers']['Dual_Constant_Force']['forces_pN'][0]
                d = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
                k = remd_toml['Replica']['Force'][f'{iforce}']
                out_toml['Tweezers']['Dual_Constant_Force']['forces_pN'] = [[k*vector[0]/d, k*vector[1]/d, k*vector[2]/d],]

            ## Remove 'Replica' field (and the children)
            out_toml.pop('Replica', None)


            inpfile = out_dir + f'/input_{irep:04d}.toml'
            f_inp = open(inpfile, 'w')
            toml.dump(out_toml, f_inp)
            f_inp.close()

            flog = open(out_dir + f'/{irep:04d}.log', 'w')
            ferr = open(out_dir + f'/{irep:04d}.err', 'w')
            subprocess.run([path_to_sis, inpfile], stdout=flog, stderr=ferr)
            flog.close()
            ferr.close()
