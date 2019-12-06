#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 16:09:56 2019

@author: luca
"""
import os

#CREATE f2py Makefile
Debug = False
if Debug:
    F90Flags = "-g -fbacktrace -ffpe-trap=invalid,zero -ffpe-summary=invalid,zero -fbounds-check"
else:
    F90Flags = "-g -O2"

Program_Dir = 'TwoMoment_Flavor_Relaxation_SQA'
# These are the modules defined in the local Program Directory
MyModules = ['InterpolationModule','ReadProfileModule','InitializationModule',\
             'InputOutputRelaxationModule', 'ImplicitSolverModule', 'IntegrationModule', \
             'OscillationsUtilsModule','OscillationsModule']
f90Compiler = 'mpif90'
PyModuleName = 'Thor2SQA'

PythonWrap_Dir = 'Executables'
Executables_dir = Program_Dir + '/' + 'Executables'
# These are the modules to be "turned" into python modules
InterfaceNames = ['ThornadoSQAInterfaceModule', 'ProgramStartEndModule']

# These are othr directories you need to specify
Thornado_Dir = os.environ['THORNADO_DIR']
WeakLib_Dir = os.environ['WEAKLIB_DIR']
HDF5DIR = '/usr/local/hdf5-1.10.5/hdf5'
HDF5_Inc = HDF5DIR + '/include'
HDF5_Lib = HDF5DIR + '/lib'
LAPACK_Inc = '/usr/include'
LAPACK_Lib = '/usr/lib/x86_64-linux-gnu/lapack'

for InterfaceName in InterfaceNames:
    os.system('cp ' + Program_Dir + '/' + InterfaceName + '.f90 ' + PythonWrap_Dir) 

# Move module files in Modules directory
for file in os.listdir(Executables_dir):
    filename = os.fsdecode(file)
    SkipNextIteration = False

    if filename.endswith(".mod"):
        
        for InterfaceName in InterfaceNames:
            if filename == InterfaceName.lower() + ".mod":
                SkipNextIteration = True
        
        if SkipNextIteration:
            continue
            
        os.system('cp ' + Executables_dir + '/' + filename + ' ' + PythonWrap_Dir)
        
    elif filename.endswith(".o"):
        
        for InterfaceName in InterfaceNames:
            if filename == InterfaceName + ".o":
                SkipNextIteration = True
                continue

        if SkipNextIteration:
            continue
            
        os.system('cp ' + Executables_dir + '/' + filename + ' ' + PythonWrap_Dir)


F90_Modules = []
''' WeakLib ObjectFiles '''
Build_Modules_File = WeakLib_Dir + r'/Distributions/Build/Makefile_WeakLib_ObjectFiles'
    
ObjectFiles_List = []
with open(Build_Modules_File,'r') as f:
    for line in f:
        if line.strip()[0:7] == 'include':
            for s,i in zip(line,range(len(line))):
                if s == '/':
                    ObjectFiles_List.append(line[i+1:-1])
                    break

for Makefile in ObjectFiles_List:
    with open(WeakLib_Dir + '/' + Makefile) as f:
        for line in f:
            if line.strip()[-2:] == '.o':
                F90_Modules.append('-I. ' + line.split()[0])
                
            if line.strip()[-1] == '\\':
                if line.split()[0][-2:] == '.o':
                    F90_Modules.append('-I. ' + line.split()[0])


''' Thornado ObjectFiles '''
Build_Modules_File = Thornado_Dir + r'/Build/Makefile_Thornado_ObjectFiles'

ObjectFiles_List = []
with open(Build_Modules_File,'r') as f:
    for line in f:
        if line.strip()[0:7] == 'include':
            for s,i in zip(line,range(len(line))):
                if s == '/':
                    ObjectFiles_List.append(line[i+1:-1])
                    break

for Makefile in ObjectFiles_List:
    with open(Thornado_Dir + '/' + Makefile) as f:
        for line in f:
            if line.strip()[-2:] == '.o':
                F90_Modules.append('-I. ' + line.split()[0])
                break
            
            if line.strip()[-1] == '\\':
                if line.split()[0][-2:] == '.o':
                    F90_Modules.append('-I. ' + line.split()[0])

''' Program-Specific ObjectFiles '''
for MyModule in MyModules:                    
    F90_Modules.append('-I. ' + MyModule + '.o \\\n')


with open('Executables/Makefile_ObjectFiles','w') as f:
    f.write('F90_MODULES = \\\n')
    
    # Include LAPACK and HDF5
    f.write('\t' + '-I' + LAPACK_Inc + ' \\\n')
    f.write('\t' + '-I' + HDF5_Inc + ' \\\n')
    
    # Include Object Files
    for f90M in F90_Modules:
        f.write('\t' + f90M + ' \\' + '\n')

    # Link LAPACK and HDF5 Libraries  
    f.write('\t' + '-L' + HDF5_Lib + ' -lhdf5_fortran -lhdf5 \\\n')
    f.write('\t' + '-L' + LAPACK_Lib + ' -llapack -lblas\n')


with open('Executables/Makefile','w') as f:
    f.write('include Makefile_ObjectFiles\n\n')
    f.write('default:\n')
    f.write('\t' + 'f2py -c --f90exec=' + f90Compiler + ' $(F90_MODULES) \\\n')
    f.write('\t' +  '--f90flags="-fopenmp ' + F90Flags + '" -lgomp \\\n')
    f.write('\t' + '-m ' + PyModuleName + ' ') 
    for InterfaceName in InterfaceNames:
        f.write(InterfaceName + '.f90 ')    
    
    f.write('\n')
    f.write('clean:\n')
    f.write('\t' + 'rm -f *.so')

#Create dictionary that detects double precision for python
with open(r'Executables/.f2py_f2cmap','w') as f:
    f.write("dict(real=dict(sp='float', dp='double'))\n")

