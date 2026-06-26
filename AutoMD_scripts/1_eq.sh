#!/bin/bash

### Addition of parameters of Job manager

### Loading runtime environment for LAMMPS
# module load mpi
# module load cuda
# module load lammps

### Loading python environment
# conda activate radonpy


####################
# RadonPy settings #
####################
## Set IDs and SMILES
export RadonPy_DBID="PE_1"
export RadonPy_SMILES="*CC*"
export RadonPy_Monomer_ID="PE"
export RadonPy_TER_ID='CH3'
#export RadonPy_TER_ID2=''

## Parallel computing settings for LAMMPS
export RadonPy_OMP=2
export RadonPy_MPI=64
export RadonPy_GPU=0
# export RadonPy_LAMMPS_INTEL='auto'
# export RadonPy_LAMMPS_OPT='auto'

## Settings for construction of initial structure and simulation model
# export RadonPy_NAtom='1000'
# export RadonPy_NChain='10'
# export RadonPy_Tacticity='atactic'
# export RadonPy_Copoly_Ratio='1'
# export RadonPy_Copoly_Type='random'
# export RadonPy_Ini_Density='0.05'
# export RadonPy_FF='GaFF2_mod'

## MD simulation conditions
# export RadonPy_Temp='300'
# export RadonPy_Press='1.0'

## Directory path and operation settings
# export RadonPy_TMP_Dir=${TMP_DIR}
# export RadonPy_Monomer_Dir=/path/to/directory/of/repeating_unit/file
# export RadonPy_Ter_Dir=/path/to/directory/of/terminating_group/file
# export RadonPy_RetryEQ=0

python3 /path/to/AutoMD_scripts/1_eq.py
