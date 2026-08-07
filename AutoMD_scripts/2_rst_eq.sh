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
## Set ID
export RadonPy_DBID="PE_1"

## Parallel computing settings for LAMMPS
export RadonPy_OMP=2
export RadonPy_MPI=64
export RadonPy_GPU=0
# export RadonPy_LAMMPS_INTEL='auto'
# export RadonPy_LAMMPS_OPT='auto'

## MD simulation conditions
# export RadonPy_Temp='300'
# export RadonPy_Press='1.0'

## Directory path and operation settings
# export RadonPy_TMP_Dir=${TMP_DIR}
# export RadonPy_JSON_File=''
# export RadonPy_Pickle_File=''
export RadonPy_RetryEQ='2'
# export RadonPy_Del_Traj='True'

python3 /path/to/AutoMD_scripts/2_rst_eq.py
