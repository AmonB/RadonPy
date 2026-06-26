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
export RadonPy_EFDP_Freq=50.0*1e+9    # unit: Hz
# export RadonPy_EFDP_EField_Axis='z'
# export RadonPy_Temp='300'
# export RadonPy_Press='1.0'
# export RadonPy_EFDP_Tuning_Max_EField='1.0'
# export RadonPy_EFDP_Tuning_Rate='2.0e-7'
# export RadonPy_EFDP_Tuning_Plot='False'
# export RadonPy_EFDP_Number_of_Waves='10'
# export RadonPy_EFDP_Ensemble='npt'

## Directory path and operation settings
# export RadonPy_TMP_Dir=${TMP_DIR}
# export RadonPy_JSON_File=''
# export RadonPy_Pickle_File=''
# export RadonPy_DP_Force='False'

python3 /path/to/AutoMD_scripts/6_ef_dp.py
