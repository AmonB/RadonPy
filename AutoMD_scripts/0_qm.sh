#!/bin/bash

### Addition of parameters of Job manager

### Loading runtime environment for LAMMPS
# module load mpi
# module load cuda
# module load lammps
# export LAMMPS_EXEC=/opt/bin/lammps/lammps

### Loading python environment
# conda activate radonpy


####################
# RadonPy settings #
####################
## Set IDs and SMILES
export RadonPy_DBID="PE_1"
export RadonPy_SMILES="*CC*"
export RadonPy_Monomer_ID="PE"

## QM calculation settings
# export RadonPy_QM_Method='wb97m-d3bj'
# export RadonPy_Charge='RESP'
# export RadonPy_FF='GAFF2_mod'
# export RadonPy_Do_TDDFT='False'
# export RadonPy_QM_TD_Method='cam-b3lyp-d3bj'

## Computational resource settings
export RadonPy_OMP_Psi4='8'
export RadonPy_MEM_Psi4='20000'    # unit: MByte
# export RadonPy_Conf_MM_OMP='1'
# export RadonPy_Conf_MM_MPI='1'
# export RadonPy_Conf_MM_GPU='0'
# export RadonPy_Conf_MM_MP='0'
# export RadonPy_Conf_Psi4_OMP='4'
# export RadonPy_Conf_Psi4_MP='0'

## Terminating group calculation settings
# export RadonPy_Do_Ter='False'
# export RadonPy_SMILES_TER='C*'
# export RadonPy_TER_ID='CH3'
# export RadonPy_SMILES_TER2=''
# export RadonPy_TER_ID2=''

## Directory path
# export RadonPy_TMP_Dir=${TMP_DIR}

python3 /path/to/AutoMD_scripts/0_qm.py
