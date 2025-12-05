#!/usr/bin/env python3

#  Copyright (c) 2025. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found inthe LICENSE file.

__version__ = '0.1.0'

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os

# Preset module for solubility parameter calculation require to install "TALLY" package in LAMMPS
# For user installed LAMMPS
# from radonpy.core import const
# const.lammps_exec = '/path/to/lammps/binary'

from radonpy.core import utils
from radonpy.sim import helper
from radonpy.sim.preset import sp


if __name__ == '__main__':
    data = {
        'DBID': os.environ.get('RadonPy_DBID'),
        **helper.get_version(),
        'preset_sp_ver': sp.__version__,
    }
    
    omp = int(os.environ.get('RadonPy_OMP', 0))
    mpi = 16 # mpi should always be 16, 64 or 128 for rerun.
    #mpi = int(os.environ.get('RadonPy_MPI', utils.cpu_count()))
    gpu = int(os.environ.get('RadonPy_GPU', 0))
    rst_pickle_file = os.environ.get('RadonPy_Pickle_File', None)
    rst_data_file = os.environ.get('RadonPy_LAMMPS_Data_File', None)
    
    work_dir = './%s' % data['DBID']
    save_dir = os.path.join(work_dir, 'analyze')
    io = helper.IO_Helper(work_dir, save_dir)
    
    # Load results.csv or input_data.csv file
    data = io.load_md_csv(data)

    # Load JSON file, pickle file, or LAMMPS data file
    mol = io.load_md_obj(rst_pickle_file=rst_pickle_file)

    # SP calculation (using rerun)
    spmd = sp.Rerun(mol, work_dir=work_dir)
    mol, sp_results = spmd.exec(temp=data['temp'], press=data['press'], mpi=mpi, omp=omp, gpu=gpu)    

    # Reload MD csv data
    data = io.load_md_csv(data)

    #tg_results is the dictionary including results of SPMD
    data.update(sp_results)

    # Data output after SPMD
    io.output_md_data(data)

