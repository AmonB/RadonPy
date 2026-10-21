#!/usr/bin/env python3

#  Copyright (c) 2026. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

__version__ = '1.0b2'

import matplotlib

matplotlib.use('Agg')

import pandas as pd
import os
import sys


from radonpy.core import const
# const.mpi_cmd = 'mpiexec -n %i'
# const.check_package_disable = True
const.lammps_exec = '/home/amon/LAMMPS/bin/lmp'
# os.environ['RadonPy_LAMMPS_INTEL'] = 'off'
# os.environ['RadonPy_LAMMPS_OPT'] = 'off'
# os.environ['RadonPy_No_Traj'] = str(os.environ.get('RadonPy_No_Traj', 'True'))

from radonpy.core import utils
from radonpy.sim import helper, lammps
from radonpy.sim.preset import saos

if __name__ == '__main__':
    data = {
        'DBID': os.environ.get('RadonPy_DBID'),
        **helper.get_version(),
        'preset_lt_ver': saos.__version__,
    }

    omp = int(os.environ.get('RadonPy_OMP', 0))
    mpi = int(os.environ.get('RadonPy_MPI', 1))
    gpu = int(os.environ.get('RadonPy_GPU', 1))
    intel = os.environ.get('RadonPy_LAMMPS_INTEL', 'auto')
    opt = os.environ.get('RadonPy_LAMMPS_OPT', 'auto')
    rst_json_file = os.environ.get('RadonPy_JSON_File', None)
    rst_pickle_file = os.environ.get('RadonPy_Pickle_File', None)
    no_traj = bool(os.environ.get('RadonPy_No_Traj', False) == 'True')
    lt_force = bool(os.environ.get('RadonPy_Lt_Force', False) == 'True')

    work_dir = './%s' % data['DBID']
    save_dir = os.path.join(work_dir, 'analyze')
    io = helper.IO_Helper(work_dir, save_dir)

    # Load results.csv or input_data.csv file
    data = io.load_md_csv(data)
    if not data['check_eq'] and not lt_force:
        print('check_eq: FALSE')
        sys.exit(0)

    # Load JSON file, pickle file, or LAMMPS data file
    # mol = io.load_md_obj(rst_json_file=rst_json_file, rst_pickle_file=rst_pickle_file)
    mol = utils.pickle_load(os.path.join(save_dir, 'eq.pickle'))


    # Lt calculation
    ltmd = saos.SAOS(mol, work_dir=work_dir, no_traj=no_traj)
    mol, lt_results = ltmd.exec(temp=data['temp'], amplitude=0.1, period=5e4, cycle=40,
                                mpi=mpi, omp=omp, gpu=gpu, intel=intel, opt=opt)

    # Reload MD csv data
    data = io.load_md_csv(data)

    #lt_results is the dictionary including results of LTMD
    data.update(lt_results)

    # Data output after LTMD
    io.output_md_data(data)
