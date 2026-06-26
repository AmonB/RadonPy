#!/usr/bin/env python3

#  Copyright (c) 2026. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

__version__ = '1.0b2'

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os
import math

# For user installed LAMMPS
# from radonpy.core import const
# const.lammps_exec = '/path/to/lammps/binary'

from radonpy.core import utils
from radonpy.sim import helper
from radonpy.sim.preset import ef_dp


if __name__ == '__main__':
    data = {
        'DBID': os.environ.get('RadonPy_DBID'),
        'temp': float(os.environ.get('RadonPy_Temp', 300.0)),
        'press': float(os.environ.get('RadonPy_Press', 1.0)),
        **helper.get_version(),
        'preset_ef_dp_ver': ef_dp.__version__,
    }

    freq = float(os.environ.get('electric_field_freq', 50.0*1e+9))    # unit: Hz
    axis = os.environ.get('electric_field_axis', 'z')

    tuning_maxe = float(os.environ.get('tuning_max_EField', 1.0))    # unit: V/angstrom
    tuning_rate = float(os.environ.get('tuning_rate', 2.0e-7))    # unit: V/angstrom/fs
    tuning_plot = os.environ.get('tuning_plot', 'False') == 'True'
    wave_num = int(os.environ.get('number_of_waves', 10))
    efield = os.environ.get('electric_field_value', 'auto')
    if efield == "auto":
        evalue = tuning_maxe
    else:
        evalue = [float(i) for i in efield.split(',')]
    ef_ensemble = os.environ.get('ef_ensemble', 'npt')

    omp = int(os.environ.get('RadonPy_OMP', 0))
    mpi = int(os.environ.get('RadonPy_MPI', utils.cpu_count()))
    gpu = int(os.environ.get('RadonPy_GPU', 0))
    intel = os.environ.get('RadonPy_LAMMPS_INTEL', 'auto')
    opt = os.environ.get('RadonPy_LAMMPS_OPT', 'auto')
    rst_pickle_file = os.environ.get('RadonPy_Pickle_File', None)
    dp_force = os.environ.get('RadonPy_DP_Force', 'False') == 'True'
    
    work_dir = './%s' % data['DBID']
    save_dir = os.path.join(work_dir, 'analyze')
    io = helper.IO_Helper(work_dir, save_dir)
    
    # Load results.csv or input_data.csv file
    data = io.load_md_csv(data)

    # Load pickle file or LAMMPS data file
    mol = io.load_md_obj(rst_pickle_file=rst_pickle_file)

    if not data['check_eq'] and not dp_force:
        print('check_eq: FALSE')
        sys.exit(0)


    ## Calculation section
    raw_ef_step = wave_num * 1e+9 / freq
    round_ef_step = round(raw_ef_step, 1-math.floor(math.log10(raw_ef_step)))   # unit = ns
    if freq >= 1e+10:
        simulation_freq = round_ef_step * 1000
    else:
        simulation_freq = round_ef_step * 200

    ef_dp_results = {}
    
    xyz_direction = [direction for direction in axis.split(',')]
    for direction in xyz_direction:
        if efield == "auto":
            tuning_time = evalue / tuning_rate
            round_tuning_time = round(tuning_time, 1-math.floor(math.log10(tuning_time)))
            tuning_step = round_tuning_time / 1000000
            tuning_freq = round_tuning_time / 5000

            # Execute tuning
            efeq = ef_dp.further_Additional(mol, work_dir=work_dir, axis=direction, evalue=evalue, tuning_rate=tuning_rate, process='tuning')
            if os.path.isfile(os.path.join(work_dir, 'ef_%s_tuning_last.data' % direction)):
                utils.radon_print('Tuning data exists.', level=1)
            else:
                efeq.exec_dp(temp=data['temp'], press=data['press'], eq_step=tuning_step,
                                 mpi=mpi, omp=omp, gpu=gpu,
                                 thermo_freq=tuning_freq, dump_freq=tuning_freq,
                                 ensemble_s=ef_ensemble, intel=intel, opt=opt)
            estimated_evalue = round(efeq.tuning(tuning_plot), 2)

            # Execute efield-MD
            efeq = ef_dp.further_Additional(mol, work_dir=work_dir, axis=direction, evalue=estimated_evalue, freq=freq, process='dp')
            if os.path.isfile(os.path.join(work_dir, 'ef_%sGHz_%sEF_%s_last.data' % (round(freq/1e+9, 2), estimated_evalue, direction))):
                utils.radon_print('Calculating permittivity (%s direction, %1.2f V/A)' %(direction, estimated_evalue), level=1)
            else:
                efeq.exec_dp(temp=data['temp'], press=data['press'], eq_step=round_ef_step,
                             mpi=mpi, omp=omp, gpu=gpu,
                             thermo_freq=simulation_freq, dump_freq=simulation_freq,
                             ensemble_s=ef_ensemble, intel=intel, opt=opt)
            results = efeq.dipole_anal(props=data, printout=True, save=True)

            results_dict = results.iloc[0, :].to_dict()
            ef_dp_results.update({
                'ef_dp_efield_%iGHz_%s' % (round(freq/1e+9, 2), direction): estimated_evalue,
                'ef_dp_tuning_max_EField_%iGHz_%s' % (round(freq/1e+9, 2), direction): tuning_maxe,
                'ef_dp_tuning_rate_%iGHz_%s' % (round(freq/1e+9, 2), direction): tuning_rate,
                'ef_dp_nwave_%iGHz_%s' % (round(freq/1e+9, 2), direction): wave_num,
                'ef_dp_D0_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['D0'],
                'ef_dp_delta_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['delta'],
                'ef_dp_dynamic_dielectric_constant_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['real part of permittivity'],
                'ef_dp_dynamic_dielectric_constant_corr_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['(+ refractive_index^2 - 1)'],
                'ef_dp_dielectric_loss_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['imaginary part of permittivity'],
                'ef_dp_loss_tangent_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['loss tangent'],
                'ef_dp_r2_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['r2'],
                'ef_dp_rmspe_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['rmspe'],
                'ef_dp_remarks_%iGHz_%s' % (round(freq/1e+9, 2), direction): results_dict['remarks'],
            })

        # else:
        #     for e in evalue:
        #         efeq = ef_dp.further_Additional(mol, work_dir=work_dir, axis=direction, evalue=e, freq=freq, process='dp')
        #         if os.path.isfile(os.path.join(work_dir, 'ef_%sGHz_%sEF_%s_last.data' % (round(freq/1e+9, 2), e, direction))):
        #             utils.radon_print('Calculating permittivity (%s direction, %1.2f V/A)' %(direction, e), level=1)
        #         else:
        #             # Execute efield-MD
        #             efeq.exec_dp(temp=data['temp'], press=data['press'], eq_step=round_ef_step,
        #                          mpi=mpi, omp=omp, gpu=gpu,
        #                          thermo_freq=simulation_freq, dump_freq=simulation_freq,
        #                          ensemble_s=ef_ensemble, intel=intel, opt=opt)
        #         results = efeq.dipole_anal(props=data, printout=True, save=True)


    # Reload MD csv data
    data = io.load_md_csv(data)

    #tg_results is the dictionary including results of TGMD
    data.update(ef_dp_results)

    # Data output after TGMD
    io.output_md_data(data)

