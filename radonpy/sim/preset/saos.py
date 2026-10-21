#!/usr/bin/env python3

#  Copyright (c) 2026. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# ******************************************************************************
# sim.preset.saos module
# ******************************************************************************
import os
import glob
import re
import datetime
import numpy as np
import pandas as pd
import math
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from functools import partial
from rdkit import Geometry as Geom
from ...core import calc, const, utils
from .. import lammps, preset
from ..md import MD
import shutil

from matplotlib import pyplot as plt

__version__ = '0.1.0'

### Plot global settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
# plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['mathtext.rm'] = 'Times New Roman'


### backup log file before run
def backup_log_file(work_dir=None, log_file=None):
    _ = os.path.join(work_dir, log_file)
    if not os.path.exists(_):
        return None

    counter = 1
    while True:
        backup_path = f"{_}.{counter}"
        if os.path.exists(backup_path):
            counter += 1
        else:
            shutil.move(_, backup_path)
            break


class EQMD(preset.Preset):
    def __init__(self, mol, prefix='', work_dir=None, save_dir=None, solver_path=None, no_traj=True, **kwargs):
        """
        preset.eq.Equilibration

        Base class of equilibration preset

        Args:
            mol: RDKit Mol object
        """
        super().__init__(mol, prefix=prefix, work_dir=work_dir, save_dir=save_dir, solver_path=solver_path, **kwargs)
        self.data = {}
        self.curdir = os.getcwd()
        self.work_dir_lt = os.path.join(self.work_dir, 'md_lt')
        if not os.path.exists(self.work_dir_lt):
            os.mkdir(self.work_dir_lt)

        self.pre_in_file = kwargs.get('pre_in_file', '%slt_pre.in' % self.prefix)
        self.pre_log_file = kwargs.get('pre_log_file', '%slt_pre.log' % self.prefix)
        self.pre_dat_file = kwargs.get('pre_dat_file', '%slt_pre.data' % self.prefix)
        self.pre_pickle_file = kwargs.get('pre_pickle_file', '%slt_pre.pickle' % self.prefix)
        self.pre_json_file = kwargs.get('pre_json_file', '%slt_pre.json' % self.prefix)

        self.eq_in_file = kwargs.get('eq_in_file', '%slt_eq.in' % self.prefix)
        self.eq_log_file = kwargs.get('eq_log_file', '%slt_eq.log' % self.prefix)
        self.eq_dat_file = kwargs.get('eq_dat_file', '%slt_eq.data' % self.prefix)
        self.eq_pickle_file = kwargs.get('eq_pickle_file', '%slt_eq.pickle' % self.prefix)
        self.eq_json_file = kwargs.get('eq_json_file', '%slt_eq.json' % self.prefix)
        self.eq_restart_file = kwargs.get('eq_restart_file', '%slt_eq.rst' % self.prefix)

        self.csv_file = kwargs.get('csv_file', '%slt.csv' % self.prefix)
        self.in_file = kwargs.get('in_file', '%slt.in' % self.prefix)
        self.log_file = kwargs.get('log_file', '%slt.log' % self.prefix)
        self.dat_file = kwargs.get('dat_file', '%slt.data' % self.prefix)
        self.pickle_file = kwargs.get('pickle_file', '%slt.pickle' % self.prefix)
        self.json_file = kwargs.get('json_file', '%slt.json' % self.prefix)
        if no_traj:
            self.dump_file = kwargs.get('dump_file')
            self.xtc_file = kwargs.get('xtc_file')
        else:
            self.dump_file = kwargs.get('dump_file', '%slt.dump' % self.prefix)
            self.xtc_file = kwargs.get('xtc_file', '%slt.xtc' % self.prefix)
        self.restart_file = kwargs.get('restart_file', '%slt.rst' % self.prefix)
        self.rg_file = kwargs.get('rg_file', '%slt.profile' % self.prefix)

    def equilibration(self, temp, press=1.0, time_step=1.0, step=5e5, p_dump=1000, **kwargs):
        md = MD()
        md.pair_style = self.pair_style
        md.cutoff_in = self.cutoff_in
        md.cutoff_out = self.cutoff_out
        md.kspace_style = self.kspace_style
        md.kspace_style_accuracy = self.kspace_style_accuracy
        md.bond_style = self.bond_style
        md.angle_style = self.angle_style
        md.dihedral_style = self.dihedral_style
        md.improper_style = self.improper_style
        md.neighbor = '%s bin' % self.neighbor_dis
        #        md.dump_freq = kwargs.get('dump_freq',   1e10)
        md.dump_file = None
        md.xtc_file = None
        md.thermo_freq = 100 if int(step / 1000) < 100 else int(step / 1000) if int(step / 1000) < 1000 else 1000
        md.log_file = kwargs.get('eq_log_file', self.eq_log_file)
        md.dat_file = kwargs.get('pre_dat_file', self.pre_dat_file)
        if kwargs.get('set_init_velocity', False):
            md.set_init_velocity = temp
        md.add_md('npt', int(step), time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                  p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('eq_dat_file', self.eq_dat_file)

        return md

    def saos(self, temp, axis='xy', amplitude='0.1', period='5e4', cycle=40, time_step=1.0, p_dump=1000, **kwargs):
        md = MD()
        step = int(period *cycle)
        md.pair_style = self.pair_style
        md.cutoff_in = self.cutoff_in
        md.cutoff_out = self.cutoff_out
        md.kspace_style = self.kspace_style
        md.kspace_style_accuracy = self.kspace_style_accuracy
        md.bond_style = self.bond_style
        md.angle_style = self.angle_style
        md.dihedral_style = self.dihedral_style
        md.improper_style = self.improper_style
        md.neighbor = '%s bin' % self.neighbor_dis
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.thermo_freq = 100 if int(step / 1000) < 100 else int(step / 1000) if int(step / 1000) < 1000 else 1000
        md.log_file = kwargs.get('log_file', self.log_file)
        md.dat_file = kwargs.get('eq_dat_file', self.eq_dat_file)
        md.add_md('nvt', step, time_step=time_step, shake=False, t_start=temp, t_stop=temp, p_dump=p_dump, **kwargs)
        md.wf[-1].add_deform(dftype='wiggle', axis=axis, amplitude=amplitude, period=period, **kwargs)
        md.wf[-1].add_timeave_wiggle()
        md.write_data = kwargs.get('dat_file', self.dat_file)

        return md

    def analyze(self, ignore_log=[], **kwargs):

        analy = SAOS_analyze(
            log_file  = os.path.join(self.work_dir_lt, kwargs.get('log_file', self.log_file)),
            ignore_log = ignore_log,
            **kwargs
        )

        return analy


class SAOS_analyze(lammps.Analyze):
    '''
    SASO_analyze class reads and analyzes dump/log files.
    '''
    def __init__(self, log_file=None, save_dir=None, ignore_log=[], **kwargs):
        super().__init__(log_file=log_file, ignore_log=ignore_log, **kwargs)
        self.save_dir = save_dir


    def saos(self, axis, Tp, last=10):
        def sine_curve(t, A, phi, C, omega):
            return A * np.sin(omega * t + phi) + C

        def sine_fitting(t, y, omega):
            M = np.column_stack([np.sin(omega*t), np.cos(omega*t), np.ones_like(t)])
            (a, b, C), *_ = np.linalg.lstsq(M, y, rcond=None)
            A, phi = np.hypot(a, b), np.arctan2(b, a)
            return A, phi, C

        def get_upper(x):
            order = math.floor(math.log10(abs(x)))
            scale = 10 ** order
            return math.ceil(x * 1.2 / scale) * scale


        df = self.read_log(self.log_file)[-1]

        # take last N cycles for curve fitting
        df_last = df.iloc[-int(last*Tp):].loc[:, ['Time', 'v_gamma', f'v_avg_{axis}']].reset_index(drop=True)
        df_last['Time'] = (df_last['Time'] - df_last['Time'].iloc[0]) / 1e3   # ps
        t = df_last['Time'].to_numpy()
        gamma = df['v_gamma'].to_numpy()
        gamma_0 = gamma.max()
        df_last[f'v_avg_{axis}'] = df_last[f'v_avg_{axis}'] * const.atm2pa / 1e6    # Mpa
        sigma_raw = df_last[f'v_avg_{axis}']

        # save process data
        df_last.to_csv(os.path.join(self.save_dir, f'os_last_{last}.csv'), index=False)

        omega = 2 * np.pi / Tp
        sigma_0, delta, C = sine_fitting(t, sigma_raw, omega)
        sigma = sine_curve(t, sigma_0, delta, C, omega)


        fig, ax1 = plt.subplots(figsize=(6, 6))
        ax2 = ax1.twinx()
        # stress
        ax1.set_xlabel('Time [ps]', fontsize=12)
        ax1.set_ylabel('Stress [MPa]', color='royalblue', fontsize=12)
        ax1.tick_params(axis='y', labelcolor='royalblue')
        ax2.spines['left'].set_color('royalblue')
        sigma_upper = get_upper(sigma_0)
        ax1.set_ylim(-sigma_upper, sigma_upper)
        line1 = ax1.scatter(t, sigma_raw, color='#7F8C8D', alpha=0.25, label='Raw Stress Data')
        line2, = ax1.plot(t, sigma, color='royalblue', linewidth=1.5, label='Fitted Stress')

        # strain

        ax2.set_ylabel('Strain', color='#179148', fontsize=12)
        ax2.tick_params(axis='y', labelcolor='#179148')
        ax2.spines['right'].set_color('#179148')
        strain_upper = get_upper(gamma_0)
        ax2.set_ylim(-strain_upper, strain_upper)
        line3, = ax2.plot(t, gamma, color='#179148', linewidth=1.5, label='Strain')

        lines = [line1, line2, line3]
        ax1.legend(lines, [l.get_label() for l in lines], loc='lower right', frameon=True,
                   prop={'family': 'Arial', 'size': 9})

        plt.tight_layout()
        fig.savefig(os.path.join(self.save_dir, f'Oscillatory_shear_{axis}.png'))
        plt.close(fig)

        data = {}
        data['sigma_0'] = sigma_0
        data['delta'] = delta
        data['gamma_0'] = gamma_0
        data['storage_modulus'] = sigma_0 * np.cos(delta) / gamma_0
        data['loss_modulus'] = sigma_0 * np.sin(delta) / gamma_0
        data['loss_tan'] = np.tan(delta)

        if data['loss_modulus'] < 0:
            return None
        else:
            return data


class SAOS(EQMD):
    '''
    SASO class
    '''
    def exec(self, temp=300.0, time_step=1.0, amplitude=0.1, period=5e4, confId=0, axis='xy',
             cycle=40, omp=1, mpi=1, gpu=0, intel='auto', opt='auto', **kwargs):
        lmp = lammps.LAMMPS(work_dir=self.work_dir_lt, solver_path=self.solver_path)


        # Equilibration algorithm for loss modulus
        if os.path.exists(os.path.join(self.save_dir, self.eq_pickle_file)):
            self.mol = utils.pickle_load(os.path.join(self.save_dir, self.eq_pickle_file))
            self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
            utils.radon_print("Equilibration (lt_eq) has already completed.", level=1)
        else:
            dt1 = datetime.datetime.now()
            utils.radon_print('Equilibration (lt_eq) by LAMMPS is running...', level=1)
            md1 = self.equilibration(temp=temp, step=5e6, **kwargs)
            lmp.make_input(md1, file_name=self.eq_in_file)
            backup_log_file(self.work_dir_lt, self.eq_log_file)
            self.mol = lmp.run(md1, mol=self.mol, confId=confId, input_file=self.eq_in_file, omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            dt2 = datetime.datetime.now()
            utils.radon_print('Complete Equilibration (lt_eq). Elapsed time = %s' % str(dt2-dt1), level=1)
            utils.MolToJSON(self.mol, os.path.join(self.save_dir, self.eq_json_file))
            utils.pickle_dump(self.mol, os.path.join(self.save_dir, self.eq_pickle_file))

        # Small amplitude oscillation shear algorithm for loss modulus
        if os.path.exists(os.path.join(self.save_dir, self.pickle_file)):
            self.mol = utils.pickle_load(os.path.join(self.save_dir, self.pickle_file))
            self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
            utils.radon_print("SAOS (lt_saos) has already completed.", level=1)
        else:
            dt1 = datetime.datetime.now()
            utils.radon_print('SAOS (lt_saos) by LAMMPS is running...', level=1)
            md2 = self.saos(temp=temp, axis=axis, amplitude=amplitude, period=period, cycle=cycle,  **kwargs)
            lmp.make_input(md2, file_name=self.in_file)
            backup_log_file(self.work_dir_lt, self.in_file)
            self.mol = lmp.run(md2, mol=self.mol, confId=confId, input_file=self.in_file, omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            dt2 = datetime.datetime.now()
            utils.radon_print('Complete SASO (lt_saos). Elapsed time = %s' % str(dt2-dt1), level=1)
            utils.MolToJSON(self.mol, os.path.join(self.save_dir, self.json_file))
            utils.pickle_dump(self.mol, os.path.join(self.save_dir, self.pickle_file))

        result = SAOS_analyze(log_file=os.path.join(self.work_dir_lt, self.log_file),
                              save_dir=self.save_dir).saos(axis=axis, Tp=period, last=10)

        if result is None:
            utils.radon_print('Loss tangent could not be calculated in SAOS_analyze.', level=3)
        else:
            self.data.update(result)
            pd.DataFrame(self.data, index=[0]).to_csv(os.path.join(self.save_dir, self.csv_file))
            utils.radon_print(f"Storage modulus G' = {result['storage_modulus']}\n "
                              f"Loss modulus G'' {result['loss_modulus']} = \n "
                              f"Loss tangent tan(δ) = {result['loss_tan']}", level=1)

        return self.mol, self.data

