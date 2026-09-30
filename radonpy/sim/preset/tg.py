#  Copyright (c) 2026. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.
#  Author: Hiroki Sugisawa @ Mitsubishi Chemical Corp.
# ******************************************************************************
# sim.preset.tg module
# ******************************************************************************

import os
import glob
import re
import datetime
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from rdkit import Geometry as Geom
from ...core import calc, const, utils
from .. import lammps, preset
from ..md import MD

from matplotlib import pyplot as plt

__version__ = '0.1.5'

### Plot global settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
# plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['mathtext.rm'] = 'Times New Roman'

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
        self.work_dir_tg = os.path.join(self.work_dir, 'md_tg')
        if not os.path.exists(self.work_dir_tg):
            os.mkdir(self.work_dir_tg)
      
        self.pre_in_file     = kwargs.get('pre_in_file',    '%stg_pre.in' % self.prefix)
        self.pre_log_file    = kwargs.get('pre_log_file',   '%stg_pre.log' % self.prefix)
        self.pre_dat_file    = kwargs.get('pre_dat_file',   '%stg_pre.data' % self.prefix)
        self.pre_pickle_file = kwargs.get('pre_pickle_file','%stg_pre.pickle' % self.prefix)
        self.pre_json_file   = kwargs.get('pre_json_file',  '%stg_pre.json' % self.prefix)
        
        self.eq_in_file     = kwargs.get('eq_in_file',     '%stg_eq.in' % self.prefix)
        self.eq_log_file    = kwargs.get('eq_log_file',    '%stg_eq.log' % self.prefix)
        self.eq_dat_file    = kwargs.get('eq_dat_file',    '%stg_eq.data' % self.prefix)
        self.eq_pickle_file = kwargs.get('eq_pickle_file', '%stg_eq.pickle' % self.prefix)
        self.eq_json_file   = kwargs.get('eq_json_file',   '%stg_eq.json' % self.prefix)
        self.eq_restart_file= kwargs.get('eq_restart_file','%stg_eq.rst' % self.prefix)
        
        self.csv_file    = kwargs.get('csv_file',    '%stg.csv' % self.prefix)
        self.in_file     = kwargs.get('in_file',     '%stg.in' % self.prefix)
        self.log_file    = kwargs.get('log_file',    '%stg.log' % self.prefix)
        self.dat_file    = kwargs.get('dat_file',    '%stg.data' % self.prefix)
        self.pickle_file = kwargs.get('pickle_file', '%stg.pickle' % self.prefix)
        self.json_file   = kwargs.get('json_file',   '%stg.json' % self.prefix)
        if no_traj:
            self.dump_file   = kwargs.get('dump_file')
            self.xtc_file    = kwargs.get('xtc_file')
        else:
            self.dump_file   = kwargs.get('dump_file',   '%stg.dump' % self.prefix)
            self.xtc_file    = kwargs.get('xtc_file',    '%stg.xtc' % self.prefix)
        self.restart_file= kwargs.get('restart_file','%stg.rst' % self.prefix)    
        self.rg_file     = kwargs.get('rg_file',     '%stg.profile' % self.prefix)
        
    def density_checker(self, t_start, t_stop, press=1.0, time_step=1.0, step=1e5, p_dump=1000, **kwargs):
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
        md.thermo_freq = 100 if int(step/1000)<100 else int(step/1000) if int(step/1000)<1000 else 1000
        md.dat_file = kwargs.get('pre_dat_file',   self.pre_dat_file)
        md.log_file = kwargs.get('pre_log_file',   self.pre_log_file)
        md.add_md('npt', step, time_step=time_step, shake=True, t_start=t_start, t_stop=t_stop,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('pre_dat_file', self.pre_dat_file)
        return md
    
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
        md.thermo_freq = 100 if int(step/1000)<100 else int(step/1000) if int(step/1000)<1000 else 1000
        md.log_file = kwargs.get('eq_log_file',   self.eq_log_file)
        md.dat_file = kwargs.get('pre_dat_file',   self.pre_dat_file)
        if kwargs.get('set_init_velocity', False):
            md.set_init_velocity = temp
        md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('eq_dat_file', self.eq_dat_file)

        return md

    def step_wise(self, max_temp, min_temp, press=1.0, time_step=1.0, p_dump=1000, cooling_rate=1e3, interval_temp=10, **kwargs):
        '''
        max_temp = 1200
        min_temp = 50
        cooling_rate = 1e3 # 10 K/ cooling_rate[step]*time_step[fs]
        '''
        step = int(cooling_rate/2)
        
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
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.thermo_freq = 100 if int(step/1000)<100 else int(step/1000) if int(step/1000)<1000 else 1000
        md.dump_freq = md.thermo_freq
        md.log_file = kwargs.get('log_file',    self.log_file)
        md.dat_file = kwargs.get('eq_dat_file', self.eq_dat_file)
        
        temp_range = np.arange(max_temp, min_temp, -interval_temp)
        n_sample = np.ones(len(temp_range)) *  cooling_rate
        for temp, n_sample in zip(temp_range, n_sample):
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp-interval_temp,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.add_md('npt', step, time_step=time_step, shake=True, t_start=min_temp, t_stop=min_temp,
            p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('dat_file', self.dat_file)
        
        return md
    
    def quench(self, max_temp, min_temp, press=1.0, time_step=1.0, p_dump=1000, cooling_rate=1e3, interval_temp=10, **kwargs):
        '''
        max_temp = 1200
        min_temp = 50
        cooling_rate = 1e3 # 10 K/ cooling_rate[step]*time_step[fs]
        '''
        step = int(cooling_rate/2)
        
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
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.thermo_freq = 100 if int(step/1000)<100 else int(step/1000) if int(step/1000)<1000 else 1000
        md.log_file = kwargs.get('log_file',    self.log_file)
        md.dat_file = kwargs.get('eq_dat_file', self.eq_dat_file)
        
        temp_range = np.arange(max_temp, min_temp, -interval_temp)
        n_sample = np.ones(len(temp_range)) *  cooling_rate
        for temp, n_sample in zip(temp_range, n_sample):
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp-interval_temp, t_stop=temp-interval_temp,
                   p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.add_md('npt', step, time_step=time_step, shake=True, t_start=min_temp, t_stop=min_temp,
            p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('dat_file', self.dat_file)
        
        return md
    
    def continuous(self, max_temp, min_temp, press=1.0, time_step=1.0, p_dump=1000, cooling_rate=1e3, interval_temp=10, **kwargs):
        
        step = int( cooling_rate*(max_temp - min_temp)/interval_temp )
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
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.thermo_freq = 100 if int(step/1000)<100 else int(step/1000) if int(step/1000)<1000 else 1000
        md.log_file = kwargs.get('log_file',    self.log_file)
        md.dat_file = kwargs.get('eq_dat_file', self.eq_dat_file)
        md.add_md('npt', step, time_step=time_step, shake=True, t_start=max_temp, t_stop=min_temp,
            p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
        md.write_data = kwargs.get('dat_file', self.dat_file)
        return md


    def analyze(self, ignore_log=[], **kwargs):

        analy = TGMD_analyze(
            log_file  = os.path.join(self.work_dir_tg, kwargs.get('log_file', self.log_file)),
            ignore_log = ignore_log,
            **kwargs
        )

        return analy


    def step_wise_new(self, temp, press=1.0, time_step=1.0, p_dump=1000, cooling_rate=1e3, interval_temp=10,
                      idx=0, is_last=False, **kwargs):
        step = int(cooling_rate / 2)
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
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.thermo_freq = 100 if int(step / 1000) < 100 else int(step / 1000) if int(step / 1000) < 1000 else 1000
        md.dump_freq = md.thermo_freq
        md.log_file = 'Tg_conf_T%04d.log' % temp
        if idx == 0:
            md.dat_file = kwargs.get('eq_dat_file', self.eq_dat_file)
        else:
            md.dat_file = 'Tg_conf_T%04d.data' % temp
        backup_log_file(self.work_dir_tg, md.log_file)

        if not self.dump_file:
            md.dump_file = 'Tg_conf_T%04d.dump' % temp
        else:
            md.dump_file = None

        if not self.xtc_file:
            md.xtc_file = 'Tg_conf_T%04d.xtc' % temp
        else:
            md.xtc_file = None

        if is_last:
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                       p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
            md.write_data = 'Tg_conf_T%04d_eq.data' % temp
        else:
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp,
                       p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
            md.write_data = 'Tg_conf_T%04d_eq.data' % temp
            md.add_md('npt', step, time_step=time_step, shake=True, t_start=temp, t_stop=temp - interval_temp,
                       p_start=press, p_stop=press, p_dump=p_dump, **kwargs)
            md.write_data = 'Tg_conf_T%04d.data' % (temp - interval_temp)
        return md

class TGMD_analyze(lammps.Analyze):
    '''
    TGMD_analyze class reads and analyzes dump/log files. 
    '''
    def __init__(self, log_file=None, save_dir=None, ignore_log=[], **kwargs):
        super().__init__(log_file=log_file, ignore_log=ignore_log, **kwargs)
        self.save_dir = save_dir
    
    def density_checker(self, index):
        def function_fit(x, a, b):
            if index == 0:
                return -np.abs(a)*np.power(np.abs(x), 3) + b
            else:
                return -np.abs(a)*np.abs(x) + b
    
        def function_fit_rev(y, a, b):
            if index == 0:
                return np.power( -(y-b)/np.abs(a), 1/3)
            else:
                return -(y-b)/np.abs(a)
        
        df = self.read_log(self.log_file)[-1]
        df_x = df['Temp'].to_numpy()
        df_y = df['Density'].to_numpy()
        param, cov = curve_fit(function_fit, df_x, df_y)
        
        data = {}
        data['tg_next_temp'] = int(function_fit_rev(0.45, param[0], param[1])/50+0.5)*50
        data['tg_init_density'] = function_fit(np.max(df_x), param[0], param[1])
        if data['tg_init_density'] < 0.5:
            data['tg_init_density_check'] = True
        elif data['tg_next_temp'] <= df_x[-1]:
            data['tg_init_density_check'] = True
        else:
            data['tg_init_density_check'] = False
        
        return data
    
    def step_wise(self, ):
        #prop_list = ['Temp', 'E_vdwl', 'Volume', 'Density']
        prop_list = ['Temp', 'Density']
        def arg_get(index):
            return [True if p == index else False for p in prop_list]
        
        def linregression(X, Y):
            X = np.array([ [x, 1] for x in X[0]])
            return np.linalg.lstsq(X, Y[0], rcond=None)[0]

        df = self.read_log(self.log_file)
        df_mean = pd.DataFrame( [d[prop_list].mean() for d in df[0::2]] ).to_numpy().T
        df_index = df_mean[arg_get('Density')] <= 0.45
        df_arg   = np.argsort(df_mean[arg_get('Temp'),df_index[0]])[0]
        df_mean = df_mean[:,df_arg:]
        ndim, ndata = df_mean.shape
        
        temp_max = np.max( df_mean[arg_get('Temp')][0] )
        temp_min = np.min( df_mean[arg_get('Temp')][0] )

        nskip_min = 5   # skip 50 K
        nskip_max = 10  # skip 100 K
        flag   = 0.02
        spacer = 0
        nignor = 0
        while nignor<10:
        #while spacer<ndata/4:
            rmse   = []
            tgs    = []
            slopes = []
            diff   = np.zeros( ndata )
            for i in range(nskip_max+spacer+nignor,ndata-nskip_min-nignor):
                temp_1 = df_mean[arg_get('Temp'),   spacer:i-nignor]
                dens_1 = df_mean[arg_get('Density'),spacer:i-nignor]
                a1, b1 = linregression( temp_1, dens_1)
                diff[spacer:i-nignor] = (a1*temp_1+b1) - (dens_1[0])
            
                temp_2 = df_mean[arg_get('Temp'),   i+nignor:]
                dens_2 = df_mean[arg_get('Density'),i+nignor:]
                a2, b2 = linregression( temp_2, dens_2 )
                diff[i+nignor:] = ((a2*temp_2+b2) - (dens_2[0]))
            
                slopes.append( [a1, b1, a2, b2] )
                rmse.append( np.sum( np.square( diff ) ) )
                tgs.append(  (b2-b1)/(a1-a2) )
            if rmse == []:
                return None
            elif ndata < 20:
                return None
            else:
                rmse_arg = np.argmin(np.array( rmse ))
                tg = tgs[rmse_arg]
                if rmse[rmse_arg] > flag:
                    #spacer += 5
                    nignor += 1
                else:
                    break
        
        if tg>temp_max or tg > 100: 
            plt.xlabel("Temperature [K]")
            plt.ylabel(r"Density [$\mathrm{g/cm^3}$]")
            plt.xlim( temp_min, temp_max )
            x1 = np.linspace(temp_max, tg-50, 20)
            y1 = slopes[rmse_arg][0]*x1 + slopes[rmse_arg][1]
            plt.plot(x1, y1, color='gray')
            x2 = np.linspace(temp_min, tg+50, 20)
            y2 = slopes[rmse_arg][2]*x2 + slopes[rmse_arg][3]
            plt.plot(x2, y2, color='gray')
            plt.scatter(df_mean[0], df_mean[1], color='black')
            plt.axvline(x=tg, ymin=0, ymax=2.0, linestyle="dashed", color='red', label='tg=%f'%(tg))
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(self.save_dir, 'tg.png'))
            
            data = {}
            data['tg']      = tg
            data['tg_rmse'] = rmse[rmse_arg]
            data['tg_thermal_expansion_coef(upper_tg)'] = slopes[rmse_arg][0]
            data['tg_thermal_expansion_intercept(upper_tg)'] = slopes[rmse_arg][1]
            data['tg_thermal_expansion_coef(below_tg)'] = slopes[rmse_arg][2]
            data['tg_thermal_expansion_intercept(below_tg)'] = slopes[rmse_arg][3]
            return data
        else:
            return None
        
    def quench():
        print('hello')
    def continuous():
        print('hello')

        
class TGMD(EQMD):
    '''
    TGMD class 
    '''
    def exec(self, pre_temp_start, pre_temp_stop=800, confId=0, min_temp=50.0, time_step=1.0, interval_temp=10,
             cooling_rate=1e3, eq_step=0.5, omp=1, mpi=1, gpu=0, intel='auto', opt='auto', **kwargs):
        lmp = lammps.LAMMPS(work_dir=self.work_dir_tg, solver_path=self.solver_path)
    
        #
        # Density check algorithms
        #

        # CHECK HAS DONE
        if os.path.exists(os.path.join(self.save_dir, self.pre_pickle_file)):
            self.mol = utils.pickle_load(os.path.join(self.save_dir, self.pre_pickle_file))
            self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
            utils.radon_print("Density Search algorithm has already completed.", level=1)

        else:
            pre_cooling_rate = 1e6
            t_start, t_stop = pre_temp_start, pre_temp_stop
            # step = (t_stop-t_start)*pre_cooling_rate
            lmp.make_dat(self.mol, file_name=self.pre_dat_file, confId=confId)
            for i in range(10):
                dt1 = datetime.datetime.now()
                utils.radon_print('Density Search algorithm (tg_pre): No%02d (%d K -> %d K)'%(i+1, t_start, t_stop), level=1)
                md1 = self.density_checker(t_start=t_start, t_stop=t_stop, step=(t_stop-t_start)*pre_cooling_rate)
                lmp.make_input(md1, file_name=self.pre_in_file)
                self.mol = lmp.run(md1, mol=self.mol, confId=confId, input_file=self.pre_in_file, omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
                dt2 = datetime.datetime.now()
                result = TGMD_analyze(log_file=os.path.join(self.work_dir_tg, self.pre_log_file), save_dir=self.save_dir).density_checker(i)
                utils.radon_print('Complete Density Search algorithm (tg_pre: %0.2f). Elapsed time = %s' %(result['tg_init_density'], str(dt2-dt1)), level=1)
                if result['tg_init_density_check']:
                    result['tg_max_temp'] = t_stop
                    break
                else:
                    t_start = t_stop
                    t_stop  = result['tg_next_temp']
            if 'tg_max_temp' not in result:
                utils.radon_print('Density Search algorithm failed. Using the last t_stop = %d K as tg_max_temp.' % t_stop, level=3)
                result['tg_max_temp'] = t_stop
            utils.MolToJSON(self.mol, os.path.join(self.save_dir, self.pre_json_file))
            utils.pickle_dump(self.mol, os.path.join(self.save_dir, self.pre_pickle_file))
            self.data.update(result)
            pd.DataFrame(self.data, index=[0]).to_csv(os.path.join(self.save_dir, self.csv_file))
        
        ####################
        #self.mol = utils.pickle_load( os.path.join(self.save_dir, self.pre_pickle_file) )
        #self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
        ####################
        
        #
        # Equilibration algorithm for Tg
        #
        if os.path.exists(os.path.join(self.save_dir, self.eq_pickle_file)):
            self.mol = utils.pickle_load(os.path.join(self.save_dir, self.eq_pickle_file))
            self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
            utils.radon_print("Equilibration (tg_eq) has already completed.", level=1)
        else:
            dt1 = datetime.datetime.now()
            utils.radon_print('Equilibration (tg_eq) by LAMMPS is running...', level=1)
            md2 = self.equilibration(temp=float(self.data['tg_max_temp']), step=5e6, **kwargs)
            lmp.make_input(md2, file_name=self.eq_in_file)
            self.mol = lmp.run(md2, mol=self.mol, confId=confId, input_file=self.eq_in_file,omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            dt2 = datetime.datetime.now()
            utils.radon_print('Complete Equilibration (tg_eq). Elapsed time = %s' % str(dt2-dt1), level=1)
            utils.MolToJSON(self.mol, os.path.join(self.save_dir, self.eq_json_file))
            utils.pickle_dump(self.mol, os.path.join(self.save_dir, self.eq_pickle_file))
        
        ####################
        #self.mol = utils.pickle_load( os.path.join(self.save_dir, self.eq_pickle_file) )
        #self.data = pd.read_csv(os.path.join(self.save_dir, self.csv_file)).iloc[0].to_dict()
        ####################
        
        
        #
        # Glass transition (Tg) algorithm
        #
        self.data['tg_min_temp'] = min_temp
        self.data['tg_cooling_rate'] = cooling_rate
        self.data['tg_interval_temp'] = interval_temp
        # previous version
        # t_step = (self.data['tg_max_temp'] - self.data['tg_min_temp'])/self.data['tg_interval_temp']*cooling_rate
        # dt1 = datetime.datetime.now()
        # utils.radon_print('Glass transition algorithm (tg) by LAMMPS is running...', level=1)
        # utils.radon_print('Total Computational Time: %d steps (= %0.2f ns)'%(t_step, t_step/1e6), level=1)
        # md3 = self.step_wise(max_temp=self.data['tg_max_temp'],
        #                      min_temp=self.data['tg_min_temp'],
        #                      cooling_rate=self.data['tg_cooling_rate'],
        #                      interval_temp=self.data['tg_interval_temp'], **kwargs)
        # lmp.make_input(md3, file_name=self.in_file)
        # self.mol = lmp.run(md3, mol=self.mol, confId=confId, input_file=self.in_file, omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
        # dt2 = datetime.datetime.now()
        # utils.radon_print('Complete glass transition algorithm (tg). Elapsed time = %s' % str(dt2-dt1), level=1)

        # new version: add checkpoint
        temp_range = np.arange(self.data['tg_max_temp'], self.data['tg_min_temp'], -self.data['tg_interval_temp'])
        save_points = [int(t) for t in temp_range] + [int(min_temp)]

        # get the checkpoint
        start_idx = 0
        for i, temp in enumerate(save_points):
            pkl = os.path.join(self.save_dir, 'Tg_conf_T%04d.pickle' % temp)
            if os.path.exists(pkl):
                start_idx = i + 1
            else:
                break

        if start_idx >= len(save_points):
            utils.radon_print('All temperature points has already done. Skip Tg MD.', level=1)

        else:
            if start_idx > 0:
                last_temp = save_points[start_idx-1]
                pkl = os.path.join(self.save_dir, 'Tg_conf_T%04d.pickle' % last_temp)
                self.mol = utils.pickle_load(pkl)
                utils.radon_print('Resume from %d K (Point %d / %d)' % (last_temp, start_idx, len(save_points)), level=1)
            else:
                utils.radon_print('Start from %d K' % save_points[0], level=1)

            dt1 = datetime.datetime.now()
            utils.radon_print('Glass transition algorithm (tg) by LAMMPS is running...', level=1)

            # run Tg MD
            for idx in range(start_idx, len(save_points)):
                temp = save_points[idx]
                is_last = (idx == len(save_points) - 1)
                md3 = self.step_wise_new(temp=temp, idx=idx, is_last=is_last,
                                         cooling_rate=self.data['tg_cooling_rate'],
                                         interval_temp=self.data['tg_interval_temp'], **kwargs)
                lmp.make_input(md3, file_name=self.in_file)
                self.mol = lmp.run(md3, mol=self.mol, confId=confId, input_file=self.in_file, omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)

                # record the start temp and save it after cooling.
                utils.MolToJSON(self.mol, os.path.join(self.save_dir, 'conf_T%04d.json' % temp))
                utils.pickle_dump(self.mol, os.path.join(self.save_dir, 'conf_T%04d.pickle' % temp))
                utils.radon_print('Complete %d K (%d/%d)' % (temp, idx + 1, len(save_points)), level=1)

            dt2 = datetime.datetime.now()
            utils.radon_print('Complete glass transition algorithm (tg). Elapsed time = %s' % str(dt2-dt1), level=1)

        # merge log for analysis
        merge_log = os.path.join(self.work_dir_tg, self.log_file)
        with open(merge_log, 'w') as f:
            for temp in save_points:
                log_path = os.path.join(self.work_dir_tg, 'Tg_conf_T%04d.log' % temp)
                with open(log_path, 'r') as fin:
                    f.write(fin.read())
                    f.write('\n\n')

        result = TGMD_analyze(log_file=merge_log, save_dir=self.save_dir).step_wise()

        if result is None:
            utils.radon_print('Tg could not be calculated in TGMD_analyze.', level=3)
        else:
            self.data.update(result)
        
        return self.mol, self.data

        
class TGMD_Additional():
    '''
    TGMD_analyze class reads and analyzes dump/log files. 
    '''
    def equilibration():
        print('hello')
    def step_wise():
        print('hello')
    def quench():
        print('hello')
    def continuous():
        print('hello')

