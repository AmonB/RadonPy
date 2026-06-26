#!/usr/bin/env python3

#  Copyright (c) 2026. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

__version__ = '1.0b2'

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os

from radonpy.core import utils, calc
from radonpy.ff.gaff import GAFF
from radonpy.ff.gaff2 import GAFF2
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.ff.dreiding import Dreiding, Dreiding_UT
from radonpy.sim import qm, helper


if __name__ == '__main__':
    data = {
        'DBID': os.environ.get('RadonPy_DBID'),
        'monomer_ID': os.environ.get('RadonPy_Monomer_ID', None),
        'smiles_list': os.environ.get('RadonPy_SMILES'),
        'smiles_ter_1': os.environ.get('RadonPy_SMILES_TER', '*C'),
        'ter_ID_1': os.environ.get('RadonPy_TER_ID', 'CH3'),
        'smiles_ter_2': os.environ.get('RadonPy_SMILES_TER2', None),
        'ter_ID_2': os.environ.get('RadonPy_TER_ID2', None),
        'qm_method': os.environ.get('RadonPy_QM_Method', 'wb97m-d3bj'),
        'charge': os.environ.get('RadonPy_Charge', 'RESP'),
        'qm_td_method': os.environ.get('RadonPy_QM_TD_Method', 'cam-b3lyp-d3bj'),
        'forcefield': str(os.environ.get('RadonPy_FF', 'GAFF2_mod')),
        'remarks': os.environ.get('RadonPy_Remarks', ''),
        **helper.get_version()
    }

    omp_psi4 = int(os.environ.get('RadonPy_OMP_Psi4', 4))
    mem_psi4 = int(os.environ.get('RadonPy_MEM_Psi4', 1000))

    conf_mm_omp = int(os.environ.get('RadonPy_Conf_MM_OMP', 1))
    conf_mm_mpi = int(os.environ.get('RadonPy_Conf_MM_MPI', utils.cpu_count()))
    conf_mm_gpu = int(os.environ.get('RadonPy_Conf_MM_GPU', 0))
    conf_mm_mp = int(os.environ.get('RadonPy_Conf_MM_MP', 0))
    conf_psi4_omp = int(os.environ.get('RadonPy_Conf_Psi4_OMP', omp_psi4))
    conf_psi4_mp = int(os.environ.get('RadonPy_Conf_Psi4_MP', 0))

    do_ter = bool(os.environ.get('RadonPy_Do_Ter', False)=='True')
    do_tddft = bool(os.environ.get('RadonPy_Do_TDDFT', False)=='True')


    work_dir = './%s' % data['DBID']
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    save_dir = os.path.join(work_dir, 'analyze')
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    tmp_dir = os.environ.get('RadonPy_TMP_Dir', work_dir)
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
        
    io = helper.IO_Helper(work_dir, save_dir)

    smi_list = data['smiles_list'].split(',')
    if data['monomer_ID']: monomer_id = data['monomer_ID'].split(',')

    if data['forcefield'] == 'GAFF':
        ff = GAFF()
    elif data['forcefield'] == 'GAFF2':
        ff = GAFF2()
    elif data['forcefield'] == 'GAFF2_mod':
        ff = GAFF2_mod()
    elif data['forcefield'] == 'Dreiding':
        ff = Dreiding()
    elif data['forcefield'] == 'Dreiding_UT':
        ff = Dreiding_UT()
    else:
        raise ValueError("Force field %s is not available." % data['forcefield'])

    for i, smi in enumerate(smi_list):
        monomer_data = {
            'smiles': smi,
            'qm_method': data['qm_method'],
            'charge': data['charge'],
            'remarks': data['remarks'],
            **helper.get_version(),
        }
        data['smiles_%i' % (i+1)] = smi

        # Conformation search and RESP charge calculation of a repeating unit
        mol = utils.mol_from_smiles(smi)
        mol, energy = qm.conformation_search(mol, ff=ff, work_dir=work_dir, tmp_dir=tmp_dir, opt_method=data['qm_method'],
            psi4_omp=conf_psi4_omp, psi4_mp=conf_psi4_mp, mpi=conf_mm_mpi, omp=conf_mm_omp, gpu=conf_mm_gpu,
            mm_mp=conf_mm_mp, log_name='monomer%i' % (i+1), memory=mem_psi4)
        qm.assign_charges(mol, charge=data['charge'], work_dir=work_dir, tmp_dir=tmp_dir, omp=omp_psi4, opt=False, log_name='monomer%i' % (i+1), memory=mem_psi4)

        # Dump pickle file
        if data['monomer_ID']:
            data['monomer_ID_%i' % (i+1)] = monomer_data['monomer_ID'] = monomer_id[i]
            utils.pickle_dump(mol, os.path.join(save_dir, 'monomer_%s.pickle' % monomer_id[i]))
            utils.MolToJSON(mol, os.path.join(save_dir, 'monomer_%s.json' % monomer_id[i]))
        else:
            utils.pickle_dump(mol, os.path.join(save_dir, 'monomer%i.pickle' % (i+1)))
            utils.MolToJSON(mol, os.path.join(save_dir, 'monomer%i.json' % (i+1)))

        # Get monomer properties
        update = {
            'mol_weight': calc.molecular_weight(mol),
            'vdw_volume': calc.vdw_volume(mol)
        }

        # Output monomer properties
        data, monomer_data = io.update_monomer_data(update, data, monomer_data, monomer_idx=i)

        # Single point calculation
        sp_data = qm.sp_prop(mol, opt=False, work_dir=work_dir, tmp_dir=tmp_dir, sp_method=data['qm_method'],
                            omp=omp_psi4, log_name='monomer%i' % (i+1), memory=mem_psi4)
        data, monomer_data = io.update_monomer_data(sp_data, data, monomer_data, monomer_idx=i)

        # Polarizability calculation
        polar_data = qm.polarizability(mol, opt=False, work_dir=work_dir, tmp_dir=tmp_dir, polar_method=data['qm_method'], omp=conf_psi4_omp, mp=conf_psi4_mp,
                                        log_name='monomer%i' % (i+1), memory=mem_psi4)
        data, monomer_data = io.update_monomer_data(polar_data, data, monomer_data, monomer_idx=i)

        if do_tddft:
            # Frequency dependent polarizability calculation
            fd_polar_data = qm.polarizability_sos(mol, wavelength=[486, 589, 656], p_state=0.003, opt=False, work_dir=work_dir, save_dir=save_dir, tmp_dir=tmp_dir,
                                            td_method=data['qm_td_method'], omp=conf_psi4_omp, mp=conf_psi4_mp, log_name='monomer%i' % (i+1), memory=mem_psi4)
            data, monomer_data = io.update_monomer_data(fd_polar_data, data, monomer_data, monomer_idx=i)


    if do_ter:
        ter1 = utils.mol_from_smiles(data['smiles_ter_1'])
        qm.assign_charges(ter1, charge=data['charge'], work_dir=work_dir, tmp_dir=tmp_dir, opt_method=data['qm_method'], omp=omp_psi4, log_name='ter1', memory=mem_psi4)
        if data['ter_ID_1']:
            utils.pickle_dump(ter1, os.path.join(save_dir, 'ter_%s.pickle' % data['ter_ID_1']))
            utils.MolToJSON(ter1, os.path.join(save_dir, 'ter_%s.json' % data['ter_ID_1']))
        else:
            utils.pickle_dump(ter1, os.path.join(save_dir, 'ter1.pickle'))
            utils.MolToJSON(ter1, os.path.join(save_dir, 'ter1.json'))

        if data['smiles_ter_2'] is not None:
            ter2 = utils.mol_from_smiles(data['smiles_ter_2'])
            qm.assign_charges(ter2, charge=data['charge'], work_dir=work_dir, tmp_dir=tmp_dir, opt_method=data['qm_method'], omp=omp_psi4, log_name='ter2', memory=mem_psi4)
            if data['ter_ID_2']:
                utils.pickle_dump(ter2, os.path.join(save_dir, 'ter_%s.pickle' % data['ter_ID_2']))
                utils.MolToJSON(ter2, os.path.join(save_dir, 'ter_%s.json' % data['ter_ID_2']))
            else:
                utils.pickle_dump(ter2, os.path.join(save_dir, 'ter2.pickle'))
                utils.pickle_dump(ter2, os.path.join(save_dir, 'ter2.json'))
