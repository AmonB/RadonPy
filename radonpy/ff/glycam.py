#  Copyright (c) 2025. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# ******************************************************************************
# ff.glycam module
# ******************************************************************************

import numpy as np
import os
import json
from itertools import permutations
from rdkit import Chem
from ..core import calc, utils
from . import ff_class, leap
import subprocess
import sys
from . import leap

__version__ = '1.0b1'


class GLYCAM_06j():
    """
    glycam.GLYCAM_06j() class

    Forcefield object with typing rules for Glycam model.

    Attributes:
        ff_name: glycam
        pair_style: lj
        bond_style: harmonic
        angle_style: harmonic
        dihedral_style: fourier
        improper_style: cvff
    """

    def __init__(self, db_file=None, work_dir=None):
        self.name = 'glycam'
        self.pair_style = 'lj'
        self.bond_style = 'harmonic'
        self.angle_style = 'harmonic'
        self.dihedral_style = 'fourier'
        self.improper_style = 'cvff'
        self.work_dir = work_dir

    def ff_assign(self, mol, charge=None, retryMDL=True, useMDL=True):
        """
        GLYCAM_06j.ff_assign

        GLYCAM_06j force field assignment for RDkit Mol object

        Args:
            mol: rdkit mol object

        Optional args:
            charge: Method of charge assignment. If None, charge assignment is skipped. 
            retryMDL: Retry assignment using MDL aromaticity model if default aromaticity model is failure (boolean)
            useMDL: Assignment using MDL aromaticity model (boolean)

        Returns: (boolean)
            True: Success assignment
            False: Failure assignment
        """

        if useMDL:
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)

        mol.SetProp('ff_name', str(self.name))

        lmp = leap.LEaP_FF_LAMMPS()
        lmp.assign(mol, self.work_dir, 'leaprc.GLYCAM_06j-1')
        
        if charge == "lmp":
            for i, q in enumerate(lmp.charges):
                mol.GetAtomWithIdx(i).SetDoubleProp('AtomicCharge', q)
        
        result = self.assign_ptypes(mol, lmp.atoms, lmp.atom_type_names,
                                    lmp.pair_coeffs)
        if result: result = self.assign_btypes(mol, lmp.bonds, lmp.bond_coeffs)
        if result: result = self.assign_atypes(mol, lmp.angles, lmp.angle_coeffs)
        if result: result = self.assign_dtypes(mol, lmp.dihedrals)
        if result: result = self.assign_itypes(mol)
        if result and charge is not None and charge != "lmp":
            result = calc.assign_charges(mol, charge=charge)

        return result


    def assign_ptypes(self, mol, atoms, atypenames, pair_coeffs):
        """
        GLYCAM_06j.assign_ptypes

        GLYCAM_06j specific particle typing rules.

        Args:
            mol: rdkit mol object

        Returns:
            boolean
        """
        result_flag = True
        mol.SetProp('pair_style', self.pair_style)
        
        for i, atom in enumerate(atoms):
            it = atoms[i][0]
            p = mol.GetAtomWithIdx(i)
            self.set_ptype(p, atypenames[it], pair_coeffs[it])

        return result_flag

    def set_ptype(self, p, pt, pair_coeff):
        p.SetProp('ff_type', pt)
        p.SetDoubleProp('ff_epsilon', pair_coeff[0])
        p.SetDoubleProp('ff_sigma', pair_coeff[1])

        return p
        
        
    def assign_btypes(self, mol, bonds, bond_coeffs):
        """
        GLYCAM_06j.assign_btypes

        GLYCAM_06j specific bond typing rules.

        Args:
            mol: rdkit mol object

        Returns:
            boolean
        """
        result_flag = True
        mol.SetProp('bond_style', self.bond_style)

        for i, b in enumerate(mol.GetBonds()):
            a1 = b.GetBeginAtom()
            a2 = b.GetEndAtom()
            ba = a1.GetProp('ff_type')
            bb = a2.GetProp('ff_type')

            ia1 = a1.GetIdx()
            ia2 = a2.GetIdx()

            for bond in bonds:
                it = bond[0]
                i1 = bond[1]
                i2 = bond[2]
                if ((ia1 == i1 and ia2 == i2) or 
                    (ia1 == i2 and ia2 == i1)):
                    coeff = bond_coeffs[it]
                    bt = '%s,%s,%f,%f' % (ba, bb, coeff[0], coeff[1])
                    self.set_btype(b, bt, coeff)
                    
        return result_flag
    

    def set_btype(self, b, bt, bond_coeff):

        b.SetProp('ff_type', bt)
        b.SetDoubleProp('ff_k', bond_coeff[0])
        b.SetDoubleProp('ff_r0', bond_coeff[1])

        return True
        

    def assign_atypes(self, mol, angles, angle_coeffs):
        """
        GLYCAM_06j.assign_atypes

        GLYCAM_06j specific angle typing rules.

        Args:
            mol: rdkit mol object

        Returns:
            boolean
        """
        result_flag = True
        mol.SetProp('angle_style', self.angle_style)
        setattr(mol, 'angles', {})
        
        for angle in angles:
            at = angle[0]
            i0 = angle[1]
            i1 = angle[2]
            i2 = angle[3]

            p0 = mol.GetAtomWithIdx(i0)
            p1 = mol.GetAtomWithIdx(i1)
            p2 = mol.GetAtomWithIdx(i2)
            
            coeff = angle_coeffs[at]

            pt0 = p0.GetProp('ff_type')
            pt1 = p1.GetProp('ff_type')
            pt2 = p2.GetProp('ff_type')
            at = '%s,%s,%s,%f,%f' % (pt0, pt1, pt2, coeff[0], coeff[1])
            
            result = self.set_atype(mol, i0, i1, i2, at, coeff)
            if not result:
                result_flag = False
            
        return result_flag

    def set_atype(self, mol, a, b, c, at, coeff):
            
        angle = utils.Angle(
            a=a, b=b, c=c,
            ff=ff_class.Angle_harmonic(
                ff_type=at,
                k=coeff[0],
                theta0=coeff[1]
            )
        )
        
        key = '%i,%i,%i' % (a, b, c)
        mol.angles[key] = angle
        
        return True

    def assign_dtypes(self, mol, dihedrals):
        """
        GLYCAM_06j.assign_dtypes

        GLYCAM_06j specific dihedral typing rules.
        
        Args:
            mol: rdkit mol object

        Returns:
            boolean
        """
        result_flag = True
        mol.SetProp('dihedral_style', self.dihedral_style)
        setattr(mol, 'dihedrals', {})

        dih_dict = {}
        
        for d in dihedrals:
            i0 = d[0]
            i1 = d[1]
            i2 = d[2]
            i3 = d[3]
            key = '%i,%i,%i,%i' % (i0, i1, i2, i3)

            if key in dih_dict:
                dih_dict[key].append(d)
            else:
                dih_dict[key] = [d]

        for key, dih_list in dih_dict.items():
            d = dih_list[0]
            i0 = d[0]
            i1 = d[1]
            i2 = d[2]
            i3 = d[3]

            p0 = mol.GetAtomWithIdx(i0)
            p1 = mol.GetAtomWithIdx(i1)
            p2 = mol.GetAtomWithIdx(i2)
            p3 = mol.GetAtomWithIdx(i3)
            
            p0t = p0.GetProp('ff_type')
            p1t = p1.GetProp('ff_type')
            p2t = p2.GetProp('ff_type')
            p3t = p3.GetProp('ff_type')
            dt = '%s,%s,%s,%s' % (p0t, p1t, p2t, p3t)
            for d in dih_list:
                dt += ',%f,%f,%f' % (d[4], d[5], d[6])
            
            result = self.set_dtype(mol, i0, i1, i2, i3, dt, dih_list)
            if not result:
                result_flag = False
                
        return result_flag


    def set_dtype(self, mol, a, b, c, d, dt, dih_list):

        key = '%i,%i,%i,%i' % (a, b, c, d)
        m  = len(dih_list)
        d0 = [ d[4] for d in dih_list ]
        k  = [ d[5] for d in dih_list ]
        n  = [ d[6] for d in dih_list ]
        
        ff=ff_class.Dihedral_fourier(
            ff_type=dt,
            k=k,
            d0=d0,
            m=m,
            n=n,
        )
        dihedral = utils.Dihedral(a=a, b=b, c=c, d=d, ff=ff)
        mol.dihedrals[key] = dihedral
        
        return True

    def assign_itypes(self, mol):
        """
        GLYCAM_06j.assign_itypes

        GLYCAM_06j specific improper typing rules.

        Args:
            mol: rdkit mol object

        Returns:
            boolean
        """
        mol.SetProp('improper_style', self.improper_style)
        setattr(mol, 'impropers', {})
        
        return True            

    def load_ff_json(self, json_file):
        with open(json_file) as f:
            j = json.loads(f.read())

        ff = self.Container()
        ff.pt = {}
        ff.bt = {}
        ff.at = {}
        ff.dt = {}
        ff.it = {}

        ff.ff_name = j.get('ff_name')
        ff.ff_class = j.get('ff_class')
        ff.pair_style = j.get('pair_style')
        ff.bond_style = j.get('bond_style')
        ff.angle_style = j.get('angle_style')
        ff.dihedral_style = j.get('dihedral_style')
        ff.improper_style = j.get('improper_style')
        
        for pt in j.get('particle_types'):
            pt_obj = self.Container()
            for key in pt.keys():
                setattr(pt_obj, key, pt[key])
            ff.pt[pt['name']] = pt_obj
        
        for bt in j.get('bond_types'):
            bt_obj = self.Container()
            for key in bt.keys():
                setattr(bt_obj, key, bt[key])
            ff.bt[bt['name']] = bt_obj
            ff.bt[bt['rname']] = bt_obj
        
        for at in j.get('angle_types'):
            at_obj = self.Container()
            for key in at.keys():
                setattr(at_obj, key, at[key])
            ff.at[at['name']] = at_obj
            ff.at[at['rname']] = at_obj
        
        for dt in j.get('dihedral_types'):
            dt_obj = self.Container()
            for key in dt.keys():
                setattr(dt_obj, key, dt[key])
            ff.dt[dt['name']] = dt_obj
            ff.dt[dt['rname']] = dt_obj
        
        for it in j.get('improper_types'):
            it_obj = self.Container()
            for key in it.keys():
                setattr(it_obj, key, it[key])
            ff.it[it['name']] = it_obj
        
        return ff
            
    
    class Container(object):
        pass


    ## Backward compatibility
    class Angle_ff():
        """
            GLYCAM_06j.Angle_ff() object
        """
        def __init__(self, ff_type=None, k=None, theta0=None):
            self.type = ff_type
            self.k = k
            self.theta0 = theta0
            self.theta0_rad = theta0*(np.pi/180)

        def to_dict(self):
            dic = {
                'ff_type': str(self.type),
                'k': float(self.k),
                'theta0': float(self.theta0),
            }
            return dic
        
        
    class Dihedral_ff():
        """
            GLYCAM_06j.Dihedral_ff() object
        """
        def __init__(self, ff_type=None, k=[], d0=[], m=None, n=[]):
            self.type = ff_type
            self.k = np.array(k)
            self.d0 = np.array(d0)
            self.d0_rad = np.array(d0)*(np.pi/180)
            self.m = m
            self.n = np.array(n)
        
        def to_dict(self):
            dic = {
                'ff_type': str(self.type),
                'k': [float(x) for x in self.k],
                'd0': [float(x) for x in self.d0],
                'm': int(self.m),
                'n': [int(x) for x in self.n],
            }
            return dic

        
    class Improper_ff():
        """
            GLYCAM_06j.Improper_ff() object
        """
        def __init__(self, ff_type=None, k=None, d0=-1, n=None):
            self.type = ff_type
            self.k = k
            self.d0 = d0
            self.n = n
        
        def to_dict(self):
            dic = {
                'ff_type': str(self.type),
                'k': float(self.k),
                'd0': float(self.d0),
                'n': int(self.n),
            }
            return dic

