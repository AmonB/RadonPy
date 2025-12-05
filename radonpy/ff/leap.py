#  Copyright (c) 2025. RadonPy developers. All rights reserved.
#  Use of this source code is governed by a BSD-3-style
#  license that can be found in the LICENSE file.

# ******************************************************************************
# ff.leap module
# ******************************************************************************

import numpy as np
import os
import json
from itertools import permutations
from rdkit import Chem
from ..core import calc, utils
from . import ff_class
import subprocess
import sys

__version__ = '1.0b1'

class LEaP_FF_LAMMPS():
    """
    leap.LEaP_FF_LAMMPS() class

    LAMMPS Forcefield data class assigned with LEaP.

    Attributes:
        atom_type_names
        pair_coeffs
        bond_coeffs
        angle_coeffs
        dihedral_coeffs
        atoms
        bonds
        angles
        dihedrals
        charges
        cmap_params
        cmaps
    """

    def __init__(self):
        self.atom_type_names = []
        self.pair_coeffs = []
        self.bond_coeffs = []
        self.angle_coeffs = []
        self.dihedral_coeffs = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.charges = []
        self.cmap_params = []
        self.cmaps = []
    

    def assign(self, mol, work_dir, leaprc):
        
        cwd = os.getcwd()
        if work_dir is not None:
            os.chdir(work_dir)

        utils.MolToPDBFile(mol, 'tleap.pdb')

        with open('tleap.in', 'w') as f:
            f.write('source %s\n' % leaprc)
            f.write('mol = loadpdb tleap.pdb\n')
            f.write('SaveAmberParm mol mol.prmtop mol.inpcrd\n')
            f.write('quit\n')
            f.flush()
            if hasattr(os, 'fdatasync'):
                os.fdatasync(f.fileno())
            else:
                os.fsync(f.fileno())
                
            cmd  = 'tleap -f tleap.in'
            utils.radon_print(cmd)

            cp = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='UTF-8')
            with open('tleap.log', 'a') as f:
                f.write(cmd+'\n')
                f.write(cp.stdout+'\n')
                f.write(cp.stderr+'\n')
                f.write('tLEaP returncode = %s \n' % (str(cp.returncode)))

        cmd  = 'intermol-convert --amb_in mol.prmtop mol.inpcrd --lammps'
        utils.radon_print(cmd)

        cp = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, encoding='UTF-8')
        with open('intermol.log', 'a') as f:
            f.write(cmd+'\n')
            f.write(cp.stdout+'\n')
            f.write(cp.stderr+'\n')
            f.write('intermol-convert returncode = %s \n' % (str(cp.returncode)))

        with open('mol.prmtop', 'r') as f:
            prm_lines = f.readlines()

        with open('mol_from_amber.top', 'r') as f:
            top_lines = f.readlines()
            
        with open('mol_converted.input', 'r') as f:
            inp_lines = f.readlines()
            
        with open('mol_converted.lmp', 'r') as f:
            lmp_lines = f.readlines()
            
        os.chdir(cwd)

        natom = 0
        nbond = 0
        nangle = 0
        ndihedral = 0
        natom_type = 0
        nbond_type = 0
        nangle_type = 0
        ndihedral_type = 0

        bond_coeffs = []
        angle_coeffs = []
        dihedral_coeffs = []

        atoms = []
        bonds = []
        angles = []
        dihedrals = []

        charges = []

        for i, line in enumerate(lmp_lines):
            line = line.strip()
            if line.find('atoms') > -1:
                natom = int(line.split()[0])
                continue
            if line.find('bonds') > -1:
                nbond = int(line.split()[0])
                continue
            if line.find('angles') > -1:
                nangle = int(line.split()[0])
                continue
            if line.find('dihedrals') > -1:
                ndihedral = int(line.split()[0])
                continue
            if line.find('atom types') > -1:
                natom_type = int(line.split()[0])
                continue
            if line.find('bond types') > -1:
                nbond_type = int(line.split()[0])
                continue
            if line.find('angle types') > -1:
                nangle_type = int(line.split()[0])
                continue
            if line.find('dihedral types') > -1:
                ndihedral_type = int(line.split()[0])
                continue
            if line.find('Bond Coeffs') == 0:
                for j in range(nbond_type):
                    words = lmp_lines[i+2+j].strip().split()
                    k  = float(words[2])
                    r0 = float(words[3])
                    bond_coeffs.append([k, r0])
                continue
            if line.find('Angle Coeffs') == 0:
                for j in range(nangle_type):
                    words = lmp_lines[i+2+j].strip().split()
                    k      = float(words[2])
                    theta0 = float(words[3])
                    angle_coeffs.append([k, theta0])
                continue
            if line.find('Atoms') == 0:
                for j in range(natom):
                    words = lmp_lines[i+2+j].strip().split()
                    atom_type = int(words[2]) - 1
                    q = float(words[3])
                    atoms.append([atom_type])
                    charges.append(q)
                continue
            if line.find('Bonds') == 0:
                for j in range(nbond):
                    words = lmp_lines[i+2+j].strip().split()
                    btype = int(words[1]) - 1
                    i0 = int(words[2]) - 1
                    i1 = int(words[3]) - 1
                    bonds.append([btype, i0, i1])
                continue
            if line.find('Angles') == 0:
                for j in range(nangle):
                    words = lmp_lines[i+2+j].strip().split()
                    atype = int(words[1]) - 1
                    i0 = int(words[2]) - 1
                    i1 = int(words[3]) - 1
                    i2 = int(words[4]) - 1
                    angles.append([atype, i0, i1, i2])
                continue
            
        pair_coeffs = {}
        for i, line in enumerate(inp_lines):
            line = line.strip()
            if line.find('pair_coeff') == 0:
                words = line.split()
                i1 = int(words[1])
                i2 = int(words[2])
                eps = float(words[3])
                sig = float(words[4])
                if i1 == i2:
                    pair_coeffs[i1-1] = [eps, sig]
                continue
        
        atypenames = []
        for i, line in enumerate(top_lines):
            if line.find('[ atomtypes ]') == 0:
                for j in  range(natom_type):
                    l = top_lines[i+2+j].strip()
                    atypenames.append(l.split()[0])
                continue
            if line.find('[ dihedrals ]') == 0:
                for j in  range(ndihedral):
                    l = top_lines[i+2+j].strip()
                    words = l.split()
                    i0 = int  (words[0]) - 1
                    i1 = int  (words[1]) - 1
                    i2 = int  (words[2]) - 1
                    i3 = int  (words[3]) - 1
                    d0 = float(words[5])
                    k  = float(words[6]) / 4.184 # KJ -> kcal
                    n  = int  (words[7])
                    dihedrals.append([i0, i1, i2, i3, d0, k, n])
                continue

        # Atoms may be reordered, so match them by atomic name
        atomnames = []
        for i, line in enumerate(prm_lines):
            if line.find('%FLAG ATOM_NAME') == 0:
                nline = (natom - 1) // 20 + 1
                l = "".join(prm_lines[i+2:i+2+nline]).replace("\n", "")
                for j in range(natom):
                    s = l[4*j:4*(j+1)]
                    atomnames.append(s.strip())
                continue

        res_idx = []
        for i, atom in enumerate(mol.GetAtoms()):
            ires = atom.GetPDBResidueInfo().GetResidueNumber() - 1
            res_idx.append(ires)

        ires_prev = -1
        res_ptr = []
        for i, ires in enumerate(res_idx):
            if ires_prev != ires:
                res_ptr.append(i)
                ires_prev = ires
        res_ptr.append(natom)
            
        fwd = [-1] * natom
        rev = [-1] * natom
        for i, name in enumerate(atomnames):
            ires = res_idx[i]
            for j in range(res_ptr[ires], res_ptr[ires+1]):
                mol_name = mol.GetAtomWithIdx(j).GetPDBResidueInfo().GetName()
                if name == mol_name.strip():
                    fwd[j] = i
                    rev[i] = j
                    break
            
        ncmap = 0
        ncmap_type = 0
        cmap_params = []
        cmaps = []
        for i, line in enumerate(prm_lines):
            if line.find('%FLAG CMAP_COUNT') == 0:
                l = prm_lines[i+2]
                ncmap      = int(l[0: 8])
                ncmap_type = int(l[8:16])
                continue
            if line.find('%FLAG CMAP_PARAMETER_') == 0:
                l = prm_lines[i+1]
                comment = l.replace('%COMMENT', '').strip()
                cmap_param = []
                for j in range(72): # 24*24/8
                    l = prm_lines[i+3+j]
                    for k in range(8):
                        s = l[k*9:(k+1)*9]
                        cmap_param.append(float(s))
                cmap_params.append([comment, cmap_param])
                continue
            if line.find('%FLAG CMAP_INDEX') == 0:
                for j in range(ncmap):
                    l = prm_lines[i+2+j]
                    cmap = []
                    for k in range(6):
                        s = l[k*8:(k+1)*8]
                        cmap.append(int(s) - 1)
                    cmaps.append([cmap[5], cmap[:5]])
                continue


        atoms     = [atoms[fwd[i]] for i in range(natom)]
        bonds     = [[b[0], rev[b[1]], rev[b[2]]] for b in bonds]
        angles    = [[a[0], rev[a[1]], rev[a[2]], rev[a[3]]] for a in angles]
        dihedrals = [[rev[d[0]], rev[d[1]], rev[d[2]], rev[d[3]], d[4], d[5], d[6]]
                     for d in dihedrals]
        charges   = [charges[fwd[i]] for i in range(natom)]
        for c in cmaps:
            c[1] = [rev[i] for i in c[1]]

        self.atom_type_names = atypenames
        self.pair_coeffs = pair_coeffs
        self.bond_coeffs = bond_coeffs
        self.angle_coeffs = angle_coeffs
        self.dihedral_coeffs = dihedral_coeffs
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.charges = charges
        self.cmap_params = cmap_params
        self.cmaps = cmaps

        return

