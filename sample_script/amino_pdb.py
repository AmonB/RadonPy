from radonpy.core import utils, poly
from radonpy.sim import qm
from radonpy.sim.preset import eq, tc
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry as Geom
import pickle
import os
import subprocess

def make_residue_pdb(res, workdir='./', save_dir='./'):

    cwd = os.getcwd()
    if work_dir is not None:
        os.chdir(work_dir)

    with open('tleap.in', 'w') as f:
        f.write('source leaprc.protein.ff19SB\n')
        f.write('mol = sequence { ACE %s NME }\n' % res)
        f.write('charge mol\n')
        f.write('savepdb mol mol.pdb\n')
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

    lines = []
    try:
        f = open('mol.pdb', 'r')
        lines = f.readlines()
        f.close()
    except:
        f.close()
        print("ERROR: Failed to read mol.pdb")
        os.chdir(cwd)
        sys.exit()

    with open('mol.prmtop', 'r') as f:
        prm_lines = f.readlines()
        
    os.chdir(cwd)

    # read charge
    charges = []
    for i, line in enumerate(prm_lines):
        if line.find('%FLAG CHARGE') == 0:
            s = ''
            for j in  range(len(prm_lines)):
                l = prm_lines[i+2+j]
                if l.find('%') == 0:
                    break
                s = s + l[:80]
            for j in range(len(s) // 16):
                q = s[j*16:(j+1)*16]
                if q.strip() == '':
                    break
                charges.append(float(q) / 18.2223)
            break

    atom_lines = []
    i = 0
    for line in lines:
        line = line.strip()
        if line[:4] == "ATOM":
            resseq = int(line[22:26])
            q = 0.0
            if resseq == 2:
                q = charges[i]
            line += '              %7.4f' % q
            atom_lines.append(line)
            i += 1
        if line[:3] == "TER":
            resseq = int(atom_lines[-1][22:26])
            if resseq == 1:
                atom_lines = []
            if resseq == 2:
                break

    new_lines = []
    try:
        for i, line in enumerate(atom_lines):
            name   = line[12:16].strip()
            resseq = int(line[22:26])
            if (resseq == 1 and name == 'C') or \
               (resseq == 3 and name == 'N'):
                line = line[:12] + ' H  ' + line[16:]
                new_lines.append(line)
            if resseq == 2:
                new_lines.append(line)

    except Exception as e:
        print(e)
        print("ERROR: Invalid syntax in mol.pdb")
        sys.exit()

    resname = res[-3:]
    resname = "%3s" % resname
    for i, line in enumerate(new_lines):
        iatom = "%5i" % (i+1)
        name = line[12:16]
        new_lines[i] = line[:6] + iatom + line[11:17] \
            + resname + line[20:22] + '   1' + line[26:]

    new_lines.append('TER   ')
    new_lines.append('END   ')

    pdbfile = os.path.join(save_dir, res + '.pdb')
    with open(pdbfile, 'w') as f:
        f.write("\n".join(new_lines))
        f.flush()
        if hasattr(os, 'fdatasync'):
            os.fdatasync(f.fileno())
        else:
            os.fsync(f.fileno())
    
if __name__ == '__main__':

    work_dir = './work_dir'
    save_dir = './pdb'

    residues = [
        # all_amino19.lib
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH",
        "GLN", "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU",
        "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
        "VAL",
        # all_aminont12.lib
        "ACE" , "NALA", "NARG", "NASN", "NASP", "NCYS", "NCYX", "NGLN",
        "NGLU", "NGLY", "NHID", "NHIE", "NHIP", "NILE", "NLEU", "NLYS",
        "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
        # all_aminoct12.lib
        "CALA", "CARG", "CASN", "CASP", "CCYS", "CCYX", "CGLN", "CGLU",
        "CGLY", "CHID", "CHIE", "CHIP", "CHYP", "CILE", "CLEU", "CLYS",
        "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL",
        "NHE" , "NME"]

    for res in residues:
        print(res)
        make_residue_pdb(res, work_dir, save_dir)
    
    sys.exit()
    

