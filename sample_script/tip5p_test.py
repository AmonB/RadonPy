from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.ff.tip import TIP5P
from radonpy.sim import qm
from radonpy.sim.preset import eq, tc
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry as Geom
import pickle
import os

work_dir = './work_dir'
ff_tip = TIP5P()
ff = GAFF2_mod()
temp = 300
omp_psi4 = 8
mpi = 8
omp = 1
gpu = 0
mem = 10000

if __name__ == '__main__':

    wat = TIP5P.mol()

    mol = utils.mol_from_smiles("CC")
    mol, energy = qm.conformation_search(mol, ff=ff, work_dir=work_dir,
                                         psi4_omp=omp_psi4, mpi=mpi, omp=omp,
                                         memory=mem, log_name='monomer1')
    qm.assign_charges(mol, charge='RESP', opt=False, work_dir=work_dir,
                      omp=omp_psi4, memory=mem, log_name='monomer1')

    #print(Chem.MolToMolBlock(mol))

    # Force field assignment
    result = ff.ff_assign(mol)
    if not result:
        print('[ERROR: Can not assign force field parameters.]') 

    # Force field assignment
    result = ff_tip.ff_assign(wat)
    if not result:
        print('[ERROR: Can not assign force field parameters.]') 

    # Generate simulation cell
    ac = poly.amorphous_cell([mol, wat], [1, 10], density=0.05)
    #ac = poly.amorphous_cell([mol], [1], density=0.05)
    #ac = poly.amorphous_cell([wat], [10], density=0.05)

    # Equilibration MD
    eqmd = eq.EQ21step(ac, work_dir=work_dir)
    ac = eqmd.exec(temp=temp, press=1.0, mpi=mpi, omp=omp, gpu=gpu)

    analy = eqmd.analyze()
    prop_data = analy.get_all_prop(temp=temp, press=1.0, save=True)
    result = analy.check_eq()

    # Additional equilibration MD
    for i in range(4):
        if result: break
        eqmd = eq.Additional(ac, work_dir=work_dir)
        ac = eqmd.exec(temp=temp, press=press, mpi=mpi, omp=omp, gpu=gpu)
        analy = eqmd.analyze()
        prop_data = analy.get_all_prop(temp=temp, press=press, save=True)
        result = analy.check_eq()

    if not result:
        print('[ERROR: Did not reach an equilibrium state.]')
        sys.exit()
        
    # Non-equilibrium MD for thermal condultivity
    nemd = tc.NEMD_MP(ac, work_dir=work_dir)
    ac = nemd.exec(decomp=True, temp=temp, mpi=mpi, omp=omp, gpu=gpu)
    nemd_analy = nemd.analyze()
    TC = nemd_analy.calc_tc(decomp=True, save=True)
    if not nemd_analy.Tgrad_data['Tgrad_check']:
        print('[ERROR: Low linearity of temperature gradient.]')

    print('Thermal conductivity: %f' % TC) 
    
    sys.exit()
    

