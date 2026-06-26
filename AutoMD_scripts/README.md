# AutoMD Scripts
These scripts were used in the data generation for our database project. Calculation conditions and parameters can be controlled by environment variables.
Run the scripts in the order listed below.
- 0_qm.py
- 1_eq.py
- 2_rst_eq.py
- 3_tc.py
- 4_tg.py
- 5_sp.py
- 6_ef_dp.py

## Descriptions of each script
- 0_qm.py
  - Description: For a repeating unit, 3D coordinates generation, conformation search, geometry optimization (DFT), atomic charge (RESP) calculation (HF), polarizability calculation (DFT).
  - Input environment variables:
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_SMILES : (Required) SMILES strings of a repeating unit in polymer for calculation target.
    - RadonPy_Monomer_ID : (Optional) Arbitrary strings as repeating unit ID when loading a pre-calculated repeating unit file in 1_eq.py.
    - RadonPy_QM_Method : (Optional) A string of functional used in DFT calculations. (default: wb97m-d3bj)
    - RadonPy_Charge : (Optional) A string of atomic charge model. (default: RESP)
    - RadonPy_FF : (Optional) A string of force field name used in conformation search by MM level. (default: GAFF2_mod, available: GAFF | GAFF2 | GaFF2_mod | Dreiding | Dreiding_UT)
    - RadonPy_Do_TDDFT : (Optional) A boolean of performing additional TD-DFT calculations. (default: False)
    - RadonPy_QM_TD_Method : (Optional) A string of functional used in TD-DFT calculations. (default: cam-b3lyp-d3bj)
    - RadonPy_OMP_Psi4 : (Optional) An int value of pallarel number of OpenMP in Psi4. (default: 4)
    - RadonPy_MEM_Psi4 : (Optional) An int value of limitation of memory amount (MB) used in Psi4. (default: 1000)
    - RadonPy_Conf_MM_OMP : (Optional) An int value of pallarel number of OpenMP in the conformation search using LAMMPS. (default: 1)
    - RadonPy_Conf_MM_MPI : (Optional) An int value of pallarel number of MPI in the conformation search using LAMMPS. (default: number of available CPU cores)
    - RadonPy_Conf_MM_GPU : (Optional) An int value of number of using GPUs in the conformation search using LAMMPS. (default: 0)
    - RadonPy_Conf_MM_MP : (Optional) An int value of pallarel number of Multi-processing of Python in the conformation search using LAMMPS. (default: 0)
    - RadonPy_Conf_Psi4_OMP : (Optional) An int value of pallarel number of OpenMP in the conformation search using Psi4. (default: RadonPy_OMP_Psi4 value)
    - RadonPy_Conf_Psi4_MP : (Optional) An int value of pallarel number of Multi-processing of Python in the conformation search using Psi4. (default: 0)
    - RadonPy_Do_Ter : (Optional) A boolean of performing QM calculation for terminating groups. (default: False)
    - RadonPy_SMILES_TER : (Optional) SMILES strings of a terminating group (head) in polymer.
    - RadonPy_TER_ID : (Optional) Arbitrary strings as terminating group ID (head) when loading a pre-calculated repeating unit file in 1_eq.py.
    - RadonPy_SMILES_TER2 : (Optional) SMILES strings of a terminating group (tail) in polymer. If it is omitted, it use RadonPy_SMILES_TER value.
    - RadonPy_TER_ID2 : (Optional) Arbitrary strings as terminating group ID (tail) when loading a pre-calculated repeating unit file in 1_eq.py. If it is omitted, it use RadonPy_TER_ID value.

- 1_eq.py
  - Description: Generation of polymer chains, assignment of force field types, generation of a simulation cell as an initial structure, equilibration MD simulation, check equilibrium state using the preset `sim.preset.eq`.
  - Input environment variables:
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_SMILES : (Required) SMILES strings of a repeating unit in polymer for calculation target.
    - RadonPy_Monomer_Dir : (Optional) Directory path of pre-calculated repeating unit files.
    - RadonPy_Monomer_ID : (Optional) Strings as repeating unit ID when loading a pre-calculated repeating unit file.
    - RadonPy_Ter_Dir : (Optional) Directory path of pre-calculated repeating unit files.
    - RadonPy_TER_ID : (Optional) Arbitrary strings as terminating group ID (head) when loading a pre-calculated repeating unit file in 1_eq.py.
    - RadonPy_TER_ID2 : (Optional) Arbitrary strings as terminating group ID (tail) when loading a pre-calculated repeating unit file in 1_eq.py. If it is omitted, it use RadonPy_TER_ID value.
    - RadonPy_NAtom : (Optional)  (default: 1000)
    - RadonPy_NChain : (Optional)  (default: 10)
    - RadonPy_Tacticity : (Optional)  (default: atactic, available: atactic | isotactic | syndiotactic)
    - RadonPy_Copoly_Ratio : (Optional) 
    - RadonPy_Copoly_Type : (Optional)  (default: random, available: random | alternating | block)
    - RadonPy_Ini_Density : (Optional)  (default: 0.05, unit: g/cm^3)
    - RadonPy_FF : (Optional) A string of force field name used in conformation search by MM level. (default: GAFF2_mod, available: GAFF | GAFF2 | GaFF2_mod | Dreiding | Dreiding_UT)
    - RadonPy_Temp : (Optional) A float value of simulation temperature. (unit: K)
    - RadonPy_Press : (Optional) A float value of simulation pressure. (unit: atm)
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_RetryEQ : (Optional)  (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)

- 2_rst_eq.py
  - Description: Restarting equilibration MD simulation if system does not reach equilibrium state after 1_eq.py.


- 3_tc.py
  - Description: Calculation of thermal conductivity by non-equilibrium MD (NEMD) simulation using the preset `sim.preset.tc`


- 4_tg.py
  - Description: Calculation of Tg by NEMD simulation using the preset `sim.preset.tg`


- 5_sp.py
  - Description: Calculation of solubility parameters using the preset `sim.preset.sp`


- 6_ef_dp.py
  - Description: Calculation of dynamic dielectric constants and loss tangent by NEMD simulation using `sim.preset.ef_dp`

