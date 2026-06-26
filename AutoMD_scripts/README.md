# AutoMD Scripts
These scripts have been used in the data generation for our database project. Calculation conditions and parameters can be controlled by environment variables, see in .sh files. Run the scripts in the order listed below.
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
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_Temp : (Optional) A float value of simulation temperature. (unit: K)
    - RadonPy_Press : (Optional) A float value of simulation pressure. (unit: atm)
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_RetryEQ : (Optional)  (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)
    - RadonPy_LAMMPS_INTEL : (Optional)  (default: auto)
    - RadonPy_LAMMPS_OPT : (Optional)  (default: auto)
    - RadonPy_JSON_File : (Optional) 
    - RadonPy_Pickle_File : (Optional) 


- 3_tc.py
  - Description: Calculation of thermal conductivity by non-equilibrium MD (NEMD) simulation using the preset `sim.preset.tc`
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_Temp : (Optional) A float value of simulation temperature. (unit: K)
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)
    - RadonPy_LAMMPS_INTEL : (Optional)  (default: auto)
    - RadonPy_LAMMPS_OPT : (Optional)  (default: auto)
    - RadonPy_JSON_File : (Optional) 
    - RadonPy_Pickle_File : (Optional) 
    - RadonPy_TC_Force : (Optional)  (default: False)


- 4_tg.py
  - Description: Calculation of Tg by NEMD simulation using the preset `sim.preset.tg`
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)
    - RadonPy_LAMMPS_INTEL : (Optional)  (default: auto)
    - RadonPy_LAMMPS_OPT : (Optional)  (default: auto)
    - RadonPy_JSON_File : (Optional) 
    - RadonPy_Pickle_File : (Optional) 


- 5_sp.py
  - Description: Calculation of solubility parameters using the preset `sim.preset.sp`
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)
    - RadonPy_LAMMPS_INTEL : (Optional)  (default: auto)
    - RadonPy_LAMMPS_OPT : (Optional)  (default: auto)
    - RadonPy_JSON_File : (Optional) 
    - RadonPy_Pickle_File : (Optional) 


- 6_ef_dp.py
  - Description: Calculation of dynamic dielectric constants and loss tangent by NEMD simulation using `sim.preset.ef_dp`
    - RadonPy_DBID : (Required) An arbitrary string as data ID, which is used in not only data ID but also the name of working directry in calculation.
    - RadonPy_EFDP_Freq : (Optional) 
    - RadonPy_EFDP_EField_Axis : (Optional) 
    - RadonPy_EFDP_Tuning_Max_EField : (Optional) 
    - RadonPy_EFDP_Tuning_Rate : (Optional) 
    - RadonPy_EFDP_Tuning_Plot : (Optional) 
    - RadonPy_EFDP_Number_of_Waves : (Optional) 
    - RadonPy_EFDP_Ensemble : (Optional) 
    - RadonPy_Temp : (Optional) A float value of simulation temperature. (unit: K)
    - RadonPy_Press : (Optional) A float value of simulation pressure. (unit: atm)
    - RadonPy_OMP : (Optional) An int value of pallarel number of OpenMP in LAMMPS. (default: 1)
    - RadonPy_MPI : (Optional) An int value of pallarel number of MPI in LAMMPS. (default: number of available CPU cores)
    - RadonPy_GPU : (Optional) An int value of number of using GPUs in LAMMPS. (default: 0)
    - RadonPy_TMP_Dir : (Optional) Directory path of temporary files (default: working directry in calculation)
    - RadonPy_LAMMPS_INTEL : (Optional)  (default: auto)
    - RadonPy_LAMMPS_OPT : (Optional)  (default: auto)
    - RadonPy_JSON_File : (Optional) 
    - RadonPy_Pickle_File : (Optional) 
    - RadonPy_DP_Force : (Optional)  (default: False)



## Descriptions of output file
| Columns | Description |
| --- | --- |
| UUID | Database ID (UUID) |
| DBID | Database ID (arbitrary ID) |
| monomer_ID | Structure ID of the repeating unit. A unique ID is assigned to each unique chemical structure. If the polymer contains multiple repeating units, such as in a copolymer, multiple IDs are listed as comma-separated values |
| ter_ID_1 | Structure ID of the terminal group on the head side. If `ter_ID_2` is not specified, the same terminal group as `ter_ID_1` is also used on the tail side |
| ter_ID_2 | Structure ID of the terminal group on the tail side |
| smiles_list | SMILES of the repeating unit. If the polymer contains multiple repeating units, such as in a copolymer, multiple SMILES strings are listed as comma-separated values |
| monomer_dir |Directory from which precomputed monomer data were loaded |
| smiles_ter_1 | SMILES of the terminal group on the head side. If `smiles_ter_2` is not specified, the same terminal group as `smiles_ter_1` is also used on the tail side |
| smiles_ter_2 | SMILES of the terminal group on the tail side |
| qm_method | Quantum chemical calculation method used |
| charge | Charge model used for charge calculation |
| copoly_ratio_list | Composition ratio of repeating units in the copolymer |
| copoly_type | Copolymerization type of the copolymer: random copolymerization: `random`; alternating copolymerization: `alternating`; block copolymerization: `block` |
| input_natom | Input setting for the number of atoms in one polymer chain |
| input_nchain | Input setting for the number of polymer chains in the simulation cell |
| ini_density | Input setting for the density of the initial structure |
| temp | Input setting for the MD simulation temperature (K) |
| press | Input setting for the MD simulation pressure (atm) |
| input_tacticity | Input setting for tacticity |
| tacticity | Tacticity of the polymer chain actually generated |
| remarks | Remarks, error messages, etc. |
| Python_ver | Python version in the execution environment |
| RadonPy_ver | RadonPy version in the execution environment |
| RDKit_ver | RDKit version in the execution environment |
| Psi4_ver | Psi4 version in the execution environment |
| LAMMPS_ver | LAMMPS version in the execution environment |
| preset_eq_ver | Version of the equilibration preset used |
| preset_tc_ver | Version of the thermal conductivity calculation preset used |
| check_eq | `True` if the system was judged to have reached equilibrium |
| check_tc | `True` if the validity check for the thermal conductivity calculation was passed |
| do_TC | `True` if the thermal conductivity calculation is allowed to be performed |
| smiles_1 | SMILES of the first repeating unit, corresponding to the order in `monomer_ID` or `smiles_list`. Columns such as `smiles_2` indicate properties of the second repeating unit |
| monomer_ID_1 | Structure ID of the first repeating unit |
| copoly_ratio_1 | Composition ratio of the first repeating unit in the copolymer |
| qm_method_monomer1 | QM calculation method for the first repeating unit |
| charge_monomer1 | Charge calculation method for the first repeating unit |
| remarks_monomer1 | Remarks on the QM calculation of the first repeating unit |
| Python_ver_monomer1 | Python version used for the QM calculation of the first repeating unit |
| RadonPy_ver_monomer1 | RadonPy version used for the QM calculation of the first repeating unit |
| RDKit_ver_monomer1 | RDKit version used for the QM calculation of the first repeating unit |
| Psi4_ver_monomer1 | Psi4 version used for the QM calculation of the first repeating unit |
| LAMMPS_ver_monomer1 | LAMMPS version used for the QM calculation of the first repeating unit |
| mol_weight_monomer1 | Molecular weight of the first repeating unit |
| vdw_volume_monomer1 | van der Waals volume of the first repeating unit (unit: angstroms^3) |
| qm_total_energy_monomer1 | Total energy obtained from the QM calculation of the first repeating unit (unit: hartree) |
| qm_homo_monomer1 | HOMO energy of the first repeating unit (unit: eV) |
| qm_lumo_monomer1 | LUMO energy of the first repeating unit (unit: eV) |
| qm_dipole_monomer1 | The dipole moment of the first repeating unit (unit: debye) |
| qm_dipole_x_monomer1 | x-component of the dipole moment of the first repeating unit (unit: debye) |
| qm_dipole_y_monomer1 | y-component of the dipole moment of the first repeating unit (unit: debye) |
| qm_dipole_z_monomer1 | z-component of the dipole moment of the first repeating unit (unit: debye) |
| qm_polarizability_monomer1 | Polarizability of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_xx_monomer1 | xx-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_yy_monomer1 | yy-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_zz_monomer1 | zz-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_xy_monomer1 | xy-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_xz_monomer1 | xz-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| qm_polarizability_yz_monomer1 | yz-component of the polarizability tensor of the first repeating unit (unit: angstroms^3) |
| DP | Degree of polymerization |
| n_mol | Number of molecules in the simulation cell |
| n_atom | Number of atoms per molecule |
| n_atom_mean | Mean number of atoms per molecule |
| n_atom_var | Variance of the number of atoms per molecule |
| mol_weight | Molecular weight of each molecule |
| Mn | Number-average molecular weight |
| Mw | Weight-average molecular weight |
| Mw/Mn | Ratio of Mw to Mn |
| density | Density (unit: g/cm^3) |
| Rg | Radius of gyration (unit: angstroms) |
| Scaled_Rg | Radius of gyration scaled by molecular weight (Rg/M^0.6) |
| self-diffusion | Self-diffusion coefficient (unit: m^2/s) |
| Cp | Specific heat capacity at constant pressure (unit: J/(kg K)) |
| Cv | Specific heat capacity at constant volume (unit: J/(kg K)) |
| compressibility | Compressibility (unit: 1/Pa) |
| isentropic_compressibility | Isentropic compressibility (unit: 1/Pa) |
| bulk_modulus | Bulk modulus (unit: Pa) |
| isentropic_bulk_modulus | Isentropic bulk modulus (unit: Pa) |
| volume_expansion | Volumetric thermal expansion coefficient (unit: 1/K) |
| linear_expansion | Linear thermal expansion coefficient (unit: 1/K) |
| r2 | Mean-square end-to-end distance (unit: angstroms^2) |
| static_dielectric_const | Static dielectric constant |
| dielectric_const_dc | Static dielectric constant - 1 + refractive index^2 |
| nematic_order_parameter | Orientational order parameter of the repeating unit or molecule |
| refractive_index | Refractive index |
| thermal_conductivity | Thermal conductivity (unit: W/(m K)) |
| thermal_diffusivity | Thermal diffusivity (unit: m^2/s) |
| TC_ke | Component decomposition of thermal conductivity: kinetic energy term (unit: W/(m K)) |
| TC_pe | Component decomposition of thermal conductivity: potential energy term (unit: W/(m K)) |
| TC_pair | Component decomposition of thermal conductivity: electrostatic, excluding long-range contribution, and vdW term (unit: W/(m K)) |
| TC_bond | Component decomposition of thermal conductivity: bond term (unit: W/(m K)) |
| TC_angle | Component decomposition of thermal conductivity: angle term (unit: W/(m K)) |
| TC_dihed | Component decomposition of thermal conductivity: dihedral term (unit: W/(m K)) |
| TC_improper | Component decomposition of thermal conductivity: improper term (unit: W/(m K)) |
| TC_kspace | Component decomposition of thermal conductivity: long-range electrostatic term (unit: W/(m K)) |
| preset_tg_ver | Version of the equilibration preset used |
| tg_next_temp |  |
| tg_init_density |  |
| tg_init_density_check |  |
| tg_max_temp |  |
| tg_min_temp |  |
| tg_cooling_rate |  |
| tg_interval_temp |  |
| tg | Glass transition temperature (unit: K) |
| tg_rmse |  |
| tg_thermal_expansion_coef(upper_tg) |  |
| tg_thermal_expansion_intercept(upper_tg) |  |
| tg_thermal_expansion_coef(below_tg) |  |
| tg_thermal_expansion_intercept(below_tg) |  |
| sp_ced | Cohesive energy density (unit: J cm^-3) |
| sp_total | Hildebrand solubility parameter (unit: MPa^0.5) |
| sp_vdw |Dispersion-force contribution to the solubility parameter (unit: MPa^0.5) |
| sp_ele | Electrostatic interaction contribution to the solubility parameter, corresponding to the sum of the polar and hydrogen-bonding contributions in the Hansen solubility parameters (unit: MPa^0.5) |
| sp_ele_short | Short-range component of the electrostatic interaction contribution to the solubility parameter, calculated as the pairwise Coulomb interaction term in the MD calculation (unit: MPa^0.5) |
| sp_ele_long | Long-range component of the electrostatic interaction contribution to the solubility parameter, calculated using PPPM in the MD calculation (unit: MPa^0.5) |



