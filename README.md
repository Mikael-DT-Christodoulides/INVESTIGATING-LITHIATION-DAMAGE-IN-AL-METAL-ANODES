# Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach #
Software, scripts and data developed during the creation of the paper "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach".

**Contents:**
_SimulationRuns_
- lithiation.py: _Python-wrapped LAMMPS simulation for the insertion of Li atoms into a defined region of an aluminium sample_
- delithiation.py: _Python-wrapped LAMMPS simulation for the removal of Li atoms from an aluminium sample_
- relax.py: _Simple minimization script to relax a sample of aluminiuim after removal of Li atoms_
- AlLi.meam, library.meam: _MEAM potentials files developed Roy et al. and hosted on OpenKIM below: https://openkim.org/id/MEAM_LAMMPS_RoyDuttaChakraborti_2021_AlLi__MO_971738391444_000_
- Al_input.zip: _Unlithiated LAMMPS input file created using Atomsk for simulation runs L1-D_
- Definitions.py, Master.py, Preprocess.py: _Python stack for analysing LAMMPS output data using the disorder map, X._
- DMGMetric_MATLAB.m: _MATLAB Script for analysing EBSD output data using the disorder map, X._
  
_Data_
- StressStrainT4.txt: _Stress-strain data for tensile sample T4_
- TensileDmg.txt: _Disorder (X) output data for tensile samples T1-T3_
- SimDmgData.zip: _Disorder (X) output data and z positions for simulation runs L1-D_


Supplementary results, including stress-strain data, IPF maps and damage heatmaps for tensile tests, are available upon reasonable request to the corresponding author: chrmik001@myuct.ac.za

All scripts in this repository are open source and freely available for use, modification, or redistribution (see LICENSE).
Where the scripts in this repository are used, citations of the paper below are greatly appreciated:

- M. Christodoulides et al., "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", [_Submission in progress_]
