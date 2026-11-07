# License: MIT 
# This script is intended for educational and research purposes, and may be freely modified and used for (non-)commercial purposes

# Artificial intelligence tools were not used to generate the script below, and it is the intellectual property of the author. 
# For any papers or publications that use this script, please cite the original paper, listed below.

# "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", M. Christodoulides et al., [In submission 2026, DOI tbc]

###################################################################################################################################################################

# Minimization Script
# Mikael D.T. Christodoulides | 6 Nov 2026
# V1.0.0

from lammps import PyLammps
from tqdm import tqdm

L = PyLammps(cmdargs=["-screen", "none", "-log", "none"])
# remove logs due to tqdm, can change if needed


########################################### Script Parameters ###########################################

L.command("variable eTol equal 5.0e-5")              # energy tolerance for minimization (heavily affects runtime)
L.command("variable fTol equal 5.0e-7")              # force tolerance for minimization (should not be reached)
L.command("variable maxIter equal 500")              # max iterations for minimization (heavily affects runtime if energy tolerance not reached)

L.command("variable midway equal 30")                # define as the midway point on the z axis (depends on box dimensions)
L.command("variable baselen equal 10")               # the width of the base layer (must be <= length of core, 10A works well)

L.command("variable temp equal 300")                 # relaxation temp, room temp

L.command("variable mdsteps equal 500")              # MD relax steps after Li insertion (heavily affects runtime)



########################################### Basic setup ###########################################

L.command("units metal")                             # sets internal conversions
L.command("atom_style atomic")                       # default atom style
L.command("boundary s s m")                          # shrink wrapped in x and y to allow for expansion, shrink wrapped with min in z to keep dimensions

L.command("read_data Al_delithiated.data extra/atom/types 2")
# read Atomsk input file, defines two extra types for Li and seeds
# types: Al=1, Li=2, seed=3

L.command("mass 2 6.940")                            # molar mass Li
L.command("mass 3 26.980")                           # molar mass Al (for seeds)

L.command("pair_style meam")                         # pairwise interactions computation (modified embedded atom method)

L.command("pair_coeff * * library.meam Al Li AlLi.meam Al Li Al")
# uses published potentials from the paper below
# ' S. Roy, A. Dutta, and N. Chakraborti, "A novel method of determining interatomic potential for Al and Al-Li alloys and studying strength of Al-Al" '


########################################### Regions/groups ###########################################

L.command("region reg_core block INF INF INF INF INF ${midway}")
L.command("region reg_surface block INF INF INF INF ${midway} INF")
L.command("region reg_base block INF INF INF INF INF ${baselen}")
# set the surface, core and base regions, infinite in x,y and limited in z
# note the base region is a tiny fixed end layer at the bottomof the core

L.command("group core region reg_core")
L.command("group surface region reg_surface")
L.command("group base region reg_base")
# groups the atoms in these regions together

L.command("variable is_surface atom gmask(surface)")
L.command("variable is_base atom gmask(base)")
# sets layer variables for visualizations in OVITO (is_base not outputted, can be changed if desired)




########################################### Minimize structure ###########################################

L.command("thermo 10")
L.command("thermo_style custom step pe fmax")


L.command("min_style cg")
L.command("minimize ${eTol} ${fTol} ${maxIter} 5000")



########################################### Outputs ###########################################
L.command("dump 1 all custom 1 Al_minimized.dump id type x y z v_is_surface")
L.command("dump_modify 1 sort id")
L.command("run 0")
L.command("write_data Al_minimized.data")

###############################################################################################


