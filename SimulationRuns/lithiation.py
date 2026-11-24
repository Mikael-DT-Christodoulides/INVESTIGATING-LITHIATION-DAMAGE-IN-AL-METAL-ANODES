# License: MIT 
# This script is intended for educational and research purposes, and may be freely modified and used for (non-)commercial purposes

# Artificial intelligence tools were not used to generate the script below, and it is the intellectual property of the author. 
# For any papers or publications that use this script, please cite the original paper, listed below.

# "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", M. Christodoulides et al., [In submission 2026, DOI tbc]

###################################################################################################################################################################

# Lithiation Script
# Mikael D.T. Christodoulides | 6 Nov 2026
# V1.0.3

from lammps import PyLammps
from tqdm import tqdm

L = PyLammps(cmdargs=["-screen", "none", "-log", "none"])
# remove logs due to tqdm, can change if needed

########################################### Script Parameters ###########################################

coding = 0                                           # set to 0, 1, 2 to run normal, skip minimization or skip Li loop (used while developing)
loopbreak = 10                                       # hard stop if loop gets stuck running infinitely (will insert cycles*loopbreak atoms)
cyclecount = 1000                                    # Li atoms per cycle (used at 1000 throughout paper)
liCeil = 3000                                       # Li atoms to add total

L.command("variable eTol equal 1.0e-5")              # energy tolerance for minimization (heavily affects runtime)
L.command("variable fTol equal 1.0e-7")              # force tolerance for minimization (should not be reached)
L.command("variable maxIter equal 500")              # max iterations for minimization (heavily affects runtime if energy tolerance not reached)

L.command("variable midway equal 30")                # define as the midway point on the z axis (depends on box dimensions)
L.command("variable baselen equal 10")               # the width of the base layer (must be <= length of core, 10A works well)

L.command("variable temp equal 300")                 # relaxation temp, room temp
L.command("variable engoffset equal 0.15")           # eV over minimum for Al to count as "high-energy" (forgiving by design, too high will limit insertion)
# furthermore, this value is not of immense importance if not directly studying the surface layer (which is advised)

L.command("variable depth equal 1.5")                # z offset for Li atom insertion
L.command("variable rsphere equal 3.0")              # radius of the insertion sphere
L.command("variable overlapdist equal 1.8")          # minimum Li-Al distance during insertion (larger values may cause issues with insertion)
L.command("variable mdsteps equal 300")              # MD relax steps after Li insertion (heavily affects runtime)



########################################### Basic atomic setup ###########################################

L.command("units metal")                            # sets internal conversions
L.command("atom_style atomic")                      # default atom style
L.command("boundary s s m")                         # shrink wrapped in x and y to allow for expansion, shrink wrapped with min in z to keep dimension

L.command("read_data Al_sample_input.lmp extra/atom/types 2")
# read the Atomsk input files, defines two extra types for Li and seeds
# types: Al=1, Li=2, seed=3

L.command("mass 2 6.940")                            # molar mass Li
L.command("mass 3 26.980")                           # molar mass Al (for seeds)


L.command("pair_style meam")                         # pairwise interactions computation (modified embedded atom method)

L.command("pair_coeff * * library.meam Al Li AlLi.meam Al Li Al")
# uses published potentials from the paper below
# ' S. Roy, A. Dutta, and N. Chakraborti, "A novel method of determining interatomic potential for Al and Al-Li alloys and studying strength of Al-Al" '



# ########################################### DELETE STRAY ATOMS ###########################################
# L.command("region keep block -50 350 -50 350 -30 120 units box")              # dimensions depend on file
# L.command("group keeping region keep")
# L.command("group strayatoms subtract all keeping")                            # complement of keeping group
# L.command("variable nstray equal count(strayatoms)")                          # count of keeping group
# L.command("delete_atoms group strayatoms")                                    # remove atoms (if deleting too many at once, may cause problems)

# L.command("group strayatoms delete")
# L.command("group keeping delete")
# print("Deleted", int(L.variables["nstray"].value), "stray atoms")
# ##########################################################################################################

# ^ this code block should be included if the inputted file contains atoms outside the bulk which may cause simulation issues



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

L.command("group Li empty")
# note this should be change if the file contains Li atoms which should be counted towards the total

L.command("group energetic empty")
L.command("group seeds empty")
# more groups for the Li atoms, high energy atoms and seeds 



########################################### Minimize structure (initial/pre-Li) ###########################################

L.command("thermo 10")
L.command("thermo_style custom step pe fmax")

if coding == 0:
    L.command("min_style cg")
    L.command("minimize ${eTol} ${fTol} ${maxIter} 5000")





########################################### Simulation Setup (freeze base, etc.) ###########################################

L.command("timestep 0.001")                          # time is default to ps so this is one femtosecond
L.command("neighbor 2.0 bin")                        # bins neighbor list with a parameter of 2 angstroms
L.command("neigh_modify every 1 delay 0 check yes")  # builds a neighbor list every step, no delay until start and check is the criterion below:
# "only build if at least one atom has moved half the skin distance or more"

L.command("fix md all nve")                                     # plain time integration for position/velocity
L.command("fix therm all langevin ${temp} ${temp} 0.1 12345")   # Langevin thermostat, 0.1ps damping, temp=RT                  
# Together, these make for a Brownian dynamics setup

L.command("fix anchor base setforce 0.0 0.0 0.0")
L.command("velocity base set 0.0 0.0 0.0")
# set F=0N and v=0m/s to all atoms in the base (fix in place)






########################################### Calculated figures ########################################### 

L.command("compute poteng all pe/atom")                     # compute PE
L.command("compute coremineng core reduce min c_poteng")    # computes the lowest PE in the core

L.command("variable nLi equal count(Li)")                   # number of Li atoms (should be 0)
L.command(f"variable nLiMax equal {liCeil}")                # number of Li atoms desired total

L.command('''variable energeticatom atom "type == 1 && gmask(surface) && c_poteng > v_engthresh"''')
# variable declaration for high-energy sites, needed below


nLi = L.variables["nLi"].value
nLiMax = L.variables["nLiMax"].value
# python vars


########################################### Lithiation (Li insertion) ###########################################

with tqdm(total=(nLiMax), desc="Adding Li atoms...") as pbar: # progress bar (this method takes awhile)

    if coding != 2:
        # each cycle follows the following pattern
        # 1. compute energy
        # 2. pick seed sites
        # 3. insert Li near seeds
        # 4. minimize
        i = 0

        while (nLi < nLiMax) and (i < loopbreak):
            i+= 1                                   # i += 1 (loop var)                

            L.command("run 0 post no")
            L.command("variable engthresh equal $(c_coremineng+v_engoffset)")
            # energy threshhold for high-energy Al atoms

            L.command("group energetic clear")
            L.command("group energetic variable energeticatom")
            # find 'high energy' atoms

            L.command("variable nEnerg equal count(energetic)")
            nCount = L.variables["nEnerg"].value

            if nCount == 0:
                print('No high-energy sites left, aborting')
                i = loopbreak
                continue

            print(f"\nCycle #{i}: candidates: {nCount}, Li thus far: {int(L.variables['nLi'].value)}")
            # logs


            for j in range(0, cyclecount):
                randseed = i * 673 + j
                L.command(f"variable randseed equal {randseed}")

                L.command("set group energetic type/subset 3 1 ${randseed}")
                # picks a random energetic atom and marks it as a seed

                L.command("group seeds clear")
                L.command("group seeds type 3")
                # updates seeds group

                L.command("variable xs equal bound(seeds,xmin)")
                L.command("variable ys equal bound(seeds,ymin)")
                L.command("variable zs equal bound(seeds,zmin)")
                # extracts position of seed

                L.command("set group seeds type 1")

                L.command("region site sphere ${xs} ${ys} $(v_zs-v_depth) ${rsphere} units box")
                # define region for insertion sphere based on seed

                before = L.lmp.get_natoms()
                L.command("create_atoms 2 random 1 ${randseed} site overlap ${overlapdist} maxtry 1000")
                # insert Li atom into a random position in this region, trying 1000 times
                
                L.command("region site delete")
                

            L.command("group Li type 2")
            L.command("group surface union surface Li")

            if coding == 0:
                L.command("minimize ${eTol} ${fTol} ${maxIter} 10000")
                print("minimized")
                L.command("run ${mdsteps}")
                print("run")
                L.command("minimize ${eTol} ${fTol} ${maxIter} 10000")
                print("minimized again")

            nLi = L.variables["nLi"].value
            nLiMax = L.variables["nLiMax"].value
            # update Python vars for loop condition

            pbar.n = int(nLi)
            pbar.refresh()
            # update tqdm bar




########################################### Outputs ###########################################
L.command("thermo_style custom step temp pe c_coremineng v_nLi")
L.command("dump 1 all custom 1 Al_lithiated.dump id type x y z v_is_surface c_poteng")
L.command("dump_modify 1 sort id")
L.command("run 0")
L.command("write_data Al_lithiated.data")
###############################################################################################
