# License: MIT 
# This script is intended for educational and research purposes, and may be freely modified and used for (non-)commercial purposes

# Artificial intelligence tools were not used to generate the script below, and it is the intellectual property of the author. 
# For any papers or publications that use this script, please cite the original paper, listed below.

# "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", M. Christodoulides et al., [In submission 2026, DOI tbc]

###################################################################################################################################################################

# Delithiation Script
# Mikael D.T. Christodoulides | 6 Nov 2026
# V1.0.1

from lammps import PyLammps

L = PyLammps(cmdargs=["-screen", "none", "-log", "none"])
print("Starting delithiation process")



########################################### Basic atomic setup ###########################################

L.command("units metal")                             # sets internal conversions
L.command("atom_style atomic")                       # default atom style
L.command("boundary s s s")                          # shrink wrapped in x and y to allow for expansion, shrink wrapped with min in z to keep dimension

L.command("read_data Al_lithiated.data")             # read lithiated input

L.command("pair_style meam")                         # pairwise interactions computation (modified embedded atom method)
L.command("pair_coeff * * library.meam Al Li AlLi.meam Al Li Al")
# uses published potentials from the paper below
# ' S. Roy, A. Dutta, and N. Chakraborti, "A novel method of determining interatomic potential for Al and Al-Li alloys and studying strength of Al-Al" '





########################################### Regions/groups ###########################################

L.command("region reg_surface block INF INF INF INF 30 INF")    # midway is 30
L.command("group surface region reg_surface")
L.command("variable is_surface atom gmask(surface)")            # surface atoms, defined once

L.command("region reg_base block INF INF INF INF INF 10")
L.command("group base region reg_base")                         # groups the atoms in base region encore pour la meme raison qu'avant

print("Groups defined")



########################################### Simulation Setup (freeze base, etc.) ###########################################

L.command("timestep 0.001")                        # time is default to ps so this is one femtosecond
L.command("fix anchor base setforce 0.0 0.0 0.0")
L.command("velocity base set 0.0 0.0 0.0")
# set F=0N and v=0m/s to all atoms in the base (fix them in place)

L.command("neighbor 2.0 bin")                        # bins neighbor list with a parameter of 2 angstroms
L.command("neigh_modify every 1 delay 0 check yes")  # builds a neighbor list every step, no delay until start and check is the criterion below:
# "only build if at least one atom has moved half the skin distance or more"





########################################### Relaxation ###########################################
L.command("group li type 2")
L.command("variable nLi equal count(li)")    
nRem = L.variables["nLi"].value
# Li atoms tracker


print("Initializing delithiation")

z = 30 # starting z position for delithiation
zInc = 1 # z increment for delithiation, starting at core/surface interface and moving up
firstPass = True



########################################### Clear atoms not in the solid mass ###########################################
L.command("region bulk block INF INF INF INF -25 95 side out") # the complement of our sample space
L.command("group bulk region bulk") # group the atoms outside of the sample space
L.command("variable nBulk equal count(bulk)") # count the number of atoms outside of the sample space
print("Deleting ", L.variables["nBulk"].value, " atoms outside of the sample space")
L.command("delete_atoms region bulk") # delete atoms outside of the sample space
#########################################################################################################################


L.command(f"region reg_delith block INF INF INF INF INF {z}")
L.command("group delith region reg_delith")

L.command("group void intersect li delith") # Li atoms in the delithiation region
L.command("variable nVoid equal count(void)") # number of such atoms

print("Beginning delithiation process...")
while nRem > 0 and z < 95:  # 95 is safety net, should be done by then but depends on sample size (obviously)
    ########################################### Delithiation logic ###########################################
    if (not firstPass): # excludes first run
        #### Clear atoms since expelled from the solid mass ####
        L.command("region bulk delete") # delete the previous bulk region
        L.command("group bulk delete") # delete the previous bulk group

        L.command("region bulk block INF INF INF INF -25 95 side out") # the complement of our sample space
        L.command("group bulk region bulk") # group the atoms outside of the sample space
        L.command("variable nBulk equal count(bulk)") # count the number of atoms outside of the sample space
        print("Deleting ", L.variables["nBulk"].value, " atoms outside of the sample space")
        L.command("delete_atoms region bulk") # delete atoms outside of the sample space
        ##################################################


    print("Void contains ", L.variables["nVoid"].value, " Li atoms, ~", L.variables["nVoid"].value%1000, "cycle(s) before z increase.")

    # delete 1000 random atoms in void
    L.command("delete_atoms random count 1000 no void NULL 43678") # surfcce set to NULL so only impose group criterion
    L.command("group mobile subtract all base") # group all atoms except the base for relaxation

    # Relax (minimization not too harsh to preserve damage sites from lithiation)
    L.command("fix md all nve")

    L.command("fix therm mobile langevin 600.0 600.0 0.1 24688")
    L.command("run 100")
    L.command("fix therm mobile langevin 300.0 300.0 0.1 24688")
    L.command("run 100")
    # <NOTE>: These runs at elevated temperature allows for atom movement, as desired, but as such cause the sample to become more amorphous. This is an almost unavoidable
    # effect (we cannot have the atoms move and relax from their strained delithiation positions while also preserving crystallinity), but it should be noted.

    L.command("unfix md")
    L.command("unfix therm")

    # quick minimize after deletion and MD runs
    L.command("minimize 1e-5 1e-7 200 10000")


    L.command("variable nVoid equal count(void)") # new count
    L.command("variable nLi equal count(li)") # new Li
    nVoid = L.variables["nVoid"].value
    nRem = L.variables["nLi"].value

    print("Remaining Li atoms: ", nRem, " remaining in void: ", nVoid)

    if (nVoid == 0): # out of Li atoms to delete, need to update z
        print("No Li atoms remaining in the delithiation region, increasing z...")
        L.command("run 0 post no")
        L.command(f"write_dump all custom Al_stg7_z{z}.dump id type x y z v_is_surface") # write output for current z
        z += zInc

        print("New delithiation region: z = ", z, (z/60)*100, "% complete")

        L.command("region reg_delith delete") # delete the previous delithiation region
        L.command(f"region reg_delith block INF INF INF INF INF {z}")
        L.command("group delith region reg_delith")

        L.command("group void intersect li delith") # Li atoms in the delithiation region
        L.command("variable nVoid equal count(void)") # number of such atoms

    if (firstPass):
        firstPass = False # now on second iteration



########################################### Output ###########################################

L.command("write_data Al_relaxed.data")
# data file since dump files produced in loop above

###############################################################################################

