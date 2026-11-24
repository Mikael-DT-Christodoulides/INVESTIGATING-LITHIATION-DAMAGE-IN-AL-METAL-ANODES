# License: MIT 
# This script is intended for educational and research purposes, and may be freely modified and used for (non-)commercial purposes

# Artificial intelligence tools were not used to generate the script below, and it is the intellectual property of the author. 
# For any papers or publications that use this script, please cite the original paper, listed below.

# "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", M. Christodoulides et al., [In submission 2026, DOI tbc]

###################################################################################################################################################################

# MASTER SCRIPT
# Mikael D.T. Christodoulides
# V1.0.1 - 16th September 2026

# IMPORTS
import numpy as np
import ast
import os
from scipy.spatial.transform import Rotation
from Preprocess import organizefile
from Definitions import O, omega, get_kernel, calc_GD
from tqdm import tqdm
import math

fname = input("Please enter the input filename (without extension): ")
alpha = float(input("Please enter an alpha value for lambda: "))
filename = f"inputs/{fname}" # we keep this var for later methods

print("Welcome to the Damage Metric 3000!\nBeginning cutting edge research...")
if (os.path.exists(filename + "_MAT.txt") == False): # check if rotation file exists (default in quaternions)
    print("Rotation matrix file not found, creating one...")
    sname = organizefile(filename)
else:
    print("Rotation file found, proceeding...")
    sname = filename + "_MAT.txt"
 
#NB: This now uses the rotation file, each atom is stored as "(posx, posy, poz) grainID [11 12 13 21 22 23 31 32 33 <rot matrix>, surface]"

file = open(sname, "r")
byte_len = os.path.getsize(sname) # for progress bar

totalCount = 0 # atom count
atoms = {} # list of all simulated atoms/unit cells indexed by unique ID (we treat atoms as unit cells for OVITO compatability simplicity)
posDict = {} # dictionary to get ID from position tuple
grainDict = {} # dictionary to get ID for a particular grain
uID = -1 # unique ID for each atom (0 indexed to match arrays)

# extract each atom and store it in our dictionaries
with tqdm(total=byte_len, unit="B", unit_scale=True, desc="Processing atoms...") as pbar:
    for i, line in enumerate(file):
        pbar.update(len(line.encode('utf-8')))  # Update progress bar based on the byte length of the line
        if (i < 3):
            continue # skip the first *three* lines of the file (headers)

        sPos, sId, sMat, sSurface = line.split('\t')
        gId = int(sId)
        posx, posy, posz = map(float, sPos[1:-1].split(":")) # position triplet
        mat = np.array(ast.literal_eval(sMat), dtype=float).reshape(3, 3) # rotation matrix
        surface = 1 if (sSurface[:-1]=="True") else 0 # surface var

        uID += 1 # unique ID

        atoms[uID] = [(posx, posy, posz), gId, mat, surface] # master dictionary for indexing by ID <NOTE>: This matches our line number, important for later
        posDict[(posx, posy, posz)] = uID # dict for pos -> ID
        grainDict[gId] = grainDict.get(gId, []) + [uID] # dict for grain -> ID list


# all atoms are now stored in the following data structure, inspired by a relational database (primary/foreign keys)
    # a master dict with each atom position, grainID and rotation matrix keyed by ID
    # translation dicts, for converting positions to IDs and grainID to a list of atom IDs for that grain

file.close()



coreatoms = {key: value for key, value in atoms.items() if value[3] == 0} # the atoms we want to include in GD, damage calculations

print(f"{len(coreatoms)} core atoms found our of {len(atoms)} total.")


print("Data organized, beginning calculations.\nComputing global distortion...")
GD = -1
if (os.path.exists(filename + "_GD.txt") == False): # check if GD file exists
    print("No stored value for global distortion, computing...")
    GD = calc_GD(coreatoms, grainDict, filename)
else:
    print("Global distortion value found, reading...")
    with open(filename + "_GD.txt", "r") as f:
        GD = float(f.readline())
        print(f"GD value of {GD} read from previous run.")


import pickle
omega_cache = filename + "_OMEGA.pkl"
table = {}

if (os.path.exists(omega_cache)):
    print("Distortion cache found, loading...")
    with open(omega_cache, "rb") as f:
        table = pickle.load(f)
else:
    print("No distortion cache found, computing...")
    for i in tqdm(coreatoms.keys(), desc="Caching omega values..."):
        for j in get_kernel(i, filename, atoms, posDict):
            key = (i, j) if i < j else (j, i)
            # Check if the omega value is already cached
            if key in table:
                continue # already computed, skip

            p1 = atoms[key[0]][2]
            p2 = atoms[key[1]][2]
            dmg, p2f = omega(p1, p2)
            table[key] = dmg  # store computed omega value in a dictionary

    with open(omega_cache, "wb") as f:
        pickle.dump(table, f)


dmgvals = {}

with tqdm(total=len(coreatoms.keys()), desc="Computing local distortion (chi)") as pbar:
    for id in atoms.keys(): # compute chi(p) for every CORE atom in the simulation

        if (id not in coreatoms): # this is a surface atom
            dmgvals[id] = -1        # we store a dummy value, so our indices are kept for our final output file (needs to match input, which contained surface atoms)
            continue

        x = 0
        x_al = 0
        N = 0
        nAliens = 0
        grain = atoms[id][1] # grain of atom
        kernel = get_kernel(id, filename, atoms, posDict)
        for i in kernel:

            ##################################### Adding lambda logic #####################################
            gId = atoms[i][1] # grains of kernel
            if (gId != grain): # alien
                akey = (id, i) if id < i else (i, id)
                x_al += table[akey]  # use cached omega value
                nAliens += 1
            ##################################### Adding lambda logic #####################################

            key = (id, i) if id < i else (i, id)
            x += table[key]  # use cached omega value
            N += 1

            # p1 = atoms[id][2] # rotation matrix 1 (should be a CORE atom)
            # p2 = atoms[i][2] # rotation matrix 2 (might be a surface atom, intentional due to the interface)
            # dmg, p2f = omega_lookup(p1, p2)
            # x += dmg
            
        if (N == 0):
            X = 0 # isolated atom, not relevant to us, but we need to keep it for indexing
            L = 0
        else:
            X = x / N # local distortion value

            ### LAMBDA CHANGE ###
            L = (alpha)*(((-1)*(nAliens)*(x-x_al)/((N)*(N-nAliens))) + (x_al/N)) # lambda value 
            #####################

        dmg = math.log(((X-L)*GD + 1), 10) # as per paper's definition
        if id in dmgvals:
            print("Critical error: overwriting a damage value!")
            exit() 
        dmgvals[id] = dmg
        pbar.update(1)


with open(f"{filename}.txt", "r") as infile, open(f"{filename}_RESULTS_{str(alpha)}.txt", "w") as outfile:
    for i, line in enumerate(infile, start=-2):
        if (i == -2): # first 2 lines headers
            outfile.write(line) # write as-is
            continue 
        if (i == -1):
            outfile.write(line[:-1] + ":Damage:R:1" + "\n") # add new var for damage
            continue 

        data = line.strip().split()
        uId = posDict[tuple(map(float, data[0:3]))] # read position of fileline, and convert to ID
        # <NOTE:> In theory uID = 1, but this was making some strange results
        # Note from the future: this was foreshadowing, the entire script didn't work because of this. It has since been fixed.

        outfile.write(line[:-1] + " " + str(dmgvals[uId]) + "\n")
