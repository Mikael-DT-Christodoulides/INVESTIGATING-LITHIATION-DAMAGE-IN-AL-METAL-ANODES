# PREPROCESSING SCRIPT
# Mikael D.T. Christodoulides 
# V1.0.1 - 16th September 2026

# IMPORTS
import numpy as np
from scipy.spatial.transform import Rotation
from Definitions import O, omega, central_pixel 
from tqdm import tqdm

def organizefile(filename):
    """Processes data input to create a rotation-matrix-based dictionary"""
    file = open(filename + ".txt", "r")
    #NB: Each atom is initially stored as [posx, posy, posz, orientw, orientx, orienty, orientz, grainID]

    atoms = {} # list of all simulated atoms/unit cells (we treat atoms as unit cells for OVITO simplicity)

    # extract each atom's line and convert it into a rotation matrix and coordinate triple
    layered = False # is the simulation split into a surface and a core? 
    for i, line in enumerate(file):
        if (i < 2):
            continue # skip the first two lines of the file (headers)

        data = line.strip().split()

        # if we have 9 rows, we need to not which are non-surface atoms (which is column 9's value)
        # if we have 8 rows, we're assuming we have pure Al
        if (i == 2): # check first row for format
            layered = (len(data) == 9)
            if (layered):
                print("Surface layer found, noting...")

        surface = False # by default, everything is a non-surface atom
        if (layered):
            surface = (float(data[8]) == 1)

        grainID = int(data[7]) # 8th entry is always grainID, <NOTE>: Don't index negatively (len changes)
        qw, qx, qy, qz = [float(x) for x in data[3:7]]

        if np.sqrt(qw**2 + qx**2 + qy**2 + qz**2) == 0: # zero quaternion
            R = np.eye(3)  # identity rotation
        else:
            R = Rotation.from_quat([qx, qy, qz, qw]).as_matrix()

        atoms[grainID] = atoms.get(grainID, []) + [[tuple(map(float, data[0:3])), R, surface]]
        # this is a master array, contains everything. BUT, when we work out GD, damage vals, we're going
        # to focus on non-surface atoms. That said, we need the surface atoms since they might be neighbours
        # at the interface

    #NB: Each atom is now stored as [(posx, posy, posz), rotation_matrix, surface] in the grains dictionary
    file.close()



    # # <NOTE>: debugging
    # file = open(filename + "checkingsomething.txt", "w")
    # file.write("8529")
    # file.write('''Lattice="67.6962947206 0.0 0.0 0.0 189.8402592401 0.0 0.0 0.0 87.9182135676" Origin="-5.6912447096 -45.6811732596 -16.9234696204" Properties=pos:R:3:Grain:I:1:v_is_surface:R:1\n''')
    # for grain in atoms:
    #     for itatom in atoms[grain]: #for each grain
    #         x, y, z = itatom[0]
    #         surface = 1 if itatom[2] else 0
    #         file.write(f"{x} {y} {z} {grain} {surface}\n")

    # print("alr theres ur file for chekcing")
    # file.close()
    # exit()
    # debugging



    atomCount = 0
    print("Completed initial parsing \nGrain ID\tAtom Count")
    for grainID, atom in atoms.items():
        print(f"{grainID}\t\t{len(atom)}")
        atomCount += len(atom)

    # rebuild formatted input file
    outStr = "Atomic Simulation Results - Matrix Version\nPosition\t\t\tGrainID\t\t\tRotation\t\t\tSurface\n________________________________________________________________________________________________________________________________________________________________________________________________________________________________\n"

    with tqdm(total=atomCount, desc="Building rotation matrix file") as pbar: # progress bar (this method takes awhile)
        for grainID, grainAtoms in atoms.items():
            for i, atom in enumerate(grainAtoms): # for every atom in the sim
                pos = "("+str(atom[0][0]) + ":" + str(atom[0][1]) + ":" + str(atom[0][2])+")"
                rotation = "[" + ",".join(map(str, atom[1].flatten())) + "]"
                surface = str(atom[2]) # boolean value for surface atoms
                outStr = outStr + pos + "\t" + str(grainID) + "\t" +  rotation + "\t" + surface + "\n"
                pbar.update(1)

    sname = filename + "_MAT.txt"

    with open(sname, "w") as f:
        f.write(outStr)

    f.close()

    print("Data processed, written to", sname)
    return sname