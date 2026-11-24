# License: MIT 
# This script is intended for educational and research purposes, and may be freely modified and used for (non-)commercial purposes

# Artificial intelligence tools were not used to generate the script below, and it is the intellectual property of the author. 
# For any papers or publications that use this script, please cite the original paper, listed below.

# "Towards Spatially Mapping the Degradation of Aluminium Anodes: A Group Theoretic Approach", M. Christodoulides et al., [In submission 2026, DOI tbc]

###################################################################################################################################################################

# DEFINITIONS & FUNCTIONS
# Mikael D.T. Christodoulides 
# V1.0.1 - 16th September 2026

import numpy as np
from scipy.spatial.transform import Rotation
from itertools import product
from ovito.io import import_file
from ovito.data import NearestNeighborFinder, CutoffNeighborFinder
from tqdm import tqdm


# CONSTANTS
O = [np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]),
    np.array([[ 1, 0, 0], [ 0, 0,-1], [ 0, 1, 0]]),
    np.array([[ 1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]),
    np.array([[ 1, 0, 0], [ 0, 0, 1], [ 0,-1, 0]]),
    np.array([[ 0, 0, 1], [ 0, 1, 0], [-1, 0, 0]]),
    np.array([[-1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]),
    np.array([[ 0, 0,-1], [ 0, 1, 0], [ 1, 0, 0]]),
    np.array([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0, 1]]),
    np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0, 1]]),
    np.array([[ 0, 1, 0], [-1, 0, 0], [ 0, 0, 1]]),
    np.array([[ 0, 0, 1], [ 1, 0, 0], [ 0, 1, 0]]),
    np.array([[ 0, 1, 0], [ 0, 0, 1], [ 1, 0, 0]]),
    np.array([[ 0, 0,-1], [ 1, 0, 0], [ 0,-1, 0]]),
    np.array([[ 0, 1, 0], [ 0, 0,-1], [-1, 0, 0]]),
    np.array([[ 0, 0, 1], [-1, 0, 0], [ 0,-1, 0]]),
    np.array([[ 0,-1, 0], [ 0, 0, 1], [-1, 0, 0]]),
    np.array([[ 0, 0,-1], [-1, 0, 0], [ 0, 1, 0]]),
    np.array([[ 0,-1, 0], [ 0, 0,-1], [ 1, 0, 0]]),
    np.array([[ 0, 1, 0], [ 1, 0, 0], [ 0, 0,-1]]),
    np.array([[ 0,-1, 0], [-1, 0, 0], [ 0, 0,-1]]),
    np.array([[ 0, 0, 1], [ 0,-1, 0], [ 1, 0, 0]]),
    np.array([[ 0, 0,-1], [ 0,-1, 0], [-1, 0, 0]]),
    np.array([[-1, 0, 0], [ 0, 0, 1], [ 0, 1, 0]]),
    np.array([[-1, 0, 0], [ 0, 0,-1], [ 0,-1, 0]])]


# the 26 offset positions for a given atom (i.e. kernel coordinate set) excluding the origin
KERNEL_OFFSETS = [o for o in product((-1, 0, 1), repeat=3) if o != (0, 0, 0)]

# DEFINITIONS
def omega(p1, p2):
    """ Returns the distortion beween p1O and p2O, the orientations of these matrices"""
    # to determine omega, we only need loop over O for one of the orientations
    
    w = 360 # initial min value 
    p2f = p2;
    for o in O:
        p2O = p2 @ o
        i = np.clip(((np.trace(p1.T @ p2O)-1) * 0.5), -1.0, 1.0) # clip rounding errors (e.g. -1.0000000000000002)
        i = np.degrees(np.arccos(i))

        if (i < w):
            w = i
            p2f = p2O # neighbour orientation

    return [round(w, 4), p2f]


def central_pixel(ids, atoms):
    """Returns some x_G in phi(p_G) for a given grain (an average orientation matrix)"""
    arbInd = ids[0] # arbitrary atom in our grain
    x = atoms[arbInd][2]  # arbitrary rotation matrix from this atom

    x_ks = [x] # list of matrices and neighbours
    n = 0

    for id in ids:
        if (id == arbInd):
            continue  # skip first one (already included)

        x_k = atoms[id][2] # rotation matrix

        # obtain neighbour matrices
        w, x_k = omega(x, x_k)
        x_ks.append(x_k)
        n += 1

    x_ave = sum(x_ks) / n
    # average 'rotation' matrix, obtained from elementwise average

    u, s, vt = np.linalg.svd(x_ave) # single value decomp. of average

    x_so3 = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt  # closest matrix to x_ave in SO(3)
    # convert to matrix in SO(3) 
    return x_so3


def calc_GD(atoms, grainDict, filename):
    wsum = 0
    count = 0

    with tqdm(total=len(grainDict.keys()), desc="Computing global distortion (GD)") as pbar: # progress bar (this method takes awhile)
        for gID in grainDict.keys(): # all grains in our sample, including surface ones
            grainIDs = [g for g in grainDict[gID] if g in atoms]  # the IDs of CORE atoms in our grain, not grain IDs
            if (len(grainIDs) == 0): 
                continue # grain is entirely cointained in the surface, exclude

            xG = central_pixel(grainIDs, atoms)

            for id in grainIDs: # this now excludes surface atoms
                x = atoms[id][2] # rotation matrix
                val, mat = omega(xG, x)
                wsum += val
                count += 1

            pbar.update(1)

    GD = wsum / count + 1

    with open(filename + "_GD.txt", "w") as f:
        f.write(str(GD))

    f.close()

    print("Global Distortion:", GD)

    return GD


# returns the kernel of a given pixel

kernel_state = {}

def init_kernel(filename, posToID, cutoff=5): # cutoff value  based on the interatomic distance of Al atoms, scaled by sqrt(3)
    data = import_file(f"{filename}.txt").compute()
    pos = [tuple(map(float, p)) for p in data.particles.positions]
    return {
        "finder": CutoffNeighborFinder(cutoff, data),
        "ovitotoID": {posToID[p]: i for i, p in enumerate(pos)}, # atom ID > OVITO ID
        "idOf": {i: posToID[p] for i, p in enumerate(pos)}, # OVITO ID > atom ID
    }



def get_kernel(atomID, filename, atoms, posToID):
    """Returns an atom's kernel set"""
    # we use the OVITO API for simplicity, which allows us to define neighbours within a shell
    
    if filename not in kernel_state:
        kernel_state[filename] = init_kernel(filename, posToID) # intialize kernel
    state = kernel_state[filename]

    return[state["idOf"][neigh.index] for neigh in state["finder"].find(state["ovitotoID"][atomID])]

    # positions array, keyed by OVITO indices (we're going to assume order is preserved, in which
    # case they'd be the same as ours. <NOTE: This is something that will cause bugs if not the case

