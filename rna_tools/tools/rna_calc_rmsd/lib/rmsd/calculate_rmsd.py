#!/usr/bin/env python

"""
Calculate RMSD between two XYZ files

by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Andersen Bratholm <larsbratholm@gmail.com>
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE

"""
import numpy as np
import re
from rna_tools.tools.extra_functions.select_fragment import is_in_selection

def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = rotate(P, Q)
    return rmsd(P, Q)


def rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def get_coordinates(filename, selection, ignore_selection, fmt, ignore_hydrogens):
    """Get coordinates from filename."""
    return get_coordinates_pdb(filename, selection, ignore_selection, ignore_hydrogens)


def get_coordinates_pdb(filename, selection, ignore_selection, ignore_hydrogens):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    V = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.
    atoms = []
    resi_set = set()
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            # hmm...
            # of models: 490
            # Error: # of atoms is not equal target (1f27_rpr.pdb):641 vs model (struc/1f27_rnakbnm_decoy0001_amb_clx_rpr.pdb):408
            #if line.startswith("TER") or line.startswith("END"):
            #    break
            if line.startswith("ATOM"):
                curr_chain_id = line[21]
                curr_resi = int(line[22:26])
                curr_atom_name = line[12:16].strip()
                if selection:
                    if curr_chain_id in selection:
                        if curr_resi in selection[curr_chain_id]:
                            # ignore if to be ingored (!)
                            #try:
                                    resi_set.add(curr_chain_id + ':' + str(curr_resi))
                                    x = line[30:38]
                                    y = line[38:46]
                                    z = line[46:54]
                                    if ignore_selection:
                                        if not is_in_selection(ignore_selection, curr_chain_id, curr_resi, curr_atom_name):
                                            V.append(np.asarray([x,y,z],dtype=float))
                                    else:
                                        V.append(np.asarray([x,y,z],dtype=float))
                else:
                                    x = line[30:38]
                                    y = line[38:46]
                                    z = line[46:54]
                                    V.append(np.asarray([x,y,z],dtype=float))

    V = np.asarray(V)
    #print filename, resi_set, len(resi_set)
    return len(V), V
