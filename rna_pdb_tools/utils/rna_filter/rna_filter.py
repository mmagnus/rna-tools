#!/usr/bin/env python
"""rna_filter"""

from rna_pdb_tools.utils.rmsd_calc.lib.rmsd.calculate_rmsd import get_coordinates
from rna_pdb_tools.utils.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment
from rna_pdb_tools.utils.simrna_trajectory.simrna_trajectory import SimRNATrajectory

import argparse
import re
import numpy as np
import sys

def parse_logic(restraints_fn, verbose):
    txt = ''
    with open(restraints_fn) as f:
        for l in f:
            if not l.startswith('#'):
                txt += l.strip()
    if verbose:
        print txt
    restraints = re.findall('\(d:(?P<start>.+?)-(?P<end>.+?)\s*(?P<operator>\>\=|\=|\<|\<\=)\s*(?P<distance>[\d\.]+)\s+(?P<weight>.+?)\)', txt)
    return restraints # [('A9', 'A41', '10.0', '1'), ('A10', 'A16', '10', '1')]

def get_distance(a,b):
    """
    """
    diff = a - b
    return np.sqrt(np.dot(diff, diff)) 

def parse_pdb(pdb_fn, selection):
    """
{'A9': {'OP1': array([ 53.031,  21.908,  40.226]), 'C6': array([ 54.594,  27.595,  41.069]), 'OP2': array([ 52.811,  24.217,  39.125]), 'N4': array([ 53.925,  30.861,  39.743]), "C1'": array([ 55.611,  26.965,  43.258]), "C3'": array([ 53.904,  25.437,  43.809]), "O5'": array([ 53.796,  24.036,  41.353]), 'C5': array([ 54.171,  28.532,  40.195]), "O4'": array([ 55.841,  25.746,  42.605]), "C5'": array([ 54.814,  23.605,  42.274]), 'P': array([ 53.57 ,  23.268,  39.971]), "C4'": array([ 55.119,  24.697,  43.283]), "C2'": array([ 54.563,  26.706,  44.341]), 'N1': array([ 55.145,  27.966,  42.27 ]), "O2'": array([ 55.208,  26.577,  45.588]), 'N3': array([ 54.831,  30.285,  41.747]), 'O2': array([ 55.76 ,  29.587,  43.719]), 'C2': array([ 55.258,  29.321,  42.618]), "O3'": array([ 53.272,  24.698,  44.789]), 'C4': array([ 54.313,  29.909,  40.572])}}
    """
    V = {}
    with open(pdb_fn) as f:
        for line in f:
            if line.startswith("ATOM"):
                curr_chain_id = line[21]
                curr_resi = int(line[22:26])
                curr_atom_name = line[12:16].strip()
                if selection:
                    if selection.has_key(curr_chain_id):
                        if curr_resi in selection[curr_chain_id]:
                            x = line[30:38]
                            y = line[38:46]
                            z = line[46:54]
                            #V.append(np.asarray([x,y,z],dtype=float))
                            if V.has_key(curr_chain_id + str(curr_resi)):
                                V[curr_chain_id + str(curr_resi)][curr_atom_name] = np.asarray([x,y,z],dtype=float)
                            else:
                                V[curr_chain_id + str(curr_resi)] = {}
                                V[curr_chain_id + str(curr_resi)][curr_atom_name] = np.asarray([x,y,z],dtype=float)
    return V
        
def check_condition(condition, wight):
    """return True/False, score"""
    pass

def get_residues(pdb_fn, restraints, verbose):
    residues = set()
    for h in restraints:
        a = h[0]
        b = h[1]
        a = a[0] + ':' + a[1:]
        residues.add(a) # A19
        b = b[0] + ':' + b[1:]
        residues.add(b)
    # set(['A:41', 'A:9', 'A:10', 'A:16'])

    selection = ','.join(residues)
    selection_parsed = select_pdb_fragment(selection, separator=",", splitting="[,:;]")

    residues = parse_pdb(pdb_fn, selection_parsed)

    # get mb
    for r in residues:
        if residues[r].has_key('N9'): # A,G
            residues[r]['mb'] = residues[r]['N9'] - ((residues[r]['N9'] - residues[r]['C6']) / 2 )
        else: # A,G
            residues[r]['mb'] = residues[r]['N1'] - ((residues[r]['N1'] - residues[r]['C4']) / 2 )
    for r in residues:
        #print 'mb for ' + str(r) + ' is ' + residues[r]['mb']
        print ' mb for ', str(r), residues[r]['mb']
    return residues

# main
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="prog [<options>] <pdb files: test_data/*>")

    parser.add_argument('-r',"--restraints_fn",
                        dest="restraints_fn",
                        required=True,
                        help="""restraints_fn:
Format:
(d:A9-A41 <  10.0  1)|(d:A41-A9 <= 10 1)
""")

    parser.add_argument("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="be verbose")

    parser.add_argument('-s', dest="structures", help='structures', nargs='+')#, type=string)

    parser.add_argument('-t', dest="trajectory")#help='structures', nargs='+')#, type=string)

    args = parser.parse_args()

    pdb_files = args.structures
    verbose = args.verbose
    restraints_fn = args.restraints_fn
    #score = 1
    #print ((True|True)|(False|False)), score

    restraints = parse_logic(restraints_fn, verbose)
    print ' restraints', restraints

    # h = ('A1', 'A2', '<', '10.0', '1')
    if args.structures:
        for pdb_fn in pdb_files:
            print '\n', pdb_fn
            residues = get_residues(pdb_fn, restraints, verbose)
            for h in restraints:
                dist = get_distance(residues[h[0]]['mb'], residues[h[1]]['mb'])
                if verbose:
                    print '  d:' + h[0] + '-' + h[1] + ' ' + str(dist)

    if args.trajectory:
        print
        f = (line for line in open(args.trajectory).xreadlines())
        c = 0
        while 1:
            try:
                header = f.next().strip()
            except StopIteration: # not nice
                break
            c += 1
            coords = f.next().strip()
            traj = SimRNATrajectory()
            traj.load_from_string(c, header + '\n' + coords)
            frame = traj.frames[0]
            print(c)
            for h in restraints:
                a = int(h[0].replace('A','')) - 1 # A1 -> 0 (indexing Python-like)
                b = int(h[1].replace('A','')) - 1 
                a_mb = frame.residues[a].get_center()
                b_mb = frame.residues[b].get_center()
                #print '  mb for A' + str(a+1), a_mb
                #print '  mb for A' + str(b+1), b_mb
                dist = get_distance(a_mb, b_mb)
                print '   d:A' + str(a+1) + "-A" + str(b+1),  dist
