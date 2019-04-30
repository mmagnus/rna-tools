#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This is a quick and dirty method of comparison two RNA structures (stored in pdb files).
It measures the distance between the relevan atoms (C4') for nucleotides defined as "x" in the
sequence alignment.

author: F. Stefaniak, modified by A. Zyla,  supervision of mmagnus
"""
from __future__ import print_function
from Bio.PDB import PDBParser
from scipy.spatial import distance
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt

# logger
logger = logging.getLogger()
handler = logging.StreamHandler()
logger.addHandler(handler)


def get_seq(alignfn, seqid):
    """Get seq from an alignment with gaps.

    Args:

       alignfn (str): a path to an alignment
       seqid   (str): seq id in an alignment

    Usage::

        >>> get_seq('test_data/ALN_OBJ1_OBJ2.fa', 'obj1')
        SeqRecord(seq=SeqRecord(seq=Seq('GUUCAG-------------------UGAC-', SingleLetterAlphabet()), id='obj1', name='obj1', description='obj1', dbxrefs=[]), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])

    Returns:
         SeqRecord
    """
    # alignment = AlignIO.read(alignfn, 'fasta')
    alignment = SeqIO.index(alignfn, 'fasta')
    # print SeqRecord(alignment[seqid])
    sequence = SeqRecord(alignment[seqid])
    return sequence


def open_pdb(pdbfn):
    """Open pdb with Biopython.

    Args:
       pdbfn1 (str): a path to a pdb structure

    Returns:
       PDB Biopython object: with a pdb structure

    """
    parser = PDBParser()
    return parser.get_structure('', pdbfn)


def find_core(seq_with_gaps1, seq_with_gaps2):
    """.

    Args:
        seq_with_gaps1 (str): a sequence 1 from the alignment
        seq_with_gaps1 (str): a sequence 2 from the alignment

    Usage::

        >>> find_core('GUUCAG-------------------UGAC-', 'CUUCGCAGCCAUUGCACUCCGGCUGCGAUG')
        'xxxxxx-------------------xxxx-'


    Returns:
        core="xxxxxx-------------------xxxx-"

    """
    core = "".join(["x" if (a != '-' and b != '-') else "-" for (a, b)
                    in zip(seq_with_gaps1, seq_with_gaps2)])
    return core


def map_coords_atom(structure):
    """.

        Args:
        structure (pdb): PDB Biopython object: with a pdb structure

        Returns:
            struct1dict: a list of coords for atoms
            structure1realNumber: a list of residues
        """

    resNumber = 0
    struct1dict = {}
    structure1realNumbers = {}

    for res in structure.get_residues():
        # print res
        for atom in res:
            name, coord = atom.name.strip(), atom.coord
            # G C5' [-15.50800037  -7.05600023  13.91800022]
            if name == atomToCompare:
                # print name, coord
                struct1dict[resNumber] = coord
                structure1realNumbers[resNumber] = res.get_full_id()[3][1]

        resNumber += 1
    return struct1dict,  structure1realNumbers


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("seqid1", help="seq1 id in the alignemnt")
    parser.add_argument("seqid2", help="seq2 id in the alignemnt")
    parser.add_argument("alignfn", help="alignemnt in the Fasta format")
    parser.add_argument("pdbfn1", help="pdb file1")
    parser.add_argument("pdbfn2", help="pdb file2")
    return parser


# main
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    args = get_parser().parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    # if not args.outfn:
        # args.outfn = args.pdbfn.replace('.pdb', '_out.pdb')

    seq_with_gaps1 = get_seq(args.alignfn, args.seqid1)
    seq_with_gaps2 = get_seq(args.alignfn, args.seqid2)
    pdb1 = open_pdb(args.pdbfn1)
    pdb2 = open_pdb(args.pdbfn2)
    core = find_core(seq_with_gaps1, seq_with_gaps2)
    atomToCompare = "C4'"
    sep = '\t'

    structure1 = pdb1
    structure2 = pdb2

    struct1dict, structure1realNumbers = map_coords_atom(pdb1)
    struct2dict, structure2realNumbers = map_coords_atom(pdb2)

    stats = []
    stats.append(["res1", "res2", "distance [A]"])

    resNumber = 0
    seq1number = -1
    seq2number = -1

    for char in core:
        # local sequence numbering
        if seq_with_gaps1[resNumber] != '-':
            seq1number += 1
            # print "seq1", seq1number, seq1[resNumber]
        if seq_with_gaps2[resNumber] != '-':
            seq2number += 1
            # print "seq2", seq2number, seq2[resNumber]

        # alignment checking (iksy)
        if char == 'x':
            vect1 = struct1dict[seq1number]
            vect2 = struct2dict[seq2number]
            stats.append([str(structure1realNumbers[seq1number]),
                          str(structure2realNumbers[seq2number]),
                          str(distance.euclidean(vect1, vect2))])
            # print vect1,vect2

        resNumber += 1

    # struc = renumber(seq_with_gaps, pdb, args.residue_index_start)
    # write_struc(struc, args.outfn)
    list_res=[]
    for i in stats:
        print(sep.join(i))
        table=(sep.join(i))
        list_res.append(i)
    #print (list_res)
    res_matrix = np.array(list_res[1:])

    '''
    Creating a plot

    '''

    new_resi = list_res[1:]
    new_resis = []

    for i in new_resi:
        #print (j)
        new_resis.append(str(i[0]) + '/' + str(i[1]))


    #print (new_resis)

    list2_matrix= new_resis

    list2_matrix1 = map(float, list(res_matrix[:,2]))
    #print (list2_matrix1)

    plt.bar(list2_matrix,list2_matrix1,facecolor='pink' )

    plt.suptitle('Distance between C4 atoms of residues')
    plt.ylabel("distance [A]")
    plt.xlabel('Nr of residue')
    plt.show()
