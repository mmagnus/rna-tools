#!/usr/bin/env python
"""renum_pdb_to_aln.py - renumber a pdb file based on the alignment.

author: A. Zyla under supervision of mmagnus

.. warning:: works only for single chain! and requires Biopython (tested with v1.68)
"""

import logging
import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Atom import PDBConstructionWarning
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

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
       pdbfn (str): a path to a pdb structure

    Returns:
       PDB Biopython object: with a pdb structure

    """
    parser = PDBParser()
    return parser.get_structure('struc', pdbfn)


def renumber(seq_with_gaps, struc, residue_index_start):
    """Renumber a pdb file.

    Args:
        seq_with_gaps (str): a target sequence extracted from the alignment
        struc         (pdb): a structure
        residue_index_start (int): starting number

    Returns:
        BioPython Structure object

    """
    new_numbering = []
    for nt in seq_with_gaps:
        if nt != '-':
            nt_num_a = [residue_index_start, nt]
            new_numbering.append(residue_index_start)
            logger.info(nt_num_a)
        residue_index_start = residue_index_start + 1

    logger.info(new_numbering)

    # works only for single chain
    for struc in pdb:
        for chain in struc:
            for residue, resi in zip(chain, new_numbering):
                residue.id = (residue.id[0], resi, residue.id[2])
    return struc


def write_struc(struc, outfn):
    """Write renumbered pdb with Biopython.

    Args:
        struc (pdb): a renumbered structure
        outfn (str): a path to a new, renumbered pdb file

    Returns:
        none: writes to a file

    """
    io = PDBIO()
    io.set_structure(struc)
    io.save(outfn)
    logger.info('Structure written to %s' % outfn)


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    parser.add_argument("--residue_index_start",
                        help="renumber starting number (default: 1)",
                        default=1, type=int)
    parser.add_argument("--outfn", help="output pdb file (default: pdbfn .pdb -> _out.pdb)")
    parser.add_argument("seqid", help="seq id in the alignemnt")
    parser.add_argument("alignfn", help="alignemnt in the Fasta format")
    parser.add_argument("pdbfn", help="pdb file")
    return parser


# main
if __name__ == '__main__':
    args = get_parser().parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if not args.outfn:
        args.outfn = args.pdbfn.replace('.pdb', '_out.pdb')

    seq_with_gaps = get_seq(args.alignfn, args.seqid)
    pdb = open_pdb(args.pdbfn)
    struc = renumber(seq_with_gaps, pdb, args.residue_index_start)
    write_struc(struc, args.outfn)
