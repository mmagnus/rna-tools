#!/usr/bin/env python

"""Extract full atom structures from a SimRNA trajectory file.

Options:

  SIMRNA_DATA_PATH has to be properly defined in ``rpt_config_local``.

"""

from simrna_trajectory import *
import argparse
import os

from rna_pdb_tools.rpt_config import SIMRNA_DATA_PATH

def get_parser():
    """Get parser of arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template', help="template PDB file used for reconstruction to full atom models", required=True)
    parser.add_argument('-f', '--trafl', help="SimRNA trafl file", required=True)
    parser.add_argument('-c', '--cleanup', action='store_true', help="Keep only *_AA.pdb files, remove *.ss_detected and *.pdb")
    return parser

def get_data():
    """Get a link to SimRNA data folder in cwd."""
    print 'getting SimRNA data folder in cwd ...'
    os.system('ln -s %s %s' % (SIMRNA_DATA_PATH, os.getcwd()))
        
def extract(template, trafl):
    """Run SimRNA_trafl2pdb to extract all full atom structures in the trajectory."""
    os.system('SimRNA_trafl2pdbs %s %s AA :' % (template, trafl))

def cleanup(trafl):
    """Create _<trafl> with all CG structures and ss_detected and in current directory 
    keep only full atom structures.
    
    e.g ::

            _1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low

            [mm] sim_a99 ls _1* | head
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000001.pdb
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000001.ss_detected
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000002.pdb
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000002.ss_detected
    """
    traflfn = trafl.replace('.trafl','')
    try:
            os.mkdir('_%s' % traflfn)
    except OSError:
            pass
    os.system('mv %s*.pdb _%s' % (traflfn, traflfn) )
    os.system('mv %s*.ss_detected _%s' % (traflfn, traflfn) )
    os.system('mv _%s/*AA.pdb .' % traflfn)

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()    

    get_data()
    extract(args.template, args.trafl)
    if args.cleanup:
            cleanup(os.path.basename(args.trafl))
