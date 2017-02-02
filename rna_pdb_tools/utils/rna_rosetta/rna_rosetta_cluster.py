#!/usr/bin/env python

"""**run_rosetta** - wrapper to ROSETTA tools for RNA modeling

Based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015.
http://www.sciencedirect.com/science/article/pii/S0076687914000524 

    rna_rosetta_cluster.py ade_min.out 20000

"""

import platform
import argparse
import os
import glob
import subprocess
import math
import logging
import shutil

logging.basicConfig(level=logging.INFO)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='ade.out')
    parser.add_argument('n', type=int, help='# of total structures')
    return parser

def get_no_structures(file):
    p = subprocess.Popen('cat ' + file + ' | grep SCORE | wc -l', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().strip()
    if stderr:
        print stderr
    return int(p.stdout.read().strip()) - 1

def get_selected(file, nc):
    """Get selected for clustering"""
    cmd = 'silent_file_sort_and_select.py ' + file +' -select 1-' + str(nc) + ' -o selected.out'
    logging.info(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().strip()
    if stderr:
        print stderr


def cluster(radius=1):
    """Internal function of cluster_loop: It removes cluster.out first."""
    try: os.remove("cluster.out")
    except OSError: pass

    if platform.system() == "Darwin":
        cmd = 'cluster.macosclangrelease'
    if platform.system() == "Linux":
        cmd = "cluster.default.linuxgccrelease"
    cmd += " -in:file:silent selected.out -in:file:fullatom -out:file:silent_struct_type binary -export_only_low false -out:file:silent cluster.out -cluster:limit_clusters 1 -cluster:radius " + str(radius)

    #'-out:prefix ade_ "#| tee clustering.out"
    logging.info(cmd)
    os.system(cmd)
    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #p.wait()
    #stderr = p.stderr.read().strip()
    #if stderr:
    #    print stderr

def extract():
    #cmd = "~/src/rna-pdb-tools/rna_pdb_tools/utils/rna_rosetta/rna_extract_lowscore_decoys.py cluster0.out 5"
    if platform.system() == "Darwin":
        cmd = 'extract_pdbs.default.macosclangrelease -in::file::silent cluster.out'
    if platform.system() == "Linux":
        cmd = 'extract_pdbs.default.linuxgccrelease -in::file::silent cluster.out'
    os.system(cmd)
    
def cluster_loop(ns):
    """Go from radius 1 to get 1/6 of structures of ns (# of selected structures)
    in the first cluster, then it stops."""
    radius = 6 # should be 1
    while 1:
        cluster(radius)
        n = get_no_structures('cluster.out')
        if n > ns * .16: # 1/6
            break
        radius += 1
    logging.info('radius %i n %i' % (radius, n))
        
def run():
    """Pipline for modeling RNA"""
    args = get_parser().parse_args()

    n = args.n #get_no_structures(args.file) # if you mini then # is the total number of structures
    nc = int(math.ceil(n * 0.005)) # nc no for clustring
    print '# total:', n, ' # selected:', nc
    get_selected(args.file, nc)
    ns = get_no_structures('selected.out')
    print '# selected:', ns

    ns = get_no_structures('selected.out')

    cluster_loop(ns) # loop over selected to get 1/6 in the biggest cluster

    print '# of structures in the biggest cluster", get_no_structures('cluster.out')

    extract()

#main
if __name__ == '__main__':
    run()
