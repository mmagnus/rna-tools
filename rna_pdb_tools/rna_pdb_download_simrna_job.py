#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Download model files for a given SimRNAweb job.

  $ rna_pdb_download_simrna_job.py d86c07d9-9871-4454-bfc6-fb2e6edf13fc # --traj #@todo: download trajectory

The names will be shorten: ``d86c07d9-9871-4454-bfc6-fb2e6edf13fc_ALL_thrs12.50A_clust01-000001_AA.pdb -> d86c07d9-thrs12.50A_clust01X.pdb``.
"""
# Ideas:  71707ff4-fe16-4b78-8340-78913312a547_ALL_thrs12.50A_clust01-000001_AA keep
#                                              // 


import argparse
import os
import urllib3
import sys
sys.tracebacklimit = 0

class SimRNAwebError(Exception):
    pass

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('job_id', help='job_id')
    parser.add_argument('-p', '--prefix', help='prefix to the name, withouth _')
    parser.add_argument('-t', '--trajectory', action='store_true', help='download also trajectory')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    job_id = args.job_id
    
    job_id = job_id.replace('genesilico.pl/SimRNAweb/jobs/', '').replace('http://','').replace('/', '') # d86c07d9-9871-4454-bfc6-fb2e6edf13fc/

    # download models, get propare names of pdb files
    # http://genesilico.pl/SimRNAweb/media/jobs/d86c07d9-9871-4454-bfc6-fb2e6edf13fc/output_PDBS/d86c07d9-9871-4454-bfc6-fb2e6edf13fc_ALL_thrs12.50A_clust01-000001_AA.pdb
    http = urllib3.PoolManager()
    url = "http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + "/output_PDBS/"
    response = http.request('GET', url)
    if not response.status == 200: raise SimRNAwebError('Job not found on the server: %s' % job_id)
    html = response.data

    for l in html.split('\n'):
        if l.find('AA.pdb') > -1 and l.find('clust') > -1:
            # find fn
            fn = l.split('"')[1]

            # shorten names
            nfn = fn.replace("-000001", '').replace('_AA','X').replace('_ALL_','-')
            parts = nfn.split('-')
            print(parts)
            nfn = '-'.join([fn[:12], ''.join(parts[-1])])

            # wget
            cmd = "wget http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + "/output_PDBS/" + fn + " -O " + nfn
            os.system(cmd)

    # trajectory link
    # http://iimcb.genesilico.pl/SimRNAweb/media/jobs/
    # rp12aawlpk-8713ed35/processing_results/rp12aawlpk-8713ed35_ALL.trafl

    if args.trajectory:
        cmd = "wget http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + \
          "/processing_results/" + job_id + "_ALL.trafl "#-O " + fn
        os.system(cmd)

    # d2b57aef_ALL-thrs8.40A_clust01X.pdb -> gba_pk_d2b57aef_ALL-thrs8.40A_clust01X.pdb
    if args.prefix:
        os.system("rename 's/^/" + args.prefix.strip() + "_/' *")
