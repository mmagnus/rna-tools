#!/usr/bin/env python
"""rna_simrnaweb_download_job.py - download model files, trajectory for a given SimRNAweb job.

Usage::

    rp17pk$ rna_pdb_download_simrna_job.py 27b5093d -m -t -x
    # download more clusters, trajectory, extract100

    cp771_pk$ rna_pdb_download_simrna_job.py -t -x -m cf8f8bb2 -p cp771_pk
    # download with a trajectory, and cluster #4 and #5, add to all pdb files
    # prefix: cp771_pk

Example::

    rna_pdb_download_simrna_job.py -t -x -m 20569fa1 -p zmp_pk

    [mm] zmp_pk ls
    20569fa1_ALL_100low.trafl
    _20569fa1-thrs7.10A_clust04
    _20569fa1-thrs7.10A_clust05
    _20569fa1_ALL_100low
    data
    rna_simrna_extract.log
    subset.png
    zmp_pk_20569fa1-thrs7.10A_clust01X.pdb
    zmp_pk_20569fa1-thrs7.10A_clust02X.pdb
    zmp_pk_20569fa1-thrs7.10A_clust03X.pdb
    zmp_pk_20569fa1-thrs7.10A_clust04X.pdb
    zmp_pk_20569fa1-thrs7.10A_clust05X.pdb

.. downloaded clusters from 1 to 5, all pdb files have added prefix `zmp_pk`.
"""

# Ideas:  71707ff4-fe16-4b78-8340-78913312a547_ALL_thrs12.50A_clust01-000001_AA keep
#                                              //


import argparse
import os
import urllib3
import shutil
sys.tracebacklimit = 0
import subprocess


class SimRNAwebError(Exception):
    pass


def extract100(args):
    if args.extract100:
        os.system('rna_simrna_lowest.py *_ALL.trafl')
        os.system('rm *_ALL.trafl*')
        os.system('rna_simrna_extract.py -t *01X.pdb -f *low.trafl -c')


def download_trajectory(args):
    if args.trajectory:
        cmd = "wget http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + \
            "/processing_results/" + job_id + "_ALL.trafl "  # -O " + fn
        os.system(cmd)


def add_prefix(args):
    # d2b57aef_ALL-thrs8.40A_clust01X.pdb -> gba_pk_d2b57aef_ALL-thrs8.40A_clust01X.pdb
    if args.prefix:
        # print('disable right now')
        os.system("rename 's/^/" + args.prefix.strip() + "_/' *pdb")
        # os.system("rename 's/^/" + args.prefix.strip() + "_/' *_ALL.tarfl*")


def more_clusters(args):
    print('more clusters')
    http = urllib3.PoolManager()
    url = "http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + "/processing_results/"
    response = http.request('GET', url)
    if not response.status == 200:
        raise SimRNAwebError('Job not found on the server: %s' % job_id)
    html = response.data

    for l in html.split('\n'):
        if 'clust04.trafl' in l or 'clust05.trafl' in l:
            fn = l.split('"')[1]
            print(fn)

            # shorten names
            nfn = fn.replace("-000001", '').replace('_AA',
                                                    'X').replace('_ALL_', '-')
            parts = nfn.split('-')
            print(parts)
            # nfn = '-'.join([fn[:12], ''.join(parts[-1])])

            # wget
            cmd = "wget http://genesilico.pl/SimRNAweb/media/jobs/" + \
                job_id + "/processing_results/" + fn + " -O " + nfn
            os.system(cmd)

            if 'clust04.trafl' in l:
                os.system(
                    'rna_simrna_extract.py -t *01X.pdb -f *04.trafl -c -n 1')
                os.remove(nfn)
                shutil.move(nfn.replace('.trafl', '-000001_AA.pdb'),
                            nfn.replace('.trafl', 'X.pdb'))
            if 'clust05.trafl' in l:
                os.system(
                    'rna_simrna_extract.py -t *01X.pdb -f *05.trafl -c -n 1')
                os.remove(nfn)
                # 27b5093d-thrs6.20A_clust05-000001_AA.pdb -> 27b5093d-thrs6.20A_clust05X.pdb
                shutil.move(nfn.replace('.trafl', '-000001_AA.pdb'),
                            nfn.replace('.trafl', 'X.pdb'))


def download_models(args):
    # download models, get propare names of pdb files
    # http://genesilico.pl/SimRNAweb/media/jobs/d86c07d9-9871-4454-bfc6-fb2e6edf13fc/output_PDBS/d86c07d9-9871-4454-bfc6-fb2e6edf13fc_ALL_thrs12.50A_clust01-000001_AA.pdb
    http = urllib3.PoolManager()
    url = "http://genesilico.pl/SimRNAweb/media/jobs/" + \
        job_id + "/output_PDBS/"
    response = http.request('GET', url)
    if not response.status == 200:
        raise SimRNAwebError('Job not found on the server: %s' % job_id)
    html = response.data

    for l in html.split('\n'):
        if l.find('AA.pdb') > -1 and l.find('clust') > -1:
            # find fn
            fn = l.split('"')[1]

            # shorten names
            nfn = fn.replace("-000001", '').replace('_AA',
                                                    'X').replace('_ALL_', '-')
            parts = nfn.split('-')
            print(parts)
            # nfn = '-'.join([fn[:12], ''.join(parts[-1])]) ### ? nfn

            # wget
            cmd = "wget http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + "/output_PDBS/" + \
                  fn + " -O " + nfn
            subprocess.check_call(cmd, shell=True)


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('job_id', help='job_id')
    parser.add_argument(
        '-p', '--prefix', help='prefix to the name, withouth _, be careful with this')
    parser.add_argument('-x', '--extract100',
                        action='store_true', help='extract 100 the lowest')
    parser.add_argument('-t', '--trajectory',
                        action='store_true', help='download also trajectory')
    parser.add_argument('-m', '--more_clusters',
                        action='store_true', help='download also cluster 4 and 5')
    return parser


# main
if __name__ == '__main__':
    # trajectory link
    # http://iimcb.genesilico.pl/SimRNAweb/media/jobs/
    # rp12aawlpk-8713ed35/processing_results/rp12aawlpk-8713ed35_ALL.trafl
    parser = get_parser()
    args = parser.parse_args()
    job_id = args.job_id
    job_id = job_id.replace('genesilico.pl/SimRNAweb/jobs/', '').replace(
        'http://', '').replace('/', '')  # d86c07d9-9871-4454-bfc6-fb2e6edf13fc/
    download_models(args)
    download_trajectory(args)
    more_clusters(args.more_clusters)
    extract100(args)
    add_prefix(args)
