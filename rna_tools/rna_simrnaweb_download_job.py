#!/usr/bin/env python
"""rna_simrnaweb_download_job.py - download model files, trajectory for a given SimRNAweb job.

Usage::

    rp17pk$ rna_pdb_download_simrna_job.py 27b5093d -m -t -x
    # download more clusters, trajectory, extract100

    cp771_pk$ rna_pdb_download_simrna_job.py -t -x -m cf8f8bb2 -p cp771_pk
    # download with a trajectory, and cluster #4 and #5, add to all pdb files
    # prefix: cp771_pk

    $ rna_simrnaweb_download_job.py --web-models rp17_well_d10_e1-a43d3ab5 --prefix tar
    # prefix added will be tar_XXXX

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
import subprocess
import sys
import logging
import time

sys.tracebacklimit = 0

SIMRNAWEB_ARCHIVE = "/home/magnus/work/simrnaweb-archive"

class SimRNAwebError(Exception):
    pass


def extract(job_id, nstruc, remove_trajectory, lowest=True):
    """using rna_simrna_extract.py"""
    if lowest:
        # Ok, this is off for now because I can download top100 and top200
        cmd = 'rna_simrna_lowest.py -n ' + str(nstruc) + ' ' + job_id + '*.trafl'
        print(cmd)
        os.system(cmd)
    if remove_trajectory:
        os.system('rm *_ALL.trafl*')
    cmd = 'rna_simrna_extract.py -t ' + job_id + '*01X.pdb -f *' + job_id + '*top' + str(nstruc) + '.trafl -n ' + str(nstruc) + ' -c'
    print(cmd)
    os.system(cmd)


def download_trajectory():
    cmd = "wget --restrict-file-names=nocontrol http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + \
        "/processing_results/" + job_id + "_ALL.trafl "  # -O " + fn
    print(cmd)
    os.system(cmd)

def download_trajectory_top100():
    cmd = "wget --restrict-file-names=nocontrol http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + \
        "/processing_results/top100/" + job_id + "_ALL_top100.trafl "  # -O " + fn
    print(cmd)
    os.system(cmd)

def download_trajectory_top200():
    cmd = "wget --restrict-file-names=nocontrol http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + \
        "/processing_results/top200/" + job_id + "_ALL_top200.trafl "  # -O " + fn
    print(cmd)
    os.system(cmd)

def run_cmd(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def copy_trajectory(args):
    """Using xfind https://github.com/mmagnus/Xfind

    Do I have to copy it? Nope, you can make a link (ln)"""
    # OK, find this file on the drive
    # cmd = "mdfind -name " +  args.job_id  + " 000001_AA.pdb"
    print('[Get Trajectory]')
    cmd = "mdfind -name " + args.job_id + " | grep .ALL.trafl"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_fn = out.split('\n')[0]
        cmd = 'ln -s ' + path_to_fn + " " + job_id + "_ALL.trafl "
        print(cmd)
        os.system(cmd)
        return

    clustername = "malibu"
    cmd = "ssh " + clustername + " xfind " + args.job_id + " | grep .ALL.trafl"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_trafl = out.split('\n')[0]
        os.system('scp ' + clustername + ':' + path_to_trafl + " " + SIMRNAWEB_ARCHIVE + '/' + job_id + "_ALL.trafl ")
        return

    raise Exception("Trajectory not found on the cluster")

def copy_template(args):
    """Copy a template for the trajectory, on cluster and locally.

    Args:
      args.job_id

    If it fails, raise an Exception"""
    # OK, find this file on the drive
    # cmd = "mdfind -name " +  args.job_id  + " 000001_AA.pdb"
    cmd = "mdfind -name " +  args.job_id  + " | grep 'clust01X.pdb$'"  #'000001_AA.pdb$'"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_fn = out.split('\n')[0]
        print('Template found %s' % path_to_fn)
        # cmd = 'cp ' + path_to_fn + " " + job_id + "-01X.pdb"
        shutil.copyfile(path_to_fn, job_id + "-01X.pdb")
        return

    cmd = "mdfind -name " +  args.job_id  + " | grep '01X.pdb$'"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_fn = out.split('\n')[0]
        print('Template found %s' % path_to_fn)
        # cmd = 'cp ' + path_to_fn + " " + job_id + "-01X.pdb"
        shutil.copyfile(path_to_fn, job_id + "-01X.pdb")
        return

    cmd = "mdfind " +  args.job_id  + " | grep '000001_AA.pdb$'"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_fn = out.split('\n')[0]
        print('Template found %s' % path_to_fn)
        # cmd = 'cp ' + path_to_fn + " " + job_id + "-01X.pdb"
        shutil.copyfile(path_to_fn, job_id + "-01X.pdb")
        return

    # find it on the cluster
    clustername = "malibu"
    cmd = "ssh " + clustername + " xfind " + args.job_id + " | grep 'clust01-000001.pdb$'"
    print(cmd)
    out, err = run_cmd(cmd)
    if out:
        path_to_fn = out.split('\n')[0]
        print('Template found %s' % path_to_fn)
        cmd = 'scp ' + clustername + ':' + path_to_fn + " " + job_id + "-01X.pdb"
        print(cmd)
        os.system(cmd)
        return

    # hack with returns, not pretty
    raise Exception("Trajectory template not found on the cluster")

def add_prefix(prefix):
    """This is very crude way to do it"""
    # d2b57aef_ALL-thrs8.40A_clust01X.pdb -> gba_pk_d2b57aef_ALL-thrs8.40A_clust01X.pdb
    # print('disable right now')
    os.system("rename 's/^/" + prefix.strip() + "_/' *pdb")
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
                cmd = 'rna_simrna_extract.py -t *01X.pdb -f *04.trafl -c -n 1'
                subprocess.check_call(cmd, shell=True)
                os.remove(nfn)
                shutil.move(nfn.replace('.trafl', '-000001_AA.pdb'),
                            nfn.replace('.trafl', 'X.pdb'))
            if 'clust05.trafl' in l:
                cmd = 'rna_simrna_extract.py -t *01X.pdb -f *05.trafl -c -n 1'
                subprocess.check_call(cmd, shell=True)
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
            cmd = "wget --restrict-file-names=nocontrol http://genesilico.pl/SimRNAweb/media/jobs/" + job_id + "/output_PDBS/" + \
                  fn + " -O " + nfn.replace('%2B', '+')  # fix
            print(cmd)
            subprocess.check_call(cmd, shell=True)


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('job_id', help='job_id')
    parser.add_argument(
        '-p', '--prefix', help='prefix to the name, withouth _, be careful with this.'
        'If you have already some files with the given folder, their names might'
        'be changed.')
    parser.add_argument('-n', '--nstruc',
                       help='extract nstruc the lowest energy, this option must go with --web', type=int, default=100)
    parser.add_argument('-e', '--extract',
                       help='extract nstruc the lowest energy, this option must go with --web', action="store_true")
    ## parser.add_argument('-t', '--trajectory',
    ##                     action='store_true', help='download also trajectory')
    parser.add_argument('-m', '--more_clusters',
                        action='store_true', help='download also cluster 4 and 5')
    parser.add_argument('-r', '--remove-trajectory',
                        action='store_true', help='remove trajectory after analysis', default=False)
    parser.add_argument('-c', '--cluster',
                        action='store_true', help='get trajectory from cluster OR local on your computer (mdfind for macOS)', default=False)
    parser.add_argument('-d', '--download-trajectory',
                        action='store_true', help='web', default=False)
    parser.add_argument('--top100', action='store_true', help='download top100 trajectory', default=False)
    parser.add_argument('--top200', action='store_true', help='download top200 trajectory', default=False)
    parser.add_argument('--web-models',
                        action='store_true', help='web models download', default=False)

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

    if args.web_models or args.top100 or args.top200 or args.download_trajectory:
        download_models(args)  # models are needed if you want to extract anything from the trajectory

    # trajectory
    if args.top100:
        download_trajectory_top100()

    if args.top200:
        download_trajectory_top200()

    if args.download_trajectory:
        download_trajectory()

    if args.cluster:
        copy_trajectory(args)

    if (args.top100 or args.top200):
        extract(args.job_id, args.nstruc, args.remove_trajectory, lowest=False)

    if (download_trajectory or args.cluster) and args.nstruc:
        extract(args.job_id, args.nstruc, args.remove_trajectory, lowest=True)

    if args.extract:
        extract(args.job_id, args.nstruc, args.remove_trajectory)

    if args.prefix:
        add_prefix(args.prefix)

    #if not args.nolog:
    if 1:
        logging.basicConfig(filename='rna_simrnaweb_download_job.log',level=logging.INFO)
        #logging.basicConfig(level=logging.INFO)
        logging.info(time.strftime("%Y-%m-%d %H:%M"))
        logging.info(args)
