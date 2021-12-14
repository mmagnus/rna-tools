#!/usr/bin/python

"""
"""
from __future__ import print_function

import os
import subprocess
import sys
import time
import sendmail
import argparse
import shutil
from datetime import datetime, timedelta

sys.path.append(os.path.dirname(__file__))
os.environ['DJANGO_SETTINGS_MODULE'] = 'web.settings'
from app import models
JOB_STATUSES = dict(map(reversed, models.JOB_STATUS))
SAXS_FORMATS = dict(map(reversed, models.SAXS_FORMAT))
from django import db
from web import settings

MAX_RUNNING = 4  # 3 jobs only running
DISK_SPACE_LIMIT = 2  # each 4 of jobs 5gb = 20g free must be


def run_job(j):
    """A function to run one job"""
    print(time.asctime())

    j.status = JOB_STATUSES['running']
    if not DEBUG: j.status = JOB_STATUSES['running']
    j.started = datetime.now()
    j.save()

    # if not DEBUG:
    if j.email: sendmail.send_email_start(j.job_id, j.email)

    # ##################### CODE ##################################
    os.chdir(settings.JOBS_PATH + os.sep + j.job_id)
    print(cmdx('pwd')[0])
    print(cmdx('ls')[0])
    print('http://localhost:8667/RNAMasonry/jobs/' + j.job_id)

    # step
    steps = ' --steps ' + str(j.nsteps) + ' '

    # freeze
    freeze = ''
    if j.residues_to_freeze:
        freeze += "--freeze='" + j.residues_to_freeze.strip() + "' "
    # remove folder
    shutil.rmtree('rms_input_fa', ignore_errors=True)

    # sax format
    sax = ''
    if j.restraints_fn and j.saxs_format == SAXS_FORMATS['foxs']:
        sax = ' --saxs-data=' + j.restraints_fn
    if j.restraints_fn and j.saxs_format == SAXS_FORMATS['crysol']:
        sax = ' --crysol --saxs-data=' + j.restraints_fn

    input = ''
    if j.pdb_fn:  # it's also cif!
        input = " -i '" + j.pdb_fn + "'"
    else:
        input = " -i input.fa "

    cmd = settings.PATH_TO_RM + '/rnamasonryweb.sh --project_dir rms_input_fa --webservice ' + input + ' ' + sax + ' ' + freeze + ' ' + steps + ' | tee log.txt # job_id:' + j.job_id
    print(cmd)
    if DEBUG:
        debug = ''  # define some testing coditions
        os.system(settings.PATH_TO_RM + '/rnamasonryweb.sh -i input.fa ' + debug + ' | tee log.txt')  # for testing
        out = ''
        err = ''
    else:
        out, err = cmdx(cmd)
        print(out)
        print(err)

    if 'Normal termination of RNA Masonry' in out:
        print('OK - daemon')
        # if not DEBUG:
        j.status = JOB_STATUSES['finished']
        j.finished = datetime.now()
        j.save()
    else:
        # is stopped, then it means that the job was killed (which would cause an error)
        # here so ignore this as an error, and just keep stopped
        # if j.status == JOB_STATUSES['stopped']:
        #    print('Daemon to kill', err)
        if err == 'Killed':  # if killed this is fine, hopefully it was killed by
            # the other deamon to stop the job
            pass
        else:
            print('Error - daemon', err)
            j.get_error()
            j.status = JOB_STATUSES['finished with errors']

        j.finished = datetime.now()
        j.save()

    if not DEBUG:
        if j.email: sendmail.send_email_done(j.job_id, j.email)


def cmdx(cmd):
    print('> ', cmd)
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=False, bufsize=0)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def check_and_run(verbose, force_run, debug_id):
    if verbose:
        print('debug_id: %s' % debug_id)

    # to be able to run debug_id as a job
    if not debug_id:
        try:
            job_in_queue = models.Job.objects.filter(status=JOB_STATUSES['waiting'])
        except IndexError:
            # print time.asctime() + ' the queue is empty!'
            return
    else:
        job_in_queue = models.Job.objects.filter(job_id=debug_id)
        job_in_queue[0].status = JOB_STATUSES['waiting']

    # check how many running
    job_running = models.Job.objects.filter(status=JOB_STATUSES['running'])
    if verbose:
        print('running jobs:', job_running)
        print('# of running jobs is %i' % len(job_running))

    for j in job_in_queue:
        if j.email in settings.POWER_USERS:
            if verbose:
                print('run as root !!!!! %s ' % j.job_id)
            run_job(j)
            return

    # if force run then simply run it, I want to use this for testing
    if job_in_queue:
        if force_run:
            run_job(job_in_queue[0])

    # unlock this to run only magnus@
    if job_in_queue:
        if len(job_running) < MAX_RUNNING:
            # run the job
            run_job(job_in_queue[0])


def check_if_running(verbose=False):
    """Check all jobs statused `running` and check if are running, if not, change status to `finished`

    Off for now."""
    job_running = models.Job.objects.filter(status=JOB_STATUSES['running'])
    if verbose:
        print('running jobs:', job_running)
        print('# of running jobs is %i' % len(job_running))
    for j in job_running:
        if verbose: print(j)
        out = commands.getoutput('ps -aux 2> /dev/null | grep %s ' % j.job_id + '')  # | grep npdock')
        if out.find('npdock.py') > -1:
            if verbose: print(out)
        else:
            j.status = JOB_STATUSES['finished']
            j.save()


def kill_all_stopped():
    """Kill stopped jobs but not killed yet! to .. hmm.. kill them, yes!"""
    job_stopped_to_kill = models.Job.objects.filter(status=JOB_STATUSES['stopped']).filter(killed=0)
    print('stopped jobs:', job_stopped_to_kill)
    print('# of stopped jobs is %i' % len(job_stopped_to_kill))
    for j in job_stopped_to_kill:
        kill_job(j.job_id)
        j.killed = True
        j.save()


def kill_job(job_id):
    """
    ps gets this kind of listing, so the script tries to get this and find pid of all process and kill in a loop all of them one
    by one. For whatever reason, if you kill only one process, the rest (even children) are happy running. Confusing.
    Now it's a brute force solution but works.

    rnamaso+  52271  52269  0 May28 pts/3    00:00:00          \_ python manage.py runserver --settings web.settings 0.0.0.0:8667
    rnamaso+ 104238  52271  5 10:41 pts/3    00:00:01              \_ /home/rnamasonry/rnamasonryweb_env/bin/python manage.py runserver --settings web.settings 0.0.0.0:8667
    rnamaso+ 104243      1  0 10:42 ?        00:00:00 /bin/bash /home/rnamasonry/rnamasonryweb_env/rnamasonry-web/daemon.sh
    rnamaso+ 104245 104243  7 10:42 ?        00:00:00  \_ python daemon.py --prod --verbose
    rnamaso+ 104254 104245  0 10:42 ?        00:00:00      \_ /bin/sh -c /home/rnamasonry/rnamasonry//rnamasonryweb.sh --project_dir rms_input_fa --webservice  -i input.fa   --saxs-data=tRNA_saxs.dat   --steps 5  | tee log.txt # job_id:tRNA_with_SAXS-800d6744
    rnamaso+ 104255 104254  0 10:42 ?        00:00:00          \_ /bin/sh /home/rnamasonry/rnamasonry//rnamasonryweb.sh --project_dir rms_input_fa --webservice -i input.fa --saxs-data=tRNA_saxs.dat --steps 5
    rnamaso+ 104268 104255 99 10:42 ?        00:00:08          |   \_ /home/rnamasonry/rnamasonry/virtualenv/cctbx_build/../bin/python2 -Qnew -u /home/rnamasonry/rnamasonry/rnamasonry.py --project_dir rms_input_fa --webservice -i input.fa --saxs-data=tRNA_saxs.dat --steps 5
    rnamaso+ 104300 104268  0 10:42 ?        00:00:00          |       \_ /home/rnamasonry/rnamasonry/virtualenv/cctbx_build/../bin/python2 -Qnew -u /home/rnamasonry/rnamasonry/rnamasonry.py --project_dir rms_input_fa --webservice -i input.fa --saxs-data=tRNA_saxs.dat --steps 5
    rnamaso+ 104256 104254  0 10:42 ?        00:00:00          \_ tee log.txt
    """
    # find the parent process, deamon
    cmd = 'ps -aef --forest | grep rnama'
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps = o.stdout.read().strip().decode()

    # ugly hack: replace this start of the set of process to tee log.txt to be easily split this
    ps = ps.replace('python daemon.py --prod --verbose', 'tee log.txt')
    for set_of_processes in ps.split('tee log.txt'):
        if job_id in set_of_processes:
            for line in set_of_processes.split('\n'):
                # kill it all !
                pid = line.split()[1]
                if 'job_id' not in pid:
                    # kill -9 job_id:tRNA_with_SAXS-800d6744
                    # kill -9 104255
                    os.system('kill -9 ' + pid)


def kill_all_finished(verbose=0):
    """To make sure that all jobs finished in the db, are really not running on the cluster. You can change status of
    a job to `finished` and then it will be killed.

    This is off right now. "this It works for jobs from the last 10 days.". Why this should work only for last 10d?
    OK, other wise this is suppppppppppppppppppper long!

    kill -9 $(ps ax | grep d82544f4-a6e8-4655-af4e-b9346b85f92c | fgrep -v grep | awk '{ print $1 }')

    Off for now.
    """
    job_finished = models.Job.objects.filter(status=JOB_STATUSES['finished'], created__gte=datetime.now() - timedelta(days=10))
    if verbose:
        print('running jobs:', job_finished)
        print('# of running jobs is %i' % len(job_finished))
    for j in job_finished:
        if verbose: print(j)
        cmd = "kill -9 $(ps ax | grep 'npdock.py  -j " + j.job_id + "' | fgrep -v grep | awk '{ print $1 }')"
        if verbose: print(cmd)


def check_disc_space(verbose=True):
    """155G, return 1 to stop"""
    # disk hack
    #########################################
    cmd = "df -h | grep " + settings.DISK_TO_TRACK + " | awk '{ print $4 }'"
    if verbose: print(cmd)
    out = commands.getoutput(cmd)
    # if verbose: print cmd
    # out = '102G'
    # print out
    out = out.replace('Gi', 'G')  # on mac it's 6.5Gi

    if out.endswith('G'):
        space = out.replace('G', '')
        if verbose: 'space', space
        if float(space) < DISK_SPACE_LIMIT:  # for debugging + 3000:
            print('Error - NOT enough space to run a job')  # >>sys.stderr,
            # run cleaning
            # print >>sys.stderr, exe("/home/rpdock/web/cleanup.sh")
            return 1
        else:
            # say nothing if OK
            # if verbose: print 'OK, enough space to run a job' # mute it, this print is shown triggers an e-mail from my crontab
            return 0
    else:
        print('Error - NOT enough space to run a job')  # >>sys.stderr,
        # print >>sys.stderr, exe("/home/rpdock/web/cleanup.sh")
        return 1


def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


# main
if __name__ == '__main__':
    parser = argparse.ArgumentParser("")
    parser.add_argument("--debug", help="Debug mode")
    parser.add_argument("--prod", help="Production mode", action='store_true')
    parser.add_argument("--verbose", help="be verbose", action='store_true')
    parser.add_argument("--force-run", help="force run", action='store_true')

    args = parser.parse_args()

    print(settings.SERVER_NAME, '=' * 60)

    if len(sys.argv) == 1:
        print(parser.format_help())  # prints help if no arguments
        sys.exit(1)

    if args.debug:
        DEBUG = True
    if args.prod or args.force_run:
        DEBUG = False

    # check_if_running(args.verbose)
    # kill_all_finished(args.verbose)
    kill_all_stopped()

    # if 10gb free
    if not check_disc_space():
        check_and_run(args.verbose, args.force_run, args.debug)

    job_running = models.Job.objects.filter(status=JOB_STATUSES['running'])
    for j in job_running:
        print(j, j.created)

    db.connections.close_all()
