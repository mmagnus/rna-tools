#!/usr/bin/python
"""
Remove jobs from VM. MAILS_TO_KEEP_DATA set up to keep the data. REAL_REMOVE

RNA models for: 260c8ff6-f24e-4eff-9760-1831407fc770/ (1kh6 revised without restraints)
RNA models for: f3c68edc-99c2-4950-9efd-0486812f15ae (3l0u ss + rst)
RNA models for: 4371aff4-cea1-4deb-88ce-fbb9dc6ae5a4 (tetraloop)

Remember that this is 20 days from a moment when a job has been created.
"""
import os
import subprocess
import sys
import time
USE_TZ = False
import datetime
import argparse
from os import listdir
import shutil
import glob
import datetime

import warnings
from django.core.wsgi import get_wsgi_application

warnings.filterwarnings(
        'error', r"DateTimeField .* received a naive datetime",
        RuntimeWarning, r'DateTimeField received a naive datetime')

sys.path.append('/home/simrna/simrnaweb_env/simrnaweb')
JOB_DIR = '/home/simrna/simrnaweb_env/simrnaweb/media/jobs/' # slash at the end must be
os.environ['DJANGO_SETTINGS_MODULE'] = 'web.settings'
application = get_wsgi_application()

from app import models
JOB_STATUSES = dict(map(reversed, models.JOB_STATUS))

# Settings
# number of days to keep the data on the server
OLD_IS_MORE = 14
# jobs to be kept
MAILS_TO_KEEP_DATA = ['iamb@genesilico.pl', 'rpluta@genesilico.pl', 'magnus@genesilico.pl','azyla@genesilico.pl']

# don't remove these jobs
JOBS_TO_KEEP = ['260c8ff6-f24e-4eff-9760-1831407fc770', 'f3c68edc-99c2-4950-9efd-0486812f15ae', '4371aff4-cea1-4deb-88ce-fbb9dc6ae5a4']

from web import settings

from django.utils import timezone

def remove_old_jobs(args):
    """time of job is calculated beased on the files! 
    If there is now file then you don't estimate the time

    args.debug vs args.prod. prod kick off REAL_REMOVE

    Keep jobs that are on JOBS_TO_KEEP list"""

    def del_job(j, REAL_REMOVE=True, not_in_db=False):
        """Low-level function of removing
        j = j.job_id
        not_in_db False, remove all files!"""
        if True:
            try:
                if True:
                        try:
                            if not_in_db: # remove all files, this job is not in the db
                                print 'not_in_db', j
                                if REAL_REMOVE: shutil.rmtree(JOB_DIR + j)
                                return 
                            else:
                                print 'full removing!!!!', j, 
                                if REAL_REMOVE: shutil.rmtree(JOB_DIR + j)
                                return
                        except OSError:
                            print 'no such file or directory -- skip (it means that there is demo job in the db , but not on the drive -- its fine! :-)', 
                else:
                    print 'removing (moving to be removed later!)', j, 
                    #shutil.move('/home/rpdock/web/media/jobs/' + j, '/home/magnus/to-remove-dir')
            except IOError:
                print 'IOError -- problem',
        else:
            print 'to be removed', j,

    if args.debug:
        REAL_REMOVE=False
    if args.prod:
        REAL_REMOVE=True

    jobs = models.Job.objects.filter()
    jobs_job_ids = []

    d = datetime.timedelta(days=OLD_IS_MORE)

    for j in jobs:
        #jobs_job_ids.append(j.job_id)
        fn = JOB_DIR + '/' + str(j) #+ '/' + 'log.txt'

        try:
            t = datetime.date.fromtimestamp(os.path.getmtime(fn))
            #tt = datetime.datetime.strptime(t)
            #print 't - d', t - d,
            today = datetime.date.today()
            delta = today - t
            if delta > d:
                print j.job_id, j.email, '-- old',    
                if j.job_id in JOBS_TO_KEEP:
                    print '-- to KEEP (JOBS_TO_KEEP)'
                    continue
                if j.email not in MAILS_TO_KEEP_DATA:
                    if j.job_id not in JOB_STATUSES:
                        print '-- to remove',
                        del_job(j.job_id, REAL_REMOVE)
                else:
                    print '-- to keep it',
                    pass
                print
            #print 'fn', fn, "last modified: %s" % time.ctime(os.path.getmtime(fn)), j.email
            pass
        
        except OSError:
            #print
            pass

#main
if __name__ == '__main__':
    parser = argparse.ArgumentParser("")
    parser.add_argument("--debug", help="Debug mode", action='store_true')
    parser.add_argument("--prod", help="Production mode", action='store_true')
    args = parser.parse_args()

    remove_old_jobs(args)
