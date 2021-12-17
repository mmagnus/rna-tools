import os
import sys

from web.settings import JOBS_PATH
from django.db import models
from django.contrib.auth.models import User
from django.utils.html import format_html

import uuid

JOB_STATUS = (
        (0, 'waiting'),
        (1, 'running'),
        (2, 'finished'),
        (3, 'finished with errors'),
        (4, 'removed'),
        (5, 'queue'),
    (6, 'stopped'),
    (7, 'urgent')
    )
JOB_STATUSES = dict(map(reversed, JOB_STATUS))

SAXS_FORMAT = (
        (0, 'foxs'),
        (1, 'crysol'),)       

SAXS_FORMATS = dict(map(reversed, SAXS_FORMAT))

class Setting(models.Model):
    """Setting"""
    run_when_usage_is_lower_than = models.FloatField()
    def __str__(self, ):
        return 'max_jobs: ' + str(self.run_when_usage_is_lower_than)
    def __repr__(self, ):
        return 'max_jobs: ' + str(self.run_when_usage_is_lower_than)
    
class Job(models.Model):
    """Job class"""
    job_id = models.CharField(max_length=50)
    job_title = models.CharField(max_length=50)

    status = models.IntegerField(max_length=1, choices=JOB_STATUS)

    email = models.EmailField(max_length=100, blank=True)

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    comment = models.CharField(max_length=200, blank=True)

    seq = models.TextField(max_length=300, blank=True)
    ss = models.TextField(max_length=300, blank=True)
    ss_raw = models.TextField(max_length=300, blank=True)

    pdb_fn = models.CharField(max_length=200, blank=True)
    restraints_fn = models.CharField(max_length=200, blank=True)

    nsteps = models.IntegerField(max_length=100, default=100)
    percentage_of_frames = models.IntegerField(blank=True, default=1, null=True)

    error_text = models.TextField(blank=True, null=True)

    residues_to_freeze =  models.CharField(max_length=200, blank=True, null=True)

    saxs_format = models.IntegerField(max_length=1, choices=SAXS_FORMAT)

    progress = models.FloatField(blank=True, null=True)
    #stopped = models.IntegerField(default=False, null=True)
    killed = models.BooleanField(default=False)

    seq_len = models.IntegerField(max_length=100, default=0)    
    ip = models.CharField(max_length=100, blank=True, null=True)

    def url(self):
        return format_html("<a target='_blank' href='/rnamasonry/jobs/" + str(self.job_id) + "/'>" + str(self.job_id) + "</a>")

    def error_text_stub(self):
        if self.error_text:
            return self.error_text[:30] + '...'
        
    def __str__(self, ):
        return self.job_id
    def get_log(self):
        """Get the content of log file"""
        log_filename = os.path.join(JOBS_PATH,self.job_id, 'log.txt')
        if os.path.exists(log_filename):
            log = open(log_filename).read()
            return log.strip().split('>')[-1][1:].replace(' ', '&nbsp')
        else:
            return None

    def get_log_reverse(self):
        """Get the content of log file"""
        log_filename = os.path.join(JOBS_PATH,self.job_id, 'log.txt')
        if os.path.exists(log_filename):
            # reverse order of showing log lines
            log = ''
            for l in open(log_filename):
                #       39% |########               | ETA:  00:00:20
                if 'ETA' in l:
                    continue
                else:
                    log += l #+ '\n'
            lines = log.split('\n') # .strip().split('>')[-1][1:].replace(' ', ' ').split('\n')
            lines.reverse()
            return '\n'.join(lines)
        else:
            return None

    def get_input(self):
        fn = os.path.join(JOBS_PATH,self.job_id, 'input.fa')
        if os.path.exists(fn):
            return open(fn).read()
        else:
            return None

    def get_restraints(self):
        fn = os.path.join(JOBS_PATH,self.job_id, self.restraints_fn)
        if os.path.exists(fn):
            return '\n'.join(open(fn).read().split('\n')[:5]) + '\n(...)'
        else:
            return None

    def get_processing(self):
        """Get the content of log file"""
        log_filename = os.path.join(JOBS_PATH,self.job_id, 'processing.txt')
        if os.path.exists(log_filename):
            return open(log_filename).read()
        else:
            return None

    def get_nframes(self,cluster=False):
        """Get the content of log file"""
        if cluster:
            log_filename = os.path.join(SIMRNA_WORKSPACE,self.job_id, 'processing.txt')
        else:
            log_filename = os.path.join(JOBS_PATH,self.job_id, 'processing.txt')            
        if os.path.exists(log_filename):
            for l in open(log_filename):
                if l.startswith('nframes'):
                   return l.split()[1]  #nframes: 322
        else:
            return None

    def get_error(self):
        """Get the content of log file"""
        log_filename = os.path.join(JOBS_PATH,self.job_id, 'log.txt')
        if os.path.exists(log_filename):
            error_box_true = False
            error_box = ''
            for l in open(log_filename):
                if l.find('ERROR') > -1:
                    return "ERROR:\n" + ''.join(open(log_filename).readlines()[-3:-1])
                if l.startswith('Traceback'):
                    error_box_true = True
                if error_box_true:
                    error_box += l
            return error_box
        else:
            return self.error_text

    def get_setup_ok(self):
        """Get the content of log file"""
        log_filename = os.path.join(JOBS_PATH,self.job_id, 'log.txt')
        if os.path.exists(log_filename):
            error_box_true = False
            error_box = ''
            for l in open(log_filename):
                if l.find('It seems that the setup is correct') > -1:
                    return True
        return False

    def get_aa_files(self):
        """Get a list like: [u'jobname_traj_0.200.pdb', u'jobname_traj_200.288.pdb']"""
        path = os.path.join(JOBS_PATH,self.job_id,'output_PDBS')
        if os.path.exists(path):        
            files = os.listdir(path)
            aa_files = []
            for f in files:
                if f.find('AA')>-1:
                    aa_files.append(f)
            aa_files.sort()
            return aa_files
        else:
            return None
                    
    def get_ss_detected_files(self):
        """Get a list like: [u'jobname_traj_0.200.pdb', u'jobname_traj_200.288.pdb']"""
        path = os.path.join(JOBS_PATH,self.job_id,'output_PDBS')

        if os.path.exists(path):        
            files = os.listdir(path)
            aa_files = []
            for f in files:
                if f.find('ss_detected')>-1:
                    aa_files.append(f)
            aa_files.sort()

        if os.path.exists(path):        
            ss_detected = []
            for f in aa_files:
                ss = open(os.path.join(JOBS_PATH,self.job_id,'output_PDBS', f)).read()
                ss_detected.append([f, ss.strip()])
            return ss_detected

    def get_status(self):
        """Get status, running/waiting etc."""
        return JOB_STATUS[self.status][1]

    def get_saxs_format(self):
        """Get status, running/waiting etc."""
        return SAXS_FORMAT[self.saxs_format][1] #  (1, 'crysol')

    def get_traj(self):
        """Get a list like: [u'jobname_traj_0.200.pdb', u'jobname_traj_200.288.pdb']"""
        files = os.listdir(os.path.join(JOBS_PATH,self.job_id))
        traj_files = []
        for f in files:
            if f.find('traj')>-1:
                traj_files.append(f)
        return traj_files

    ##   owner = models.ForeignKey(User, null=True, blank=True)
    def check_files_exits(self):
        """Return False if folder is remove, otherwise true"""
        try:
            os.listdir(os.path.join(JOBS_PATH,self.job_id))
        except OSError:
            self.status = JOB_STATUSES['removed']
            self.save()
            return False
        return True

    def get_progress(self, cluster=False):
        curr_nframes = self.get_nframes(cluster)
        if curr_nframes:
            progress = str(round((int(curr_nframes)-1) / (float(self.nsteps) * 80)*100,2)) # 90 %
        else:
            progress = 0
        self.progress = progress
        self.save()
        return self.progress
