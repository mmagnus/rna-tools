#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pip install django-ipware

https://stackoverflow.com/questions/4581789/how-do-i-get-user-ip-address-in-django
"""

APP_PATH = '/home/rnamasonry/rnamasonryweb_env/rnamasonry-web/'
LIMIT_PER_USER = 5  # 2
POWER_USERS = ['iamb@genesilico.pl', 'azyla@genesilico.pl', 'magnus@genesilico.pl']

from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext
from django.utils.datastructures import MultiValueDictKeyError
from web import settings
from os import sep, makedirs, listdir, umask

from django.db.models import Q

from ConfigParser import ConfigParser
import re
import uuid
import subprocess
import sys
import tempfile
import shutil
import json
import os,re,sys
import zipfile
from cStringIO import StringIO

try:
    from app.models import Job, JOB_STATUSES, SAXS_FORMATS
except:
    print("App is missing")

#from lib.yapdb_parser.pdb_parser_lib import *
#from lib.rna_pdb_edit_occupancy import rna_pdb_edit_occupancy
#from lib.rna_pdb_tools.rna_pdb_tools.utils.rna_convert_pseudoknot_formats.rna_ss_pk_to_simrna import get_multiple_lines, is_pk

try:
    from ipware.ip import get_ip
except:
    print('Install ipware')

import string

CORES = 80

intervals = (
    ('weeks', 604800),  # 60 * 60 * 24 * 7
    ('days', 86400),    # 60 * 60 * 24
    ('hours', 3600),    # 60 * 60
    ('minutes', 60),
    ('seconds', 1),
)


def display_time(seconds, granularity=2):
    result = []

    for name, count in intervals:
        value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip('s')
            result.append("{} {}".format(value, name))
    return ', '.join(result[:granularity])


def home(request):
    print 'request.path', request.path
    error = ''
    return render_to_response('home.html', RequestContext(request, {
        'load': ''
    }))


def stop(request, job_id):
    """Stop job based on job_id /stop/<job_id>. Get the job, change status to stopped
    and re-direct page to /job/,job_id>
    """
    try:
        j = Job.objects.get(job_id=job_id.replace('/', ''))
    except:  # DoesNotExist:  @hack
        return render_to_response('dont_exits.html', RequestContext(request, {
        }))

    # status stopped but not yet stopped
    # this is signal to demon to stop it and when processed killed
    # then set True .stopped
    # see daemon.py for more
    j.stopped = False 
    j.status = JOB_STATUSES['stopped']
    j.save()
    # shell way to kill it # tRNA_with_SAXS-d5a37d86
    return HttpResponseRedirect('http://genesilico.pl/rnamasonry/jobs/' + j.job_id + '/')  # @hack


def job(request, job_id):
    job_dir = settings.JOBS_PATH + sep + job_id
    load = ''
    
    try:
        j = Job.objects.get(job_id=job_id.replace('/', ''))
    except:  # DoesNotExist:  @hack
        return render_to_response('dont_exits.html', RequestContext(request, {
        }))

    error = j.get_error()
    if error:
        j.error_text = j.get_error()
        j.status = JOB_STATUSES['finished with errors']
        j.save()
        return render_to_response('errors.html', RequestContext(request, {
            'j': j
        }))

    if j.status == JOB_STATUSES['finished with errors']:
        return render_to_response('errors.html', RequestContext(request, {
            'j': j
        }))

    # if j.get_setup_ok() and j.status == JOB_STATUSES['waiting']:
    #      j.status = JOB_STATUSES['running']
    #      j.save()

    nframes = 0
    curr_nframes = j.get_nframes()
    progress = j.get_progress()

    if j.status == JOB_STATUSES['finished'] or j.status == JOB_STATUSES['stopped']:
        try:
            files = listdir(job_dir)
        except OSError:
            return render_to_response('deleted.html', RequestContext(request, {
                'job_id': job_id,
                'job_dir': job_dir,
                'j': j,
            }))

        return render_to_response('result.html', RequestContext(request, {
            'job_id': job_id,
            'job_dir': job_dir,
            'files': files,
            'nframes': nframes,
            'curr_nframes': curr_nframes,
            'progress': progress,
            'j': j,
        }))
    else:
        return render_to_response('progress.html', RequestContext(request, {
            'j': j,
            'progress': progress,
            'nframes': nframes,
            'curr_nframes': curr_nframes,
            'load': load,
            'progress_float': float(progress),
        }))


# -----------------------------------------------------------------------------

def download_project_dir(request, job_id):


    job_dir = settings.JOBS_PATH + sep + job_id
    fname="%s.zip" % job_id

    response = HttpResponse(content_type='application/zip')
    response['Content-Disposition'] = 'filename=%s'%fname
    all_files = []


    try:
        job_status = Job.objects.get(job_id=job_id.replace('/', ''))
    except:
        job_status = None

    for root, dirs, files in os.walk(job_dir):
        for fn in files:
            abs_fn = os.path.join(root,fn)
            #print os.path.relpath(abs_fn, job_dir)
            with open(abs_fn, 'rb') as ifile:
                all_files.append( (os.path.relpath(abs_fn, job_dir), ifile.read()) )

    buffer = StringIO()
    zip = zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED)
    for name, f in all_files:
        zip.writestr( os.path.join(job_id, name), f)
    zip.close()
    buffer.flush()
    #the import detail - we return the content of the buffer
    ret_zip = buffer.getvalue()
    buffer.close()
    response.write(ret_zip)



    return response


# -----------------------------------------------------------------------------

def ajax_job_status(request, job_id):

    job_dir = settings.JOBS_PATH + sep + job_id
    load = ''

    response_dict = {'reload':False}

    try:
        j = Job.objects.get(job_id=job_id.replace('/', ''))

        #if j.status == JOB_STATUSES['running']:
        #    response_dict['reload']=False
        if j.status == JOB_STATUSES['finished'] or j.status == JOB_STATUSES['stopped']:
            response_dict['reload']=True

        log_filename = os.path.join(settings.JOBS_PATH,job_id,'log.txt')
        with open(log_filename, 'r') as ifile:
            log = ifile.read()
            log = re.sub(r"[\n]", "</br>", log)

            # --> Clustering
            #log = re.sub(r"[\-]+> Clustering[\w\s]+\d+\%[\s\|#]+ETA:\s+[(\d\-)\:]+\r", "", log)
            # --> Annealing
            #log = re.sub(r"[-]+> Annealing[\w\s]+\d+\%[\s\|#]+ETA:[\d\s\:\-]+\r", "", log)
            # --> Preparing data
            #log = re.sub(r"[\s]{5,}\d+\%[\s\|#]+ETA:[\d\s\:\-]+\r", "", log)

            log = re.sub(r"[\-]+>[\w\s]+\d+\%[\s\|#]+ETA:\s+[(\d\-)\:]+\r[^\Z]", "", log)
            log = re.sub(r"[\s]{4,}\d+\%[\s\|#]+ETA:\s+[(\d\-)\:]+\r[^\Z]", "", log)


            response_dict['log']=log.replace(' ', '&nbsp')

    except:
        response_dict['log'] = ""


    return HttpResponse(json.dumps(response_dict), "application/json")


# -----------------------------------------------------------------------------


def about(request):
    return render_to_response('about.html', RequestContext(request, {}))


def help(request):
    return render_to_response('help.html', RequestContext(request, {}))


def contact(request):
    return render_to_response('contact.html', RequestContext(request, {}))


def submit(request):
    """Get submit page or if POST then get processing of input files and redirect to /job/XXX"""
    # only for not power users, check if more than 2
    # if request.POST['email'].strip() not in POWER_USERS
    load = ''
    if request.method == 'POST':
        ## if request.POST['email'] not in POWER_USERS and request.POST['email']:
        ##     ip = get_ip(request)
        ##     jobs_per_ip = len(Job.objects.filter((Q(status=JOB_STATUSES['running']) | Q(
        ##         status=JOB_STATUSES['queue']) | Q(status=JOB_STATUSES['waiting'])) & Q(ip=ip)))
        ##     if jobs_per_ip > LIMIT_PER_USER:
        ##         block_the_user = 'true'
        ##     else:
        ##         block_the_user = 'false'
        ## else:
        ##     block_the_user = 'false'
        ##     jobs_per_ip = 1000 # 1
        ## #######################################################
        ## if block_the_user == 'false':  # off this
        ##     return render_to_response('submit.html', RequestContext(request, {
        ##         'error': 'Warning! You are running ==>' + str(jobs_per_ip) + '<==.\nYou have reached the limit of jobs run per user (2 runs per user). \nPlease wait untill all your jobs are done. \nIf we have blocked you by mistake, contact us: simrnaweb@genesilico.pl',
        ##     }))

        if request.POST['captcha'] != 'RNA':
            return render_to_response('submit.html', RequestContext(request, {
                'error': 'If you are not a robot, type RNA as an answer to `What do you want to model`?',
            }))


        j = Job()

        # save ip
        ip = get_ip(request)
        j.ip = str(ip)

        # get job_id till it does not give you job_id that is already in the db
        while 1:
            job_title = request.POST['job_title'].strip()
            if job_title:
                s = job_title.replace(' ', '_')
                # source http://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
                job_id = "".join(x for x in s if x in list(string.ascii_letters) +
                                 list('+_-.') + list(string.digits))[:20]

                job_id = job_id + '-' + str(uuid.uuid4()).split('-')[0]  # name plus hash

                hits = Job.objects.filter(job_id=job_id)
                # gla_x1 gla_x2
                j.job_id = job_id + '_x' + str(len(hits))
            else:
                # e462d8a5-7079-41df-b1bb-25edcb065cca -> e462d8a5
                job_id = str(uuid.uuid4()).split('-')[0]

            if not Job.objects.filter(job_id=job_id):
                j.job_id = job_id
                break

        error = ''

        # create folder
        try:
            JOB_PATH = settings.JOBS_PATH + sep + j.job_id
            umask(0002)
            makedirs(JOB_PATH)
        except OSError:
            pass

        # job
        jt = request.POST['job_title']
        j.job_title = "".join(x for x in s if x in list(string.ascii_letters) +
                              list('#&/\|+_-,.()*> \'":;@%^~?') + list(string.digits))

        j.status = 0
        j.ss_raw = request.POST['ss']  # it's OK! empty string and empty file!
        j.ss = request.POST['ss']

        #format
        j.saxs_format = SAXS_FORMATS[request.POST['saxs_format']]

        #if is_pk(j.ss):
        #    j.ss = get_multiple_lines(j.ss)

        # pdb_fn
        try:
            # mój_ulubiony_struc_!.pdb -> moj_ulubiony_struc_.pdb
            pdb_fn = request.FILES['pdb_fn']
            pdb_fn_name = pdb_fn.name.replace(' ', '_').lower()
            j.pdb_fn = pdb_fn_name
            j.pdb_fn = "".join(x for x in pdb_fn_name if x.isalnum() or x == '_')
            j.pdb_fn = j.pdb_fn.replace('pdb', '.pdb')
            j.pdb_fn = j.pdb_fn.replace('cif', '.cif')
            with open(JOB_PATH + sep + j.pdb_fn, 'w') as pf:
                for chunk in pdb_fn.chunks():
                    pf.write(chunk)

        except MultiValueDictKeyError:
            j.pdb_fn = ''

        # restraints
        # try:
        if request.POST['demo'] == 'restraints':
            j.restraints_fn = "tRNA_saxs.dat"
            shutil.copyfile(
                APP_PATH + '/app/static/app/demo/tRNA_saxs.dat', JOB_PATH + sep + j.restraints_fn)
        else:
            try:
                restraints_fn = request.FILES['saxs']
                restraints_fn_name = restraints_fn.name.replace(' ', '_').lower()
                j.restraints_fn = restraints_fn_name
                j.restraints_fn = "".join(
                    x for x in restraints_fn_name if x.isalnum() or x == '_' or x == '.')
                with open(JOB_PATH + sep + j.restraints_fn, 'w') as pf:
                    for chunk in restraints_fn.chunks():
                        pf.write(chunk)
            except MultiValueDictKeyError:
                j.restraints_fn = ''

        j.interpret_occupancy = 1  # @todo
        j.residues_to_freeze = request.POST['residues_to_freeze'].replace(' ', '')  # remove ' '
        j.email = request.POST['email']
        j.nsteps = request.POST['number_steps']

        j.seq = request.POST['seq']
        j.seq_len = len(j.seq)

        print('j.seq', j.seq)
        
        j.save()

        ## seq = ''
        ## for l in request.POST['seq'].split('\n'):
        ##     if l.find('>') > -1:
        ##         error = "You must an RNA sequence only (not a fasta file, please convert your input into an RNA sequence)!"
        ##         j.status = JOB_STATUSES['finished with errors']
        ##         j.error_text = error
        ##         j.save()
        ##         return render_to_response('submit.html', RequestContext(request, {
        ##             'error': error}))
        ##     seq += l.strip()

        if j.seq:
            with open(JOB_PATH + sep + 'seq.fa', 'w') as f:
                f.write(j.seq)

        if j.ss:
            with open(JOB_PATH + sep + 'seq.ss', 'w') as f:
                for l in j.ss.split('\n'):
                    f.write(l.strip() + '\n')

        with open(JOB_PATH + sep + 'input.fa', 'w') as f:
            f.write('> ' + j.job_title + '\n')
            f.write(j.seq + '\n')
            for l in j.ss.split('\n'):
                    f.write(l.strip() + '\n')

        cfg = ConfigParser()
        cfg.add_section('Job')
        pattern = re.compile(r'\s+')

        cfg.set('Job', 'seq', j.seq)
        if j.ss:
            cfg.set('Job', 'ss', 'yes-there-is-ss-save-to-seq.ss')  # @hack!!!!
        else:
            cfg.set('Job', 'ss', '')  # @hack!!!!
        cfg.set('Job', 'restraints_fn', j.restraints_fn)

        cfg.set('Job', 'percentage_of_frames', j.percentage_of_frames)
        cfg.set('Job', 'interpret_occupancy', j.interpret_occupancy)

        # no seq and no pdb
        if (j.seq and j.ss) or j.pdb_fn:
            pass  # ok
        else:
            error = "You must provide submit an RNA sequence and RNA secondary structure OR an RNA structure!"
            j.status = JOB_STATUSES['finished with errors']
            j.error_text = error
            j.save()
            return render_to_response('submit.html', RequestContext(request, {
                'error': error}))

        # if len(j.seq) > 220:
        if len(j.seq) > 300:
            error = 'Sequence is too long.'
            j.status = JOB_STATUSES['finished with errors']
            j.error_text = error
            j.save()
            return render_to_response('submit.html', RequestContext(request, {
                'error': error}))

        ## # process pdb
        ## if j.pdb_fn:
        ##     pdb_fn = JOB_PATH + sep + j.pdb_fn
        ##     shutil.copy(pdb_fn, pdb_fn + '_raw.pdb')

        ##     pdb_fn_clean = JOB_PATH + sep + j.pdb_fn
        ##     with open(pdb_fn_clean) as f:
        ##         f.write(
        ##     cfg.set('Job', 'pdb_fn', j.pdb_fn)
        ##     j.save()

        ##     if error:
        ##         j.status = JOB_STATUSES['finished with errors']
        ##         j.error_text = error
        ##         j.save()

        ##         return render_to_response('submit.html', RequestContext(request, {
        ##             'error': error}))

        ## else:
        ##     cfg.set('Job', 'pdb_fn', '')

        f = open(JOB_PATH + sep + 'job.cfg', 'w')
        cfg.write(f)
        f.close()

        # @hack
        return HttpResponseRedirect('http://genesilico.pl/rnamasonry/jobs/' + j.job_id + '/') # 
    else:
        error = ''
        return render_to_response('submit.html', RequestContext(request, {
            'error': error,
        }))
