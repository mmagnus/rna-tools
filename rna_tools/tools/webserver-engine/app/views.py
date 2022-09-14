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

from configparser import ConfigParser
import re
import uuid
import subprocess
import sys
import tempfile
import shutil
import json
import os,re,sys
import zipfile
from io import StringIO
from django.http import JsonResponse

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

try:
    from app.models import Job, JOB_STATUSES, SAXS_FORMATS
except:
    print("App is missing")

#from lib.yapdb_parser.pdb_parser_lib import *
#from lib.rna_pdb_edit_occupancy import rna_pdb_edit_occupancy
#from lib.rna_pdb_tools.rna_pdb_tools.utils.rna_convert_pseudoknot_formats.rna_ss_pk_to_simrna import get_multiple_lines, is_pk

import string
from urllib.parse import unquote
import glob
import os

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
    print('request.path', request.path)
    error = ''
    return render_to_response('home.html', RequestContext(request, {
        'load': ''
    }))

def tools(request):
    return render_to_response('tools.html', RequestContext(request, {
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

    from io import BytesIO            
    buffer = BytesIO()
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


def about(request):
    return render_to_response('about.html', RequestContext(request, {}))


def help(request):
    return render_to_response('help.html', RequestContext(request, {}))

def notes(request, fn):
    note = open('/home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/snippets/' + fn + '.txt', encoding="utf-8").read()
    return render_to_response('notes.html', RequestContext(request, {
        'note' : note}))

def qr(request, fn):
    f = open('qr/' + fn, "rb")
    return HttpResponse(f.read(), content_type="image/jpeg")

def image(request, fn):
    f = open('/home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/images/' + fn, "rb")
    return HttpResponse(f.read(), content_type="image/jpeg")

def ssl(request, fn):
    f = open('/home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/ssl/' + fn, "r")
    return HttpResponse(f.read())# , content_type="image/jpeg")


def contact(request):
    return render_to_response('contact.html', RequestContext(request, {}))



def demo(request, tool, job_id):
    print(job_id)
    try:
        j = Job.objects.get(job_id=job_id.replace('/', ''))
    except:  # DoesNotExist:  @hack
        return render_to_response('dont_exits.html', RequestContext(request, {
        }))
    import os
    job_dir = settings.JOBS_PATH + sep + job_id + '/'
    p = settings.PATH + '/app/static/app/demo/'
    #demo
    if tool == 'topdb':
        f = '1xjr.cif'
        shutil.copyfile(p + f, job_dir + f)        

    if tool == 'tocif':
        f = '1xjr.pdb'
        shutil.copyfile(p + f, job_dir + f)        

    if tool in ['qrnas']:
        f = 'MissingAtomsAdded_rpr.pdb' # tetraloop_mdr.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool in ['min']:
        f = 'tetraloop_mdr.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool == 'h2a':
        f = 'gtp.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool == 'analysis':
        f = '5k7c_clean_onechain_renumber_as_puzzle_srr.pdb' #5k7c_clean_onechain_renumber_as_puzzle_srr_min.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool in ['rpr', 'mdr', 'mutate']:
        f = 'missing_op1op2_r43OK.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool in ['seq', 'ss', 'calc-rmsd', 'calc-inf', 'assess', 'clarna']:
        fs = ['21_3dRNA_1_rpr.pdb', '21_Adamiak_1_rpr.pdb', '21_ChenHighLig_1_rpr.pdb',
              '21_Das_1_rpr.pdb', '21_solution_2_rpr.pdb']
        for f in fs:
            shutil.copyfile(p + f, job_dir + f)
    if tool in ['extract', 'delete', 'edit', 'swap']:
        f = 'yC_5lj3_Exon_Intron_rpr.pdb' # yC_5lj3_Cwc2_Exon_Intron.pdb'
        shutil.copyfile(p + f, job_dir + f)        
    if tool in ['cat']:
        fs = ['yC_5lj3_CWC2_C.pdb',
              'yC_5lj3_Exon_C.pdb',
              'yC_5lj3_Intron_C.pdb']
        for f in fs:
            shutil.copyfile(p + f, job_dir + f)

    if tool in ['diffpdb']:
        fs = ['5k7c_clean_onechain_renumber_as_puzzle_srr_min.pdb',
              'pistol_thrs0.50A_clust272-000001_AA_min.pdb']
        for f in fs:
            shutil.copyfile(p + f, job_dir + f)

    if tool in ['rpl']:
        fs = ['triple-ACA.pdb',
              'to-replace.pdb']
        for f in fs:
            shutil.copyfile(p + f, job_dir + f)

    return JsonResponse({})
    
def run(request, tool, job_id):
    print(job_id)
    try:
        j = Job.objects.get(job_id=job_id.replace('/', ''))
    except:  # DoesNotExist:  @hack
        return render_to_response('dont_exits.html', RequestContext(request, {
        }))

    import os
    job_dir = settings.JOBS_PATH + sep + job_id

    job_id = job_id.replace('/', '')

    if tool == 'clear': #clear
        files = glob.glob(job_dir + "/*")
        for f in files:
            os.remove(f)
        return ############## !!!!!!!!!!!!!

    if tool == 'cat':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('rm ' + job_id + '.pdb\n')
             f.write('cat *.pdb > ' + job_id + '.pdb\n')#&> log.txt \n')
             
    if tool == 'seq':
        fasta = request.GET['fasta'].strip()
        if fasta == 'true':
            with open(job_dir + '/run.sh', 'w') as f:
                 f.write('rna_pdb_tools.py --fasta --get-seq *.pdb &> log.txt\n')
        else:
            with open(job_dir + '/run.sh', 'w') as f:
                 f.write('rna_pdb_tools.py --get-seq *.pdb &> log.txt\n')
                 #f.write('echo "ANOTHER FORMAT:" >> log.txt\n')
                 #f.write("rna_pdb_tools.py --get-seq --uniq '[:10]' --compact  *.pdb | sort &>> log.txt\n") # --chain-first#

    if tool == 'rpl':
        files = glob.glob(job_dir + "/*pdb")
        files.sort(key=os.path.getmtime)
        files = [os.path.basename(f) for f in files]
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
rna_pdb_replace.py %s %s &> log.txt\n
""" % (files[0], ' '.join(files[1:])))
             f.write("head *.pdb &> log.txt\n")

    if tool == 'ss':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('rna_pdb_tools.py --get-ss *.pdb &> log.txt\n')

    if tool == 'tocif':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('rna_pdb_tools.py --pdb2cif *.pdb > log.txt\n\n')
             
    if tool == 'topdb':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('rna_pdb_tools.py --cif2pdb *.cif > log.txt\n\n')

    if tool == 'analysis':
        with open(job_dir + '/run.sh', 'w') as f:
             # f.write('rna_x3dna.py -l *.pdb &>> log.txt\n')
            f.write('rna_x3dna.py -l *.pdb > log.txt\n')

    if tool == 'minmd':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("echo 'OpenMM does not give progress for minimization, just wait...' > log.txt \n")
             f.write('time python /home/ubuntu/rna-tools/rna_tools/tools/md/rna_minimize.py -sp *.pdb &>> log.txt\n')

    if tool == 'qrnas':
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('time rna_refinement.py -i *.pdb 2>&1 | tee log.txt\n')
             
    if tool == 'assess':
        with open(job_dir + '/run.sh', 'w') as f:
             #f.write('rna_mq_collect.py -t RASP *.pdb -m 0 -f -o mq.csv | tee log.txt\n')
             f.write('rna_mq_rasp.py *.pdb 2>&1\n')# | tee log.txt\n')
             f.write('rna_csv_sort.py --col rasp_all rasp.csv\n')
             f.write('cat rasp_sorted_rasp_all.csv > log.txt\n')
             f.write('rna_mq_dfire.py *.pdb 2>&1\n')# | tee -a log.txt\n')
             f.write('rna_csv_sort.py --col dfire dfire.csv\n')
             f.write('cat dfire_sorted_dfire.csv >> log.txt\n')
             f.write('rna_merge_two_dfs.py dfire.csv rasp.csv fn --force-writing-output >> log.txt\n')
            
    if tool == 'seq-search':
        print('run, seq,' + job_id)
        seq = request.GET['seq'].strip()       
        with open(job_dir + '/run.sh', 'w') as f:
             f.write('rna_pdb_tools.py --get-seq *.pdb &> log.txt\n')

    if tool == "rpr":
        r = request.GET['renumber'].strip()
        opt = ''
        if r == 'true':
            opt = ' --renumber-residues '
        with open(job_dir + '/run.sh', 'w') as f:
             #  2> log.txt
             f.write("for i in *.pdb; do rna_pdb_tools.py " + opt + " --get-rnapuzzle-ready $i &> ${i/.pdb/_rpr.pdb}; done\n")
             #f.write("echo '== _rpr.pdb files created ==' >> log.txt \n")
             f.write("grep 'REMARK 250  - ' *_rpr.pdb &>> log.txt\n")
             f.write("head *_rpr.pdb &> log.txt\n")

    if tool == "h2a":
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("for i in *.pdb; do rna_pdb_tools.py --replace-htm $i &> ${i/.pdb/_h2a.pdb}; done\n")
             f.write("head *.pdb >> log.txt\n") #for i in *.pdb; do echo -e '\n$i\n'; cat $i; done )             

    if tool == "mdr":
        with open(job_dir + '/run.sh', 'w') as f:
            
            f.write("for i in *.pdb; do rna_pdb_tools.py --mdr $i > ${i/.pdb/_mdr.pdb}; done\n")
            #f.write("echo '== _mdrpr.pdb files created ==' >> log.txt \n")

    if tool == 'calc-rmsd': #rmsd
        allvsall = request.GET['allvsall'].strip()
        if allvsall == 'true':
            
            with open(job_dir + '/run.sh', 'w') as f:
                 #f.write("""rna_calc_rmsd_all_vs_all.py -i ../%s -o rmsds.csv -m all-atom &> log.txt\n""" % (job_id))
                 f.write("""rna_calc_rmsd_all_vs_all.py -i . -o rmsds.csv -m all-atom &> log.txt\n""")# % (job_id))
                 #f.write('ls *.pdb >> log.txt\n\n')
                 #f.write('cat rmsds.csv >> log.txt\n\n')
        else:

            target = request.GET['target'].strip()
            targetselection = request.GET['targetselection'].strip()
            modelselection = request.GET['modelselection'].strip()
            targetignoreselection = request.GET['targetignoreselection'].strip()
            modelignoreselection = request.GET['modelignoreselection'].strip()

            if targetselection:
                targetselection = " --target-selection " + targetselection
                
            if modelselection:
                modelselection = " --model-selection " + modelselection

            if targetignoreselection:
                targetignoreselection = " --target-ignore-selection " + targetignoreselection

            if modelignoreselection:
                modelignoreselection = " --model-ignore-selection " + modelignoreselection

            files = glob.glob(job_dir + "/*pdb")
            files.sort(key=os.path.getmtime)
            files = [os.path.basename(f) for f in files]

            if not len(files):
                with open(job_dir + '/run.sh', 'w') as f:
                    f.write('echo "Empty folder, add files" >> log.txt\n')
            else:
                if not target:
                    target = files[0]
                    files = files[1:]
                    with open(job_dir + '/run.sh', 'w') as f:
                             f.write("""
                rna_calc_rmsd.py -sr -t """ + target + targetselection + modelselection + targetignoreselection + modelignoreselection + """ %s &> log.txt\n
                """ % ' '.join(files))
                             
                else:
                    try:
                        files.remove(target)
                    except:
                        with open(job_dir + '/run.sh', 'w') as f:
                            f.write('echo "Wrong name for target. Is %s OK?" >> log.txt\n' % target)
                    else:
                        with open(job_dir + '/run.sh', 'w') as f:
                             f.write("""
                rna_calc_rmsd.py -sr -t """ + target + targetselection + modelselection + targetignoreselection + modelignoreselection + """ %s &> log.txt\n
                """ % ' '.join(files))
                     #f.write('ls *.pdb >> log.txt\n\n')
                     #f.write('cat rmsds.csv >> log.txt\n\n')
             
    if tool == 'calc-inf':
        files = glob.glob(job_dir + "/*pdb")
        files.sort(key=os.path.getmtime)
        files = [os.path.basename(f) for f in files]
        target = request.GET['target'].strip()

        if not target:
            target = files[0]
            files = files[1:]                
        else:
            files.remove(target)

        #print("\n".join(files))
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
rna_calc_inf.py -t %s %s &> log.txt\n
""" % (target, ' '.join(files)))
             #f.write('ls *.pdb >> log.txt\n\n')
             f.write('column -s, -t < inf.csv >> log.txt\n\n')

    if tool == 'clarna':
        files = glob.glob(job_dir + "/*pdb")
        files.sort(key=os.path.getmtime)
        files = [os.path.basename(f) for f in files]
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
for i in *pdb; do echo $i; rna_clarna_run.py -ipdb $i; done >> log.txt;
""")
             #f.write('ls *.pdb >> log.txt\n\n')
             #f.write('column -s, -t < inf.csv >> log.txt\n\n')

    if tool == 'diffpdb':
        files = glob.glob(job_dir + "/*pdb")
        files.sort(key=os.path.getmtime)
        files = [os.path.basename(f) for f in files]
        r = request.GET['names'].strip()
        opt = ''
        if r == 'true':
            opt = ' --names '
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
diffpdb.py """ + opt + """ --method diff %s %s &> log.txt \n
""" % (files[0], files[1]))
             
    if tool == 'extract':
        opt = request.GET['extract'].strip()
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""for i in *.pdb; do rna_pdb_tools.py --extract '%s' $i > ${i/.pdb/_extract.pdb}; done;""" % opt)
             f.write("rna_pdb_tools.py --get-seq *.pdb >> log.txt")

    if tool == 'delete':
        opt = request.GET['del'].strip()
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
for i in *.pdb; do rna_pdb_tools.py --delete '%s' $i > ${i/.pdb/_delete.pdb}; done;
""" % opt)

    if tool == 'edit':
        opt = request.GET['edit'].strip()
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
for i in *.pdb; do rna_pdb_tools.py --edit '%s' $i > ${i/.pdb/_edit.pdb}; done;
""" % opt)
             f.write("rna_pdb_tools.py --get-seq *.pdb >> log.txt\n")

    if tool == 'swap':
        opt = request.GET['swap'].strip()
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("for i in *.pdb; do rna_pdb_tools.py --swap-chains '%s' $i > ${i/.pdb/_swap.pdb}; done\n" % opt)
             f.write("rna_pdb_tools.py --get-seq *.pdb >> log.txt\n")

    if tool == 'mutate':
        opt = request.GET['mutate'].strip()
        with open(job_dir + '/run.sh', 'w') as f:
             f.write("""
for i in *.pdb; do rna_pdb_tools.py --mutate '%s' $i > ${i/.pdb/_mutate.pdb}; done;
""" % opt)
             f.write('ls *.pdb > log.txt\n')

             # ; rm ".done"; rm log.txt; 

    with open(job_dir + '/run.sh', 'a') as f:
        #f.write("echo 'http://rna-tools.online/tools/%s/%s' >> log.txt\n" % (tool, job_id))
        f.write("echo 'DONE!' >> log.txt\n")
    # add conda at the beginning
    run = open(job_dir + '/run.sh').read()
    with open(job_dir + '/run.sh', 'w') as f:
        #f.write('# <a href="http://rna-tools.online/tools/%s/%s">http://rna-tools.online/tools/%s/%s</a> >> log.txt \n' % (tool, job_id, tool, job_id))
        f.write("# http://rna-tools.online/tools/%s/%s\n" % (tool, job_id)) 
        f.write('source ~/.env\n')
        f.write('rm log.txt\n')        
        # f.write('echo `date` >> log.txt\n') # does not work
        f.write(run)
    #os.system('cd %s && chmod +x run.sh && /bin/zsh && /bin/zsh source ~/.zshr && /bin/zsh run.sh && touch ".done" &' % job_dir)
    # sh: 1: source: not found
    #os.system('cd %s && chmod +x run.sh && /bin/zsh && source ~/.zshrc && /bin/zsh run.sh && touch ".done" &' % job_dir)
    os.system('cd %s && chmod +x run.sh && /bin/zsh run.sh && touch ".done" &' % job_dir)
    j.status = JOB_STATUSES['running']
    j.save()

    return JsonResponse({'post':'false'})

def file_upload(request, job_id):
    if request.method == 'POST':
        my_file=request.FILES.get('file')
        fn = re.sub(r'[\\/*?:"<>|]',"", str(my_file))
        with open(settings.JOBS_PATH + sep + job_id + '/' + fn , 'wb+') as destination:
            for chunk in my_file.chunks():
                destination.write(chunk)
    return JsonResponse({'post':'false'})

def fetch(request, job_id):
    fetch = request.GET['fetch'].strip()
    src = settings.JOBS_PATH + sep + fetch
    dst = settings.JOBS_PATH + sep + job_id
    #f = settings.JOBS_PATH + sep + job_id
    cmd = 'cp -v %s/*.pdb %s' % (src, dst)
    os.system(cmd)
    return JsonResponse({'post':'false'})


def ajax_rm_file(request, job_id_fn):
    f = settings.JOBS_PATH + sep + job_id_fn
    os.remove(f)
        
def ajax_job_status(request, job_id, tool=''):
    job_dir = settings.JOBS_PATH + sep + job_id
    job_id = job_id.replace('/', '')

    print('ajax', job_id, tool)
    #if os.path.exists(job_dir + '/.done'):
    #    #return JsonResponse({'post':'false'})

    import enum
    # Enum for size units
    class SIZE_UNIT(enum.Enum):
       BYTES = 1
       KB = 2
       MB = 3
       GB = 4
    def convert_unit(size_in_bytes, unit):
       """ Convert the size from bytes to other units like KB, MB or GB"""
       if unit == SIZE_UNIT.KB:
           return size_in_bytes/1024
       elif unit == SIZE_UNIT.MB:
           return size_in_bytes/(1024*1024)
       elif unit == SIZE_UNIT.GB:
           return size_in_bytes/(1024*1024*1024)
       else:
           return size_in_bytes
    def get_file_size(file_name, size_type = SIZE_UNIT.BYTES ):
       """ Get file in size in given unit like KB, MB or GB"""
       size = os.path.getsize(file_name)
       return round(convert_unit(size, size_type), 2)

    #try:
    #    tool = request.GET['tool']
    #except:
    #    tool = ''    
    response_dict = {'reload': False}

    if 1:
        j = Job.objects.get(job_id=job_id.replace('/', ''))

        #if j.status == JOB_STATUSES['running']:
        #    response_dict['reload']=False
        #if j.status == JOB_STATUSES['finished'] or j.status == JOB_STATUSES['stopped']:
        #    response_dict['reload'] = False
        #    return JsonResponse({'post':'false'})

        files = glob.glob(job_dir + "/*")
        files.sort(key=os.path.getmtime)
# LOG
        log = ''
        # http://rna-tools.online/tools/calc-rmsd/
        log = '<title>%s</title>JOB ID %s <a href="http://rna-tools.online/tools/%s/%s">http://rna-tools.online/tools/%s/%s</a><br>' % (j,j, tool, j, tool, j) + log
        if files:
           # FILES
           log += "<br>"
           for f in files:
               bf = os.path.basename(f)
               # The raw output files for each step of the pipeline can be found <a href="{{ {{ j.job_id }}.zip">here</a>
               # target="_blank" # for PDB files it open an empty page
               # this is not needed
               size = str(get_file_size(f, SIZE_UNIT.KB)) + ' KB'
               bfl = size.rjust(10).replace(' ', '&nbsp;') + ' <a href="/media/jobs/' + job_id + '/' + bf +'">' + bf + '</a>' # .ljust(60).replace(' ', '&nbsp') 
               log += '<a class="icon-remove-circle" href="#" onclick="del(\'' + job_id + '/' + bf + '\');"></a> ' +  bfl + '<br>'
        # log += "<br>== FILES END ==<br>"

        try:
            # SCRIPT
            with open(os.path.join(settings.JOBS_PATH, job_id, 'run.sh')) as f:
                 log += "<br>" + f.read().replace('\n', "<br>")# + "<br>== SCRIPT END ==<br>"
        except FileNotFoundError:
            pass
        
        try:
            log_filename = os.path.join(settings.JOBS_PATH, job_id, 'log.txt')
            with open(log_filename, 'r') as ifile:
                l = ifile.read()
                log += re.sub(r"[\n]", "<br>", l)
                # http://rna-tools.online/tools/calc-rmsd/
                # http://rna-tools.online/tools/calc-rmsd/
                # '<title>%s</title><a href="%s">%s</a><br>' % (j,j,j) + log
                log = log.replace('source ~/.env<br>', '')
                log = log.replace('cat rmsds.csv', '')
                log = log.replace('****************************************************************************', '<hr>')
                log = log.replace('-----------------------------------------------------------------------------', '<hr>')
                log = log.replace('PyMOL not running, entering library mode (experimental)', '')
                log = log.replace("echo 'DONE!'", "")
                log = log.replace("rm log.txt", "")

                nl = ''
                for l in log.split('<br>'):
                     if l.startswith('echo'):
                         continue
                     if 'Performing minimization step:' in l:
                         pass
                     elif 'Writing PDB file:' in l:
                         pass
                     elif 'Internal Warning:' in l:
                         pass
                     elif 'residue renamed' in l:
                         pass
                     elif 'Missing atom added' in l:
                         pass
                     else:
                         nl += l + '<br>'
                log = nl
                #log = log.replace('>>', '>&gt;') # for fasta
                #log = '<a href="#tetraloop_helix.pdb">tetraloop_helix.pdb</a>' + log
                all = re.findall('# .*? #', log) # # f #
                for a in all:
                    # # f.pdb #
                    inside = a.replace('#', '').strip()
                    na = a.replace('# ', '<h3 id="#' + inside + '">')
                    na = na.replace(' #', '</h3>')
                    log = log.replace(a, na)
                    #log = re.sub(r"# ", "<h3>", log)
                    #log = re.sub(r"#<br>", "</h3>", log)
                # collect H3 and make 
                print(log)
                
                #log += re.sub(r"^<br>%", "", log)
                # --> Clustering
                #
                # --> Annealing
                #log = re.sub(r"[-]+> Annealing[\w\s]+\d+\%[\s\|#]+ETA:[\d\s\:\-]+\r", "", log)
                # --> Preparing data
                #log = re.sub(r"[\s]{5,}\d+\%[\s\|#]+ETA:[\d\s\:\-]+\r", "", log)
                #log = re.sub(r"[\-]+>[\w\s]+\d+\%[\s\|#]+ETA:\s+[(\d\-)\:]+\r[^\Z]", "", log)
                #log = re.sub(r"[\s]{4,}\d+\%[\s\|#]+ETA:\s+[(\d\-)\:]+\r[^\Z]", "", log)

                if tool == 'cat': #log
                    log += '</div><pre>RESULTS<br>'
                    log += '<a href="/media/jobs/' + job_id + '/' + job_id  + '.pdb">' + job_id  + '.pdb</a><br>'
                    log += '</pre>'

                if tool in ['extract', 'delete', 'rpr', 'mutate', 'mdr', 'min', 'h2a', 'rpl']:
                    files = glob.glob(job_dir + "/*_" + tool + ".pdb")
                    files = [os.path.basename(f) for f in files]
                    log += '</div>RESULTS<br>' # <pre
                    if files:
                        for f in files:
                            log += '<a href="/media/jobs/' + job_id + '/' + f + '">' + f + '</a><br>'
                    log += '</pre>'

## <script type="text/javascript">
##     Info.script=`load {{ MEDIA_URL }}jobs/{{ j.job_id }}/21_ChenHighLig_1_rpr.pdb; color [x00b159]; cartoons only; \
##     load APPEND {{ MEDIA_URL }}jobs/{{ j.job_id }}/21_Adamiak_1_rpr.pdb;select 2.1; color [xd11141]; cartoons only; frame *;`
##   jmolApplet0 = Jmol.getApplet("jmolApplet0", Info)
## </script>
## <p><b>Target [green] 21_ChenHighLig_1_rpr.pdb vs Model [red] 21_Adamiak_1_rpr</b></p>
## </center>


## """
                    # log += '<center><br><a href="/tools/%s/%s">Click get JSmol view</a></center>' % (tool, j)

                    #log += j.get_status()
                
                if 'DONE!' in log:
                    if j.get_status() != 'finished' and j.get_status() != 'finished with errors': # of other then just leave it
                        import time
                        time.sleep(1)
                        j = Job.objects.get(job_id=job_id.replace('/', ''))
                        if 'Error' in log or 'Wrong' in log:
                            j.status = JOB_STATUSES['finished with errors']
                        else:
                            j.status = JOB_STATUSES['finished']
                        j.save()
                        response_dict['reload'] = True # False # True # fix it

        except FileNotFoundError:
            pass

    # clea up log
    log = log.replace('REMARK 250 ', '')
    log = log.replace(' &>> log.txt', '')
    log = log.replace(' &> log.txt', '')    
    log = log.replace(' >> log.txt', '')    
    log = log.replace(' > log.txt', '')
    log = log.replace("echo 'DONE'", '')     #log
    log = log.replace("DONE!", '')     #log

    log = log.strip()

    #if j.get_status() != 'finished':
    #    log += '</br>(...)'

    #if os.path.exists(job_dir + '/.done'):
    #   log += '<span class="label label-success">DONE</span>'
    response_dict['log'] = log
    #if j.comment != log: # so this is different, reload
    if 1:
        j.comment = log
        j.save()
        return HttpResponse(json.dumps(response_dict), "application/json") # update only 
        # log is different
    #else:
    #     return HttpResponse({}) #json.dumps(response_dict), "application/json")
    return JsonResponse({'post':'false'})


def tool(request, tool, job_id):
    if not job_id:
        print('create a job')
        id = str(uuid.uuid4()).split('-')[0]  # name plus hash
        j = Job()
        j.job_id = id
        j.status = 0
        print('tools, make:', j)
        j.save()
        # create folder
        try:
            job_dir = settings.JOBS_PATH + sep + j.job_id
            umask(0o002)
            makedirs(job_dir)
        except OSError:
            pass
        log = ''
    else:
        try:
            j = Job.objects.get(job_id=job_id.replace('/', ''))
        except:  # DoesNotExist:  @hack
            return render_to_response('dont_exits.html', RequestContext(request, {
            }))
        job_dir = settings.JOBS_PATH + sep + job_id
        try:
            with open(job_dir + '/_xxlog.txt') as f:
                log = f.read()
        except FileNotFoundError:
            log = ''

    jsmol = ''
    target = ''
    topmodel = ''
    out = ''
    headers = ''
    caption = ''
    if 1:
    ## try:
    ##     import csv
    ##     csv_fp = open(job_dir + '/rmsds.csv', 'r')
    ##     reader = csv.DictReader(csv_fp)
    ##     headers = [col for col in reader.fieldnames]
    ##     headers = ['Filename', 'RMSD']
    ##     out = []
    ##     for row in reader:
    ##         out.append((row['fn'], row['rmsd_all']))
    ##     print(out)
        # get target

        # target: 21_Das_1_rpr.pdb
        try:
            for l in open(job_dir + '/log.txt', 'r'):
                if '# target: ' in l:
                    target = l.replace('# target: ', '')
        except: # FileNotFoundError:
            pass

        try:
            for l in open(job_dir + '/rmsds.csv', 'r'):
                if '-1.0' in l:
                    continue
                if '.pdb' in l:
                    topmodel = l.split(',')[0] # .replace('# target: ', '')
                    break
        except: # FileNotFoundError:
            pass

        
        #target = '21_Adamiak_1_rpr.pdb'
        #topmodel = '21_ChenHighLig_1_rpr.pdb'
        caption =  '<b>Target [green] ' + target + ' vs Top Model [red] ' + topmodel + '</b>'
##         jsmol = """
## <script type="text/javascript">
##     Info.script=`load {{ MEDIA_URL }}jobs/{{ j.job_id }}/21_ChenHighLig_1_rpr.pdb; color [x00b159]; cartoons only; \
##     load APPEND {{ MEDIA_URL }}jobs/{{ j.job_id }}/21_Adamiak_1_rpr.pdb;select 2.1; color [xd11141]; cartoons only; frame *;`
##   jmolApplet0 = Jmol.getApplet("jmolApplet0", Info)
## </script>
## <p><b>Target [green] 21_ChenHighLig_1_rpr.pdb vs Model [red] 21_Adamiak_1_rpr</b></p>
## </center>
## """
##         jsmol = ''
    ## except:
    ##      out = ''
    ##      headers = ''

    ref = ''
    model = ''
    rows = ''
    headers = ''
    import time

    if tool == 'calc-inf':
        time.sleep(1)
        f = job_dir + '/inf.csv'
        if os.path.exists(f):
            import csv
            csv_file = open(f, 'r')
            csv_reader = csv.reader(csv_file, delimiter=',')
            c = 0
            rows = []
            for row in csv_reader:
                if c == 0:            
                    headers = row[1:] # skip target
                else:
                    rows.append(row[1:])
                c += 1
            ic(rows)

    if tool == 'calc-rmsd':
        time.sleep(1)
        f = job_dir + '/rmsds.csv'
        if os.path.exists(f):
            try:
                for l in open(job_dir + '/log.txt', 'r'):
                    if 'target: ' in l:
                        target = l.replace('target: ', '')
            except: # FileNotFoundError:
                pass

            import csv
            csv_file = open(f, 'r')
            csv_reader = csv.reader(csv_file, delimiter=',')
            c = 0
            rows = []
            for row in csv_reader:
                if c == 0:            
                    headers = row[:] # skip target
                else:
                    rows.append(row[:])
                c += 1

    if tool == 'assess':
            time.sleep(1)
            f = job_dir + '/merged.csv'
            if os.path.exists(f):
                import csv
                csv_file = open(f, 'r')
                csv_reader = csv.reader(csv_file, delimiter=',')
                c = 0
                rows = []
                for row in csv_reader:
                    if c == 0:            
                        headers = row[:] # skip target
                    else:
                        rows.append(row[:])
                    c += 1

    if tool in ['mutate', 'delete', 'edit', 'qrnas', 'min']:
        time.sleep(1)
        files = glob.glob(job_dir + "/*pdb")
        files.sort(key=os.path.getmtime)
        if len(files) > 1:
            files = [os.path.basename(f) for f in files]
            ref = files[0] # job_dir + os.sep + 
            model = files[1] # job_dir + os.sep + 

    struc = False
    if tool == 'cat':
        if os.path.exists(job_dir + '/' + job_id + '.pdb'):
            struc = True
        ic(struc)

    return render_to_response(tool + '.html', RequestContext(request, {
        'j': j,
        'log' : log,
        'tool' : tool,
        'data' : out,

        'headers' : headers,
        'rows' : rows,
        
        'jsmol': jsmol,
        'target': target,
        'topmodel' : topmodel,
        'caption' : caption,
        'struc' : struc,
        
        'ref' : ref,
        'model' : model,
        }))
