import re,os,sys,time,itertools,tempfile
from datetime import datetime
from xml.dom import minidom
from fabric.api import *
from fabric.contrib import files
import hashlib
import imp
# sys.path.append('../clarna')
from utils import PDBObject,FileNamesObject

# env['hosts']=['twalen@rnacontacts-vm']

# zgodnie z http://stackoverflow.com/questions/3077281/pythons-fabric-connect-to-a-host-listed-ssh-config
# env.use_ssh_config = True 

INSTALLATIONS = {
    'localhost': {
            'host_string': 'localhost',
            'DATA_DIR': os.path.join(os.path.dirname(__file__),'gc-data'),
            'SCRIPTS_DIR': os.path.dirname(__file__),
            'DJANGO_DIR': os.path.join(os.path.dirname(__file__),'web')
            
    },
    'ouzo': {
            'host_string': 'twalen@ouzo',
            'DATA_DIR': '/home/twalen/work/clarna/gc-data',
            'SCRIPTS_DIR': "/home/twalen/work/clarna",
            'DJANGO_DIR': "/home/twalen/work/clarna/web"
    },
    'peyote2': {
            'host_string': 'twalen@peyote2',
            'DATA_DIR': '/home/twalen/gc-data',
            'SCRIPTS_DIR': "/home/twalen/genesilico/clarna",
            'DJANGO_DIR': "/home/twalen/genesilico/clarna/web"
    },
    'peyote2-im': {
            'host_string': 'twalen@peyote2',
            'DATA_DIR': '/home/twalen/gc-data',
            'SCRIPTS_DIR': "/home/twalen/genesilico/clarna",
            'DJANGO_DIR': "/home/twalen/genesilico/clarna/web",
            'QSUB_QUEUE': "immediate.q",
            "QSUB_LIMIT": 200
    },
    'peyote2-aux': {
            'host_string': 'twalen@peyote2',
            'DATA_DIR': '/home/twalen/gc-data',
            'SCRIPTS_DIR': "/home/twalen/genesilico/clarna",
            'DJANGO_DIR': "/home/twalen/genesilico/clarna/web",
            'QSUB_QUEUE': "aux.q",
            "QSUB_LIMIT": 50
    },
    'rnacontacts-vm': {
            'host_string': 'twalen@rnacontacts-vm',
            'DATA_DIR': '/home/twalen/gc-data',
            'SCRIPTS_DIR': "/var/www.rna-contacts/clarna",
            'DJANGO_DIR': "/var/www.rna-contacts/clarna/web"
    }
}

def write_file(fn, contents):
    f = open(fn, "w")
    f.write(contents)
    f.close()

def read_file(fn):
    f = open(fn)
    res = f.read()
    f.close()
    return res

def remote_write_file(fn, contents):
    if env.INS in ['localhost']:
        write_file(fn, contents)
    else:
        f=tempfile.NamedTemporaryFile(delete=False)
        f.write(contents)
        f.close()
        tmp_fn = f.name
        with hide('running','stdout','stderr'):
            dir = os.path.dirname(fn)
            if dir!='':
                run('mkdir -p "%s"'%dir)
            put(tmp_fn, fn)
        os.unlink(tmp_fn)

def remote_read_file(fn):
    if env.INS in ['localhost']:
        return read_file(fn)
    else:
        output = ""
        with hide('running'):
            output = run("cat '%s'" % fn, pty=True, shell=False)
        return output

def remote_mkdir(dirname):
    if env.INS in ['localhost']:
        local("mkdir -p '%s'"%dirname)
    else:
        run("mkdir -p '%s'"%dirname)


def settings_set_ins(ins,force=False):
    if not force and hasattr(env,'INS'):
        return settings()
    ins_conf = INSTALLATIONS[ins]
    return settings(host_string=ins_conf['host_string'],
            INS=ins,
            DATA_DIR=ins_conf['DATA_DIR'],
            SCRIPTS_DIR=ins_conf['SCRIPTS_DIR'],
            DJANGO_DIR=ins_conf['DJANGO_DIR'])

# set_ins task
@task
def set_ins(insname,force=False):
    if not force and hasattr(env,'INS'):
        return
    ins_conf = INSTALLATIONS[insname]
    env.INS = insname
    env.host_string = ins_conf['host_string']
    env.DATA_DIR = ins_conf['DATA_DIR']
    env.SCRIPTS_DIR = ins_conf['SCRIPTS_DIR']
    env.DJANGO_DIR = ins_conf['DJANGO_DIR']
    env.QSUB_QUEUE = ins_conf.get('QSUB_QUEUE')
    env.QSUB_LIMIT = ins_conf.get('QSUB_LIMIT')

# ins decorator
def ins(insname,force=False):
    def _ins_func(f):
        def new_f(*args,**kwargs):
            with settings_set_ins(insname,force):
                return f(*args,**kwargs)
        new_f.func_name = f.func_name
        if hasattr(f,'__doc__'):
            new_f.__doc__ = f.__doc__
        return new_f
    return _ins_func

QSUB_LOG_DIR = "/home/twalen/qsub-logs"
N_TYPES = [n1+n2 for n1,n2 in itertools.product(["A","C","U","G"],repeat=2)]
CL_CAT = {
    "bp":[
        "HH_cis",
        "HH_tran",
        "HS_cis",
        "HS_tran",
        "HW_cis",
        "HW_tran",
        "SH_cis",
        "SH_tran",
        "SS_cis",
        "SS_tran",
        "SW_cis",
        "SW_tran",
        "WH_cis",
        "WH_tran",
        "WS_cis",
        "WS_tran",
        "WW_cis",
        "WW_tran"
    ],
    "stacking":[
        "<<",">>","<>","><"
    ],
    "base-phosphate":[
        "H_0BPh",
        "SW_2BPh",
        "S_1BPh",
        "W_345BPh",
        "W_6BPh",
        "H_789BPh",
    ],
    "base-ribose":[
        "H_0BR",
        "S_1BR",
        "SW_2BR",
        "W_345BR",
        #"W_3BR",
        #"W_4BR",
        #"W_5BR",
        #
        "W_6BR",
        "H_789BR",
        #"H_7BR",
        #"H_8BR",
        #"H_9BR",
        #
    ],
    "other":[
        "diagonal-c",
        "diagonal-nc-ww",
        "long-stacking-c",
    ],
    "other2":[
        "base-ribose-stacking",
    ],
    "other3":[
        "phosphate-stacking1",
        "phosphate-stacking2",
    ]
}

def format_timedelta(diff):
    sec = diff.seconds
    return "%02d:%02d:%02d" % (sec/3600,((sec/60)%60),sec%60)

def format_time(d):
    return d.strftime("%H:%M:%S")

def format_date(d):
    return d.strftime("%y-%m-%d %H:%M:%S")

def _safe(s):
    return s.replace("<","l").replace(">","g")

def cl_tuples():
    return [(sc,desc,n_type) \
        for sc,desc in sum([zip([k]*len(v),v) for k,v in CL_CAT.items()],[]) \
        for n_type in N_TYPES
    ]

class AbstractJob(object):

    def __init__(self,id,cmd=None,num=None,mem=None,np=None,dir=None):
        self.id = id
        self.num = num
        self.cmd = cmd
        self.mem = mem
        self.np  = np
        self.state = '?'
        self.dir = dir
        self.output = None

    def get_execution_result(self):
        return "?"
        
    def run(self):
        raise Exception("Abstract")

    def __str__(self):
        return "Job(%s)" % self.id
        
    def _parse_job_output(self):
        last_line = self.output.strip().split("\n")[-1]
        if last_line in ['OK','ERROR']:
            return last_line
        else:
            return '?'
    
    @staticmethod
    def create(*args,**kwargs):
        if env.host_string=='localhost':
            return LocalJob(*args,**kwargs)
        elif re.match('^.*@peyote2',env.host_string):
            kwargs['queue']=env.get('QSUB_QUEUE')
            return QSubJob(*args,**kwargs)
        elif re.match('^.*@',env.host_string):
            return RemoteJob(*args,**kwargs)
        else:
            raise Exception('unknown host_string: %s' % env.host_string)

class LocalJob(AbstractJob):

    def run(self):
        with hide('running','stdout','stderr'):
            self.output = local("cd '%s' && ( %s ) && echo 'OK' || echo 'ERROR'" % (self.dir, self.cmd), capture=True)
        self.state = 'FINISHED'
        self.result = self._parse_job_output()

    def __str__(self):
        return "LocalJob(%s)" % self.id


class RemoteJob(AbstractJob):

    def run(self):
        with hide('running','stdout','stderr'):
            self.output = run("cd '%s' && ( %s ) && echo 'OK' || echo 'ERROR'" % (self.dir, self.cmd),combine_stderr=True)
        self.state = 'FINISHED'
        self.result = self._parse_job_output()

    def __str__(self):
        return "RemoteJob(%s)" % self.id

class QSubJob(AbstractJob):

    def __init__(self,id,cmd=None,num=None,mem=None,dir=None,np=None,queue=None):
        if mem is None:
            mem = "2GB"
        super(QSubJob, self).__init__(id,cmd=cmd,num=num,mem=mem,np=np,dir=dir)
        self.queue = queue

    def __str__(self):
        return "QSubJob(%s)" % self.id

    def run(self):
        with hide('running'):
            full_cmd = ". /home/oge/default/common/settings.sh ;"
            if self.dir:
                full_cmd += "cd '%s' ;" % (self.dir)
            env_cmd = "export PYTHONPATH=$HOME/libs"
            full_cmd += "/bin/echo '( %s ; %s ) && echo OK || echo ERROR'" % (env_cmd,self.cmd)
            full_cmd += "| qsub -cwd -o %s -e %s -N '%s' -l h_vmem=%s" % (QSUB_LOG_DIR,QSUB_LOG_DIR,self.id,self.mem)
            if self.queue is not None:
                full_cmd += " -q %s" % self.queue
            if self.np:
                full_cmd += ' -pe mpi %s' % self.np
            output = run(full_cmd, pty=True, shell=False)
            m = re.match("^Your job (\d+) .* has been submitted", output)
            if m:
                self.num = m.group(1)

    def stdout_fn(self):
        if self.num is None:
            return None
        return os.path.join(QSUB_LOG_DIR,self.id+".o"+self.num)

    def get_execution_result(self):
        with hide("running","stdout","stderr"):
            stdout = run("/usr/bin/tail -1 '%s'" % (self.stdout_fn()), pty=True, shell=False)
            stdout = stdout.strip()
            if stdout in ['OK','ERROR']:
                return stdout
            else:
                return 'ERROR'
        return "?"

class AbstractJobsCollection(object):

    @staticmethod
    def create(*args,**kwargs):
        if env.host_string=='localhost':
            return LocalJobsCollection(*args,**kwargs)
        elif re.match('^.*@peyote2',env.host_string):
            kwargs['qsub_limit']=env.get('QSUB_LIMIT')
            return QSubJobsCollection(*args,**kwargs)
        elif re.match('^.*@',env.host_string):
            return RemoteJobsCollection(*args,**kwargs)
        else:
            raise Exception('unknown host_string: %s' % env.host_string)

    def __init__(self,jobs,title="Jobs"):
        self.jobs = jobs
        self.title = title

    def running(self):
        return [j for j in self.jobs if j.state=='r']

    def waiting(self):
        return [j for j in self.jobs if j.state not in ['r','FINISHED']]

    def finished(self,result=None):
        res = [j for j in self.jobs if j.state=='FINISHED']
        if result is not None:
            res = [j for j in res if j.result==result]
        return res
        
    def unfinished(self):
        return [j for j in self.jobs if j.state!='FINISHED']
        
    def run(self):
        
        t0 = datetime.now()
        while len(self.unfinished()) > 0:
            self.update_state()
            running_jobs = self.running()
            waiting_jobs = self.waiting()
            finished_jobs = self.finished()
            finished_jobs_error = self.finished("ERROR")
                
            diff = datetime.now()-t0
            msg = ""
            msg += "[%s,start=%s,elapsed=%s]" % (
                self.title,
                format_time(t0),
                format_timedelta(diff)
            )
            msg += " "
            msg += "running: %d, waiting: %d, finished: %d" % (
                len(running_jobs), 
                len(waiting_jobs),
                len(finished_jobs)
            )
            if len(finished_jobs_error)>0:
                msg += ", errors: %d" % len(finished_jobs_error)
                msg += " ("
                msg += "%s" % [str(f) for f in finished_jobs_error[0:5]]
                if len(finished_jobs_error)>5:
                    msg += ",.."
                msg += ")"
            print msg
            sys.stdout.flush()
            if len(running_jobs)==0 and len(waiting_jobs)==0:
                break
            self.wait()

        t1 = datetime.now()
        diff = t1-t0
        finished_jobs = self.finished()
        finished_jobs_ok = self.finished("OK")
        finished_jobs_error = self.finished("ERROR")
        msg = ""
        if len(finished_jobs)==len(finished_jobs_ok):
            msg = "DONE"
        else:
            msg = "DONE_WITH_ERRORS!"
        msg += " [%s,start=%s,finish=%s,elapsed=%s]" % (
                    self.title,
                    format_time(t0),
                    format_time(t1),
                    format_timedelta(diff)
              )
        print msg
        if len(finished_jobs_error) > 0:
            print " - failed jobs[%d]: %s" % (len(finished_jobs_error),[j.id for j in finished_jobs_error])
        sys.stdout.flush()
        if len(finished_jobs)==len(finished_jobs_ok):
            return True
        else:
            return False
    
    def update_state(self):
        pass

    def wait(self):
        pass

class LocalJobsCollection(AbstractJobsCollection):

    def __init__(self,jobs,title="LocalJobs",threads=1):
        super(LocalJobsCollection, self).__init__(jobs,title)
        self.threads = threads

    def update_state(self):
        r = self.waiting()
        if len(r)>0:
            r0 = r[0]
            print "** running job: %s" % r0
            r0.run()

class RemoteJobsCollection(LocalJobsCollection):

    pass

class QSubJobsCollection(AbstractJobsCollection):

    def __init__(self,jobs,title="Jobs",qsub_limit=None):
        super(QSubJobsCollection, self).__init__(jobs,title)
        self.qsub_limit = qsub_limit

    def wait(self):
        time.sleep(10)

    def waiting(self):
        return [j for j in self.jobs if j.state not in ['r','FINISHED'] and j.num is not None]

    def not_in_queue(self):
        return [j for j in self.jobs if j.num is None]

    def finished(self,result=None):
        res = [j for j in self.jobs if j.num is not None and j.state=='FINISHED']
        if result is not None:
            res = [j for j in res if j.num is not None and j.result==result]
        return res
        
    def unfinished(self):
        return [j for j in self.jobs if j.num is None or j.state!='FINISHED']

    def update_state(self):
        not_in_queue = self.not_in_queue()
        if len(not_in_queue):
            if self.qsub_limit is not None:
                how_many = max(0,self.qsub_limit-len(self.running())-len(self.waiting()))
            else:
                how_many = len(not_in_queue)
            if how_many>0:
                print "** submitting jobs to queue"
                for j in not_in_queue[0:how_many]:
                    j.run()
        qstat_info = qstat()
        qstat_info_dict = dict((row['JB_job_number'], row) for row in qstat_info)
        jobs_in_queue = set(qstat_info_dict.keys())
        for j in self.jobs:
            if j.num in jobs_in_queue:
                j.state = qstat_info_dict[j.num]['state']
            elif j.num is not None and j.state!='FINISHED':
                j.result = j.get_execution_result()
                j.state = 'FINISHED'

_my_files_cache = {}

class my_files(object):

    @staticmethod
    def exists(fn):
        global _my_files_cache
        d = os.path.dirname(fn)
        f = os.path.basename(fn)
        if not _my_files_cache.has_key(d):
            stdout = ""
            with hide('running','stdout','stderr'):
                if env.host_string=='localhost':
                    stdout = local("ls -1 '%s'" % d,capture=True)
                else:
                    stdout = run("ls -1 '%s'" % d, pty=True, shell=False)
            _my_files_cache[d] = set(stdout.replace("\r","").strip().split("\n"))
        return f in _my_files_cache[d]

    def clean_cache(d):
        global _my_files_cache
        if _my_files_cache.has_key(d):
            del _my_files_cache[d]

###########

@task
@hosts("twalen@rnacontacts-vm")
def update_vm(insdir="/var/www.rna-contacts/clarna/web"):
    """update installation on rnacontacts-vm"""
    with cd(insdir):
        run("git pull")
        run("./m.sh migrate")
        sudo("/etc/init.d/apache2 restart")
        sudo("/etc/init.d/django-celery.sh restart")

@task
def pull_all():
    """perform git pull for all rnacontacts installations"""
    local("git pull")
    for ins in INSTALLATIONS.keys():
        if ins!="localhost":
            with settings_set_ins(ins,force=True):
                with cd(env.SCRIPTS_DIR):
                    run("git pull",shell=False)

@task
def status_all():
    """show status for all rnacontacts installations"""
    local("git status")
    for ins in INSTALLATIONS.keys():
        if ins!="localhost":
            with settings_set_ins(ins,force=True):
                with cd(env.SCRIPTS_DIR):
                    run("git status",shell=False)


@task
@hosts("twalen@peyote2")
def qstat():
    with hide('running','stdout','stderr'):
        stdout = run(". /home/oge/default/common/settings.sh ; qstat -xml", pty=True, shell=False)
        res = []
        running_jobs = []
        if "Disk quotas" in stdout:
            tmp = stdout.split("\n")[4:]
            stdout = "\n".join(tmp)
        xml = minidom.parseString(stdout)
        for n in xml.getElementsByTagName('job_list'):
            row = {}
            for nn in n.childNodes:
                if nn.nodeType == nn.ELEMENT_NODE:
                    if nn.firstChild is None:
                        continue
                    row[nn.nodeName] = nn.firstChild.data
            res.append(row)
        return res
    return None

####################################

@task
@ins("peyote2")
def prepare_pdb_files(pdbids,force=False,force_mo=False,force_fr=False,force_mc=False,force_agg=False,force_residues=False,force_close_doublets=False):
    if isinstance(pdbids,str):
        pdbids = pdbids.split(",")

    prepare_peyote2()

    jobs = []
    graphs_fns = {}
    for pdbid in pdbids:
        p = pdbid.lower()
        assert re.match("^[a-z0-9_]{4,64}$",p)
        
        graphs_fns[p] = {}
        
        pdb_fn = PDBObject.pdb_fn(p,"pdb",data_dir=env.DATA_DIR)
        map_fn = PDBObject.pdb_fn(p,"res_map",data_dir=env.DATA_DIR)
        res_dict_fn = PDBObject.pdb_fn(p,"residues",data_dir=env.DATA_DIR)
        rev_desc_fn = PDBObject.pdb_fn(p,"rev_desc",data_dir=env.DATA_DIR)
        cl_d_fn = PDBObject.pdb_fn(p,"close_doublets",data_dir=env.DATA_DIR)
        rv_d_fn = PDBObject.pdb_fn(p,"contacts_RV",data_dir=env.DATA_DIR)
        mc_d_fn = PDBObject.pdb_fn(p,"contacts_MC",data_dir=env.DATA_DIR)
        fr_d_fn = PDBObject.pdb_fn(p,"contacts_FR",data_dir=env.DATA_DIR)
        mo_d_fn = PDBObject.pdb_fn(p,"contacts_MO",data_dir=env.DATA_DIR)
        
        assert my_files.exists(pdb_fn)
        if not my_files.exists(res_dict_fn) or force_residues or force:
            cmd = """./make-residue-dict.py -i '%(pdb_fn)s' -o '%(res_dict_fn)s' --with-backbone --with-neighbours""" % locals()
            jobs.append(AbstractJob.create(cmd=cmd,id="res-dict-%s"%p,mem="2GB",dir=env.SCRIPTS_DIR))
        if not my_files.exists(rev_desc_fn) or force:
            cmd = """./make-reverse-descriptions-dict.py --input='%(pdb_fn)s' --output-json='%(rev_desc_fn)s'""" % locals()
            jobs.append(AbstractJob.create(cmd=cmd,id="rev-desc-%s"%p,mem="2GB",dir=env.SCRIPTS_DIR))
        for (prg_id,out_fn,prg_options) in [
                ("cl",cl_d_fn,"--close-doublets --max-distance=4.0"),
                ("rv",rv_d_fn,"--use-rnaview --extract-all"),
                ("mc",mc_d_fn,"--use-mc-annotate --extract-all"),
                ("mo",mo_d_fn,"--use-moderna --extract-all"),
                ("fr",fr_d_fn,"--use-fr3d --extract-all"),
            ]:
            graphs_fns[p][prg_id] = out_fn
            action_req = False
            if force or not my_files.exists(out_fn):
                action_req = True
            if prg_id=='cl' and force_close_doublets:
                action_req = True
            if prg_id=='mo' and force_mo:
                action_req = True
            if prg_id=='mc' and force_mc:
                action_req = True
            if prg_id=='fr' and force_fr:
                action_req = True
            mem = "2GB"
            if prg_id=="fr":
                mem = "2GB"
            if action_req:
                cmd = ("""./extract-contacts.py -i '%(pdb_fn)s' --output-graph='%(out_fn)s' --no-pdb --dont-normalize --skip-invalid-doublets """+prg_options)%locals()
                jobs.append(AbstractJob.create(cmd=cmd,id="extract-%s-%s"%(prg_id,p),mem=mem,dir=env.SCRIPTS_DIR))
    res1 = AbstractJobsCollection.create(jobs,title="prepare-pdbs1").run()
    # assert res1==True
    
    jobs = []
    for pdbid in pdbids:
        p = pdbid.lower()
        pdb_fn = PDBObject.pdb_fn(p,"pdb",data_dir=env.DATA_DIR)
        groups_fn = PDBObject.pdb_fn(p,"groups",data_dir=env.DATA_DIR)

        gr_opts = ""
        gr_opts += " --close-doublets-graph=%s" % graphs_fns[p]['cl']
        gr_opts += " --rnaview-graph=%s" % graphs_fns[p]['rv']
        gr_opts += " --mc-annotate-graph=%s" % graphs_fns[p]['mc']
        gr_opts += " --moderna-graph=%s" % graphs_fns[p]['mo']
        gr_opts += " --fr3d-graph=%s" % graphs_fns[p]['fr']

        action_req = False
        if force or force_agg or not my_files.exists(groups_fn):
            action_req = True
        if action_req:
            cmd = "./aggregate-doublets-dict.py --pdb-id=%(pdbid)s %(gr_opts)s --output-json=%(groups_fn)s" % locals()
            jobs.append(AbstractJob.create(cmd=cmd,id="agg-gr-%s"%(p),mem="4GB",dir=env.SCRIPTS_DIR))
    res2 = AbstractJobsCollection.create(jobs,title="prepare-pdbs2").run()
    # assert res2==True

@task
@ins("peyote2")
def prepare_group_files():
    """compute rna.groups.json.gz and rna.reduced_groups.json.gz files"""
    jobs = []
    for setname,s in [('data-training.txt','training'),('data-bench.txt','bench'),('data-all.txt','all')]:
        pdbids = open(setname).read().strip().split("\n")
        merge_pdb = ",".join(pdbids)
        inp_dir = os.path.join(env.DATA_DIR,"pdb_files")
        out_fn = FileNamesObject.groups_fn(setname=s,reduced=False,data_dir=env.DATA_DIR)
        out2_fn = FileNamesObject.groups_fn(setname=s,reduced=True,data_dir=env.DATA_DIR)
        mem = "8GB"
        if setname=='data-all.txt':
            mem = "16GB"
        jobs.append(AbstractJob.create(cmd="./aggregate-doublets-dict.py --input-dir=%(inp_dir)s --output-json=%(out_fn)s --merge=%(merge_pdb)s" % locals(),
            id="gen-merged-groups-%s"%setname,mem=mem,dir=env.SCRIPTS_DIR))
        jobs.append(AbstractJob.create(cmd="./aggregate-doublets-dict.py --input-dir=%(inp_dir)s --output-json=%(out2_fn)s --merge=%(merge_pdb)s --merge-skip-groups=\"(^by-pdb|^descriptions|.*/all$|.*/Pu-Pu$|.*/Pu-Py$|.*/Py-Py$)\"" % locals(),
            id="gen-reduced-groups-%s"%setname,mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="merge-groups").run()
    assert res==True

@task
@ins("peyote2")
def prepare_pdb_files_all(setname="data-all.txt"):
    prepare_peyote2()
    pdbids = open(setname).read().strip().split("\n")
    prepare_pdb_files(pdbids)
    if setname=='data-all.txt':
        prepare_group_files()

@task
@ins("peyote2")
def prepare_pdb_files_aggregate_all():
    prepare_peyote2()
    pdbids = open("data-all.txt").read().strip().split("\n")
    prepare_pdb_files(pdbids,force_agg=True)
    prepare_group_files()

@task
@runs_once
@ins("twalen@peyote2")
def prepare_peyote2():
    if env.INS in ["peyote2","peyote2-im","peyote2-aux"]:
        with cd(env.SCRIPTS_DIR):
            run("git pull",shell=False)

@task
@ins("peyote2")
def compute_groups(sc=None,desc=None,n_type=None):
    prepare_peyote2()
    
    jobs = []
    for (_sc,_desc,_n_type) in cl_tuples():
        if sc is not None and sc!=_sc:
            continue
        if desc is not None and desc!=_desc:
            continue
        if n_type is not None and n_type!=_n_type:
            continue
        jobs.append(_cl_groups_info(_sc,_desc,_n_type))
    res1 = AbstractJobsCollection.create(jobs,title="compute-groups").run()
    assert res1==True

@ins("peyote2")
def _run_classifier_for(pdbid="1rna",extra_options="",extra_options_combine=""):
    prepare_peyote2()

    job_id = "cl-run-%s" % pdbid
    
    libs = ",".join([os.path.join(env.DATA_DIR,"classifier.%s.json.gz")%k for k in CL_CAT.keys()])
    
    inp_fn = PDBObject.pdb_fn(pdbid,"pdb",data_dir=env.DATA_DIR)
    out_fn = PDBObject.pdb_fn(pdbid,"contacts_CL",data_dir=env.DATA_DIR)

    cl_d_fn = PDBObject.pdb_fn(pdbid,"close_doublets",data_dir=env.DATA_DIR)

    rv_gr_fn = PDBObject.pdb_fn(pdbid,"contacts_RV",data_dir=env.DATA_DIR)
    mc_gr_fn = PDBObject.pdb_fn(pdbid,"contacts_MC",data_dir=env.DATA_DIR)
    mo_gr_fn = PDBObject.pdb_fn(pdbid,"contacts_MO",data_dir=env.DATA_DIR)
    fr_gr_fn = PDBObject.pdb_fn(pdbid,"contacts_FR",data_dir=env.DATA_DIR)
    groups_fn = PDBObject.pdb_fn(pdbid,"groups",data_dir=env.DATA_DIR)
    eval_fn = PDBObject.pdb_fn(pdbid,"eval1",data_dir=env.DATA_DIR)
    stats_fn = PDBObject.pdb_fn(pdbid,"eval2",data_dir=env.DATA_DIR)
    stats_new_fn = PDBObject.pdb_fn(pdbid,"eval3",data_dir=env.DATA_DIR)
    corr_fn = PDBObject.pdb_fn(pdbid,"corr",data_dir=env.DATA_DIR)
    
    cmd = ""
    if extra_options=="fr_test":
        cmd += "./extract-contacts.py -i '%(inp_fn)s' --output-graph='%(out_fn)s' --no-pdb --dont-normalize --skip-invalid-doublets --use-fr3d --extract-all" % locals()
    elif extra_options=="only_combine":
        cmd += "echo only_combine"
    else:
        cmd += "./clarna.py -i %(inp_fn)s --lib=%(libs)s --normalize-graph --save-graph=%(out_fn)s --compare-with=%(groups_fn)s --ignore-bad-doublets %(extra_options)s" % locals()

    data_dir = env.DATA_DIR
    graphs_options = "--classifier-graph=%(out_fn)s --rnaview-graph=%(rv_gr_fn)s --mc-annotate-graph=%(mc_gr_fn)s --fr3d-graph=%(fr_gr_fn)s --moderna-graph=%(mo_gr_fn)s" % locals()
    cmd += " && ./combine-contact-graphs.py --pdb-id=%(pdbid)s --classifier-graph=%(out_fn)s --groups=%(groups_fn)s -o %(eval_fn)s %(extra_options_combine)s" % locals()
    cmd += " && ./combine-contact-graphs.py --eval-mode --pdb-id=%(pdbid)s %(graphs_options)s -o %(stats_fn)s %(extra_options_combine)s" % locals()
    cmd += " && ./combine-contact-graphs.py --new-eval-mode --pdb-id=%(pdbid)s %(graphs_options)s --groups=%(groups_fn)s --data-dir=%(data_dir)s -o %(stats_new_fn)s %(extra_options_combine)s" % locals()
    cmd += " && ./combine-contact-graphs.py --correlation --pdb-id=%(pdbid)s %(graphs_options)s -o %(corr_fn)s %(extra_options_combine)s" % locals()

    assert my_files.exists(os.path.join(env.DATA_DIR,"classifier.bp.json.gz"))
    assert my_files.exists(inp_fn)
    
    return AbstractJob.create(cmd=cmd,id=job_id,mem="4GB",dir=env.SCRIPTS_DIR)

@task
@ins("peyote2")
def run_classifier_for(pdbid="1rna",extra_options="",extra_options_combine=""):
    jobs = [_run_classifier_for(pdbid,extra_options,extra_options_combine)]
    res = AbstractJobsCollection.create(jobs,title="run-cl-%s"%pdbid).run()
    assert res==True

@task
@ins("peyote2")
def import_classifier_eval(setname=None):
    jobs = []
    for s in ["small","medium","training","bench"]:
        cur_setname = "data-%s.txt" % s
        if setname is not None and setname!=cur_setname:
            continue

        merge_eval_fn = FileNamesObject.eval_fn(setname=s,t="eval1",data_dir=env.DATA_DIR)
        merge_stats_fn = FileNamesObject.eval_fn(setname=s,t="eval2",data_dir=env.DATA_DIR)
        merge_stats_new_fn = FileNamesObject.eval_fn(setname=s,t="eval3",data_dir=env.DATA_DIR)
        if my_files.exists(merge_eval_fn):
            jobs.append(AbstractJob.create(cmd="./m.sh computed_data_put --input-json=%(merge_eval_fn)s --key=new-classifier-evaluation-%(s)s" % locals(), 
                             id="import-cl-result1-%s"%s,dir=env.DJANGO_DIR,mem="8GB"))
        if my_files.exists(merge_stats_fn):
            jobs.append(AbstractJob.create(cmd="./m.sh computed_data_put --input-json=%(merge_stats_fn)s --key=new-classifier-full-stats-%(s)s" % locals(), 
                             id="import-cl-result2-%s"%s,dir=env.DJANGO_DIR,mem="16GB"))
        if my_files.exists(merge_stats_new_fn):
            jobs.append(AbstractJob.create(cmd="./m.sh computed_data_put --input-json=%(merge_stats_new_fn)s --key=new-classifier-full-stats2-%(s)s" % locals(), 
                             id="import-cl-result3-%s"%s,dir=env.DJANGO_DIR,mem="16GB"))
    res = AbstractJobsCollection.create(jobs,title="import-to-django").run()
    assert res==True

@task
@ins("peyote2")
def prepare_classifier_eval(setname=None):
    jobs = []
    for s in ['small','medium','training','bench','all']:
        cur_setname = "data-%s.txt" % s
        if setname is not None and setname!=cur_setname:
            continue
        pdbids = open(cur_setname).read().strip().split("\n")

        for t in ['eval1','eval2','eval3','corr']:
            input_fn = FileNamesObject.eval_fn(setname=s,t=t+"_in",data_dir=env.DATA_DIR)
            merge_fn = FileNamesObject.eval_fn(setname=s,t=t,data_dir=env.DATA_DIR)

            list_fns = "\n".join([PDBObject.pdb_fn(p,t,data_dir=env.DATA_DIR) for p in pdbids])
            remote_write_file(input_fn, list_fns)

            jobs.append(AbstractJob.create(cmd="./merge-dicts.py --input-from-file=%(input_fn)s -o %(merge_fn)s" % locals(), id="merge-%s-%s"%(t,s),mem="4GB",dir=env.SCRIPTS_DIR))

    res1 = AbstractJobsCollection.create(jobs,title="merging").run()
    assert res1==True

    if env.INS in ['peyote2','peyote2-im','rnacontacts-vm','ouzo']:
        import_classifier_eval(setname)
    else:
        print "skipping import to django"

@task
@ins("peyote2")
def run_classifier(setname="data-small.txt",extra_options="",extra_options_combine=""):
    if env.INS in ['peyote2','peyote2-im']:
        prepare_peyote2()

    pdbids = open(setname).read().strip().split("\n")
    jobs = [_run_classifier_for(pdbid,extra_options,extra_options_combine=extra_options_combine) for pdbid in pdbids]
    res = AbstractJobsCollection.create(jobs,title="run-classifier").run()
    assert res==True
    
    prepare_classifier_eval(setname)

@task
@ins("peyote2")
def tmp_compute_classifier_new1(gr=""):
    prepare_peyote2()
    
    cross_validation_num = None
    if gr!="":
        cross_validation_num = gr.replace("cross","").replace("v","").replace("t","")
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,cross_validation_num=cross_validation_num,data_dir=env.DATA_DIR)
    
    jobs = []
    for sc in ['bp']:
        for n_type in N_TYPES:
            mem = "4GB"
            if n_type=='CC':
                mem="8GB"
            jobs.append(AbstractJob.create(cmd="./compute-classifier.py --groups='%(groups_fn)s' --sub-category='%(sc)s' --n-type='%(n_type)s'" % locals(), id="comp-cl-%s-%s"%(sc,n_type),mem=mem,dir=env.SCRIPTS_DIR))
            
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

@task
@ins("peyote2")
def tmp_compute_classifier_new2(only_sc="",only_desc="",gr=""):
    sc_list = ['bp','stacking','base-ribose','base-phosphate','other','other2','other3']
    if only_sc != "":
        sc_list = [x.strip() for x in only_sc.split(',')]
    desc_list = None
    if only_desc != "":
        desc_list = set([x.strip() for x in only_desc.split(',')])

    cross_validation_num = None
    if gr!="":
        cross_validation_num = gr.replace("cross","").replace("v","").replace("t","")
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,cross_validation_num=cross_validation_num,data_dir=env.DATA_DIR)

    prepare_peyote2()
    
    jobs = []
    for _sc in sc_list:
        for n_type in N_TYPES:
            for desc in CL_CAT[_sc]:
                if desc_list is not None and desc not in desc_list:
                    continue
                mem = "4GB"
                if _sc in ['bp','stacking'] and n_type in ['CG','GC']:
                    mem = "24GB"
                if _sc in ['base-ribose','base-phosphate']:
                    mem = "4GB"
                safe_desc = desc.replace(">","g").replace("<","l")
                jobs.append(AbstractJob.create(cmd="""./compute-classifier.py --groups='%(groups_fn)s' --sub-category='%(_sc)s' --n-type='%(n_type)s' --desc="%(desc)s" """ % locals(), id="comp-cl-%s-%s-%s"%(_sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

@task
@ins("peyote2")
def tmp_compute_classifier_new3(only_sc="",gr=""):
    prepare_peyote2()
    
    cross_validation_num = None
    if gr!="":
        cross_validation_num = gr.replace("cross","").replace("v","").replace("t","")
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,cross_validation_num=cross_validation_num,data_dir=env.DATA_DIR)
    
    sc_list = ['bp','stacking','base-ribose','base-phosphate','other','other2','other3']
    if only_sc != "":
        sc_list = [x.strip() for x in only_sc.split(',')]
    
    jobs = []
    for sc in sc_list:
        jobs.append(AbstractJob.create(cmd="./compute-classifier.py --groups='%(groups_fn)s' --sub-category='%(sc)s' --combine" % locals(), id="combine-cl-%s"%(sc),mem="4GB",dir=env.SCRIPTS_DIR))
            
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

@task
@ins("peyote2")
def gen_classifier_pdb():
    prepare_peyote2()

    jobs = []
    for sc in ['bp','stacking','base-ribose','base-phosphate']:
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace("<","l").replace(">","g")
                mem = "4GB"
                input_json = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
                filter_keys = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
                output = os.path.join(env.DATA_DIR,"classifier",sc,n_type,"ref-%s.pdb.gz" % safe_desc)
                if sqsubc=='base-ribose':
                    type = 'base-ribose'
                elif sc=='base-phosphate':
                    type = 'base-phosphate'
                else:
                    type = 'bp'
                jobs.append(AbstractJob.create(cmd="""./gen-pdb.py --type='%(type)s' --limit=500 --input-json='%(input_json)s' --filter-keys="%(filter_keys)s" --save-json -o '%(output)s'""" % locals(), id="gen-pdb-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="gen_classifier_pdb").run()

@task
@ins("peyote2")
def gen_bench_pdb():
    prepare_peyote2()

    jobs = []
    for sc in ['bp','stacking','base-ribose','base-phosphate']:
        if sc=='bp':
            continue
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace("<","l").replace(">","g")
                mem = "4GB"
                input_json = FileNamesObject.eval_fn(setname="training",t="eval1",data_dir=env.DATA_DIR)
                filter_keys = "evaluation/ref-all/%(desc)s/%(n_type)s" % locals()
                output = os.path.join(env.DATA_DIR,"bench_pdb",sc,n_type,"ref-%s.pdb.gz" % safe_desc)
                if sc=='base-ribose':
                    type = 'base-ribose'
                elif sc=='base-phosphate':
                    type = 'base-phosphate'
                else:
                    type = 'bp'
                jobs.append(AbstractJob.create(cmd="""./gen-pdb.py --type='%(type)s' --limit=500 --input-json='%(input_json)s' --filter-keys="%(filter_keys)s" --save-json -o '%(output)s'""" % locals(), id="gen-pdb-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="gen_classifier_pdb").run()

@task
@ins("peyote2")
def tw_20130108_verify_reference_doublets():
    groups_fn = FileNamesObject.groups_fn(setname="all",reduced=True,data_dir=env.DATA_DIR)
    mem="8GB"
    jobs = []
    for n_type in N_TYPES:
        cmd = """./verify-reference-doublets.py --n-type=%(n_type)s --update-expert --groups=%(groups_fn)s""" % locals()
        jobs.append(AbstractJob.create(cmd=cmd, id="verify-%s"%n_type,mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="verify").run()

@task
@ins("peyote2")
def jb_rpt_20130122_gen_pdb_for_bad_bp_classes():
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    RES_DIR = "/home/twalen/jb_rpt_20130122"
    run("mkdir -p '%(RES_DIR)s'" % locals())
    jobs = []
    mem = "8GB"
    for n_type,desc in [('UG','WW_tran'),
        ('UC','HW_cis'),('UG','HW_cis'),
        ('AG','HW_tran'),('UA','HW_tran'),('UC','HW_tran'),
        ('AC','WH_cis'),
        ('GC','SW_cis'),('UU','SW_cis'),
        ('UA','SW_tran'),
        ('AC','WS_tran'),
        ('CG','HH_cis'),
        ('GC','SH_tran'),
    ]:
        for fn,group in (('detected','ref-ok'),('fuzzy','ref-ok-fuzzy'),('undetected','ref-undetected'),('training','ref-all')):
            if group=='ref-all':
                json = groups_fn
                key = "classifier/bp/%(desc)s/%(n_type)s" % locals()
            else:
                json = FileNamesObject.eval_fn(setname="bench",t="eval1",data_dir=env.DATA_DIR)
                key = "evaluation/%(group)s/%(desc)s/%(n_type)s" % locals()
            pdb = os.path.join(RES_DIR,"%(desc)s_%(n_type)s_%(fn)s.pdb.gz" % locals())
            cmd = """./gen-pdb.py --limit 500 --type=bp --input-json=%(json)s --filter-keys="%(key)s" -o "%(pdb)s" """ % locals()
            jobs.append(AbstractJob.create(cmd=cmd, id="gen-pdb-%s-%s-%s"%(desc,n_type,fn),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="create-pdbs").run()

@task
@ins("peyote2")
def tw_20130123_prepare_zz_files():
    prepare_peyote2()
    prepare_pdb_files(pdbids='zz01,zz02,zz03,zz04,zz05,zz06,zz07,zz08,zz09,zz10,zz11,zz12',force=True)
    prepare_group_files()


@task
@ins("peyote2")
def tw_rpt_20130128_gen_pdb_for_bad_bp_classes():
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    if env.INS in ['localhost']:
        RES_DIR = "/tmp/tw_rpt_20130128"
        local("mkdir -p '%(RES_DIR)s'"%locals())
    else:
        RES_DIR = "/home/twalen/tw_rpt_20130128"
        run("mkdir -p '%(RES_DIR)s'" % locals())
    jobs = []
    mem = "8GB"
    for sc,n_type,desc in [
        ('bp','AC','HW_cis'), 
        ('bp','CA','SW_tran'), 
        ('bp','AG','HS_cis'), 
        ('bp','GA','SH_cis'),
        ('bp','CU','HW_tran'), 
        ('bp','CU','WH_tran'),
        ('bp','CU','SW_cis'),
        ('bp','UC','HW_tran'),
        ('bp','UU','WS_cis'),
        ('bp','UU','HS_cis'),
    ]:
        for fn,group in (('detected','ref-ok'),('fuzzy','ref-ok-fuzzy'),('undetected','ref-undetected'),('training','ref-all'),('wrong-desc','ref-wrong-desc')):
            if group=='ref-all':
                json = {"":groups_fn}
                key = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
            else:
                json = {"_training": FileNamesObject.eval_fn(setname="training",t="eval1",data_dir=env.DATA_DIR),
                        "_bench": FileNamesObject.eval_fn(setname="bench",t="eval1",data_dir=env.DATA_DIR)}
                key = "evaluation/%(group)s/%(desc)s/%(n_type)s" % locals()
            for json_suffix,json_fn in json.items():
                pdb = os.path.join(RES_DIR,"%(desc)s_%(n_type)s_%(fn)s%(json_suffix)s.pdb.gz" % locals())
                cmd = """./gen-pdb.py --limit 500 --type=bp --input-json=%(json_fn)s --filter-keys="%(key)s" -o "%(pdb)s" --save-json """ % locals()
                jobs.append(AbstractJob.create(cmd=cmd, id="gen-pdb-%s-%s-%s%s"%(desc,n_type,fn,json_suffix),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="create-pdbs").run()

@task
@ins("peyote2")
def tw_20130129_check_residues_conformations():
    prepare_peyote2()
    jobs = []
    pdbids = open("data-all.txt").read().strip().split("\n")
    for pdbid in pdbids:
        mem = "2GB"
        jobs.append(AbstractJob.create(cmd="./_tmp_20130129_check_res_conf.py --pdb-id=%(pdbid)s" % locals(),
            id="check-res-conf-%s"%pdbid,mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="check-residues-conf").run()
    assert res==True

@task
@ins("peyote2")
def tw_rpt_20130130_gen_ref_pdb():
    if env.INS in ['localhost']:
        RES_DIR = "/tmp/tw_rpt_20130130"
        local("mkdir -p '%(RES_DIR)s'"%locals())
    else:
        RES_DIR = "/home/twalen/tw_rpt_20130130"
        run("mkdir -p '%(RES_DIR)s'" % locals())
    jobs = []
    mem = "4GB"
        
    for sc in ['bp','stacking']:
        for n_type in ['CG','GC','AU','UA','CU','UC','CC']:
            if sc in ['bp']:
                for orient in ['trans','cis']:
                    for edge in ['W','S','H']:
                        remote_write_file(os.path.join(RES_DIR,"%(sc)s_%(n_type)s_%(edge)s_%(orient)s.pml"%locals()), """
        cmd.load("ref_bp_H%(edge)s_%(orient)s_%(n_type)s.pdb.gz")
        cmd.load("ref_bp_S%(edge)s_%(orient)s_%(n_type)s.pdb.gz")
        cmd.load("ref_bp_W%(edge)s_%(orient)s_%(n_type)s.pdb.gz")
        cmd.color(4,"ref_bp_H%(edge)s_%(orient)s_%(n_type)s")
        cmd.color(2,"ref_bp_S%(edge)s_%(orient)s_%(n_type)s")
        cmd.color(3,"ref_bp_W%(edge)s_%(orient)s_%(n_type)s")
        set all_states,1
        orient animate=0""" % locals())
            continue
        
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace(">","g").replace("<","l")
                groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
                key = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
                pdb = os.path.join(RES_DIR,"ref_%(sc)s_%(safe_desc)s_%(n_type)s.pdb.gz" % locals())
                cmd = """./gen-pdb.py --limit 1000 --type=bp --input-json=%(groups_fn)s --filter-keys="%(key)s" -o "%(pdb)s" --save-json """ % locals()
                jobs.append(AbstractJob.create(cmd=cmd, id="gen-ref-pdb-%s-%s-%s"%(sc,safe_desc,n_type),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="create-pdbs").run()
    assert res==True
    
@task
@ins("peyote2")
def tw_20130131_compute_cl_for_selected_classes():
    prepare_peyote2()
    
    jobs = []
    for sc in ['bp']:
        for n_type in ['AG','GA']:
            for desc in ['HS_cis','SH_cis']:
                mem = "4GB"
                if sc in ['bp','stacking'] and n_type in ['CG','GC']:
                    mem = "12GB"
                safe_desc = desc.replace(">","g").replace("<","l")
                jobs.append(AbstractJob.create(cmd="""./compute-classifier.py --sub-category='%(sc)s' --n-type='%(n_type)s' --desc="%(desc)s" """ % locals(), id="comp-cl-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

    jobs = []
    for sc in ['bp']:
        jobs.append(AbstractJob.create(cmd="./compute-classifier.py --sub-category='%(sc)s' --combine" % locals(), id="combine-cl-%s"%(sc),mem="4GB",dir=env.SCRIPTS_DIR))
            
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True
    
@task
@ins("peyote2")
def tw_rpt_20130201_gen_pdb_for_bad_bp_classes():
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    if env.INS in ['localhost']:
        RES_DIR = "/tmp/tw_rpt_20130201"
        local("mkdir -p '%(RES_DIR)s'"%locals())
    else:
        RES_DIR = "/home/twalen/tw_rpt_20130201"
        run("mkdir -p '%(RES_DIR)s'" % locals())
    jobs = []
    mem = "8GB"
    for sc,n_type,desc in [
        ('bp','CA','WH_cis'), 
    ]:
        for fn,group in (('detected','ref-ok'),('fuzzy','ref-ok-fuzzy'),('undetected','ref-undetected'),('training','ref-all'),('wrong-desc','ref-wrong-desc')):
            if group=='ref-all':
                json = {"":groups_fn}
                key = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
            else:
                json = {"_training": FileNamesObject.eval_fn(setname="training",t="eval1",data_dir=env.DATA_DIR),
                        "_bench": FileNamesObject.eval_fn(setname=s,t="bench",data_dir=env.DATA_DIR)}
                key = "evaluation/%(group)s/%(desc)s/%(n_type)s" % locals()
            for json_suffix,json_fn in json.items():
                pdb = os.path.join(RES_DIR,"%(desc)s_%(n_type)s_%(fn)s%(json_suffix)s.pdb.gz" % locals())
                cmd = """./gen-pdb.py --limit 500 --type=bp --input-json=%(json_fn)s --filter-keys="%(key)s" -o "%(pdb)s" --save-json """ % locals()
                jobs.append(AbstractJob.create(cmd=cmd, id="gen-pdb-%s-%s-%s%s"%(desc,n_type,fn,json_suffix),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="create-pdbs").run()

@task
@ins("peyote2")
def tw_20130204_compute_cl_for_selected_classes():
    prepare_peyote2()
    
    jobs = []
    for sc in ['bp']:
        for n_type in ['AC','CA']:
            for desc in ['HW_cis','WH_cis']:
                mem = "4GB"
                if sc in ['bp','stacking'] and n_type in ['CG','GC']:
                    mem = "12GB"
                safe_desc = desc.replace(">","g").replace("<","l")
                jobs.append(AbstractJob.create(cmd="""./compute-classifier.py --sub-category='%(sc)s' --n-type='%(n_type)s' --desc="%(desc)s" """ % locals(), id="comp-cl-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

    jobs = []
    for sc in ['bp']:
        jobs.append(AbstractJob.create(cmd="./compute-classifier.py --sub-category='%(sc)s' --combine" % locals(), id="combine-cl-%s"%(sc),mem="4GB",dir=env.SCRIPTS_DIR))
            
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

@task
@ins("peyote2")
def tw_rpt_20130204_gen_pdb_for_bad_bp_classes():
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    if env.INS in ['localhost']:
        RES_DIR = "/tmp/tw_rpt_20130204"
        local("mkdir -p '%(RES_DIR)s'"%locals())
    else:
        RES_DIR = "/home/twalen/tw_rpt_20130204"
        run("mkdir -p '%(RES_DIR)s'" % locals())
    jobs = []
    mem = "8GB"
    for sc in ['bp']:
        for n_type in ['AG']:
            for desc in CL_CAT[sc]:
                for fn,group in (('detected','ref-ok'),('fuzzy','ref-ok-fuzzy'),('undetected','ref-undetected'),('training','ref-all'),('wrong-desc','ref-wrong-desc'),('prev-undet','prev-undetected')):
                    if group=='ref-all':
                        json = {"":groups_fn}
                        key = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
                    else:
                        json = {"_training": FileNamesObject.eval_fn(setname="training",t="eval1",data_dir=env.DATA_DIR),
                                "_bench": FileNamesObject.eval_fn(setname="bench",t="eval1",data_dir=env.DATA_DIR)}
                        key = "evaluation/%(group)s/%(desc)s/%(n_type)s" % locals()
                    for json_suffix,json_fn in json.items():
                        pdb = os.path.join(RES_DIR,"%(desc)s_%(n_type)s_%(fn)s%(json_suffix)s.pdb.gz" % locals())
                        cmd = """./gen-pdb.py --limit 500 --type=bp --input-json=%(json_fn)s --filter-keys="%(key)s" -o "%(pdb)s" --save-json """ % locals()
                        jobs.append(AbstractJob.create(cmd=cmd, id="gen-pdb-%s-%s-%s%s"%(desc,n_type,fn,json_suffix),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="create-pdbs").run()

@task
@ins("peyote2")
def tw_20130207_test_stacking_classifier():
    prepare_peyote2()
    
    jobs = []
    for sc in ['stacking']:
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                mem = "4GB"
                if sc in ['bp','stacking'] and n_type in ['CG','GC']:
                    mem = "12GB"
                safe_desc = desc.replace(">","g").replace("<","l")
                jobs.append(AbstractJob.create(cmd="""./compute-classifier.py --test-classifier --sub-category='%(sc)s' --n-type='%(n_type)s' --desc="%(desc)s" """ % locals(), id="test-cl-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="test-classifier").run()
    assert res==True

@task
@ins("localhost")
def nar_images_venn():
    import urllib2
    SERVER = "http://iimcb.genesilico.pl"
    PREFIX = "/clarna"
    RES_DIR = "/tmp/nar_images_venn"
    if not os.path.isdir(RES_DIR):
        os.mkdir(RES_DIR)
    stats_fn = os.path.join(RES_DIR,"stats.html")
    if not os.path.isfile(stats_fn):
        print "downloading stats"
        html = urllib2.urlopen(SERVER+PREFIX+"/alg/new-classifier-full-stats/html/bench?pr=1").read()
        write_file(stats_fn, html)
    else:
        html = read_file(stats_fn)
    img_fn = {}
    for line in html.split("\n"):
        m = re.match('<img id="venn-([^"]*)" src="([^"]*)"',line)
        if m:
            id = m.group(1)
            img_url = m.group(2).replace("/png/","/pdf/")
            print "id=%s img_url=%s" % (id,img_url)
            img_fn[id] = os.path.join(RES_DIR,id+".pdf")
            if not os.path.isfile(img_fn[id]):
                print "downloading %s (url=%s)" % (img_fn[id], img_url)
                img_data = urllib2.urlopen(SERVER+img_url).read()
                write_file(img_fn[id], img_data)
    
    o_fns = []
    for i,row in enumerate([("bp-classic","bp-non-classic","stacking",), ("base-phosphate","base-ribose")]):
        if i==0:
            point_size = 32
        else:
            point_size = 48
        files = " ".join(["-pointsize %d"%point_size+" -label '%s'"%x + " "+ img_fn[x] for x in row])
        o_fn = os.path.join(RES_DIR,"tmp%d.pdf" % i)
        o_fns.append(o_fn)
        if i==0:
            opts = " -density 300 "
        else:
            opts = " -density 200 "
        cmd = "montage %(opts)s %(files)s -mode concatenate -tile x1 %(o_fn)s" % locals()
        print "running %s" % cmd
        os.system(cmd)
    
    files = " ".join(o_fns)
    file1 = o_fns[0]
    file2 = o_fns[1]
    output_fn = os.path.join("../clarna/doc/","nar-figure2-classifier-benchmarks-rgb.pdf")
    cmd = "convert -density 300 %(file1)s -density 200 %(file2)s -gravity South -bordercolor '#FFFFFF' -border 0x50 -append %(output_fn)s" % locals()
    print "running %s" % cmd
    os.system(cmd)
    # info about combine: http://gotofritz.net/blog/geekery/combining-images-imagemagick/

    output_cmyk_fn = os.path.join("../clarna/doc/","nar-figure2-classifier-benchmarks-cmyk.pdf")
    cmd = _cmd_pdf_rgb_to_cmyk(output_fn,output_cmyk_fn)
    print "running %s" % cmd
    os.system(cmd)

    
@task
@ins("peyote2")
def tw_20130215_gen_histograms_for_bp_params():
    prepare_peyote2()

    if env.INS=='localhost':
        output_dir = "/Users/tomek/Desktop/bp_params_histograms.local"
    else:
        output_dir = "/home/twalen/bp_params_histograms"
    jobs = []
    for sc in ['bp']:
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                mem = "4GB"
                jobs.append(AbstractJob.create(cmd="""./_tmp_20130215_check_bp_params.py --n-type='%(n_type)s' --desc="%(desc)s" --output-dir="%(output_dir)s" """ % locals(), id="gen-histograms-%s-%s"%(n_type,desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="gen-histograms").run()
    assert res==True

@task
@ins("peyote2")
def tw_20130218_gen_bp_class_limits():
    prepare_peyote2()

    if env.INS=='localhost':
        output_dir = "/Users/tomek/Desktop/bp_class_limits.local"
    else:
        output_dir = "/home/twalen/bp_class_limits"
    jobs = []
    for sc in ['bp']:
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                mem = "4GB"
                if n_type in ['CG','GC','GG']:
                    mem = "8GB"
                jobs.append(AbstractJob.create(cmd="""./compute-bp-limits.py --n-type='%(n_type)s' --desc="%(desc)s" --output-dir="%(output_dir)s" """ % locals(), id="gen-limits-%s-%s"%(n_type,desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="gen-limits").run()
    assert res==True

    jobs = []
    mem = "1GB"
    jobs.append(AbstractJob.create(cmd=""" ( echo "limits = {}" ; cat %(output_dir)s/*.txt ) | sort > %(output_dir)s/bp_limits.py """ % locals(), id="merge",mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="merge").run()

@task
@ins("peyote2")
def tw_20130219_compute_cl_for_selected_classes():
    prepare_peyote2()
    
    jobs = []
    for sc in ['bp']:
        for n_type in ['AU','UA']:
            for desc in ['WW_cis']:
                mem = "4GB"
                if sc in ['bp','stacking'] and n_type in ['CG','GC']:
                    mem = "12GB"
                safe_desc = desc.replace(">","g").replace("<","l")
                jobs.append(AbstractJob.create(cmd="""./compute-classifier.py --sub-category='%(sc)s' --n-type='%(n_type)s' --desc="%(desc)s" """ % locals(), id="comp-cl-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

    jobs = []
    for sc in ['bp']:
        jobs.append(AbstractJob.create(cmd="./compute-classifier.py --sub-category='%(sc)s' --combine" % locals(), id="combine-cl-%s"%(sc),mem="4GB",dir=env.SCRIPTS_DIR))
            
    res = AbstractJobsCollection.create(jobs,title="compute-classifier").run()
    assert res==True

@task
@ins("peyote2")
def gen_supp_data_pdbs():
    prepare_peyote2()

    mem = "4GB"
    jobs = []
    input_json = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    input2_json = FileNamesObject.eval_fn(setname="bench",t="eval3",data_dir=env.DATA_DIR)

    count_json = os.path.join(env.DATA_DIR,"supp_data","ref","ref-counts.json")
    jobs.append(AbstractJob.create(cmd="""./json-tool.py -i '%(input_json)s' --apply-map-op=len --only-keys='classifier/' --rename-keys='classifier/,ref/' -o '%(count_json)s'""" % locals(), id="gen-ref-counts",mem=mem,dir=env.SCRIPTS_DIR))
    
    EVAL_GROUPS = ['tp','fn','fp-cis-vs-trans','fp-consistent-with-single-cl','fp-new-base-ribose','fp-not-recognized-by-others-cl','fp-others','fp-wh_cis-vs-sh_cis']
    for p in EVAL_GROUPS:
        count2_json = os.path.join(env.DATA_DIR,"supp_data","eval",p,"%s-counts.json"%p)
        jobs.append(AbstractJob.create(cmd="""./json-tool.py -i '%(input2_json)s' --apply-map-op=copy --only-keys='evaluation/CL/%(p)s/' --rename-keys='evaluation/CL/%(p)s/,%(p)s/' -o '%(count2_json)s'""" % locals(), id="gen-%s-counts"%p,mem=mem,dir=env.SCRIPTS_DIR))

    res = AbstractJobsCollection.create(jobs,title="gen_supp_data_pdb_counts").run()
    assert res==True

    jobs = []
    for sc in ['bp','stacking','base-ribose','base-phosphate']:
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace("<","l").replace(">","g")
                filter_keys = "classifier/%(sc)s/%(desc)s/%(n_type)s" % locals()
                base_fn = "%(sc)s-%(safe_desc)s-%(n_type)s.pdb" % locals()
                base_full_fn = "%(sc)s-%(safe_desc)s-%(n_type)s-full.pdb" % locals()
                output = os.path.join(env.DATA_DIR,"supp_data","ref",base_fn)
                output_full = os.path.join(env.DATA_DIR,"supp_data","ref",base_full_fn)
                if sc=='base-ribose':
                    type = 'base-ribose'
                elif sc=='base-phosphate':
                    type = 'base-phosphate'
                else:
                    type = 'bp'
                cmd = "( "
                cmd += """./gen-pdb.py --type='%(type)s' --limit=50 --input-json='%(input_json)s' --filter-keys="%(filter_keys)s" --save-json --skip-empty -o '%(output)s' """ % locals()
                cmd += " ; "
                cmd += """./gen-pdb.py --type='full' --limit=50 --input-json='%(input_json)s' --filter-keys="%(filter_keys)s" --save-json --skip-empty -o '%(output_full)s' """ % locals()

                for p in EVAL_GROUPS:
                    filter2_keys = "evaluation-doublets/CL/%(p)s/%(sc)s/%(desc)s/%(n_type)s" % locals()
                    pdb1_fn = os.path.join(env.DATA_DIR,"supp_data","eval",p,base_fn)
                    pdb2_fn = os.path.join(env.DATA_DIR,"supp_data","eval",p,base_full_fn)
                    cmd += " ; "
                    cmd += """./gen-pdb.py --type='%(type)s' --limit=50 --input-json='%(input2_json)s' --filter-keys="%(filter2_keys)s" --save-json --skip-empty -o '%(pdb1_fn)s' """ % locals()
                    cmd += " ; "
                    cmd += """./gen-pdb.py --type='full' --limit=50 --input-json='%(input2_json)s' --filter-keys="%(filter2_keys)s" --save-json --skip-empty -o '%(pdb2_fn)s' """ % locals()
                cmd += " )"

                jobs.append(AbstractJob.create(cmd=cmd, id="gen-pdb-%s-%s-%s"%(sc,n_type,safe_desc),mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,title="gen_supp_data_pbds").run()
    assert res==True

@task
@ins("peyote2")
def gen_supp_data_cl_doublets(only_sc=""):
    sc_list = ['bp','stacking','base-ribose','base-phosphate']
    if only_sc != "":
        sc_list = [x.strip() for x in only_sc.split(',')]

    prepare_peyote2()

    s_dir = os.path.join(env.DATA_DIR,"supp_data","cl-doublets")
    run("mkdir -p '%(s_dir)s'" % locals())
    mem = "2GB"
    
    jobs = []
    for sc in sc_list:
        cl_comp_dir = os.path.join(env.DATA_DIR,"comp_cl",sc)
        cmd = ""
        cmd = ""
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace("<","l").replace(">","g")
                for ext in ['json.gz','pdb.gz']:
                    ext2 = ext.replace(".gz","")
                    in_fn = os.path.join(cl_comp_dir,n_type,safe_desc,"doublets-cl."+ext)
                    out_fn = os.path.join(s_dir,"%(sc)s-%(safe_desc)s-%(n_type)s.%(ext2)s" % locals())
                    cmd += "if [ -f '%(in_fn)s' ] ; then zcat '%(in_fn)s' > '%(out_fn)s' ; fi ; " % locals()
        jobs.append(AbstractJob.create(cmd=cmd, id="copy-"+sc,mem=mem,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"copy-cl-doublets").run()
    
    jobs = []
    # render png files
    for sc in sc_list:
        for n_type in sorted(N_TYPES):
            for desc in sorted(CL_CAT[sc]):
                safe_desc = desc.replace("<","l").replace(">","g")
                json_fn = os.path.join(s_dir,"%(sc)s-%(safe_desc)s-%(n_type)s.json" % locals())
                pdb_fn = os.path.join(s_dir,"%(sc)s-%(safe_desc)s-%(n_type)s.pdb" % locals())
                png_fn = os.path.join(s_dir,"%(sc)s-%(safe_desc)s-%(n_type)s.png" % locals())
                png_all_fn = os.path.join(s_dir,"%(sc)s-%(safe_desc)s-%(n_type)s-all.png" % locals())
                w = 1000; h = 1000
                pymol_script = """
set max_threads, 1;
load %(pdb_fn)s;
remove (hydro);
viewport 1000,1000;
set_view (\
     0.825723529,   -0.558349788,    0.080163434,\
     0.561096549,    0.827611864,   -0.015140276,\
    -0.057890646,    0.057481110,    0.996666729,\
     0.000000000,    0.000000000,  -63.575225830,\
     3.354323387,    2.634973526,    0.323200226,\
    50.123191833,   77.027259827,  -20.000000000 )
bg_color white;
ray %(w)d,%(h)d;
png %(png_fn)s,%(w)d;
set cache_frames=0;mpng /tmp/frames_;set all_states, on;
ray $(w)d,$(h)d;
png %(png_all_fn)s,%(w)d;
quit
""" % locals()
                fix = """./gen-pdb.py --input-json="%(json_fn)s" --type=full -o "%(pdb_fn)s" """ % locals()
                cmd = """if [ -f "%(pdb_fn)s" ] ; then %(fix)s ; pymol -c -d "%(pymol_script)s" ; fi ; """ % locals()
                jobs.append(AbstractJob.create(cmd=cmd, id="gen-png-%s-%s-%s"%(sc,safe_desc,n_type),mem=mem,dir=env.SCRIPTS_DIR))
    
    res = AbstractJobsCollection.create(jobs,"render-png-files").run()
    assert res==True

def pdb_to_png_job(pdb_fn,png_fn,width=1000,height=1000,extra_script=""):
    pymol_script=""
    if isinstance(pdb_fn,list):
        pymol_script+=";".join(["load %s"%x for x in pdb_fn])
    else:
        pymol_script+="load %s"%pdb_fn
    pymol_script+=";remove (hydro)"
    pymol_script+=";show lines;orient"
    if extra_script != "":
        pymol_script += ";"+extra_script
    pymol_script+=";set cache_frames=0;mpng /tmp/frames_;set all_states, on;"
    pymol_script+=";set ray_opaque_background, off"
    pymol_script+=";bg_color white"
    pymol_script+=";ray %(width)d,%(height)d;png %(png_fn)s,%(width)s"%locals()
    pymol_script+=";quit"
    h = hashlib.md5(pymol_script).hexdigest()[0:8]
    cmd = 'pymol -c -d "%s"'%pymol_script
    return AbstractJob.create(cmd=cmd, id="pymol-%s"%h,mem="4GB",dir=env.SCRIPTS_DIR)

def pdb_to_png(*args,**kwargs):
    res = AbstractJobsCollection.create([pdb_to_png_job(*args,**kwargs)],"pdb_to_png").run()
    assert res==True

@task
@ins("localhost")
def nar_images_novel_interactions():
    DOUBLETS = [
            ('1FFK:0804:0803','1FFK:0781:0781'),
            ('1F7Y:B2:B53','1F7Y:B3:B54'),
            ('1ASY:R44:R45','1ASY:R8:R8'),
            ('1FFK:01287:01944',None)
    ]
    VIEWS = [
        (\
         0.065687850,    0.810122848,    0.582570136,\
        -0.628331065,    0.487128288,   -0.606553078,\
        -0.775167823,   -0.326203585,    0.541024566,\
         0.000000000,    0.000000000,  -49.521183014,\
       116.192375183,  162.173751831,   94.251525879,\
      -2130.936523438, 2229.979003906,  -20.000000000 ),
        (\
        -0.368283361,   -0.838740706,    0.401100576,\
        -0.911063850,    0.239604533,   -0.335485578,\
         0.185279310,   -0.488982618,   -0.852387905,\
         0.000000000,    0.000000000,  -58.813762665,\
        -7.967187405,   36.238452911,   27.544078827,\
        46.369216919,   71.258308411,  -20.000000000 ),
         (\
        -0.188687697,   -0.398291856,    0.897642136,\
        -0.827519417,    0.556665003,    0.073048852,\
        -0.528779924,   -0.729032218,   -0.434628636,\
         0.000000000,   -0.000000000,  -39.754032135,\
       113.881446838,   45.569709778,  -19.274129868,\
        31.342380524,   48.165683746,  -20.000000000 ),
         (\
         0.372202545,   -0.616765261,   -0.693589926,\
        -0.648490012,    0.361815035,   -0.669740021,\
         0.664021909,    0.699064612,   -0.265297264,\
         0.000000000,    0.000000000,  -49.664604187,\
        59.883518219,  159.712799072,   98.839210510,\
        39.155952454,   60.173255920,  -20.000000000 )
    ]
    LABELS = [
        """
label 'd1' and chain B and name O6, 'B';
label 'd1' and chain A and name N6, 'A';
label 'o1' and chain A and name N3, 'C';
        """,
        """
label 'd2' and chain B and name O6, 'B';
label 'd2' and chain A and name O6, 'A';
label 'o2' and chain A and name N3, 'C';
label 'o2' and chain B and name N3, 'D';
        """,
        """
label 'd3' and chain A and name O6, 'A';
label 'd3' and chain B and name N6, 'B';
label 'o3' and chain A and name N3, 'C';
        """,
        """
label 'd4' and chain A and name N6, 'A';
label 'd4' and chain B and name O3', 'B';
        """,
    ]

    tmpdir = "/tmp/a"
    png_out_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure3-novel-interactions-rgb.png")
    tiff_cmyk_out_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure3-novel-interactions-cmyk.png")
    png_tmpout_fn = os.path.join(tmpdir,"out.png")
    remote_mkdir(tmpdir)
    
    jobs1 = []
    jobs2 = []
    jobs3 = []
    
    montage_cmd = "montage -font Arial -pointsize 60 "
    
    for i,(d_id,other_id) in enumerate(DOUBLETS,start=1):
        tmp1_fn = os.path.join(tmpdir,"d%d.pdb"%i)
        tmp2_fn = os.path.join(tmpdir,"o%d.pdb"%i)
        png_fn = os.path.join(tmpdir,"fig%d.png"%i)
        png_c_fn = os.path.join(tmpdir,"fig%dc.png"%i)
        
        cmd = """./gen-pdb.py --doublet-id='%(d_id)s' -o '%(tmp1_fn)s' --type=full --dont-normalize"""%locals()
        jobs1.append(AbstractJob.create(cmd=cmd,id="gen-%d-d"%i,dir=env.SCRIPTS_DIR))
        if other_id is not None:
            cmd = """./gen-pdb.py --doublet-id='%(other_id)s' -o '%(tmp2_fn)s' --type=bp --dont-normalize"""%locals()
            jobs1.append(AbstractJob.create(cmd=cmd,id="gen-%d-o"%i,dir=env.SCRIPTS_DIR))
        else:
            remote_write_file(tmp2_fn,"")

        border_x = 20
        w = 1500
        h = 1000
        label_size = 40
        if i==1:
            border_x = 100
        elif i==2:
            w = 4000
            h = 3000
            label_size = 110
        elif i==3:
            border_x = 300
            label_size = 60
        elif i==4:
            label_size = 50

        view = VIEWS[i-1]
        labels_cmd = LABELS[i-1]
        pymol_script = """
set max_threads, 1;
load %(tmp1_fn)s;
load %(tmp2_fn)s;
remove (hydro);
viewport %(w)d,%(h)d;
bg_color white;
cmd.color(3,'d%(i)d');
cmd.color('grey40','o%(i)d');

create d_obj, 'd%(i)d';
create o_obj, 'o%(i)d';
show sticks, d_obj;
show sticks, o_obj;
set stick_radius, 0.15, d_obj;
set stick_radius, 0.08, o_obj;

set label_size, %(label_size)d;
set label_color, red, d%(i)d;
set label_color, red, o%(i)d;
set label_position, (0,1.2,0);
%(labels_cmd)s

set_view %(view)s;
ray %(w)d,%(h)d;
png %(png_fn)s,%(w)d;
quit
""" % locals()
        remote_write_file(os.path.join(tmpdir,"script%d.pml"%i),pymol_script)
        cmd = """pymol -c -d "%(pymol_script)s" """ % locals()
        cmd += " ; "
        cmd += """convert "%(png_fn)s" -trim +repage -resize 1000x1000 -bordercolor white -border %(border_x)dx20 "%(png_c_fn)s" """ % locals()
        jobs2.append(AbstractJob.create(cmd=cmd,id="gen-png-%d"%i,dir=env.SCRIPTS_DIR))
        
        montage_cmd += " -label '(%(i)d)' '%(png_c_fn)s'" % locals()

    montage_cmd = """
montage -font Arial -pointsize 60 \
    \( '%(tmpdir)s/fig1c.png' -resize 1000x1000 -gravity NorthWest -annotate 0 'A' \) \
    \( '%(tmpdir)s/fig2c.png' -resize 1400x1000 -gravity center -extent 1400x700 -gravity NorthWest -annotate 0 'B' \) \
    -gravity South -mode concatenate \
    -tile 2x \
    '%(tmpdir)s/out1.png' ;
montage -font Arial -pointsize 60 \
    \( '%(tmpdir)s/fig3c.png' -resize 1200x1200 -gravity NorthWest -annotate 0 'C' \) \
    \( '%(tmpdir)s/fig4c.png' -resize 1050x1500 -gravity center -extent x900 -gravity NorthWest -annotate 0 'D' \) \
    -gravity South -mode concatenate \
    -tile 2x \
    '%(tmpdir)s/out2.png' ;
montage '%(tmpdir)s/out1.png' '%(tmpdir)s/out2.png' \
    -mode concatenate \
    -tile x2 -geometry +0+50 \
    '%(tmpdir)s/out.png' ;
convert '%(tmpdir)s/out.png' -trim +repage '%(png_out_fn)s' ;
""" % locals()
    montage_cmd += _cmd_png_rgb_to_cmyk(png_out_fn,tiff_cmyk_out_fn)
    print "MONTAGE CMD: %s" % montage_cmd
    remote_write_file(os.path.join(tmpdir,"montage.sh"),"#!/bin/bash\n"+montage_cmd)
    jobs3.append(AbstractJob.create(cmd=montage_cmd,id="montage",dir=env.SCRIPTS_DIR))

    res = AbstractJobsCollection.create(jobs1,"gen-pdb-files").run()
    assert res==True
    res = AbstractJobsCollection.create(jobs2,"gen-png-files").run()
    assert res==True
    res = AbstractJobsCollection.create(jobs3,"montage").run()
    assert res==True


@task
@ins("peyote2")
def classifier_cross_validation():
    prepare_peyote2()
    
    SET1 = (1,2,3,4,5)
    SET2 = [chr(ord('A')+i) for i in xrange(13)]

    for num in SET2:
        gr = "cross"+str(num)+"t"
        c_fn = FileNamesObject.groups_fn(setname="bench",cross_validation_num=num,reduced=True,data_dir=env.DATA_DIR)

        tmp_compute_classifier_new1(gr=gr)
        tmp_compute_classifier_new2(gr=gr)
        tmp_compute_classifier_new3(gr=gr)
        
        run_classifier("data-training.txt",extra_options_combine="--only-doublets-from=%s" % c_fn)
        
        run("cd '%s' && tar cvzf cross-%s.tgz classifier.*.json.gz eval/new-classifier-eval-training.json.gz eval/new-classifier-full-stats-training.json.gz eval/new-classifier-full-stats2-training.json.gz"%(env.DATA_DIR,num))

@task
@ins("peyote2")
def tw_20130610_classifier_cross_validation():
    prepare_peyote2()

    for num in (1,2,3,4,5):
        gr = "cross"+str(num)+"t"
        
        c_fn = FileNamesObject.groups_fn(setname="bench",reduced=True,cross_validation_num=num,data_dir=env.DATA_DIR)

        run_classifier("data-training.txt",extra_options_combine="--only-doublets-from=%s" % c_fn)
        
        run("cd '%s' && tar cvzf cross-%s.tgz classifier.*.json.gz eval/new-classifier-eval-training.json.gz eval/new-classifier-full-stats-training.json.gz eval/new-classifier-full-stats2-training.json.gz"%(env.DATA_DIR,num))


@task
@ins("localhost")
def nar_images_benchmarks():
    N_VALUES = [100,200,400,800]
    CLARNA_THREADS = 4
    PDB_ID = "1vq6"
    PROGRAMS = ['mc_annotate','rnaview','fr3d','cl','cl_old']
    PDB_FN = PDBObject.pdb_fn(PDB_ID,"pdb",data_dir=env.DATA_DIR)
    CL_LIBS = ",".join([os.path.join(env.DATA_DIR,"classifier.%s.json.gz"%x) for x in CL_CAT.keys()])
    tmpdir = "/tmp/nar_img_bench"
    remote_mkdir(tmpdir)
    jobs1 = []
    jobs2 = []
    for n in N_VALUES:
        pdb_fn = os.path.join(tmpdir,"t%d.pdb"%n)
        if not my_files.exists(pdb_fn):
            cmd = "./normalize_pdb.py -i '%(PDB_FN)s' --residues-limit=%(n)d -o '%(pdb_fn)s'" % locals()
            jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb-%d"%n,dir=env.SCRIPTS_DIR))
        else:
            print "reusing old pdb file (%d)" % (n)

        for prg in PROGRAMS:
            time_fn = os.path.join(tmpdir,"time_%d_%s.txt"%(n,prg))
            if not my_files.exists(time_fn):
                if prg=="cl_old":
                    cmd = """python -m timeit -v -n 1 -r 1 'import os; os.system("./clarna.py --threads=%(CLARNA_THREADS)d -i %(pdb_fn)s --lib=%(CL_LIBS)s")' --use-old-params""" % locals()
                elif prg=="cl":
                    cmd = """python -m timeit -v -n 1 -r 1 'import os; os.system("./clarna.py --threads=%(CLARNA_THREADS)d -i %(pdb_fn)s --lib=%(CL_LIBS)s ")' """ % locals()
                else:
                    cmd = """python -m timeit -v -n 1 -r 1 'from utils import run_%(prg)s; run_%(prg)s("%(pdb_fn)s")' """ % locals()
                cmd += "| grep raw > %(time_fn)s""" % locals()
                jobs2.append(AbstractJob.create(cmd=cmd,id="timeit-%d-%s"%(n,prg),dir=env.SCRIPTS_DIR))
            else:
                print "reusing old results (%d,%s)" % (n,prg)
    res = AbstractJobsCollection.create(jobs1,"gen-pdbs").run()
    assert res==True
    res = AbstractJobsCollection.create(jobs2,"timeit").run()
    assert res==True

    data_fn = os.path.join(tmpdir,"datafile.dat")

    jobs3 = []
    programs_str = "["+",".join(['"'+x+'"' for x in PROGRAMS])+"]"
    cmd = """python -c 'import re;from utils import read_file,write_file;
res = "";
res += "\\t".join(["n"]+%(programs_str)s)+"\\n";
for n in %(N_VALUES)s:
    res += str(n)
    for prg in %(programs_str)s:
        res += "\\t"+re.sub(r"^.*raw times:\\s*([\\d\\.]+).*$","\\\\1",read_file("%(tmpdir)s/time_"+str(n)+"_"+prg+".txt"),flags=re.MULTILINE|re.DOTALL)
    res += "\\n"
write_file("%(data_fn)s",res);' """ % locals()
    jobs3.append(AbstractJob.create(cmd=cmd,id="datafile",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs3,"datafile").run()
    assert res==True


    tex1_fn = os.path.join(tmpdir,"a.tex")
    tex2_fn = os.path.join(tmpdir,"b.tex")
    tex3_fn = os.path.join(tmpdir,"c.tex")
    remote_write_file(tex1_fn,r"""
\documentclass[a4paper]{article}
\usepackage{tikz}
\usepackage{pgfplots}
\pagestyle{empty}
\begin{document}
\begin{tikzpicture}
\begin{axis}[
enlargelimits=true,
xtick={"""+",".join([str(x) for x in N_VALUES if x>50])+r"""},
xlabel=Number of residues,
ylabel=Time (in sec.)]
\addplot table[x=n,y=cl] {datafile.dat};
\end{axis}
\end{tikzpicture}
\end{document}
""")
    remote_write_file(tex2_fn,r"""
\documentclass[a4paper]{article}
\usepackage{tikz}
\usepackage{pgfplots}
\pagestyle{empty}
\begin{document}
\begin{tikzpicture}
\begin{axis}[
legend pos=outer north east,
enlargelimits=true,
ymode=log,
log basis y=2,
xtick={"""+",".join([str(x) for x in N_VALUES if x>50])+r"""},
xlabel=Number of residues,
ylabel=Time (in sec.)]
\addplot table[x=n,y=cl] {datafile.dat};
\addplot table[x=n,y=rnaview] {datafile.dat};
\addplot table[x=n,y=mc_annotate] {datafile.dat};
\addplot table[x=n,y=fr3d] {datafile.dat};
\legend{ClaRNA,RNA-View,MC-Annotate,FR3D}
\end{axis}
\end{tikzpicture}
\end{document}
""")

    remote_write_file(tex3_fn,r"""
\documentclass[a4paper]{article}
\usepackage{tikz}
\usepackage{pgfplots}
\pagestyle{empty}
\begin{document}
\begin{tikzpicture}
\begin{axis}[
legend pos=outer north east,
enlargelimits=true,
ymode=log,
log basis y=2,
xtick={"""+",".join([str(x) for x in N_VALUES if x>50])+r"""},
xlabel=Number of residues,
ylabel=Time (in sec.)]
\addplot table[x=n,y=cl] {datafile.dat};
\addplot table[x=n,y=cl_old] {datafile.dat};
\addplot table[x=n,y=rnaview] {datafile.dat};
\addplot table[x=n,y=mc_annotate] {datafile.dat};
\addplot table[x=n,y=fr3d] {datafile.dat};
\legend{ClaRNA,ClaRNA-noOPT,RNA-View,MC-Annotate,FR3D}
\end{axis}
\end{tikzpicture}
\end{document}
""")

    jobs4 = []
    out_fn1 = os.path.join(env.SCRIPTS_DIR,"doc","nar-fig-clarna-time1.png")
    cmd = """pdflatex -interaction=batchmode a.tex && pdfcrop a.pdf && convert -density 400 a-crop.pdf %(out_fn1)s""" % locals()
    jobs4.append(AbstractJob.create(cmd=cmd,id="latex-a",dir=tmpdir))
    out_fn2 = os.path.join(env.SCRIPTS_DIR,"doc","nar-fig-clarna-time2.png")
    cmd = """pdflatex -interaction=batchmode b.tex && pdfcrop b.pdf && convert -density 400 b-crop.pdf %(out_fn2)s""" % locals()
    jobs4.append(AbstractJob.create(cmd=cmd,id="latex-b",dir=tmpdir))
    out_fn3 = os.path.join(env.SCRIPTS_DIR,"doc","nar-fig-clarna-time3.png")
    cmd = """pdflatex -interaction=batchmode c.tex && pdfcrop c.pdf && convert -density 400 c-crop.pdf %(out_fn3)s""" % locals()
    jobs4.append(AbstractJob.create(cmd=cmd,id="latex-c",dir=tmpdir))
    res = AbstractJobsCollection.create(jobs4,"latex").run()
    assert res==True

@task
@ins("localhost")
def nar_images_cis_trans():
    """generate nar-fig-cis-trans.png - example with problems with cis/trans orientation"""
    png_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure6-cis-trans-rgb.png")
    tiff_cmyk_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure6-cis-trans-cmyk.tiff")

    # 1N8R:A2728:A2707 - REV:HwHw_pairing_parallel_trans_one_hbond_101
    # 1Q86:A2728:A2707 - REV:HwHw_pairing_parallel_cis_one_hbond_101
    d1 = "1N8R:A2728:A2707"
    d1_org = "1N8R:A2889:A2868"
    d2 = "1Q86:A2728:A2707"
    d2_org = "1Q86:A2889:A2868"
    w = 1500
    h = 1000
    label_size = 40

    tmpdir = "/tmp/nar_img_cis_trans"
    remote_mkdir(tmpdir)
    jobs1 = []
    tmp1_fn = os.path.join(tmpdir,"a1.pdb")
    tmp2_fn = os.path.join(tmpdir,"a2.pdb")
    tmp_png_fn = os.path.join(tmpdir,"a.png")
    cmd = "./gen-pdb.py --doublet-ids='%(d1)s' --type=full -o '%(tmp1_fn)s'" % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb",dir=env.SCRIPTS_DIR))
    cmd = "./gen-pdb.py --doublet-ids='%(d2)s' --type=full -o '%(tmp2_fn)s'" % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs1,"gen-pdbs").run()
    assert res==True

    jobs2 = []

    pymol_script = """
set max_threads, 1;
load %(tmp1_fn)s;
load %(tmp2_fn)s;
remove (hydro);
viewport %(w)d,%(h)d;
bg_color white;

color red, a1;
color green, a2;

show sticks, a1;
show sticks, a2;
hide lines, a1;
hide lines, a2;
set stick_radius, 0.15;
set stick_transparency, 0.1

set label_size, %(label_size)d;
set label_color, red, a1;
set label_color, green, a2;
set label_position, (-0.5,2,0);

label 'a1' and chain B and name N4, '%(d1_org)s';
label 'a2' and chain B and name O2, '%(d2_org)s';

set all_states,1;
orient animate=0;

ray %(w)d,%(h)d;
png %(tmp_png_fn)s,%(w)d;
""" % locals()
    remote_write_file(os.path.join(tmpdir,"script.pml"),pymol_script)
    cmd = """pymol -c -d "%(pymol_script)s" """ % locals()
    cmd += """ ; convert "%(tmp_png_fn)s" -trim +repage "%(png_fn)s" """ % locals()
    cmd += " ; " + _cmd_png_rgb_to_cmyk(png_fn,tiff_cmyk_fn)

    jobs2.append(AbstractJob.create(cmd=cmd,id="pymol",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs2,"pymol").run()
    assert res==True

@task
@ins("localhost")
def nar_images_bph():
    """generate nar-fig-bph.png - example with W_345BPh

Example of three doublets with different base-phosphate interactions detected by FR3D: 
4ABR:A566:A560 (W_3BPh), 3UYD:A566:A560 (W_4BPh) and 3UZM:A566:A560 (W_5BPh).    
"""
    png_fn = os.path.join(env.SCRIPTS_DIR,"doc","supp-figure-2-bph.png")

    d1 = "2QBK:B2684:B2693" # W_3BPh
    d1_org = "2QBK:B2747:B2756"
    d2 = "2QBC:B2684:B2693" # W_4BPh
    d2_org = "2QBC:B2747:B2756"
    
    d1 = "4ABR:A544:A538" # W_3BPh
    d1_label = "4ABR:A566:A560 (W_3BPh)"
    d2 = "3UYD:A544:A538" # W_4BPh
    d2_label = "3UYD:A566:A560 (W_4BPh)"
    d3 = "3UZM:A544:A538" # W_5BPh
    d3_label = "3UZM:A566:A560 (W_5BPh)"
    
    w = 1500
    h = 1000
    label_size = 30

    tmpdir = "/tmp/nar_img_bph"
    remote_mkdir(tmpdir)
    jobs1 = []
    tmp1_fn = os.path.join(tmpdir,"a1.pdb")
    tmp2_fn = os.path.join(tmpdir,"a2.pdb")
    tmp3_fn = os.path.join(tmpdir,"a3.pdb")
    tmp_png_fn = os.path.join(tmpdir,"a.png")
    cmd = "./gen-pdb.py --doublet-ids='%(d1)s' --type=base-phosphate -o '%(tmp1_fn)s'" % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb1",dir=env.SCRIPTS_DIR))
    cmd = "./gen-pdb.py --doublet-ids='%(d2)s' --type=base-phosphate -o '%(tmp2_fn)s'" % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb2",dir=env.SCRIPTS_DIR))
    cmd = "./gen-pdb.py --doublet-ids='%(d3)s' --type=base-phosphate -o '%(tmp3_fn)s'" % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb3",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs1,"gen-pdbs").run()
    assert res==True

    jobs2 = []

    pymol_script = """
set max_threads, 1;
load %(tmp1_fn)s;
load %(tmp2_fn)s;
load %(tmp3_fn)s;
remove (hydro);
viewport %(w)d,%(h)d;
bg_color white;

color red, a1;
color green, a2;
color blue, a3;
color gray, 'a1' and chain A;
color gray, 'a2' and chain A;
color gray, 'a3' and chain A;

show sticks, a1;
show sticks, a2;
show sticks, a3;
hide lines, a1;
hide lines, a2;
hide lines, a3;
set stick_radius, 0.15;
set stick_transparency, 0.1

set label_size, %(label_size)d;
set label_color, red, a1;
set label_color, green, a2;
set label_color, blue, a3;

set label_position, (-0.2,-1.7,0), a1;
set label_position, (-0.15,-1.5,0), a2;
set label_position, (-0.5,1.5,0), a3;

label 'a1' and chain B and name O5', '%(d1_label)s';
label 'a2' and chain B and name O3', '%(d2_label)s';
label 'a3' and chain B and name OP1, '%(d3_label)s';

set all_states,1;
orient animate=0;

ray %(w)d,%(h)d;
png %(tmp_png_fn)s,%(w)d;
""" % locals()
    remote_write_file(os.path.join(tmpdir,"script.pml"),pymol_script)
    cmd = """pymol -c -d "%(pymol_script)s" """ % locals()
    cmd += """ ; convert "%(tmp_png_fn)s" -trim +repage "%(png_fn)s" """ % locals()

    jobs2.append(AbstractJob.create(cmd=cmd,id="pymol",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs2,"pymol").run()
    assert res==True

@task
@ins("localhost")
def nar_images_loose_class():
    png_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-fig-loose-class-rgb.png")

    w = 1500
    h = 1000
    label_size = 40

    tmpdir = "/tmp/nar_img_loose_classes"
    remote_mkdir(tmpdir)
    jobs1 = []
    tmp_fn = os.path.join(tmpdir,"a.pdb")
    tmp_png_fn = os.path.join(tmpdir,"a.png")
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    cmd = """./gen-pdb.py --type=full --input-json="%(groups_fn)s" --filter-keys=classifier/bp/HS_cis/AG -o "%(tmp_fn)s" """ % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs1,"gen-pdbs").run()
    assert res==True

    jobs2 = []

    pymol_script = """
set max_threads, 1;
load %(tmp_fn)s;
remove (hydro);
viewport %(w)d,%(h)d;
bg_color white;

color green, a;
color red, ('a' and chain A)

set_view (\
    -0.844317377,   -0.213440791,    0.491497248,\
    -0.165031895,    0.976236105,    0.140447438,\
    -0.509796441,    0.037470378,   -0.859476984,\
     0.000000000,    0.000000000,  -62.454067230,\
    -3.066821098,    0.976489067,   -0.834476471,\
    49.975852966,   74.932296753,  -20.000000000 )

set all_states,1;

ray %(w)d,%(h)d;
png %(tmp_png_fn)s,%(w)d;
""" % locals()
    remote_write_file(os.path.join(tmpdir,"script.pml"),pymol_script)
    cmd = """pymol -c -d "%(pymol_script)s" """ % locals()
    cmd += """ ; convert "%(tmp_png_fn)s" -trim +repage "%(png_fn)s" """ % locals()

    jobs2.append(AbstractJob.create(cmd=cmd,id="pymol",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs2,"pymol").run()
    assert res==True

@task
@ins("localhost")
def nar_images_bad_mc():
    png_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure7-bad-mc-rgb.png")
    tiff_cmyk_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure7-bad-mc-cmyk.tiff")
    
    # classifier returns SH_cis
    # MC returns         WH_cis
    doublet_ids = ["3UZL:A362:A364", "3UYE:A1054:A1055", "3OI3:A1009:A1010", "4A19:11203:11204", "3UZ8:A1049:A1050"]
    doublet_ids = ["3UZL:A362:A364", "3UYE:A1054:A1055", "3OI3:A1009:A1010", "4A19:11203:11204", "3UZ8:A1049:A1050", "3I8F:A1363:A1365", "2UU9:A1122:A1123", "3OI3:A1269:A1271", "1N33:A486:A487", "4DR5:A486:A487", "3KIY:A1047:A1048", "3UZM:A486:A487", "2WDG:A1122:A1123", "3HUX:A1020:A1021", "1XNQ:A486:A487", "2XZM:A435:A437", "3V29:A1032:A1033", "2Y16:A486:A487", "3UZG:A486:A487", "3JYV:A426:A428", "2Y17:A1360:A1362", "3NKB:B31:B33", "3KIW:A1047:A1048", "3UZK:A1054:A1055", "2WRL:A1048:A1049", "2Y13:A1048:A1049", "3G4S:01064:01065", "3UYE:A1366:A1368", "2XZN:A435:A437", "3UZ2:A1051:A1052", "3TVH:A1049:A1050", "2WRJ:A1048:A1049", "3UZ3:A486:A487", "2WDH:A1122:A1123", "3UZH:A1049:A1050", "3OI5:A1009:A1010", "3UZL:A486:A487", "2Y10:A1105:A1107", "4A1D:11203:11204", "3UZN:A1049:A1050", "3KCR:81319:81321", "3OI5:A1269:A1271"]

    w = 1500
    h = 1000
    label_size = 40

    tmpdir = "/tmp/nar_img_bad_mc"
    remote_mkdir(tmpdir)
    jobs1 = []
    tmp_fn = os.path.join(tmpdir,"a.pdb")
    tmp_png_fn = os.path.join(tmpdir,"a.png")
    doublet_ids_str = ",".join(doublet_ids)
    cmd = """./gen-pdb.py --type=full --doublet-ids=%(doublet_ids_str)s -o "%(tmp_fn)s" """ % locals()
    jobs1.append(AbstractJob.create(cmd=cmd,id="gen-pdb",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs1,"gen-pdbs").run()
    assert res==True

    jobs2 = []

    pymol_script = """
set max_threads, 1;
load %(tmp_fn)s;
remove (hydro);
viewport %(w)d,%(h)d;
bg_color white;

color green, a;
color red, ('a' and chain A)

orient, animate=0;

set all_states,1;

ray %(w)d,%(h)d;
png %(tmp_png_fn)s,%(w)d;
""" % locals()
    remote_write_file(os.path.join(tmpdir,"script.pml"),pymol_script)
    cmd = """pymol -c -d "%(pymol_script)s" """ % locals()
    cmd += """ ; convert "%(tmp_png_fn)s" -trim +repage "%(png_fn)s" """ % locals()
    cmd += " ; " + _cmd_png_rgb_to_cmyk(png_fn,tiff_cmyk_fn)

    jobs2.append(AbstractJob.create(cmd=cmd,id="pymol",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs2,"pymol").run()
    assert res==True

@task
@ins("localhost")
def nar_dataset_info():
    DATA_SETS = ["data-small.txt", "data-medium.txt", "data-training.txt", "data-bench.txt"]

    jobs = []
    res_fns = []
    for ds in DATA_SETS:
        res_fn = "/tmp/info-"+ds
        res_fns.append(res_fn)
        pdb_list = read_file(ds)
        pdb_list_str = "["+",".join('"'+x+'"' for x in pdb_list.split("\n") if x!="")+"]"
        cmd = """python -c '
from utils import load_pdb; import re;
n=0; rn=0; 
nn=0; rrn=0;
for pdbid in %(pdb_list_str)s:
    c = len(list(load_pdb("gc-data/pdb_files/"+pdbid+"_rna.pdb.gz").get_residues()))
    if re.match("^zz",pdbid):
        nn+=1
        rrn+=c
    else:
        n+=1
        rn+=c
print "%(ds)s, structures: "+str(n)+" total resiudes: "+str(rn)+" artificial structures: "+str(nn)+" artificial residues: "+str(rrn)
' > %(res_fn)s""" % locals()

        jobs.append(AbstractJob.create(cmd=cmd,id="dataset_info_%s"%ds,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"dataset_info").run()
    assert res==True
    for res_fn in res_fns:
        print remote_read_file(res_fn)

@task
@ins("peyote2")
def nar_fig_cl_correlation(force=True):
    prepare_peyote2()
    
    DATA_SETS = ["data-small.txt", "data-medium.txt", "data-training.txt", "data-bench.txt", "data-all.txt"]
    jobs = []
    for ds in DATA_SETS:
        b = ds.replace("data-","").replace(".txt","")
        inp_fn = os.path.join(env.DJANGO_DIR, ds)
        out_png_fn = os.path.join(env.DATA_DIR, "cl-correlation-%s.png" % b)
        out_png2_fn = os.path.join(env.DATA_DIR, "cl-eval-%s.png" % b)
        out_json_fn = os.path.join(env.DATA_DIR, "cl-correlation-%s.json" % b)
        data_dir = env.DATA_DIR
        if (not force or force=="False") and my_files.exists(out_json_fn):
            extra = " --input-json=%(out_json_fn)s " % locals()
        else:
            extra = ""
        cmd = """./_tmp_20130429_classifiers_correlation.py --data-dir=%(data_dir)s -i %(inp_fn)s %(extra)s --output-png=%(out_png_fn)s --output-png2=%(out_png2_fn)s --output-json=%(out_json_fn)s """ % locals()
        cmd += """ && convert "%(out_png_fn)s" -trim +repage "%(out_png_fn)s" && convert "%(out_png2_fn)s" -trim +repage "%(out_png2_fn)s" """ % locals()
        jobs.append(AbstractJob.create(cmd=cmd,id="nar_fig_cl_correlation_%s"%b,dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"nar_fig_cl_cor").run()
    assert res==True

def _cmd_pdf_rgb_to_cmyk(in_fn,out_fn,dpi=300):
    return "convert -density %(dpi)d '%(in_fn)s' -colorspace CMYK '%(out_fn)s'" % locals()
    # to niestety nie dziala!
    return "gs -dSAFER -dBATCH -dNOPAUSE -dNOCACHE -sDEVICE=pdfwrite -sColorConversionStrategy=CMYK -dProcessColorModel=/DeviceCMYK -sOutputFile=%(out_fn)s %(in_fn)s" % locals()

def _cmd_png_rgb_to_cmyk(in_fn,out_fn):
    return "convert '%(in_fn)s' -colorspace CMYK -compress lzw '%(out_fn)s'" % locals()

@task
@ins("localhost")
def nar_fig_roc_curve():
    jobs = []
    inp_fn = FileNamesObject.eval_fn(setname="bench",t="eval3",data_dir=env.DATA_DIR)
    out_pdf_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure1-roc-rgb.pdf")
    out_pdf_cmyk_fn = os.path.join(env.SCRIPTS_DIR,"doc","nar-figure1-roc-cmyk.pdf")
    cmd = """./eval-graphs.py -i %(inp_fn)s --output-pdf=%(out_pdf_fn)s """ % locals()
    cmd += " && "
    cmd += _cmd_pdf_rgb_to_cmyk(out_pdf_fn,out_pdf_cmyk_fn)
    # cmd += """ && convert "%(out_png_fn)s" -trim +repage "%(out_png_fn)s" """ % locals()
    jobs.append(AbstractJob.create(cmd=cmd,id="nar_fig_roc",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"nar_fig_roc").run()
    assert res==True
    
@task
@ins("peyote2")
def tw_20130515_prepare_new_zz_files():
    pdbids = open("data-new-zz.txt").read().strip().split("\n")
    prepare_pdb_files(pdbids,force=True)

@task
@ins("peyote2")
def tw_20130603_prepare_new_zz_files():
    pdbids = ["zz%02x"%(112+x) for x in range(17)]
    prepare_peyote2()
    prepare_pdb_files(pdbids,force=True)


@task
def upload_supp_materials():
    django_dir = INSTALLATIONS['rnacontacts-vm'].get('DJANGO_DIR')
    local("rsync -avP ../clarna/doc/nar-supp/supp-* rnacontacts-vm:%(django_dir)s/supp_data/" % locals())
    
    
@task
@ins("peyote2")
def tw_20130615_compute_phosphate_stacking_cl():
    tmp_compute_classifier_new2(only_sc='other3')
    tmp_compute_classifier_new3(only_sc='other3')

@task
@ins("peyote2")
def tw_20130617_classifier_cross_validation():
    prepare_peyote2()

    for num in [5]:
        gr = "cross"+str(num)+"t"
        c_fn = FileNamesObject.groups_fn(setname="bench",cross_validation_num=num,reduced=True,data_dir=env.DATA_DIR)

        tmp_compute_classifier_new1(gr=gr)
        tmp_compute_classifier_new2(gr=gr)
        tmp_compute_classifier_new3(gr=gr)
        
        run_classifier("data-training.txt",extra_options_combine="--only-doublets-from=%s" % c_fn)
        
        run("cd '%s' && tar cvzf cross-%s.tgz classifier.*.json.gz eval/new-classifier-eval-training.json.gz eval/new-classifier-full-stats-training.json.gz eval/new-classifier-full-stats2-training.json.gz"%(env.DATA_DIR,num))

@task
@ins("peyote2")
def tw_20130619_rpt_phosphate_stackings():
    prepare_peyote2()

    rpt_dir = os.path.join(env.DATA_DIR,"tmp","rpt_ph_stackings")
    json = os.path.join(rpt_dir,"all.json")
    inp1 = FileNamesObject.eval_fn(setname="training",t="eval1",data_dir=env.DATA_DIR)
    inp2 = FileNamesObject.eval_fn(setname="bench",t="eval1",data_dir=env.DATA_DIR)
    cmd = "./merge-dicts.py --only-keys='evaluation/prev-undetected/phosphate-stacking.*' --input='%(inp1)s,%(inp2)s' -o '%(json)s'" % locals()
    jobs = []
    jobs.append(AbstractJob.create(cmd=cmd,id="prepare_json",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"prepare_json").run()
    assert res==True

    fn_counts = os.path.join(rpt_dir,"counts.txt")
    cmd = "./json-tool.py --only-keys='evaluation/prev-undetected/phosphate-stacking.*' -i '%(json)s' --apply-map-op=print-lengths | tee '%(fn_counts)s'" % locals()
    jobs = []
    jobs.append(AbstractJob.create(cmd=cmd,id="count_elements",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"count_elements").run()
    assert res==True
    
    jobs = []
    for n_type in N_TYPES:
        for desc in ['phosphate-stacking1','phosphate-stacking2']:
            group_pdb = os.path.join(rpt_dir,n_type+"_"+desc[-1]+".pdb")
            group_json = os.path.join(rpt_dir,n_type+"_"+desc[-1]+".json")
            k = 'evaluation/prev-undetected/%(desc)s/%(n_type)s' % locals()
            cmd = "./gen-pdb.py --type=base-phosphate --input-json='%(json)s' -o '%(group_pdb)s' --save-json --translate-doublet-ids --filter-keys='%(k)s'" % locals()
            jobs.append(AbstractJob.create(cmd=cmd,id="gen_%s_%s"%(n_type,desc[-1]),dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"gen_pdbs_and_json").run()
    assert res==True

    # ./json-tool.py --apply-map-op=print-lengths --only-keys='evaluation/prev-undetected/phosphate-stacking1/.*' -i gc-data/eval/new-classifier-eval-small.json.gz 

@task
@ins("peyote2")
def tw_20130619_nar_supp_fig_merged_bph_classes():
    prepare_peyote2()

    rpt_dir = os.path.join(env.DATA_DIR,"tmp","merged_bph_classes")
    pdbids = open("data-all.txt").read().strip().split("\n")
    # pdbids = pdbids[0:100]
    DATA_DIR = env.DATA_DIR
    
    GROUP1 = [('W_%dBPh'%x,'GA') for x in (3,4,5)]
    GROUP2 = [('H_%dBPh'%x,'CA') for x in (7,8,9)]
    jobs = []
    out_list = {}
    for pdbid in pdbids:
        cmd = ""
        contacts_fr = PDBObject.pdb_fn(pdbid,"contacts_FR",data_dir=env.DATA_DIR)
        for desc,n_type in GROUP1+GROUP2:
            out_fn = os.path.join(rpt_dir,pdbid+"."+desc+"."+n_type+".json.gz")
            key = desc+"_"+n_type
            if not out_list.has_key(key):
                out_list[key] = []
            out_list[key].append(out_fn)
            if cmd!="":
                cmd += " && "
            cmd += "./extract-doublets-from-graph.py -i '%(contacts_fr)s' --filter-full-desc='%(desc)s' --filter-n-type='%(n_type)s' --pdb-id='%(pdbid)s' -o %(out_fn)s" % locals()
        jobs.append(AbstractJob.create(cmd=cmd,id="gen_bph_%s"%(pdbid),dir=env.SCRIPTS_DIR))
    #res = AbstractJobsCollection.create(jobs,"extract_bph_doublets").run()
    #assert res==True

    # merging
    jobs = []
    jobs2 = []
    p345 = []
    p789 = []
    for key,out_fns in sorted(out_list.items()):
        input_fn = os.path.join(rpt_dir,key+".in")
        merge_fn = os.path.join(rpt_dir,key+".json.gz")
        pdb_fn = os.path.join(rpt_dir,key+".pdb.gz")
        if re.match("^W_[345]",key):
            p345.append(pdb_fn)
        if re.match("^H_[789]",key):
            p789.append(pdb_fn)
        remote_write_file(input_fn, "\n".join(out_fns))

        jobs.append(AbstractJob.create(cmd="./merge-dicts.py --input-from-file=%(input_fn)s -o %(merge_fn)s" % locals(), id="merge-%s"%(key),mem="4GB",dir=env.SCRIPTS_DIR))
        jobs2.append(AbstractJob.create(cmd="./gen-pdb.py --add-normalized-base --only-chain-b --type=base-phosphate --limit=500 --input-json=%(merge_fn)s -o %(pdb_fn)s" % locals(), id="gen-pdb-%s"%(key),mem="4GB",dir=env.SCRIPTS_DIR))
    res2 = AbstractJobsCollection.create(jobs,"merge").run()
    assert res2==True

    res3 = AbstractJobsCollection.create(jobs2,"gen-pdb").run()
    assert res3==True

    e1 = 'color green, %s;' % os.path.basename(p345[0]).split(".")[0]
    e1 += ';color red, %s' % os.path.basename(p345[1]).split(".")[0]
    e1 += ';color blue, %s' % os.path.basename(p345[2]).split(".")[0]
    e2 = 'color green, %s;' % os.path.basename(p789[0]).split(".")[0]
    e2 += ';color red, %s' % os.path.basename(p789[1]).split(".")[0]
    e2 += ';color blue, %s' % os.path.basename(p789[2]).split(".")[0]
    view = """set_view (\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000,  -56.317768097,\
     0.633658409,    2.570053101,   -0.102146149,\
    44.401359558,   68.234176636,  -20.000000000 )"""
    e1 += ";"+view
    e2 += ";"+view

    jobs = [
                pdb_to_png_job(p345,os.path.join(rpt_dir,"W_345BPh.png"),width=2000,height=2000,extra_script=e1),
                pdb_to_png_job(p789,os.path.join(rpt_dir,"H_789BPh.png"),width=2000,height=2000,extra_script=e2),
            ]
    res3 = AbstractJobsCollection.create(jobs,"gen-png").run()
    assert res3==True

    cmd = """
convert -verbose "W_345BPh.png" -trim +repage \
   -background none -gravity center -extent "100%"  +repage \
   -resize 1000x1000 +repage \
   +repage -gravity center -extent 1000x1000 \
   "a1.png"

convert -verbose "H_789BPh.png" -trim +repage \
   -background none -gravity center -extent "130%"  +repage \
   -resize 1000x1000 +repage \
   +repage -gravity center -extent 1000x1000 \
   "a2.png"

montage -background white -label "W_3/4/5BPh GA" a1.png -label "H_7/8/9BPh CA" a2.png -tile 'x1' -geometry 1000x1000 -pointsize 45 nar_supp_fig_merged_bph_classes.png
"""
    jobs = [AbstractJob.create(cmd=cmd,id="im",dir=rpt_dir)]
    res = AbstractJobsCollection.create(jobs,"im").run()
    assert res==True

    out_png_fn = os.path.join("doc","nar_supp_fig_merged_bph_classes.png")
    get(os.path.join(rpt_dir,"nar_supp_fig_merged_bph_classes.png"),out_png_fn)
    
@task
@ins("peyote2")
def nar_supp_data_ref_set_zip():
    prepare_peyote2()
    
    rpt_dir = os.path.join(env.DATA_DIR,"tmp","nar_supp_ref_set")
    jobs = []
    groups_fn = FileNamesObject.groups_fn(setname="training",reduced=True,data_dir=env.DATA_DIR)
    for sc in CL_CAT.keys():
        if not sc in ['bp','stacking','base-phosphate','base-ribose']:
            continue
        for n_type in N_TYPES:
            for desc in CL_CAT[sc]:
                safe_desc = desc.replace("<","l").replace(">","r")
                # generate json/zip for sc/desc/n_type
                base_id = "%s_%s" % (safe_desc,n_type)
                json_fn = os.path.join(rpt_dir,sc,base_id+".json")
                pdb_fn = os.path.join(rpt_dir,sc,base_id+".pdb")
                cmd = "( ./json-tool.py -i %(groups_fn)s -o %(json_fn)s --only-keys=\"classifier/%(sc)s/%(desc)s/%(n_type)s\" --apply-map-op=copy-values" % locals()
                cmd += " && ./translate-doublets-ids.py -i %(json_fn)s -o %(json_fn)s )" % locals()
                jobs.append(AbstractJob.create(cmd=cmd,id="g-json-%s-%s-%s"%(sc,safe_desc,n_type),mem="4GB",dir=env.SCRIPTS_DIR))

                cmd = "./gen-pdb.py --input-json=%(groups_fn)s -o %(pdb_fn)s --type=%(sc)s --filter-keys=\"classifier/%(sc)s/%(desc)s/%(n_type)s\" --randomize --limit=500" % locals()
                jobs.append(AbstractJob.create(cmd=cmd,id="g-pdb-%s-%s-%s"%(sc,safe_desc,n_type),mem="4GB",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"prepare_ref_set_elems").run()
    assert res==True

@task
@ins("peyote2")
def tmp_20131021_close_bp_ribose(setname="data-training.txt"):
    # ponowne przeliczanie tym razem z podzialem na <3.0A, 3.0..4.0 4.0..5.0 w dniu 2013-11-15
    prepare_peyote2()

    TASK_TMP_DIR = os.path.join(env.DATA_DIR,"tmp","close_bp_ribose")

    pdbids = open(setname).read().strip().split("\n")

    jobs = []
    list_fns = []
    for pdbid in pdbids:
	    inp_fn = PDBObject.pdb_fn(pdbid,"pdb",data_dir=env.DATA_DIR)
	    groups_fn = PDBObject.pdb_fn(pdbid,"groups",data_dir=env.DATA_DIR)
	    out_fn = os.path.join(TASK_TMP_DIR,pdbid+"_rna.contacts.graph.json.gz")
	    eval_fn = os.path.join(TASK_TMP_DIR,pdbid+"_rna.eval.json.gz")
	    list_fns.append(eval_fn)
	    
	    libs = ",".join([os.path.join(env.SCRIPTS_DIR,"lib",k) for k in ["close-bp-ribose.json"]])

	    cmd = "./clarna.py -i %(inp_fn)s --lib=%(libs)s --normalize-graph --save-graph=%(out_fn)s --ignore-bad-doublets" % locals()
	    cmd += " && ./combine-contact-graphs.py --pdb-id=%(pdbid)s --classifier-graph=%(out_fn)s --groups=%(groups_fn)s -o %(eval_fn)s " % locals()

	    jobs.append(AbstractJob.create(cmd=cmd,id="cl-run-%s"%(pdbid),mem="4GB",dir=env.SCRIPTS_DIR))
    res = AbstractJobsCollection.create(jobs,"run_classifier_bp_ribose").run()
    assert res==True

    # merge
    input_fn = os.path.join(TASK_TMP_DIR,"merge.in")
    merge_fn = os.path.join(TASK_TMP_DIR,"eval.json.gz")

    remote_write_file(input_fn, "\n".join(list_fns))

    jobs = []
    jobs.append(AbstractJob.create(cmd="./merge-dicts.py --only-keys='evaluation/prev-undetected/.*' --input-from-file=%(input_fn)s -o %(merge_fn)s" % locals(), id="merge",mem="4GB",dir=env.SCRIPTS_DIR))

    res = AbstractJobsCollection.create(jobs,"merge").run()
    assert res==True

