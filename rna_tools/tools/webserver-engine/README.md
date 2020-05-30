webserver-engine
===================================

The light-weight, complete webserver engine used by me for servers: NPDock (RNA/DNA-protein docking method, http://genesilico.pl/NPDock/), SimRNAweb (RNA 3D structure prediction method, http://iimcb.genesilico.pl/SimRNAweb/), mqapRNA (RNA 3D quality control, http://iimcb.genesilico.pl/mqapRNA/), RNAMasonry (building RNA models from recurrent 3D motifs, http://genesilico.pl/rnamasonry).

The engine handels all commont tasks for a simple webserver:

- accept a job and saves the data into an internal database,
- controals a Python deamon for job execution on the server,
- handle progress of a job,
- present the results,
- sends mail for a job at the start and end,
- sends reports to admins on server usage,
- with full-fledged Django admin to control jobs.

The server is easily customizable for new applications. The server by design and very light-weight and can be easily encapsulated into a single, simple virtual environment.

Let me know if you need any help to set it up for yourself (--@mmagnus).

![](docs/demo.png)

Install
-------------------------------------------------------------------------------

    (py27) [mx] src$ pip install virtualenv

    (py27) [mx] src$ python -m virtualenv rnamasonry_env
    created virtual environment CPython2.7.16.final.0-64 in 385ms
      creator CPython2Posix(dest=/Users/magnus/work/src/rnamasonry_env, clear=False, global=False)
      seeder FromAppData(download=False, pip=latest, setuptools=latest, wheel=latest, via=copy, app_data_dir=/Users/magnus/Library/Application Support/virtualenv/seed-app-data/v1.0.1)
      activators PythonActivator,CShellActivator,FishActivator,PowerShellActivator,BashActivator

    (rnamasonry_env) (py27) [mx] rnamasonry$ git:(master) pip install -r install.txt
    DEPRECATION: Python 2.7 reached the end of its life on January 1st, 2020. Please upgrade your Python as Python 2.7 is no longer maintained. A future version of pip will drop support for Python 2.7. More details about Python 2 support in pip, can be found at https://pip.pypa.io/en/latest/development/release-process/#python-2-support
    Requirement already satisfied: Django==1.8 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 1)) (1.8)
    Requirement already satisfied: argparse==1.2.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 2)) (1.2.1)
    Requirement already satisfied: decorator==4.0.6 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 4)) (4.0.6)
    Requirement already satisfied: django-extensions==1.6.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 5)) (1.6.1)
    Requirement already satisfied: django-ipware==1.1.6 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 6)) (1.1.6)
    Requirement already satisfied: ipdb==0.8.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 7)) (0.8.1)
    Requirement already satisfied: ipython==4.0.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 8)) (4.0.1)
    Requirement already satisfied: ipython-genutils==0.1.0 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 9)) (0.1.0)
    Requirement already satisfied: numpy==1.10.2 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 10)) (1.10.2)
    Requirement already satisfied: path.py==8.1.2 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 11)) (8.1.2)
    Requirement already satisfied: pexpect==4.0.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 12)) (4.0.1)
    Requirement already satisfied: pickleshare==0.5 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 13)) (0.5)
    Requirement already satisfied: ptyprocess==0.5 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 14)) (0.5)
    Requirement already satisfied: simplegeneric==0.8.1 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 15)) (0.8.1)
    Requirement already satisfied: six==1.10.0 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 16)) (1.10.0)
    Requirement already satisfied: traitlets==4.0.0 in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from -r install.txt (line 17)) (4.0.0)
    Requirement already satisfied: wsgiref==0.1.2 in /Users/magnus/miniconda2/envs/py27/lib/python2.7 (from -r install.txt (line 18)) (0.1.2)
    Requirement already satisfied: appnope; sys_platform == "darwin" in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from ipython==4.0.1->-r install.txt (line 8)) (0.1.0)
    Requirement already satisfied: gnureadline; sys_platform == "darwin" and platform_python_implementation == "CPython" in /Users/magnus/work/src/rnamasonry_env/lib/python2.7/site-packages (from ipython==4.0.1->-r install.txt (line 8)) (8.0.0)

Run
-------------------------------------------------------------------------------
    
    (py27) [mx] rnamasonry$ git:(master) ✗ source ../bin/activate
    
    (rnamasonry_env) (py27) [mx] rnamasonry$ git:(master) ✗ python manage.py runserver --settings web.settings 0.0.0.0:8667
    Performing system checks...

    System check identified some issues:

    WARNINGS:
    app.Job.interpret_occupancy: (fields.W122) 'max_length' is ignored when used with IntegerField
        HINT: Remove 'max_length' from field
    app.Job.nsteps: (fields.W122) 'max_length' is ignored when used with IntegerField
        HINT: Remove 'max_length' from field
    app.Job.seq_len: (fields.W122) 'max_length' is ignored when used with IntegerField
        HINT: Remove 'max_length' from field
    app.Job.status: (fields.W122) 'max_length' is ignored when used with IntegerField
        HINT: Remove 'max_length' from field

    System check identified 4 issues (0 silenced).
    April 16, 2020 - 18:20:58
    Django version 1.8, using settings 'web.settings'
    Starting development server at http://0.0.0.0:8667/
    Quit the server with CONTROL-C.

Configure
-------------------------------------------------------------------------------

Go to web.settings for full configuration:

    SERVER_NAME = 'RNAMasonry'
    ADMIN_JOBS_URL = 'http://0.0.0.0:8667/admin/app/job/'
    DISK_TO_TRACK = '/dev/mapper/rnamasonry--vm--vg-root'
    import os; 
    PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)) 
    JOBS_PATH = PATH + os.sep + 'media' + os.sep + 'jobs'
    URL_JOBS = "http://iimcb.genesilico.pl/rnamsonry/jobs/"  # required /
    URL = "http://iimcb.genesilico.pl/rnamasonry" # / is not needed"
    DEBUG = True
    TEMPLATE_DEBUG = False
    PATH_TO_RM = ''
    POWER_USERS = ['magnus@genesilico.pl']
    SEND_EACH_MAIL_TO_ADMIN = False
    SERVER_REPLY_MAIL = 'rnamasonry@genesilico.pl'
    ADMIN_MAIL = 'magnus@genesilico.pl' # used for sendmail
    ADMINS = (
         ('magnus', 'mag_dex@o2.pl'),
         #('chojnowski', 'gchojnowski@gmail.com')
    )
