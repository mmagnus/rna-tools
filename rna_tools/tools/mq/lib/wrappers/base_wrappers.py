#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Module with classes that are parents for all wrappers.
"""

import os
import sys
DIRNAME = os.path.dirname(__file__)
PARENT_DIRNAME = os.path.abspath(os.path.join(DIRNAME, os.path.pardir))
sys.path.append(PARENT_DIRNAME)

import subprocess
import signal
import tempfile
import logging
import shutil
import rna_tools.rna_tools_config as Config

# ogger = logging.getLogger(__name__)
DIRNAME = os.path.dirname(__file__)
LOG_DIRECTORY = Config.LOG_DIRECTORY


class WrapperError(Exception):
    """Base exception for all errors."""
    pass


class BadSequenceLength(WrapperError):
    """Exception to rise if sequence is too long for a particular method."""
    pass


class Wrapper(object):
    """General wrapper class for secondary structure calculation services. 
    This class does nothing on its own and should be used as parent class for 
    classes for more specific wrappers.
    """
    program_name = 'wrapper'
    metric = 'kcal/mol'
    max_seq_len = 100000000
    start_dir = None
    path = None
    output_fn = None
    error_fn = None
    pdb_fixes = []


    def __init__(self, sequence='', job_id=None):
        self.temp_files = []
        self.output_files = []
        self.sequence = sequence
        self.seq_len = len(self.sequence)
        self.result = None
        self.pdb_fixes = []
        self.job_id = job_id
        if job_id is not None:
            self.job_id_str = str(job_id) + ': '
        else:
            self.job_id_str = ''
        # set up logging
        self.logger = logging.getLogger('wrappers.' + self.__class__.__name__)
        if not self.logger.handlers:
            self.log_handler = logging.FileHandler(
                    LOG_DIRECTORY + os.sep + self.__class__.__name__ + '.log'
                    )
            self.logger.addHandler(self.log_handler)
        else:
            self.log_handler = self.logger.handlers[0]
        #self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
                )
        self.log_handler.setFormatter(formatter)
    
    def _prepare_sandbox(self):
        self.sandbox_dir = tempfile.mkdtemp()
        return self.sandbox_dir

    def _error_on_bad_len(self, full_seq):
        if len(full_seq) > self.max_seq_len:
            raise BadSequenceLength('Sequence must be shorter than %i' %
                self.max_seq_len)

    def _check_seq_length(self, sequence):
        self._error_on_bad_len(sequence)

    def __repr__(self):
        return "<%s wrapper>" % self.program_name

    def _pre_run(self):
        os.chdir(self.sandbox_dir)

    def _post_run(self):
        os.chdir(self.start_dir)

    def _prepare_stderr_stdout(self):
        # create output file
        self.output_file = os.path.join(self.path, self.output_fn)
        self.stdout = open(self.output_file, 'w')
        # create error file
        self.error_file = os.path.join(self.path, self.error_fn)
        self.stderr = open(self.error_file, 'w')

    def run(self):
        """The most important method in this class: run the wrapped program.
        """
        raise NotImplementedError

    def get_result(self):
        """Extract and return result from previously run program.
        """
        return self.result

    def get_fixes(self):
        """Return list of fixes applied to PDB file
        """
        return self.pdb_fixes
    
    def log(self, message, level='info', verbose=False):
        message = self.job_id_str + message
        if verbose:
            print(message)
        if level == 'info':
            self.logger.info(message)
        elif level == 'debug':
            self.logger.debug(message)
        elif level == 'warning':
            self.logger.warning(message)
        elif level == 'error':
            self.logger.error(message)
        elif level == 'critical':
            self.logger.critical(message)
        else:
            raise Exception('Wrong logging level')
        


class ProgramWrapper(Wrapper):
    """General wrapper class for calculation programs. It should be used as 
    a parent class for all programs that can be run locally.

    To get a working wrapper, you need to overload "run" method. You also need 
    to set "program_name" and "executable" as strings with program name and 
    executable file name (or lists if you want to run more than one program 
    in this wrapper). All executables are expected to be in BIN_PATH which 
    is set in Config.

    In some cases, you might need to overload also "get_result" or "cleanup" 
    and maybe other methods, but do it carefully and call "super()."
    """
    program_name = None
    executable = None
    flags = []
    input_fn = 'input_sequence'
    output_fn = 'output'
    error_fn = 'errors'
    fasta_input = False
    clustal_input = False
    max_seq_len = 10000

    def __init__(self, sequence='', seq_name='query_sequence',
        path=Config.BIN_PATH, flags=[], job_id=None):
        """
        Specify input file path (path+filename) and the
        output directory for results of all local programs.
        """
        super(ProgramWrapper, self).__init__(sequence, job_id)
        self._check_seq_length(self.sequence)
        self.path = path
        self.start_dir = Config.TMP_PATH
        self.seq_name = seq_name
        self.stdin = None
        self.stdout = None
        self.stderr = None
        self.sandbox_dir = None
        self.add_flags(*flags)
        self._prepare_sandbox()
        self.sandbox()
        self.path = self.sandbox_dir
        self._prepare_files()
        self._prepare_args()

    def add_flags(self, *flags):
        """Add command line options for program.
        """
        for flag in flags:
            self.flags.append(flag)

    def _prepare_files(self):
        # create input file
        self.input_file = os.path.join(self.path, self.input_fn)
        in_f = open(self.input_file, 'w')
        in_f.write(str(self.sequence))
        in_f.close()
        self._prepare_stderr_stdout()

    def _prepare_args(self):
        if self.executable is None:
            return
        if type(self.executable) != type([]):
            args = [os.path.join(self.path, self.executable)]
            args.extend(self.flags)
            self.flags = args
        else:
            args = []
            self.flags = []
            for e in self.executable:
                args.append([os.path.join(self.path, e)].extend(self.flags))
                self.flags.append(args)

    def sandbox(self):
        """Create a sandbox in temporary folder for running program.
        """
        return 
        # dont copy programs!
        ## if self.executable is not None:
        ##     if type(self.executable) != type([]):
        ##         binary = os.path.join(self.path, self.executable)
        ##         shutil.copy(binary, self.sandbox_dir)
        ##     else:
        ##         for e in self.executable:
        ##             binary = os.path.join(self.path, e) # path e rasp_fd
        ##             shutil.copy(binary, self.sandbox_dir)

    def _check_r_code(self, r_code):
        """Check return code from program and log errors if it's not 0.
        """
        if r_code != 0:
            self.logger.error('Wrong exit code: ' + str(r_code))
            stdout = open(self.output_file)
            stderr = open(self.error_file)
            try:
                self.logger.debug('STDOUT: ' + stdout.read())
                self.logger.debug('STDERR: ' + stderr.read())
            except:
                self.logger.debug('Couldn\'t get stdout and stderr')
            stdout.close()
            stderr.close()
            raise WrapperError('Error during processing command')

    def run(self):
        """Run program with previously set flags. This method should work if
        program wants some command some options and nothing more, otherwise
        you should override it.
        """
        self.logger.info('Running program')
        self._pre_run()
        self.proc = subprocess.Popen(args=self.flags, stdin=self.stdin,
                stdout=self.stdout, stderr=self.stderr)
        r_code = self.proc.wait()
        self._post_run()
        if self.stdout is not None:
            self.stdout.close()
        if self.stderr is not None:
            self.stderr.close()
        self._check_r_code(r_code)
        self.logger.info('Run finished')
        return self.get_result()
    
    def kill(self):
        os.killpg(self.proc.pid, signal.SIGKILL)

    def cleanup(self):
        """Run this method when you are done with this wrapper. It deletes
        sandbox, some of it's childrem may do some other cleanup here.
        """
        shutil.rmtree(self.sandbox_dir)

