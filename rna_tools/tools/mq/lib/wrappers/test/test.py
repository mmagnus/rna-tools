#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test unit for lib.wrappers"""

import os, sys
from shutil import copyfile
from unittest import TestCase, main


DIRNAME = os.path.dirname(__file__)

sys.path.append(os.path.abspath(os.path.join(DIRNAME, os.path.pardir)))
from base_wrappers import ProgramWrapper
from SubprocessUtils import run_command


TEST_COMMAND_PATH = 'echo'


class DummyWrapper(ProgramWrapper):
    program_name = 'dummy'
    max_seq_len = 2000
    executable = 'echo'  
    
    def __init__(self, sequence, seq_name='query_sequence'):
        super(DummyWrapper, self).__init__(sequence, seq_name)

    def run(self, pdb_file):    # pdb_file is an argument here only because it is like that in real wrappers
        run_command(self.sandbox_dir + os.sep + self.executable, ['foobar'])
        self.result = 42


class WrapperlibTests(TestCase):
    """ Tests for formatlib"""
    def set_up(self):
        """Wrapperlib: Sets up test suite."""
        pass

    def test_0(self):
        """Wrapperlib: test dummy wrapper"""
        wrapper = DummyWrapper('')
        wrapper.run('')
        wrapper.cleanup()
        self.assertEqual(42, wrapper.get_result())


if __name__ == '__main__':
    main()
