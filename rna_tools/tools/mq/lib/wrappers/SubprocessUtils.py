#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module wraps functions from Python's standard library
subprocess module for easier use.
"""

from time import sleep
from subprocess import Popen, PIPE


def run_command(cmd, args, cmd_input='', env=None,
        stdout_file=None, stderr_file=None, timeout=None, verbose=False):
    """Run command using subprocess, optionally send input or
    set environment variables

    In order to make timeout work, suprocess is started form a seperate thread.
    It's probably the only way it can work if this function is run outside
    the main thread.

    Arguments:
      * cmd = command
      * args = list of arguments
      * cmd_input = input from standard input
      * env = optional, dictionary of environment variables; if not set
        command is run within an empty environment
      * stdout_file = optional, filename for saving output
      * stderr_file = optional, filename for saving standard error output
      * timeout = time after which subprocess is killed; optional, if not
        set process can run infinitelly

    Output:
      * return code for command
    """
    # TODO: check if kwargs are not a better way to write this function
    extended_cmd = timeout is None and [cmd] or ['timeout', str(timeout), cmd]

    if verbose: print(' '.join(extended_cmd), ' '.join(args))
    command = Popen(extended_cmd + args,
                stdout=PIPE, stdin=PIPE, stderr=PIPE,
                shell=False, env=env)
    stdout, stderr = command.communicate(cmd_input.encode())
    if verbose:
        print(stdout.decode())
        print(stderr.decode())
    # make it possible to put output into a file
    if stdout_file is not None:
        f = open(stdout_file, 'w')
        f.write(stdout.decode())
        f.close()
    if stderr_file is not None:
        f = open(stderr_file, 'w')
        f.write(stderr.decode())
        f.close()
    return command.returncode


def run_parallel(commands, num_cores):
    '''Run commands from a list in parallel.

    Arguments:
      * commands = list of commands (as strings)
      * num_cores = maximum number of commands that are allowed to running
        at the same time
    '''
    running = []
    while commands:
        while len(running) < num_cores and commands:
            running.append(Popen(commands.pop(),
                                 stdout=PIPE, stdin=PIPE, stderr=PIPE,
                                 shell=False))
        running = [i for i in running if i.poll() is None]
        sleep(0.01)  # without this delay, appending is infinite
