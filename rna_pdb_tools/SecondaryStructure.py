#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Seq and secondary structure prediction"""

import os
import tempfile
import shutil
import subprocess
from rpt_config import *

def draw_ss(title, seq, ss, img_out, resolution=2, verbose=False):
    """Draw Secondary Structure using VARNA (you need correct configuration for this).

    If everything is OK, return None, if an error (=exception) return stderr.

    Can be used with http://geekbook.readthedocs.io/en/latest/rna.html"""
    curr = os.getcwd()
    os.chdir(VARNA_PATH)#VARNAv3-93-src')
    if verbose: print VARNA_PATH
    t = tempfile.NamedTemporaryFile(delete=False)
    t.name += '.png'

    cmd = 'java -cp ' + VARNA_JAR_NAME + ' fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN ' + seq + " -structureDBN '" + ss + "' -o " + t.name + " -title '" + title + "' -resolution '" + str(resolution) + "'"
    if verbose: print cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out = p.stderr.read().strip()
    os.chdir(curr)
    if out.find('Exception') > -1:
        return stderr
    else:
        if verbose: print t.name
        shutil.move(t.name, img_out)

    
if __name__ == '__main__':
    seq = 'AAAAAAA'
    ss =  '((...))'
    img_out = 'demo.png'
    draw_ss('rna', seq, ss, img_out)
    print 'Made %s' % img_out
