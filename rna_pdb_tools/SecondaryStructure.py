#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Seq and secondary structure prediction"""

import os
import tempfile
import shutil

VARNA_PATH  = '/Users/magnus/skills/rnax/varna_tut/'

def draw_ss(title,seq, ss, img_out):
    """"""
    curr = os.getcwd()
    os.chdir(VARNA_PATH)#VARNAv3-93-src')
    print os.getcwd()
    t = tempfile.NamedTemporaryFile(delete=False)
    t.name += '.png'
    os.system('java -cp VARNA.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN ' + seq + " -structureDBN '" + ss + "' -o " + t.name + " -title " + title + " -resolution '2.0'")
    os.chdir(curr)
    print img_out
    shutil.move(t.name, img_out)
    
if __name__ == '__main__':
    seq = 'AAAAAAA'
    ss =  '((...))'
    img_out = 'out.png'
    draw_ss('rna', seq, ss, img_out)
