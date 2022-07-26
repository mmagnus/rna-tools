#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 ScenePNG: wrote 2052x1350 pixel image to file "out.png".
(base) [mx] structures-sessions$ pymol -c /Users/magnus/work/cwc15/structures-sessions/cwc15_conserv.pse  -d 'ray; png out.png; quit'
"""
import sys
import argparse
import tempfile
import os
import subprocess
import platform
from rna_tools.rna_tools_config import BIN_PATH

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-b", "--black",
                        action="store_true", help="black bg")
    parser.add_argument("-d", "--detailed", default = 200,
                        help="show lines/sticks for files with # of lines slower than, default 200")
    parser.add_argument("-r", "--rainbow",
                        action="store_true", help="rainbow")
    parser.add_argument("files", help="", default="", nargs='+')
    return parser

if 'Darwin' == platform.system():
    is_mac = True
    #BIN = '/usr/local/bin' # TODO: fix it
    BIN = '/opt/homebrew/bin/' #BIN_PATH # "/opt/homebrew/bin/"
else:
    is_mac = False
    BIN = '/usr/bin'
    
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    print(BIN)
    
    for file in args.files:
        tf = tempfile.NamedTemporaryFile(delete=False)
        f = tf.name + '.png'
        file = file.replace(" ", "\\ ")
        rainbow = ''
        if args.rainbow:
            rainbow = 'util.chainbow;'

        # pdb mode
        if file.endswith('.pdb'):
            def exe(cmd):
                o = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out = o.stdout.read().strip().decode()
                err = o.stderr.read().strip().decode()
                return out, err
            cmd = BIN + '/pymol -c ' + file + " -d 'print(len([x for x in cmd.get_model().atom]))'"
            out, err = exe(cmd)
            # ugly way to get n of atoms
            """
            PyMOL>run ~/work/src/rna-tools/rna_tools/tools/PyMOL4RNA/bucket/mlk4.py
            65
            """
            n = int(out.split('\n')[-1])
            # if something short
            # aaa this will not work for pse files (!)
            #if len(open(file).readlines()) < args.detailed:
            sh = ''
            if n <  300: # atoms
                sh = ' show lines; '#sticks; ' # lines;
            os.system(BIN + '/pymol -c ' + file + " -d 'set ray_opaque_background, off;" + rainbow + " show cartoon; " + sh + "; rr; save " + f + "; quit'") #  ray 300,300,renderer=0 ray 800, 800;
        else:  # pse mode  for pse do nothing! # set ray_opaque_background, off; 
            os.system(BIN + '/pymol -c ' + file + " -d 'save " + f + "; quit'")

        if args.verbose:
            print(f)

        if is_mac:
            fcrop = f.replace('.png', '32.png')
            cmd = BIN + "/convert " + f + " -gravity center -crop 3:3 +repage " + fcrop
            os.system(cmd)
            cmd = 'unset PYTHONPATH && ' + BIN + 'fileicon set ' + file + ' ' + fcrop
        else:
            cmd = 'gvfs-set-attribute ' + file + ' metadata::custom-icon file://' + f  # f so the file before convert
        os.system(cmd)
