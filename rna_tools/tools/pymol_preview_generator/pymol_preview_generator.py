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
    
def exe(cmd):
                o = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out = o.stdout.read().strip().decode()
                err = o.stderr.read().strip().decode()
                return out, err

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    print(BIN)
    
    for f in args.files:
        if 0:
            tf = tempfile.NamedTemporaryFile(delete=False)
            fcover = tf.name + '.png'
            f = f.replace(" ", "\\ ")
        else:
            fcover = f.replace('.pdb', '.png')
        rainbow = ''
        if args.rainbow:
            rainbow = 'util.chainbow;'

        # pdb mode
        if f.endswith('.pdb'):
            cmd = BIN + '/pymol -c ' + f + " -d 'print(len([x for x in cmd.get_model().atom]))'"
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
                #sh = ' show lines; '#sticks; ' # lines;
                sh = ' show mesh; ' #set ray_opaque_background, off;
            #  bg_color black; 
            os.system(BIN + '/pymol -c ' + f + " -d 'set ray_opaque_background, on; hide cartoon;" + sh + "; save " + fcover + "; quit'") #ray 1000,1000; #  ,renderer=0            #os.system(BIN + '/pymol -c ' + f + " -d 'set ray_opaque_background, off;" + rainbow + " show cartoon; " + sh + "; rr; save " + fcover + "; quit'") 
        else:  # pse mode  for pse do nothing! # set ray_opaque_background, off; 
            os.system(BIN + '/pymol -c ' + f + " -d 'save " + fcover + "; quit'")

        if args.verbose:
            print(f)

        if is_mac:
            print('MAC')
            fcrop = fcover.replace('.png', '32.png')
            cmd = BIN + "/convert " + fcover + " -gravity center -crop 3:3 +repage " + fcrop
            os.system(cmd)
            print(cmd)
            cmd = 'unset PYTHONPATH && ' + BIN + 'fileicon set ' + f + ' ' + fcrop
            print(cmd)
            os.system(cmd)
            os.remove(fcover)
            os.remove(fcrop)
        else:
            cmd = 'gvfs-set-attribute ' + file + ' metadata::custom-icon file://' + f  # f so the file before convert
            os.system(cmd)
