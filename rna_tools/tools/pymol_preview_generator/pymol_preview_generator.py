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

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-b", "--black",
                        action="store_true", help="black bg")
    parser.add_argument("-r", "--rainbow",
                        action="store_true", help="rainbow")
    parser.add_argument("files", help="", default="", nargs='+')
    return parser

if 'Darwin' == platform.system():
    is_mac = True
    BIN = '/usr/local/bin' 
else:
    is_mac = False
    BIN = '/usr/bin'
    
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for file in args.files:
        tf = tempfile.NamedTemporaryFile(delete=False)
        f = tf.name + '.png'
        file = file.replace(" ", "\\ ")
        rainbow = ''
        if args.rainbow:
            rainbow = 'util.chainbow;'

        os.system(BIN + '/pymol -c ' + file + " -d 'set ray_opaque_background, off;" + rainbow + " show cartoon; save " + f + "; quit'") #  ray 300,300,renderer=0 ray 800, 800;
        if args.verbose:
            print(f)

        if is_mac:
            fcrop = f.replace('.png', '32.png')
            cmd = BIN + "/convert " + f + " -gravity center -crop 3:3 +repage " + fcrop
            os.system(cmd)
            cmd = 'unset PYTHONPATH && /usr/local/bin/fileicon set ' + file + ' ' + fcrop
        else:
            cmd = 'gvfs-set-attribute ' + file + ' metadata::custom-icon file://' + f  # f so the file before convert
        os.system(cmd)
