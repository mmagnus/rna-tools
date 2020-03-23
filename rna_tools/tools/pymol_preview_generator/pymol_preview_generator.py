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


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-b", "--black",
                        action="store_true", help="black bg")
    parser.add_argument("files", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for file in args.files:
        tf = tempfile.NamedTemporaryFile(delete=False)
        f = tf.name + '.png'
        file = file.replace(" ", "\\ ")
        os.system('/usr/local/bin/pymol -c ' + file + " -d 'set ray_opaque_background, off; save " + f + "; quit'") #  ray 300,300,renderer=0 ray 800, 800;
        if args.verbose:
            print(f)
        # crop
        fcrop = f.replace('.png', '32.png')
        cmd = "/usr/local/bin/convert " + f + " -gravity center -crop 3:3 +repage " + fcrop
        os.system(cmd)
        #cmd = 'source activate base && /usr/local/bin/fileicon set ' + args.file + ' ' + f
        cmd = 'unset PYTHONPATH && /usr/local/bin/fileicon set ' + file + ' ' + fcrop
        os.system(cmd)
