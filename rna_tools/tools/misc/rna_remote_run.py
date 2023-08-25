#!/usr/bin/env python
# -*- coding: utf-8 -*-

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    import os
    
    t = ''
    prev_t = ''
    while 1:
        os.system('./test.sh 2>&1 | tee test_tmp.compile')
        with open('test_tmp.compile') as f:
            t = f.read()
            
        t = t.replace('/net/holy-nfsisilon/ifs/rc_labs/eddy_lab/users/mmagnus/', '~/mnt/odx/')

        if t != prev_t:
            with open('test.compile', 'w') as f:
                f.write(t)
                prev_t = t
            
        
