#!/usr/bin/python

import subprocess

cmd="/home/oge/bin/lx24-amd64/qstat -u '*' | grep -v 'hqw' "
out = subprocess.getoutput(cmd)
cc = 0
for l in out.split('\n'):
    if l:
        c = l.split()[-1]
        try:
            cc += int(c)
        except ValueError:
            pass
print((cc, cc/1000.0, -(1000 - cc)))
