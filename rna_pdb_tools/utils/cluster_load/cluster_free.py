#!/usr/bin/python

import commands

cmd="/home/oge/bin/lx24-amd64/qstat -u '*'"
out = commands.getoutput(cmd).strip()
cc = 0
for l in out.split('\n'):
    if l.strip():
        if l.find('hqw') > -1:
            continue
        c = l.split()[-1]
        try:
            cc += int(c)
        except ValueError:
            pass
print 'jobs:', cc, 'load (1k max): ', cc/1000.0, ' to use:', 1000 - cc

cmd="/home/oge/bin/lx24-amd64/qstat "# -u '*'"
out = commands.getoutput(cmd)
cc = 0
for l in out.split('\n'):
    if l:
        c = l.split()[-1]
        try:
            cc += int(c)
        except ValueError:
            pass
print 'jobs:', cc, 'load (1k max): ', cc/1000.0, ' -- magnus'
