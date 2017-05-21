#!/usr/bin/env python
"""A super simple script to get some statistics of who is running at a cluster

Set MAX_JOBS to calc % of usage, it's an approximation of max number of jobs, e.g. peyote ~1k"""
import commands
MAX_JOBS = 800

print 'MAX_JOBS:', MAX_JOBS

def stats_for_cluster():
    """get stats (#jobs) per cluster"""
    
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
    return cc

def stats_for_user():
    """get stats (#jobs) per user"""
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
    return cc

def per_user():
    """get stats (#cpus) per user"""
    # {'deepak': 160, 'azyla': 8, 'magnus': 755}
    cmd="/home/oge/bin/lx24-amd64/qstat -u '*'"
    out = commands.getoutput(cmd).strip()
    per_user = {}
    for l in out.split('\n'):
        if l.startswith('job-ID') or l.startswith('---'):
            pass
        else:
            ll = l.split()
            user = ll[3]
            cpus = int(ll[-1])
            if per_user.has_key(user):
                per_user[user] += cpus
            else:
                per_user[user] = cpus            
    return per_user

if __name__ == '__main__':
    cc = stats_for_cluster()
    print '#jobs cluster', cc, 'load: ', cc/float(MAX_JOBS), ' to use:', MAX_JOBS - cc
    cc = stats_for_user()
    print '#jobs you    ', cc, 'load: ', cc/float(MAX_JOBS), ' to use:', MAX_JOBS - cc
    print per_user()
