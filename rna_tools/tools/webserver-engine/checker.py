#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
00 18 * * * /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/checker.sh

2204bf9d
862a9fae
5cb972ad
6fc6390f
check.py

drwxrwxr-x 2 ubuntu ubuntu   4096 May 18 04:26 2204bf9d
drwxrwxr-x 2 ubuntu ubuntu   4096 May 18 05:08 862a9fae
drwxrwxr-x 2 ubuntu ubuntu   4096 May 18 05:35 5cb972ad
drwxrwxr-x 2 ubuntu ubuntu   4096 May 18 08:49 6fc6390f
-rw-rw-r-- 1 ubuntu ubuntu    968 May 18 09:07 check.py

Add to crontab

* 18 * * * /home/rnamasonry/rnamasonryweb_env/rnamasonry-web/checker.sh
"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

import os
import subprocess
import smtplib

from sendmail_secret import USERNAME, PASSWORD
SERVER_NAME = 'rna-tools.online'


def send_mail_to(mail, txt):
    fromaddr = SERVER_NAME + ' report <magnus@genesilico.pl>'
    subject = SERVER_NAME + ' report'
    toaddrs = mail
    msg_text = txt
    msg = ("""From: %s\r\nTo: %s\r\nSubject: %s\r\nMIME-Version: 1.0\r\nContent-Type: text/html\r\nContent-Disposition: inline\r\n<html>\r\n<body>\r\n<pre style="font: monospace; font-size: 11px;">\r\n\r\n%s\r\n""" % (fromaddr, toaddrs, subject, msg_text))

    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(USERNAME, PASSWORD)
    server.sendmail(fromaddr, mail, msg)
    server.quit()


def _get_jobs_():  # not used !
    """time of job is calculated beased on the files!
    If there is now file then you don't estimate the time

    Keep jobs that are on JOBS_TO_KEEP list"""
    jobs = models.Job.objects.filter().order_by("-id")[:200]
    text = SERVER_NAME + '- checker - scripts shows 100 last jobs!\n\n'
    if True:
        for j in jobs:
            status = j.get_status()
            if status == 'finished with errors':
                status = '!!!!!!!!'
            text += str(j.created) + " <b>" + status.ljust(10) + "</b> " + j.email + ' ' + j.job_title + " " \
              +URL_JOBS + " " + j.job_id + ' ' + ' '
            if j.error_text:
                text += '\n' + j.error_text
            text += '\n'
    else:
        for j in jobs:
            text += "- " + j.get_status() + " " + "-" * 80 + "\n" + j.email + "\n" + j.job_title + '\n'
            text +=URL_JOBS + j.job_id + '\n'
            text += str(j.created) + '\n'
            text += '\n'
    return text


def run_cmd(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def get_space():
    """Set the correct disk to track"""
    cmd = "df -h | grep " + DISK_TO_TRACK
    out, err = run_cmd(cmd)
    return out


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    #parser.add_argument("file", help="", default="") # nargs='+')
    return parser


def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    import glob
    import os

    os.chdir('/home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/media/jobs')
    
    files = glob.glob("*")
    files.sort(key=os.path.getmtime)
    files = [os.path.basename(f) for f in files]

    #os.system('ls -l | grep ' + files[-100])
    txt = ''
    for d in files[-100:]: ## lenght set up here !!!!!!!!!
        if d == 'check.py':
            continue
        txt += '= ' + d + ' ' + '=' * 60 + ' \n'
        try:
            with open(d + '/run.sh') as f:
                t, err = exe('ls -l | grep ' + d)
                t += '\n'
                t += f.read()
                txt += t + '\n'
            with open(d + '/log.txt') as f:
                t = f.read()

            for c, l in enumerate(t.split('\n')[-40:]):
                txt += l + '\n'

        except FileNotFoundError:
            #print('empty')
            pass
            
    mail = '\n\n'.join([SERVER_NAME, 'ADMIN_JOBS_URL'])
    #mail += '\n\n' + get_space() + '\n\n'
    mail += txt
    print(txt)
    for admin in ['m.magnus@imol.institute']:
        send_mail_to(admin, txt)
        pass
