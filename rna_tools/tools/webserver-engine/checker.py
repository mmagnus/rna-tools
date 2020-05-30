#!/usr/bin/python
"""
Add to crontab

* 18 * * * /home/rnamasonry/rnamasonryweb_env/rnamasonry-web/checker.sh
"""
import os
import subprocess
import smtplib

from sendmail_secret import USERNAME, PASSWORD
from django.core.wsgi import get_wsgi_application
os.environ['DJANGO_SETTINGS_MODULE'] = 'web.settings'
application = get_wsgi_application()
from app import models
from web import settings
USE_TZ = False


def send_mail_to(mail, txt):
    fromaddr = settings.SERVER_NAME + ' report <magnus@genesilico.pl>'
    subject = settings.SERVER_NAME + ' report'
    toaddrs = mail
    msg_text = txt
    msg = ("""From: %s\r\nTo: %s\r\nSubject: %s\r\nMIME-Version: 1.0\r\nContent-Type: text/html\r\nContent-Disposition: inline\r\n<html>\r\n<body>\r\n<pre style="font: monospace">\r\n\r\n%s\r\n""" % (fromaddr, toaddrs, subject, msg_text))

    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(USERNAME, PASSWORD)
    server.sendmail(fromaddr, mail, msg)
    server.quit()


def get_jobs():
    """time of job is calculated beased on the files!
    If there is now file then you don't estimate the time

    Keep jobs that are on JOBS_TO_KEEP list"""
    jobs = models.Job.objects.filter().order_by("-id")[:100]
    text = settings.SERVER_NAME + '- checker - scripts shows 100 last jobs!\n\n'
    if True:
        for j in jobs:
            status = j.get_status()
            if status == 'finished with errors':
                status = '!!!!!!!!'
            text += str(j.created) + " <b>" + status.ljust(10) + "</b> " + j.email + ' ' + j.job_title + " " \
              + settings.URL_JOBS + " " + j.job_id + ' ' + ' '
            if j.error_text:
                text += '\n' + j.error_text
            text += '\n'
    else:
        for j in jobs:
            text += "- " + j.get_status() + " " + "-" * 80 + "\n" + j.email + "\n" + j.job_title + '\n'
            text += settings.URL_JOBS + j.job_id + '\n'
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
    cmd = "df -h | grep " + settings.DISK_TO_TRACK
    out, err = run_cmd(cmd)
    return out


if __name__ == '__main__':
    txt = '\n\n'.join([settings.SERVER_NAME, settings.ADMIN_JOBS_URL])
    txt += '\n\n' + get_space() + '\n\n'
    txt += get_jobs()
    for admin in settings.ADMINS:
        send_mail_to(admin[1], txt)

