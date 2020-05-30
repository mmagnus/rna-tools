#!/usr/bin/python

import smtplib
import os
from django.core.wsgi import get_wsgi_application
os.environ['DJANGO_SETTINGS_MODULE'] = 'web.settings'
application = get_wsgi_application()
from sendmail_secret import USERNAME, PASSWORD
from web.settings import URL, SEND_EACH_MAIL_TO_ADMIN, ADMIN_MAIL


def send_email_start(job_id, email):
    """
    job_id
    email email to
    """
    # Credentials (if needed)
    fromaddr = SERVER_NAME + ' <' + SERVER_REPLY_MAIL + '>'
    toaddrs  = email

    subject = '%s - RNAMasonry start' % job_id
    msg_text = """Dear User,

Your job has just started (%s) on the RNAMasonry server.

%s/jobs/%s""" %  (job_id, URL, job_id)

    msg_text += """
It might take between 3 and 6 hours. You will get an e-mail when your job is done.

The results are kept on the server for 4 days.

Thanks for using the RNAMasonry server.

This is an automated message. Please don't reply to this email. To contact RNAMasonry (%s) administrators, please send an email to: rnamasonry@genesilico.pl""" % (URL)

    msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s\r\n"
           % (fromaddr, toaddrs, subject, msg_text))

    # The actual mail send
    #print 'sending...'
    #print msg
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(USERNAME, PASSWORD)
    server.sendmail(fromaddr, toaddrs, msg)
    if SEND_EACH_MAIL_TO_ADMIN:
        server.sendmail(fromaddr, ADMIN_MAIL, msg)
    server.quit()


def send_email_done(job_id, email):
    """
    job_id
    email email to
    """
    # Credentials (if needed)
    fromaddr = 'RNAMasonry <rnamasonry@genesilico.pl>'
    toaddrs  = email

    subject = '%s - RNAMasonry end' % job_id
    msg_text = """Dear User,

Your job with job id %s has been finished on the RNAMasonry server.

%s/jobs/%s""" %  (job_id, URL, job_id)

    msg_text += """

The results are kept on the server for 4 days.

Thanks for using the RNAMasonry server.

This is an automated message. Please don't reply to this email. To contact RNAMasonry (%s) administrators, please send an email to: rnamasonry@genesilico.pl""" % (URL)

    msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s\r\n"
           % (fromaddr, toaddrs, subject, msg_text))

    # The actual mail send
    #print 'sending...'
    #print msg
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(USERNAME, PASSWORD)
    server.sendmail(fromaddr, toaddrs, msg)
    if SEND_EACH_MAIL_TO_ADMIN:
        server.sendmail(fromaddr, ADMIN_MAIL, msg)
    server.quit()

if __name__ == '__main__':
    job_id = 'tRNA_with_SAXS-5ad01dbb'
    send_email_start(job_id, 'mag_dex@o2.pl')
    send_email_done(job_id, 'mag_dex@o2.pl')
