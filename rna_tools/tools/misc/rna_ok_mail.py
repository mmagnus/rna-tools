#!/usr/bin/python
import argparse
import smtplib
import os
from rna_tools.rna_tools_config import MAIL_USERNAME, MAIL_PASSWORD


def send_email_done(job_id, email, msg_text, subject):
    """
    job_id
    email email to
    """
    # Credentials (if needed)
    fromaddr = 'rna-tools'
    toaddrs  = email

    subject = 'rna-tools done' 
    msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s\r\n"
           % (fromaddr, toaddrs, subject, msg_text))

    # The actual mail send
    #print 'sending...'
    #print msg
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(MAIL_USERNAME, MAIL_PASSWORD)
    server.sendmail(fromaddr, toaddrs, msg)
    server.quit()

def get_parser():
    parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                            action="store_true", help="be verbose")
    parser.add_argument("-s", "--subject", 
                           help="subject", default = '')
    parser.add_argument("-m", "--msg", 
                           help="msg", default = '')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    job_id = 'rna-tools'
    print('sending done...')
    send_email_done(job_id, 'mag_dex@o2.pl', args.msg, args.subject)
