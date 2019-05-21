#!/usr/bin/python

from optparse import OptionParser, OptionGroup
import smtplib
import sys
import os
import re

class Mail():
    def __init__(self):
        self.txt = ''

    def add_txt(self, txt):
        self.txt += txt.strip() + '\n'
    
    def send_email(self, toaddrs, subject):
        """
        job_id
        email email to
        """
        # Credentials (if needed)  
        fromaddr = 'simRNAserver01@genesilico.pl'
        toaddrs  = toaddrs

        subject = subject

        msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s\r\n" 
               % (fromaddr, toaddrs, subject, self.txt))
        # The actual mail send  
        #print 'sending...'
        #print msg
        server = smtplib.SMTP('mailhub')
        server.sendmail(fromaddr, toaddrs, msg)  
        server.quit()
        print 'Sent from %s to %s msg with subject %s ' % (fromaddr, toaddrs, subject)

class Model():
    def __init__(self, id, f):
        self.id = id
        self.txt = open(f).read().strip()

    def get(self):
        txt = 'MODEL %s\n' % self.id
        txt += self.txt
        return txt

def sort_nicely( l ):
   """ Sort the given list in the way that humans expect.

   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l

if __name__ == '__main__':
    print 'rnapuzzle_sender'
    print

    parser = OptionParser('rnapuzzle_sender.py')

    parser.add_option("-d", "--dir",
                     action="store", type="string", dest="dir", help="")

    parser.add_option("-s", "--email_subject",
                      action="store", type="string", dest="email_subject", help="email subject")

    (opt, arguments) = parser.parse_args()
    
    if opt.dir == None or opt.email_subject == None:
        parser.print_help()
        sys.exit(1)

    # sort files
    files = os.listdir(opt.dir)
    files = sort_nicely(files)
    print 'files: ', ' '.join([f for f in files])
    c = 1
    print

    mail = Mail()
    for f in files:
        print opt.dir + os.sep + f
        m = Model(c, opt.dir + os.sep + f)
        c += 1
        mail.add_txt(m.get())
        del m
    print

    mail.send_email('magnus@genesilico.pl', opt.email_subject)
    mail.send_email('mboni@genesilico.pl', opt.email_subject)
    mail.send_email('iamb@genesilico.pl', opt.email_subject)
    #mail.send_email('ibmc.cnrs@gmail.com', opt.email_subject)
