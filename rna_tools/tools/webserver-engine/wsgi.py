"""
WSGI config for rnamasonry project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

import os,sys
root = os.path.join(os.path.dirname(__file__), '')
sys.path.insert(0, root)
sys.path.insert(0, os.path.join(root, '../lib/python2.7/site-packages'))

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.settings")

#import django.core.handlers.wsgi
#application = django.core.handlers.wsgi.WSGIHandler()

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
