#!/home/ubuntu/miniconda3/bin/python
#/usr/bin/env python
"""
WSGI config for the project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/


User ubuntu
Group ubuntu

WSGIScriptAlias / /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/wsgi.py
WSGIPythonPath /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/


Alias /static/ /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/app/static/
<Directory /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/app/static/>
SetHandler None
Order allow,deny
Allow from all
</Directory>


"""

import os,sys
print('this')
root = os.path.join(os.path.dirname(__file__), '')
sys.path.insert(0, root)
sys.path.insert(0, '/home/ubuntu/.local/lib/python3.8/site-packages/')

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.settings")

#import django.core.handlers.wsgi
#application = django.core.handlers.wsgi.WSGIHandler()

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
