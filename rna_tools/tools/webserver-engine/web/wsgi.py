"""dont remove, this is used for dev server"""
import os,sys
root = os.path.join(os.path.dirname(__file__), '')
sys.path.insert(0, root)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "web.settings")

#import django.core.handlers.wsgi
#application = django.core.handlers.wsgi.WSGIHandler()
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
