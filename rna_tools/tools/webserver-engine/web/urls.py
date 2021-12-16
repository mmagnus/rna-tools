from django.conf.urls import patterns, include, url
from django.conf.urls.static import static
from web import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'submit', 'app.views.submit'),
    url(r'submit/', 'app.views.submit'),
    url(r'doc', 'app.views.help'),        
    url(r'contact', 'app.views.contact'),
    url(r'about', 'app.views.about'),
    url(r'jobs/(?P<job_id>.*)$','app.views.job'),
    url(r'jobstatus/(?P<job_id>.*)$','app.views.ajax_job_status', name='ajax_job_status'),
    url(r'download_project_dir/(?P<job_id>.*)$','app.views.download_project_dir', name='download_project_dir'),
    url(r'stop/(?P<job_id>.*)$','app.views.stop'),

    url(r'%s/jobs/(?P<job_id>.*)$' %  settings.server_name_py,'app.views.job'),

    url(r'%s/download_project_dir/(?P<job_id>.*)$' % settings.server_name_py,'app.views.download_project_dir', name='download_project_dir'),
    url(r'%s/jobstatus/(?P<job_id>.*)$' %  settings.server_name_py,'app.views.ajax_job_status', name='ajax_job_status'),
    url(r'%s/stop/(?P<job_id>.*)$' %  settings.server_name_py,'app.views.stop'),
    url(r'%s/*$' %  settings.server_name_py, 'app.views.home', name='home'),
    url(r'%s' %  settings.server_name_py, 'app.views.home', name='home'),
    url(r'^$', 'app.views.home', name='home'),
    url(r'admin/', include(admin.site.urls)),
    url('upload/', 'app.views.file_upload'),
    # url(r'^web/', include('web.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:

) + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True) + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT) 

