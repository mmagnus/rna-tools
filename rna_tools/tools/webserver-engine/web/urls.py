from django.conf.urls import patterns, include, url
from django.conf.urls.static import static
from web import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    url(r'.well-known/pki-validation/(?P<fn>.*)$', 'app.views.ssl'),
    url(r'2coejc8r22i2/(?P<fn>.*)$', 'app.views.notes'),
    url(r'4j6scj6p82zw400/(?P<fn>.*)$', 'app.views.qr'),
    url(r'^images/(?P<fn>.*)$', 'app.views.image'),

    url(r'^fetch/(?P<job_id>.*)$','app.views.fetch', name='fetch'),
    url(r'^del/(?P<job_id_fn>.*)$','app.views.ajax_rm_file', name='ajax_rm_file'),
    url(r'^jobstatus/(?P<job_id>.*)/(?P<tool>.*)$', 'app.views.ajax_job_status', name='ajax_job_status'),
    url(r'^jobstatus/(?P<job_id>.*)$','app.views.ajax_job_status', name='ajax_job_status'),

    url(r'^download_project_dir/(?P<job_id>.*)$','app.views.download_project_dir', name='download_project_dir'),
    url(r'^stop/(?P<job_id>.*)$','app.views.stop'),

    url(r'^download_project_dir/(?P<job_id>.*)$', 'app.views.download_project_dir', name='download_project_dir'),
    url(r'^stop/(?P<job_id>.*)$', 'app.views.stop'),

    url(r'^tools/(?P<tool>.*)/(?P<job_id>.*)$','app.views.tool'),
    url(r'^demo/(?P<tool>.*)/(?P<job_id>.*)$','app.views.demo'),
    url(r'^tools','app.views.tools'),
                       
    url(r'^doc', 'app.views.help'),        
    url(r'^contact', 'app.views.contact'),
    url(r'^about', 'app.views.about'),

    url(r'%s/*$' %  settings.server_name_py, 'app.views.home', name='home'),
    url(r'%s' %  settings.server_name_py, 'app.views.home', name='home'),
    url(r'^$', 'app.views.home', name='home'),
    url(r'admin/', include(admin.site.urls)),
    url('upload/(?P<job_id>.*)$', 'app.views.file_upload'),
    url('^run/(?P<tool>.*)/(?P<job_id>.*)$', 'app.views.run'),

    # url(r'^web/', include('web.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:

) + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True) + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT) 

