from django.contrib import admin
from app.models import Job, Setting


class SettingAdmin(admin.ModelAdmin):
    pass


class JobAdmin(admin.ModelAdmin):
    list_display = ('id', 'url', 'job_title', 'email', 'ip', 'seq_len', 'status', 'nsteps', 'created', 'error_text_stub', 'comment', 'progress', 'killed')
    save_on_top = True
    search_fields = ['job_id', 'comment', 'email', 'status', 'job_title']
    list_per_page = 40

admin.site.register(Job, JobAdmin)
admin.site.register(Setting, SettingAdmin)
