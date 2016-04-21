from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    # upload_mask
    url(r'^upload_mask/$', views.upload_file, name='upload_mask'),
    url(r'^upload_mask/upload_success/$', views.upload_success, name='success')
]
