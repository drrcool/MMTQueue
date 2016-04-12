from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render
from .forms import MaskUploadForm
from .models import Semester

def upload_file(request):
    if request.method == 'POST':
        form = MaskUploadForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return HttpResponseRedirect('upload_success/')
    else:
         form = MaskUploadForm()
    return render(request, 'catalogs/upload.html', {'form': form})

def upload_success(request):
    return HttpResponse("Your mask was successfully uploaded.")

def index(request):
    return HttpResponse("At the index page.")
