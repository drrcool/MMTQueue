from django.shortcuts import render
from django.http import HttpResponse
from .models import Schedule

# Create your views here.
def display_schedule(request):
    sched = Schedule()
    for x in sched.objlist:
        print(x)
    return HttpResponse("This is a test")
