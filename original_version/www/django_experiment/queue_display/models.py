from django.db import models
import pandas as pd
from . import queue_tools
import time


class Observing_block(models.Model):
    """Create the class that holds a the parameters for one field."""

    # Parameters for this observation
    duration = models.IntegerField()
    n_visits = models.IntegerField()
    objID = models.CharField(max_length=50)
    start_time = models.CharField(max_length=50)
    end_time = models.CharField(max_length=50)

    # Object parameters from FLD file
    ra = models.CharField(max_length=50)
    dec = models.CharField(max_length=50)
    pmra = models.DecimalField(max_digits=5, decimal_places=3)
    pmdec = models.DecimalField(max_digits=4, decimal_places=2)
    pa = models.DecimalField(max_digits=6, decimal_places=2)
    moon = models.CharField(max_length=6)
    exptime = models.DecimalField(max_digits=9, decimal_places=2)
    nexp = models.IntegerField()
    repeats = models.IntegerField()
    oneper = models.BooleanField()
    priority = models.DecimalField(max_digits=5, decimal_places=1)
    grism = models.CharField(max_length=20)
    filters = models.CharField(max_length=20)
    epoch = models.CharField(max_length=20)
    dithersize = models.DecimalField(max_digits=5, decimal_places=1)
    seeing = models.DecimalField(max_digits=3, decimal_places=1)
    photometric = models.BooleanField()
    readtab = models.CharField(max_length=20)
    gain = models.CharField(max_length=20)
    obstype = models.CharField(max_length=20)
    mag = models.DecimalField(max_digits=5, decimal_places=2)
    mask = models.CharField(max_length=20)

    def __str__(self):
        """The string representation of this block."""
        return self.objID + ' at ' + self.start_time

    def fill_from_schedule(self, line, fldPar):
        """Give a CSV line from the schedule output and fill."""
        self.duration = line['duration']
        self.start_time = line['start_time']
        self.n_visits = line['n_visits_scheduled']
        self.objID = line['objid']
        self.end_time = line['end_time']

        # Find the object in fldPar that matches
        objpar = fldPar[fldPar['objid'] == self.objID]
        self.ra = objpar['ra']
        self.dec = objpar['dec']
        self.pmra = objpar['pmra']
        self.pmdec = objpar['pmdec']
        self.pa = objpar['pa']
        self.moon = objpar['moon']
        self.exptime = objpar['exptime']
        self.nexp = objpar['nexp']
        self.repeats = objpar['repeats']
        self.oneper = objpar['oneper']
        self.priority = objpar['priority']
        self.grism = objpar['grism']
        self.filters = objpar['filter']
        self.epoch = objpar['epoch']
        self.dithersize = objpar['dithersize']
        self.seeing = objpar['seeing']
        self.photometric = objpar['photometric']
        self.readtab = objpar['readtab']
        self.gain = objpar['gain']
        self.obstype = objpar['obstype']
        self.mag = objpar['mag']
        self.mask = objpar['mask']


class Schedule(models.Model):
    """Read in the schedule file and create a list of blocks."""

    time0 = time.time()
    schedule_file = 'schedule.csv'
    schedule_path = '/Users/rcool/MMTQueue/experiment/'
    sched = pd.read_csv(schedule_path+schedule_file)

    # Set the run to None and trigger to read it
    run_name = None
    objlist = []

    # Parse schedule lines
    for ii in range(len(sched)):
        if run_name is None:
            run_name = sched.loc[0, 'run']

        line = sched.loc[ii, :]
        # Read the fldpars
        fldpars = queue_tools.read_all_fld_files(run_name)

        obj = Observing_block()
        obj.fill_from_schedule(line, fldpars)
        objlist.append(obj)

    print(time.time()-time0)
