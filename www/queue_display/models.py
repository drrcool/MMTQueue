from django.db import models
from django.utils import timezone
import queue_tools


class Observing_block(models.Model):
    """Create the class that holds a the parameters for one field."""

    # Parameters for this observation
    duration = models.IntegerField()
    n_visits = models.IntegerField()
    objID = models.CharField(max_length=50)
    start_time = models.CharField(max_length=50)

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
        _, self.duration, self.n_visits, self.objID, _, self.start_time = \
            line.strip().split()

        # Find the object in fldPar that matches
        objpar = fldPar[fldPar['objid'] == self.objID]
        self.ra =






class schedule(models.Models):
    """Read in the schedule file and create a list of blocks."""

    schedule_file = 'schedule.csv'
    schedule_path = '/Users/rcool/MMTQueue/experimenta/'


    f = open(schedule_path + schedule_file, 'r')
    f.readline()  # Read in the header

    # Set the run to None and trigger to read it
    run_name = None

    # Parse schedule lines
    for line in f.readlines():
        if run_name = None:
            run_name = line.strip().split()[4]

            # Read the fldpars
            fldpars = queue_tools.read_all_fld_files(run_name)

        obj = Observing_block()
        obj.fill_from_schedule(line, fldpar)
