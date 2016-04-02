from django.db import models
from django.conf import settings
import os

class Mask(models.Model):
    """Create the mask model.

    Note that most of this *wont* be supplied by the user and we'll need to
    write code that parses the header of the mask file. This is probably good
    as we want to do something that checks the quality of the mask.
    """
    mask = models.CharField(max_length=50)
    ra = models.CharField(max_length=20)
    dec = models.CharField(max_length=20)
    position_angle = models.DecimalField(max_digits=7, decimal_places=3)
    PI = models.CharField(max_length=20)
    trimester = models.CharField(max_length=50)

    def __str__(self):
        return self.mask

class Target(models.Model):
    """Create an observation request.

    Targets:
    --------

    I think we don't care if the targets get repeated.  We need to figure out
    how to deal with this when we generate an observing catalog, but that's
    step N
    """
    object_name = models.CharField(max_length=20)
    PI = models.CharField(max_length=20)
    run_name = models.CharField(max_length=10)
    program = models.CharField(max_length=20)
    ra = models.CharField(max_length=20)
    dec = models.CharField(max_length=20)
    ra_propermotion = models.DecimalField(max_digits=5, decimal_places=2)
    dec_propermotion = models.DecimalField(max_digits=5, decimal_places=2)
    position_angle = models.DecimalField(max_digits=7, decimal_places=3)
    moon_phase = models.CharField(max_length=6)
    exposure_time = models.IntegerField()
    positions_in_dither_pattern = models.IntegerField()
    number_of_dither_repeats = models.IntegerField()
    oneper = models.BooleanField()
    priority = models.IntegerField()
    grism = models.CharField(max_length=20)
    filters = models.CharField(max_length=20)
    required_seeing = models.DecimalField(max_digits=3, decimal_places=1)
    need_photometric = models.BooleanField()
    readtab = models.CharField(max_length=20)
    gain = models.CharField(max_length=20)
    obstype = models.CharField(max_length=20)
    mag = models.DecimalField(max_digits=5, decimal_places=2)
    mask = models.CharField(max_length=20)


class Semester(models.Model):
    # Create a class to contain the semester info

    rootdir = settings.ROOT_PATH
    run_config_dir = os.path.join(rootdir, 'run_config')
    pi_list_dir = os.path.join(run_config_dir, 'pi_per_run')

    def get_trimester(self):

        current_trimester_file = os.path.join(self.run_config_dir, 'current_semester')
        f = open(current_trimester_file, 'r')
        trimester = f.readline().strip()
        f.close()
        return trimester

    def get_allocated_pi(self):
        # Create a dictionary of all the PIs given time for each semester


        trimester = self.get_trimester()
        trimester_list_file = os.path.join(self.pi_list_dir, trimester)
        f = open(trimester_list_file, 'r')
        pi_list = []
        prog_list = []
        for line in f.readlines():
            if line.strip() != '':
                split = line.strip().split()
                pi_list.append(split[1])
                prog_list.append(split[0])
        f.close()

        return [(i, x) for i, x in enumerate(sorted(pi_list))]
