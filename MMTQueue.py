import ipdb
""" Queue software for observations at the MMT observatory.

This version is an overhaul of the first system to make it more object
oriented and remove a lot of useless recalculation.

See the ancillary README and documentation for file formats and instructions.
"""

from os import walk
import ephem as pyEphem
import datetime
import numpy as np
import pandas as pd
import os
import sys
import re
from random import randint
import matplotlib.pyplot as plt


class Target:

    """Define the Target class."""

    def __init__(self, ra, dec, position_angle):
        """Intialize the target object."""

        self.ra = ra
        self.dec = dec
        self.position_angle = float(position_angle)
        self.ephem = pyEphem.FixedBody()
        self.ephem._ra = self.ra
        self.ephem._dec = self.dec
        self.ephem._epoch = pyEphem.J2000
        self.MMT = MMT()

    def isObservable(self, timestamp):
        """Is the target observable at the given datetime."""

        # Check the airmass. Checking for < 1.0 is for numerical issues
        airmass_cutoff = 1.8
        airmass = self.airmass(timestamp)  # Using it twice, only calc once
        if (airmass < 1.0) or (airmass > airmass_cutoff):
            return 0

        # Check the rotator limit
        # TODO Add longslit check for +180 and 0
        rotator_limits = [-180, 164]
        rotAngle = self.rotator_angle(timestamp)
        if rotAngle < rotator_limits[0] or \
           rotAngle > rotator_limits[1]:
            return 0

        return 1

    def AltAz(self, timestamp):
        """Calculate the Alt/Az position of a target.

        Inputs:
            ra : Right Ascension
            dec : Declination
            timestamp : timestamp timestamp
            observatory : pyEphem.Observer object
        """
        self.MMT.MMTEphem.date = timestamp
        self.ephem.compute(self.MMT.MMTEphem)

        # Return the alt and az calculated
        return self.ephem.alt * 180.0 / np.pi, \
            self.ephem.az * 180 / np.pi

    def airmass(self, timestamp):
        """Given a targets position and time, calculate the airmass."""

        # Calculate the alt az
        alt, az = self.AltAz(timestamp)
        zenith_angle = 90.0 - alt
        za_radians = zenith_angle / 180.0 * np.pi
        return 1.0 / np.cos(za_radians)

    def parallactic_angle(self, timestamp):
        """Calculate the parallactic angle at a given observation time."""
        self.MMT.MMTEphem.date = timestamp
        self.ephem.compute(self.MMT.MMTEphem)
        return self.ephem.parallactic_angle()

    def rotator_angle(self, timestamp):
        """Calculate the rotator angle needed for the PA and time."""
        parAngle = self.parallactic_angle(timestamp) * 180.0 / np.pi

        rotAngle = parAngle - self.position_angle

        # Account for the fact that 355 and -5 are the same angle
        if rotAngle < -180:
            rotAngle += 360.0

        return rotAngle

    def separation(self, Target):
        """Calculate the distance between this target and another"""
        # Convert coordinates to radians
        return AngSep(self.ephem._ra, self.ephem._dec, Target.ra, Target.dec)

    def lunar_distance(self, timestamp):
        """Calculate the distance to the moon."""
        moon_ra, moon_dec = self.MMT.moonPosition(timestamp)
        return AngSep(self.ephem._ra, self.ephem._dec, moon_ra, moon_dec)


class MMT:
    """Define an object to hold details about the MMT."""

    def __init__(self):
        """Initialize the MMT object."""
        self.MMTEphem = pyEphem.Observer()
        self.MMTEphem.pressure = 0
        self.MMTEphem.lat = "31:41:19.6"
        self.MMTEphem.lon = "-110:53:04.4"
        self.MMTEphem.elevation = 2600

        # This holds the definition of twilight
        self.twilight_horizon = "-12"
        self.MMTEphem.horizon = self.twilight_horizon

    def moon_age(self, timestamp):
        """Return the age of the moon at the given time."""
        d1 = pyEphem.next_new_moon(timestamp).datetime()
        d2 = pyEphem.previous_new_moon(timestamp).datetime()

        # Check format of provided time. If it wasn't a timestamp, make it one.
        if isinstance(timestamp, datetime.datetime) is False:
            timestamp = timestamp.datetime()

        # Find the total time since new moon and then convert to days
        if (d1 - timestamp) < (timestamp - d2):
            return (timestamp - d1).total_seconds() / 3600. / 24.
        else:
            return (timestamp - d2).total_seconds() / 3600. / 24.

    def moonPosition(self, timestamp):
        """Return the position of the moon at the given time."""
        # Construct the lunar ephemeris and find ra an dec given time
        j = pyEphem.Moon()
        j.compute(timestamp)
        return j.ra, j.dec

    def string_timestamp_to_noonMST(self, timestamp):
        """Convert a string timestamp to reference noon MST"""
        return timestamp.split()[0] + ' 19:00'

    def evening_twilight(self, timestamp):
        """Calculate the next evening twilight."""
        # Check type of timestamp
        if type(timestamp) == str:
            timestamp = self.string_timestamp_to_noonMST(timestamp)

        # Set the calculation to the specified date
        self.MMTEphem.date = timestamp
        return self.MMTEphem.next_setting(pyEphem.Sun()).datetime()

    def morning_twilight(self, timestamp):
        """Calculate the next morning twilight."""
        # Check the type of timestamp
        if type(timestamp) == str:
            timestamp = self.string_timestamp_to_noonMST(timestamp)

        # Set the calculation to be done at the specified date
        self.MMTEphem.date = timestamp
        return self.MMTEphem.next_rising(pyEphem.Sun()).datetime()

    def is_moon_up(self, timestamp):
        """Calculate if the moon is up at the given timestamp."""

        # Calculate the last and next moonrise
        self.MMTEphem.date = timestamp
        self.MMTEphem.horizon = '-0:34'
        prev_moonset = self.MMTEphem.previous_setting(pyEphem.Moon())
        prev_moonrise = self.MMTEphem.previous_rising(pyEphem.Moon())
        self.MMTEphem.horizon = self.twilight_horizon  # Restore default

        if prev_moonset > prev_moonrise:
            return 0
        else:
            return 1


def get_cmap_tuple(ii):
    """Create colormap."""
    # rgbt = plt.cm.Pastel1(ii)
    rgbt = plt.cm.Set3(ii)
    out_cmap = (int(rgbt[0]*255),
                int(rgbt[1]*255),
                int(rgbt[2]*255))
    return out_cmap


def tuple_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def hms2dec(string):
    """Convert a string hms to a float value.

    Inputs : string in formation (+/-)HH:MM:SS
    Output : float value.

    Note, ra will return decimal *hours*, so the conversion factor
    of 15 should be applied after using this function.
    """
    # We need to treat the sign in a special way (-3:30 isn't -2.5, it's -3.5)
    if string[0] == '-':
        string = string[1:]
        sign = -1.0
    else:
        sign = 1.0

    h, m, s = map(float, string.split(':'))

    return sign*(h+m/60.0+s/3600.0)


def AngSep(ra1, dec1, ra2, dec2):
    """Calcualte the angular separation between two points on the sky.

    All inputs are Ephem Angles (so decimal radians)

    Output is in decimal degrees.
    """

    y = np.cos(dec1) * np.cos(dec2)
    z = np.sin(dec1) * np.sin(dec2)
    x = np.cos(ra1 - ra2)

    rad = np.arccos(z+y*x)

    # For small separations, use Euclidean distance
    if (rad < 0.000004848):
        sep = np.sqrt((np.cos(dec1)*(ra1-ra2))**2 + (dec1-dec2)**2)
    else:
        sep = rad

    return sep*180.0 / np.pi


def mmirs_overheads(fldPar):
    """Return the expected overhead time (in seconds) for the observation.

    For now, we assume the overhead is constant per configuration. It's
    quite possible that this assumption will need to be re-evaluated (for
    example to account for checking alignment every few hours on longer
    exposures).

    These numbers are also fairly pessimistic. As observer efficiency increases
    we can decrease these to match what we're seeing in operations.
    """

    obstype = fldPar['obstype'].values[0]

    if obstype == 'mask':
        return 2700.0
    elif obstype == 'longslit':
        return 1800.0
    elif obstype == 'imaging':
        return 120.0
    else:
        # If none of these were given, we need to throw an error
        raise AssertionError("Unexpected value of OBSTYPE in " +
                             fldPar['objid'])


def parse_mask_position_angle(mask, runname):

    """Parse the .msk file for a mask to get it's position angle.

    For the March 2016 run, the position angle was not written to the FLD
    files for any mask observations. This makes checking the rotator limits
    impossible. This code parses the .msk file to get the needed position
    angle to add to the fldPar.
    """
    # Read the mask file
    maskfile = 'mmirs_catalogs/' + runname + '/' + mask + '.msk'
    f = open(maskfile, 'r')

    # Check each line in the maskfile and find the line that starts with 'pa'
    mask_position_angle = False
    for line in f.readlines():
        sline = line.strip().split()
        if len(sline) > 1 and sline[0] == 'pa':
            mask_position_angle = float(sline[1])
    f.close()

    return mask_position_angle


def read_allocated_time(runname):
    """Read the allocated time input file for the given run.

    Inputs:
        runname : name given to dates belonging to one run. If multiple runs
            are being summed to give one queue session, we all should have the
            same runname.  This means we could separate by trimester (but
            can also do finer divisions if needed).
    """
    # TODO: Document how to updated the allocated_time file
    filename = "AllocatedTime.dat"
    f = open(filename, 'r')

    # Initialize the dictionary to hold this time
    allocated_time = {}
    for line in f.readlines():
        if line[0] == '#':
            # Comment line, so skip
            continue

        date, PI, runID = line.strip().split()

        # The run for this night doesn't match the specified run, so skip
        if runID != runname:
            continue

        mmt = MMT()
        date = pyEphem.date(date).datetime()

        # Calculate the night length and covert to hours
        night_length = mmt.morning_twilight(date) - mmt.evening_twilight(date)
        night_length = night_length.total_seconds() / 3600.

        # Update the allocated time dictionary
        if PI in allocated_time:
            allocated_time[PI] += night_length
        else:
            allocated_time[PI] = night_length

    f.close()
    return allocated_time


def read_single_fld_file(filename, runname):
    """Read a FLD file and return a dictionary with the contained data."""
    # Initialize the output dictionary
    obspars = {}

    f = open(filename, 'r')
    # Get the PI:
    _, PI_name = f.readline().strip().split()
    obspars['PI'] = PI_name

    # Get the program ID
    _, prog_ID = f.readline().strip().split()
    obspars['progID'] = prog_ID

    # Add the filename for bookkeeping
    obspars['fldfile'] = filename

    # Parse the remaining column names
    keywords = f.readline().strip().split()
    f.readline()  # Remove the line of "------"
    values = f.readline().strip().split()

    # Fill the dictionary
    for key, val in zip(keywords, values):
        obspars[key] = val

    # Fill in the position angle if this is a mask
    if obspars['obstype'] == 'mask':
        obspars['pa'] = parse_mask_position_angle(
            obspars['mask'], runname)
    f.close()

    # Include the object for target
    obspars['ephem'] = Target(obspars['ra'], obspars['dec'],
                              obspars['pa'])
    return obspars


def read_all_fld_files(runname):
    """Read all of the fld files for a run and output a dataframe."""

    # Set the path
    path = 'mmirs_catalogs/' + runname

    # Get the list of files that end in .fld in the specified path
    filelist = []
    for (dirpath, dirnames, filenames) in walk(path):
        filelist.extend([dirpath+'/' + f
                         for f in filenames if f[-4:] == '.fld'])

    # Create a list of each of the dictionaries from the individual files
    fld_list = []
    for file in filelist:
        fld = read_single_fld_file(file, runname)
        fld_list.append(fld)

    # Convert the list to a dataframe and output
    return pd.DataFrame(fld_list)


def create_done_mask(obspars, runname):
    """Return the donepar dataframe for this run.

    If this is the first time the queue has been run,
    we opt to use a blank file
    otherwise, read the existing donefile.
    """

    blank_donepar = create_blank_done_mask(obspars, runname)
    donefile = 'mmirs_catalogs/' + runname + '/donefile.dat'
    if os.path.isfile(donefile):
        print("Found Existing Set of Finished Observations:")
        print("THIS IS NOT IMPLEMENTED, STOP NOW!!!!!")
    else:
        print("No existing donefile.csv for run %s, initializing..." % runname)
        # Write this for future iterations

    return blank_donepar


def create_blank_done_mask(obspars, runname):
    """Create a blank structure to hold information about progress on fields.

    I'm not a big fan of this format. So this could change.

    Inputs:
        obspars: dataframe containing each of the requested observations
                 for this run
        runname: signifier for the given run
                 Used when determining allocated time
    """
    # Read the allocated time
    allocated_time = read_allocated_time(runname)

    # Create the blank format
    blank_dict = {}
    blank_dict['complete'] = 0
    blank_dict['done_visits'] = 0
    blank_dict['allocated_time'] = 0.0
    blank_dict['time_for_PI'] = 0.0
    blank_dict['previous_weight'] = 1.0
    blank_dict['current_weight'] = 0.0
    blank_dict['PI'] = ''
    blank_dict['objid'] = ''

    dict_list = []
    for ii in range(len(obspars)):
        copy_dict = blank_dict.copy()
        copy_dict['objid'] = obspars.loc[ii, 'objid']
        copy_dict['PI'] = obspars.loc[ii, "PI"]
        if copy_dict['PI'] in allocated_time:
            copy_dict['allocated_time'] = allocated_time[copy_dict['PI']]
        else:
            # To avoid divide by zero errors, set anyone without
            # allocated time to 1/1000 of an hour.
            copy_dict['allocated_time'] = 0.001
        dict_list.append(copy_dict)

    return pd.DataFrame(dict_list)


def moon_flag(fldpar, start_time, end_time, moon_up_time, moon_up_array):
    """Calculate a flag for a given field based on the lunar conditions.

    Inputs:
        fldpar -- field parameter entry for a single field
        startTime/endTime -- starting and ending time for block

    Output:
        0/1 flag marking if the field suffers from any lunar issues.
    """
    # Parse the target ephemeris to get lunar position and age
    # moon_up_at_start = fldpar.iloc[0]['ephem'].MMT.is_moon_up(start_time)
    # moon_up_at_end = fldpar.iloc[0]['ephem'].MMT.is_moon_up(end_time)
    mintime_to_start = min([abs(x - start_time) for x in moon_up_time])
    mintime_to_end = min([abs(x - end_time) for x in moon_up_time])

    for ii, itime in enumerate(moon_up_time):
        if abs(itime-start_time) == mintime_to_start:
            moon_up_at_start = moon_up_array[ii]
        if abs(itime-end_time) == mintime_to_end:
            moon_up_at_end = moon_up_array[ii]

    # This is set if either start or end is set
    moon_up = max(moon_up_at_end, moon_up_at_start)

    moon_age = fldpar.iloc[0]['ephem'].MMT.moon_age(start_time)
    lunar_distance = fldpar.iloc[0]['ephem'].lunar_distance(start_time)
    moon_requirement = fldpar.iloc[0]['moon']

    # Check the brightness compared to requirement
    if moon_requirement == 'bright' or moon_up == 0:
        # Either the moon is down or bright was requested, so we're good
        illum_flag = 1
    elif (moon_requirement == 'grey') and \
         (abs(moon_age) < 9) and \
         (lunar_distance < 90):
        illum_flag = 1
    elif (moon_requirement == 'dark') and \
         (abs(moon_age) < 4.5) and \
         (lunar_distance >= 90):
        illum_flag = 1
    else:
        illum_flag = 0

    # Finally, check to be sure we aren't just plain too close
    if lunar_distance < 10:
        illum_flag = 0

    return illum_flag


def calculate_observation_duration(exp_per_visit, n_repeats, fldpar):
    """Calcualte the duration of an observation.

    This happens a number of places, so removed it to not repeatcode.
    """
    total_target_time = exp_per_visit * n_repeats + \
        mmirs_overheads(fldpar)
    return datetime.timedelta(seconds=total_target_time)


def does_field_fit(fldpar, start_time, donepar):
    """Calculate a weight if a field fits.

    Here, the weights are either 0 or 1 (does the observation fit).  We also
    calculate the fraction of the observation that fits in this window.

    Inputs:
        fldpar -- parameters for the field being considered
        start_time -- time to begin the search
        donepar -- tracker of all the finished fields (contains weights)
    """

    # This only works if fldpar has one entry
    if len(fldpar) > 1:
        raise AssertionError("There are more tahn 1 fields with name %s" %
                             fldpar['objid'])

    # First, get the end of the night time

    night_end = fldpar.iloc[0]['ephem'].MMT.morning_twilight(start_time)
    time_remaining = (night_end - start_time).total_seconds()  # seconds

    # Is the target observable at the start_time?
    if fldpar.iloc[0]['ephem'].isObservable(start_time) == 0:
        return 0.0, start_time, 0

    # Get the exposure parameters for this object
    # Exposure time is stored in minutes
    exptime = float(fldpar.iloc[0]['exptime'])*60.0
    n_repeats = float(fldpar.iloc[0]['repeats']) - \
        float(donepar.iloc[0]['done_visits'])
    nexp_per_visit = float(fldpar.iloc[0]['nexp'])

    # Calculate the number of visits that fit (keep it whole -- partial visits
    # aren't useful).
    exptime_per_visit = exptime * nexp_per_visit
    possible_visits = np.floor(
        (time_remaining - mmirs_overheads(fldpar)) / exptime_per_visit)

    # Can we fit all the requested repeats in?
    if possible_visits >= n_repeats:

        duration = \
            calculate_observation_duration(exptime_per_visit,
                                           n_repeats, fldpar)

        # Is the target still observable at the end_time?
        if fldpar.iloc[0]['ephem'].isObservable(start_time+duration) == 1:
            return 1.0, start_time + duration, n_repeats
        else:
            possible_visits = n_repeats - 1

    # At this point, we can't fit the whole observation, so let's figure out
    # how many repeats we can fit in
    nrepeats_observable = possible_visits

    while nrepeats_observable >= 1:
        duration = calculate_observation_duration(exptime_per_visit,
                                                  nrepeats_observable, fldpar)
        if fldpar.iloc[0]['ephem'].isObservable(start_time+duration) == 1:
            return 1.0, start_time + duration, nrepeats_observable
        else:
            nrepeats_observable -= 1  # Decrement and try again

    # If we get here, we didn't find a combo that works, so return 0
    return 0, start_time, 0


def calc_same_target_flag(fldpar, prev_target):
    """Return a weight to upweight a target that we are already pointing at.

     This upweights the chances of observing a field we're already looking at
     rather than paying the overhead fee multiple times.
     """
    dist = fldpar.separation(prev_target)
    dist_weight = 1000.0  # Increase chances of observing nearby target
    if dist < 10./3600.:
        return dist_weight
    else:
        return 1.0


def calc_single_weight(fldpar, obj_donepar, start_time,
                       moon_up_time, moon_up_array, prev_target=None):
    """Calcualte the weight for a single object."""
    # Intialize a dictionary to hold our weights
    obs_weight = {}
    obs_weight['objid'] = fldpar.iloc[0]['objid']

    # Determine if this field fits
    fit_weight, end_time, obs_visits = \
        does_field_fit(fldpar, start_time, obj_donepar)

    # Check the moon conditions
    moon_weight = moon_flag(fldpar, start_time, end_time, moon_up_time,
                            moon_up_array)

    # Check to see if the previous target we looked at was in this field
    if prev_target is not None:
        dist_weight = calc_same_target_flag(prev_target)
    else:
        dist_weight = 1.0

    obs_weight['end_time'] = end_time
    obs_weight['duration'] = (end_time-start_time).total_seconds()
    obs_weight['n_visits_scheduled'] = obs_visits
    obs_weight['target'] = fldpar.iloc[0]['ephem']

    # TODO : Implement priority scheduling both between programs and in PI
    # Calculate the TAC weight
    tac_weight = 1.0 * obj_donepar.iloc[0]['time_for_PI'] / \
        obj_donepar.iloc[0]['allocated_time']
    # Don't allow this to be exactly zero (divide by zero errors)
    if tac_weight <= 0:
        tac_weight = 1e-5

    # Extract the previlous weight
    prev_weight = obj_donepar.iloc[0]['previous_weight']

    # Calculate the final weight
    total_weight = dist_weight / tac_weight / prev_weight * \
        fit_weight * moon_weight
    obs_weight['total_weight'] = total_weight

    return obs_weight


def calc_field_weights(obspar, donepar, start_time, moon_up_time,
                       moon_up_array, prev_target=None):
    """Determine the weight for every object at this start time.

    Using this weight we will select the best target to observe here
    """
    # Loop through all of the fields and calculate a weight
    weight_list = []
    append = weight_list.append

    for objID in obspar['objid']:

        fldpar = obspar[obspar['objid'] == objID]

        # Have we already observed this target?
        obj_donepar = donepar[donepar['objid'] == objID]
        if obj_donepar.iloc[0]['complete'] == 1:
            continue

        obs_weight = calc_single_weight(fldpar,
                                        obj_donepar,
                                        start_time,
                                        moon_up_time,
                                        moon_up_array,
                                        prev_target=prev_target)
        append(obs_weight)

    return pd.DataFrame(weight_list)


def UpdateRow(obspar, donepar, start_time, moon_up_time, moon_up,
              prev_target=None, runname="unspecified"):
    """Coordinate the weight calculation and target selection."""

    obs_weight = calc_field_weights(obspar, donepar,
                                    start_time, moon_up_time, moon_up,
                                    prev_target=prev_target)

    # Find were the weight is the maximum

    max_weight = max(obs_weight['total_weight'])
    # If the largest weight is 0, no target was selected

    if max_weight == 0:
        return None, None

    max_index = []  # Will contain locations of max weights

    for ii, val in enumerate(obs_weight['total_weight']):
        if val == max_weight:
            max_index.append(ii)

    # Do the tie breaking.
    # TODO: Add tie breaking based on priority or previous observations
    selected_index = max_index[randint(0, len(max_index)-1)]  # randomly choose
    selected_object = obs_weight.iloc[selected_index]

    # Create the final schedule entry
    schedule = {}   # Initialize
    schedule['n_visits_scheduled'] = selected_object['n_visits_scheduled']
    schedule['duration'] = selected_object['duration']
    schedule['objid'] = selected_object['objid']
    schedule['end_time'] = selected_object['end_time']
    schedule['start_time'] = start_time
    schedule['run'] = runname  # This is not elegant
    prev_target = selected_object['target']

    # Update the donepar
    index = [i for i, x in
             enumerate(donepar['objid'] == schedule['objid']) if x]
    index = index[0]
    PI = donepar.iloc[index]['PI']
    for ii in range(len(obspar)):
        if donepar.loc[ii, 'PI'] == PI:
            donepar.loc[ii, 'time_for_PI'] += \
                selected_object['duration'] / 3600.0
            donepar.loc[ii, 'current_weight'] = \
                donepar.loc[ii, 'time_for_PI'] / \
                donepar.loc[ii, 'allocated_time']
    donepar.loc[index, 'done_visits'] += \
        selected_object['n_visits_scheduled']

    # Check to see if this object is now done
    requested = int(obspar[obspar['objid'] == schedule['objid']]['repeats'])
    if donepar.loc[index, 'done_visits'] >= requested:
        donepar.loc[index, 'complete'] = 1

    return schedule, prev_target


def obsOneNight(obspar, donepar, date, runname):
    """Fully schedule one night."""
    mmt = MMT()

    # Start at twilight and add observations. Each time, increment the
    # current time by the previous duration. If nothing add 20 minutes and try
    # again.

    # Cache if the moon is up
    current_time = mmt.evening_twilight(date)
    tdelta = datetime.timedelta(minutes=20)
    time_array = [current_time]
    moon_up = [mmt.is_moon_up(current_time)]
    while (time_array[-1] < mmt.morning_twilight(date)):
        time_array.append(time_array[-1] + tdelta)
        moon_up.append(mmt.is_moon_up(time_array[-1]))

    schedule = []
    # Set some flags
    all_done = False
    prev_target = None

    while (current_time < mmt.morning_twilight(date)) and (all_done is False):

        new_sched, new_target = \
            UpdateRow(obspar, donepar, current_time, time_array, moon_up,
                      prev_target=prev_target, runname=runname)

        # Was anything observed?
        if new_sched is None:
            current_time += datetime.timedelta(minutes=20)
        else:
            # Append new entry to schedule
            schedule.append(new_sched)

            current_time += \
                datetime.timedelta(seconds=new_sched['duration'])
            if min(donepar['complete']) == 1:
                all_done = True
    return schedule


def obsAllNights(obspar, donepar, all_dates, iter_number, runname):
    """Iterate through nights and schedule."""
    full_schedule = []
    for date in all_dates:
        sys.stdout.write("\r Working on Date %s of iteration %d" %
                         (date, iter_number+1))
        schedule = obsOneNight(obspar, donepar, date, runname)
        for line in schedule:
            full_schedule.append(line)
    return pd.DataFrame(full_schedule)


def read_fitdates():
    """Read the list of dates to fit. This is fitdates.dat.

    This definitely needs to change to become more general and
    requiring less hardcoding.
    """
    date_file = 'fitdates.dat'
    f = open(date_file, 'r')
    all_dates = []
    for line in f.readlines():
        if line[0] != '#' and line.strip() != '':
            all_dates.append(line.strip())
    f.close()

    return all_dates


def schedule_to_json(schedule, obspars, outfile='schedule.json'):
    """Create a JSON file containing all the information to be parsed.

    This is in the format of the JavaScript FullCalendar module.

    The API for Event Objects is found at
    http://fullcalendar.io/docs/event_data/Event_Object/
    """
    json_schedule_list = []

    # Get a color Table
    color = [tuple_to_hex(get_cmap_tuple(x/len(obspars)))
             for x in range(len(obspars))]
    obspars['color'] = color

    for ii in range(len(schedule)):
        entry = schedule.loc[ii, :]

        obs = obspars[obspars['objid'] == entry['objid']]

        # Create a blank dictionary to hold the output
        json_template = {}
        json_template['objid'] = entry['objid']
        # Times need to be in YYYY-MM-DDTHH:MM:SS

        json_template['backgroundColor'] = obs.iloc[0]['color']
        json_template['borderColor'] = obs.iloc[0]['color']
        json_template['start'] = \
            str(entry['start_time']).replace(' ', 'T')[0:19]
        json_template['end'] = \
            str(entry['end_time']).replace(' ', 'T')[0:19]
        json_template['url'] = 'fields/' + entry['objid']

        # Append the keys needed for the tooltip
        copy_columns = ['PI', 'dec', 'dithersize', 'exptime',
                        'filter', 'gain', 'grism',
                        'mag', 'mask', 'moon', 'obstype',
                        'pa', 'photometric', 'ra', 'readtab',
                        'seeing', 'repeats']
        json_template['n_visits_scheduled'] = entry['n_visits_scheduled']
        for col in copy_columns:
            json_template[col] = obs.iloc[0][col]

        # Parse the field title
        reString = "^[a-zA-Z]+-[A-Za-z0-9]+_(.*)$"
        m = re.search(reString, entry['objid'])

        json_template['title'] = m.group(1)
        json_schedule_list.append(json_template)

    outframe = pd.DataFrame.from_dict(json_schedule_list)
    outframe.to_json(outfile, orient='records')


def main(args):
    """Main processing function."""
    if len(args) < 2:
        raise Exception("Must specify a run name")
    else:
        runname = args[1]

    # Read in the objects for this run
    obspars = read_all_fld_files(runname)

    # Create a blank donepar
    orig_donepar = create_done_mask(obspars, runname)
    donepar = orig_donepar.copy()   # This is the working copy

    # Read in the dates to be fit
    all_dates = read_fitdates()

    # Iterate
    finished_flag = False
    number_of_iterations = 5
    iter_number = 0
    while (iter_number < number_of_iterations) and finished_flag is False:
        schedule = obsAllNights(obspars, donepar, all_dates,
                                iter_number, runname)

        # Create a copy of the donepar
        new_donepar = donepar.copy()

        # Restore the original values from the donefiles
        new_donepar.loc[:, 'complete'] = orig_donepar['complete']
        new_donepar.loc[:, 'time_for_PI'] = orig_donepar['time_for_PI']
        new_donepar.loc[:, 'done_visits'] = orig_donepar['done_visits']

        # Check to see if all fields are done
        done_dict = {}
        for ii in range(len(donepar)):
            pi = donepar.loc[ii, 'PI']
            if pi not in done_dict:
                completed = donepar[donepar['PI'] == pi]['complete'].values
                if min(completed) == 1:
                    done_dict[pi] = True
                else:
                    done_dict[pi] = False
            new_donepar.loc[ii, 'previous_weight'] = \
                donepar.loc[ii, 'current_weight']

        # Check to see if we're done
        if min(done_dict.values()) == 1:
            finished_flag = True
        new_donepar.loc[:, 'current_weight'] = 0.0
        donepar = new_donepar
        iter_number += 1

    # Write out the schedule
    outfile = 'schedule.csv'
    schedule.to_csv(outfile, index_label='index')

    # Now, parse the schedule to JSON for plotting
    jsonfile = 'schedule.json'
    schedule_to_json(schedule, obspars, outfile=jsonfile)


if __name__ == "__main__":
    main(sys.argv)
