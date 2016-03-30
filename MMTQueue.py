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


class Target:

    """Define the Target class."""

    def __init__(self, ra, dec, position_angle):
        """Intialize the target object."""

        self.ra = ra
        self.dec = dec
        self.position_angle = position_angle
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
        alt, az = self.AltAz(timestamp, self.MMT.MMTEphem)
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
        return AngSep(self.ra, self.dec, Target.ra, Target.dec)

    def lunar_distance(self, timestamp):
        """Calculate the distance to the moon."""
        moon_ra, moon_dec = self.MMT.moonPosition(timestamp)
        return AngSep(self.ra, self.dec, moon_ra, moon_dec)


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

        # Calculate the last and next moonrise at this time.
        self.MMTEphem.date = timestamp
        self.MMTEphem.horizon = '-0:34'
        prev_moonset = self.MMTEphem.previous_setting(pyEphem.Moon())
        prev_moonset = prev_moonset.datetime()
        prev_moonrise = self.MMTEphem.previous_rising(pyEphem.Moon())
        prev_moonrise = prev_moonrise.datetime()
        self.MMTEphem.horizon = self.twilight_horizon  # Restore default

        if prev_moonset > prev_moonrise:
            return 0
        else:
            return 1


if __name__ == "__main__":
    filename = 'mmirs_catalogs/2016a/UAO-S143-DStark/UAO-S143_gdnel31_J.fld'
    runname = '2016a'
    fldpar = read_all_fld_files('2016a')
    print(fldpar.head())
