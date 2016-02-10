"""Create MMTEphem Class that contains the ephemeris for a given night.

Input : Local date
"""
import ephem
import math
import datetime
import matplotlib.pyplot as plt


def targetObservability(time, airmass):
    """Return 0/1 logical bit that defines if a target is observable.

    time = single time (datetime format)
    airmass = airmass target is at at time
    targetRise = time target rises
    targetSet = time target sets

    Returns : 0 if not observable 1 if observable
    """
    airmassCutoff = 2.5
    # When the target sets the "airmass" goes negative. That's bad
    if (airmass > 1.0) and (airmass < airmassCutoff):
        return 1
    else:
        return 0


def moonAge(time):
    """Return the age of the moon at the given time."""
    d1 = ephem.next_new_moon(time).datetime()
    d2 = ephem.previous_new_moon(time).datetime()

    time = time.datetime()
    print(d1-time)
    print(d2-time)
    if (d1-time) < (time-d2):
        return abs(d1-time)
    else:
        return(time-d2)


def moonPosition(datetime):
    """Return the position of the moon at given time."""
    j = ephem.Moon()
    j.compute(datetime)

    ra, dec = j.ra, j.dec
    return ra, dec


def mmtObserver():
    """Return an mmt XEphem observer for given night."""
    mmt = ephem.Observer()
    mmt.pressure = 0

    # Set the USNO definition of sunset (to match the almanac)
    # This will get overwritten when doing sunrise/sunset
    mmt.horizon = "-0:34"
    mmt.lat = "31:41:19.6"
    mmt.lon = "-110:53:04.4"
    mmt.elevation = 2600

    return mmt


def mmtTwilightDefinition():
    """Return the default value for twilight.

    Returns the degrees the sun must be below the horizon
    to define twilight.

    Nautical twilight : -12 degrees
    Astronomical twilight : -18 degrees

    For now, I'm setting to nautical, but this may need to be
    a switch for Binospec.
    """
    twilight = "-12"

    return twilight


def mmtAltAz(mmt, ra, dec, dateTime):
    """Return the Alt/AZ for a given RA/DEC at specified time."""
    star = ephem.FixedBody()
    star._ra = ra
    star._dec = dec
    star._epoch = ephem.J2000
    mmt.mmtObserver.date = dateTime
    star.compute(mmt.mmtObserver)

    return star.alt*180.0/math.pi, star.az*180/math.pi


def alt2airmass(alt):
    """Given an altitude angle (in degrees), calcualte the airmass."""
    zenithAngle = 90.0 - alt
    za_radians = zenithAngle/180.0*math.pi
    airmass = 1.0/math.cos(za_radians)

    return airmass


def airmassCurve(mmt, ra, dec):
    """Return the run of airmass as a function of time."""
    startTime = mmt.sunset.datetime()
    endTime = mmt.sunrise.datetime()

    timedelta = datetime.timedelta(minutes=5)

    timeArray = [startTime]
    while (timeArray[-1] < endTime):
        timeArray.append(timeArray[-1]+timedelta)

    airmass = []
    for time in timeArray:
        alt, az = mmtAltAz(mmt, ra, dec, time)
        airmass.append(alt2airmass(alt))

    return timeArray, airmass


class MMTEphem(object):
    """Return object with important times at MMT telescope."""

    def __init__(self, dateObs):
        """Initialize object."""
        # Append 19:00 to make sure the date references
        # noon MST.
        self.dateObs = dateObs + " 19:00"

        # Calculate sunrise and sunset
        mmt = mmtObserver()
        mmt.date = self.dateObs
        mmt.horizon = "-0:34"
        self.sunset = mmt.next_setting(ephem.Sun())
        self.sunrise = mmt.next_rising(ephem.Sun())

        # Calculate morning and evening twilight
        mmt.horizon = mmtTwilightDefinition()
        self.eveningTwilight = mmt.next_setting(ephem.Sun())
        self.morningTwilight = mmt.next_rising(ephem.Sun())

        # Calculate Moon rise and moon set
        self.moonrise = mmt.next_rising(ephem.Moon())
        self.moonset = mmt.next_setting(ephem.Moon())

        # Save the observer in case other calculations need it.
        self.mmtObserver = mmt


class ObjEphem(object):
    """Return an objects with observability parameters for a given target."""

    def __init__(self, ra, dec, dateObs):
        """Initialize the object."""
        # Get the MMTEphemeris
        mmtEphem = MMTEphem(dateObs)

        # Calculate when this target rises and sets
        self.dateObs = mmtEphem.dateObs
        mmtEphem.mmtObserver.date = mmtEphem.dateObs

        # Set a hard limit of observablity to be 2.5 airmasses
        timeArray, airmass = airmassCurve(mmtEphem, ra, dec)

        observability = []
        for ii in range(len(airmass)):
            obs = targetObservability(timeArray[ii], airmass[ii])
            observability.append(obs)

        self.observable = observability
        self.ra = ra
        self.dec = dec
        self.time = timeArray
        self.airmass = airmass


if __name__ == "__main__":
    mmt = ObjEphem("8:00:00", "30:00:00", "2016/02/09")

    date = ephem.Date("2016/02/10 19:40")
    ra, dec = moonPosition(date)
    print(ra)
    print(dec)
