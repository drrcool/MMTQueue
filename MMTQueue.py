"""Queue software for observations at the MMT Observatory.

This is the main module for the MMT queue sytem.  A number of
ancillary files will be utilized as well to keep some of the
software modular.

$ MMTQueue.py startDate endDate


Input : startDate -- first night in the observing block
        endDate -- last night in the observing block
Output : Schedule for each night starting with startDate and
         ending with endDate (inclusive) assuming good conditions.

More to come as code is fleshed out.
"""

# Imports
import sys
import os.path
import numpy as np
import queueTools
import MMTEphem
import math
import datetime
import ephem as pyEphem


def angSep(ra1, dec1, ra2, dec2):
    """Calculate the angular separation between two points on the sky.

    ALL INPUTS ARE Ephem Angles (so decimal ***RADIANS***).

    OUTPUT IS IN DECIMAL DEGREES
    """
    y = math.cos(dec1) * math.cos(dec2)
    z = math.sin(dec1) * math.sin(dec2)
    x = math.cos(ra1-ra2)

    rad = np.arccos(z+y*x)

    if (rad < 0.000004848):
        sep = math.sqrt((math.cos(dec1)*(ra1-ra2))**2 +
                        (dec1-dec2)**2)
    else:
        sep = rad

    return sep * 180 / math.pi


def mmirsOverhead(fld):
    """Return the expected overhead time (in seconds) for MMIRS observation."""

    obstype = fld['obstype'].values[0]

    if obstype == 'mask':
        return 600.0
    elif obstype == 'longslit':
        return 300.0
    elif obstype == 'imaging':
        return 120.0
    else:
        # Make sure we return one of the modes, otherwise throw a fit
        raise AssertionError("Unexpected value of OBSTYPE in " + fld["objid"])

def lunarDistance(fld, startTime, endTime):
    """Calculate the distance to the moon at starttime and end time."""
    # Position of the target cast into pyeEphem Angle
    tarRa = pyEphem.hours(fld['ra'].values[0])
    tarDec = pyEphem.degrees(fld['dec'].values[0])

    # Moon Position
    moonRa1, moonDec1 = MMTEphem.moonPosition(startTime)
    moonRa2, moonDec2 = MMTEphem.moonPosition(endTime)

    dist1 = angSep(moonRa1, moonDec1, tarRa, tarDec)
    dist2 = angSep(moonRa2, moonDec2, tarRa, tarDec)

    return (dist1+dist2)/2.0



def moonUpDuringObs(startTime, endTime, mmt):
    """Return a 0/1 if moon is down/up at any point during observation."""
    # Is the moon up during the window?
    obs = mmt.mmtObserver

    # Do a check at startTime
    obs.date = startTime
    obs.horizon = "-0:34"
    prevMoonSet = obs.previous_setting(pyEphem.Moon()).datetime()
    prevMoonRise = obs.previous_rising(pyEphem.Moon()).datetime()

    # By definition, startTime is between prev and next, so
    # we need to check and see if the previous rise or set was most
    # recent
    if (prevMoonSet > prevMoonRise):
        # This is darktime since the most recent event was the
        # moon setting.
        moonAtStart = False
    else:
        moonAtStart = True

    # Do a check at endTime
    obs.date = endTime
    prevMoonSet = obs.previous_setting(pyEphem.Moon()).datetime()
    prevMoonRise = obs.previous_rising(pyEphem.Moon()).datetime()

    # Check again at end time.  Is the most recent event still
    # the moon setting? If yes, the moon is still down.
    if (prevMoonSet > prevMoonRise):
        # This is darktime
        moonAtEnd = False
    else:
        moonAtEnd = True

    isMoonUp = (moonAtStart | moonAtEnd)  # Logical OR, only care if both are 0
    return int(isMoonUp)


def calcMoonFlag(fld, startTime, endTime, mmt):
    """Calculate the IGNORE_FLAG based on lunar brightness and position.

    Inputs:
        fld -- field parameter entry from obsPars
        startTime -- datetime formatted starting time
        endTime -- datetime formatted ending time

    Output:
        flag -- 0/1 flag marking if the field is too close to the moon or
                the moon is brighter than specified.
    """
    moonUp = moonUpDuringObs(startTime, endTime, mmt)
    moonAge = MMTEphem.moonAge(startTime)

    # Get the distance to the moon
    moonDist = lunarDistance(fld, startTime, endTime)

    # Now do brightness flag
    moonReq = fld['moon'].values[0]
    if (moonReq == 'bright') | (moonUp == 0):
        # Anything works, either we were asked for bright
        # time or the moon isn't up during
        # the entirety of the observation
        illumFlag = 0
    elif (moonReq == 'grey') & (abs(moonAge) < 9) & (moonDist < 90):
        # We were asked for grey time and it's grey time
        illumFlag = 0
    elif (moonReq == 'dark') & (abs(moonAge) < 4.5) & (moonDist > 90.0):
        # We were asked for dark time and it's dark time
        illumFlag = 0
    else:
        illumFlag = 1e6   # This isn't going to work here!

    return illumFlag


def willItFitWeight(fld, startTime, mmt):
    """Calculate a weight that measures the fraction of a field that fits.

    Here, smaller weights are better, so we look at the fraction of the
    observation that will *not* fit (so 0 means it fits, so it has a better
    chance of being observed).

    Inputs:
        fld : single entry from the FLDpar DataFrame
        startTime : datetime formatted time the observation window starts

        mmt : MMTEphem information for the given night.

    Outputs:
        weight : weight from 0 to 1 with fraction of field not able to be
                 observed
    """

    # Throw an error if there are more than one
    if len(fld) > 1:
        raise AssertionError("There are more than 1 fields with name " +
                             fld['objid'])

    # First, get the end of the night time
    nightEnd = mmt.morningTwilight
    timeRemaining = (nightEnd - startTime).total_seconds()  # In seconds

    # This loop may be specific for MMIRS as BINOSPEC
    # may not have the nvisit/nexposure definitions.

    # Exposure time is stored in minutes
    exptime = float(fld['exptime'].values[0]) * 60.0
    repeats = float(fld['repeats'].values[0])
    nexp = float(fld['nexp'].values[0])

    # Calculate the number of visits that fit (keep it whole)
    # We don't want to worry about partial visits
    expPerVisit = exptime * nexp
    possVisits = math.floor((timeRemaining - mmirsOverhead(fld)) / expPerVisit)

    if possVisits > repeats:

        totalTargetTime = expPerVisit*repeats + mmirsOverhead(fld)
        tdelta = datetime.timedelta(seconds=totalTargetTime)

        # The number of possible visits is larger than requested so
        # just return the starting time incremented by total length
        return 0, startTime + tdelta
    else:

        totalTargetTime = expPerVisit*possVisits + mmirsOverhead(fld)
        tdelta = datetime.timedelta(seconds=totalTargetTime)
        return 1.0 * possVisits / repeats, startTime+tdelta


def obsUpdateRow(fldPar, donePar, startTime, mmt):
    """Calculate observing weight for given observation.

    Input:
        fldpar : DataFrame for a single observation. Comes from
                 readAllFLDFiles (i.e. not a single dictionary)
        donepar : Dataframe with stats pertaining to a fields doneness
        starTime : datetime string for beginning time for observation

    Output:
        weight : float (0 to 1e6) representing how well this fits into the
                 observation (lower = better). Observations that get
                 IGNORE_FLAG are set to 1e6.
    """
    # The idea here is to loop through all of the fields and calculate
    # the weight for this startTime.
    for objID in fldPar["objid"]:
        fld = fldPar[fldPar["objid"] == objID]
        fitWeight, endTime = willItFitWeight(fld, startTime, mmt)
        moonFlag = calcMoonFlag(fld, startTime, endTime, mmt)
        #
        print(objID, startTime, endTime, moonFlag, mmt.moonset)




def main(args):
    """Main module where the bulk of the work is completed."""
    # Get all of the observations for this dataset
    obsPars = queueTools.readAllFLDfiles()

    # Create the blank done file if it doesn't exist
    donefile = 'mmirs_catalogs/donelist.txt'
    if os.path.isfile(donefile):
        # Read in the existing file
        pass
    else:
        donePar = queueTools.createBlankDoneMask(obsPars)

    # Run one call of obsUpdateRow as a test
    date = "2016/02/16"
    mmt = MMTEphem.ephem(date)
    startTime = mmt.eveningTwilight

    obsUpdateRow(obsPars, donePar, startTime, mmt)


if __name__ == "__main__":
    main(sys.argv)
