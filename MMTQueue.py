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
import pandas as pd
import queueTools
import MMTEphem
import math
import datetime
import plotSchedule
import ephem as pyEphem
from random import randint
from pprint import pprint


def isObservable(objEphem, endTime):
    """Return a binary bit if target is observable at given time."""
    tarray = [abs(x-endTime).total_seconds() for x in objEphem.time]
    val, idx = min((val, idx) for (idx, val) in enumerate(tarray))

    return objEphem.observable[idx]


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
        return 3600.0
    elif obstype == 'longslit':
        return 1800.0
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
        illumFlag = 1
    elif (moonReq == 'grey') & (abs(moonAge) < 9) & (moonDist < 90):
        # We were asked for grey time and it's grey time
        illumFlag = 1
    elif (moonReq == 'dark') & (abs(moonAge) < 4.5) & (moonDist > 90.0):
        # We were asked for dark time and it's dark time
        illumFlag = 1
    else:
        illumFlag = 0   # This isn't going to work here!

    return illumFlag


def willItFitWeight(fld, startTime, mmt, objEphem, donePar):
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
        endTime : time this observation would finish. If not observable
                  will return the starttime. If only partially observable
                  will return the time at end of full visit.
        nVisits : the number of visits observed.
    """

    # Throw an error if there are more than one
    if len(fld) > 1:
        raise AssertionError("There are more than 1 fields with name " +
                             fld['objid'])

    # First, get the end of the night time
    nightEnd = mmt.morningTwilight
    timeRemaining = (nightEnd - startTime).total_seconds()  # In seconds

    # Is the target observable at the StartTime?
    if isObservable(objEphem, startTime) == 0:
        return 0.0, startTime, 0

    # This loop may be specific for MMIRS as BINOSPEC
    # may not have the nvisit/nexposure definitionss

    # Exposure time is stored in minutes
    exptime = float(fld['exptime'].values[0]) * 60.0
    repeats = float(fld['repeats'].values[0]) - \
        float(donePar['doneVisit'].values[0])
    nexp = float(fld['nexp'].values[0])

    # Calculate the number of visits that fit (keep it whole)
    # We don't want to worry about partial visits
    expPerVisit = exptime * nexp
    possVisits = math.floor((timeRemaining - mmirsOverhead(fld)) / expPerVisit)

    if possVisits > repeats:

        totalTargetTime = expPerVisit*repeats + mmirsOverhead(fld)
        tdelta = datetime.timedelta(seconds=totalTargetTime)

        # Now, let's see if the target is still observable at the end of
        # the block
        endTime = startTime + tdelta
        if isObservable(objEphem, endTime) == 1:
            return 1.0, endTime, repeats
        else:
            possVisits = repeats

    # Ok, we couldn't fit all the observations in. Let's figure out
    # how much we can
    nVisitsObservable = possVisits

    while nVisitsObservable >= 1:
        totalTargetTime = expPerVisit*nVisitsObservable + mmirsOverhead(fld)
        tdelta = datetime.timedelta(seconds=totalTargetTime)
        endTime = startTime + tdelta

        # Check to see if it's observable at the end of the window
        if isObservable(objEphem, endTime) == 1:
            return 1.0 * nVisitsObservable/repeats, endTime, nVisitsObservable
        else:
            nVisitsObservable -= 1

    # If we get here, the object is never observable, return bad weight
    return 0, startTime, 0


def obsUpdateRow(fldPar, donePar, startTime, mmt):
    """Calculate observing weight for given observation.

    Input:
        fldpar : DataFrame for a single observation. Comes from
                 readAllFLDFiles (i.e. not a single dictionary)
        donepar : Dataframe with stats pertaining to a fields doneness
        starTime : datetime string for beginning time for observation

    """
    # The idea here is to loop through all of the fields and calculate
    # the weight for this startTime.
    weightList = []  # Will store array of dicts with weight info

    for objID in fldPar["objid"]:
        # Is this already observed?
        donefld = donePar[donePar['objid'] == objID]
        fld = fldPar[fldPar["objid"] == objID]

        if donefld['complete'].values[0] == 1:
            continue

        # Initialize the weight dictionary
        obsWeight = {}
        obsWeight['objid'] = objID

        objEphem = MMTEphem.ObjEphem(fld['ra'].values[0],
                                     fld['dec'].values[0],
                                     startTime, mmt)
        fitWeight, endTime, fitVisits = willItFitWeight(fld,
                                                        startTime, mmt,
                                                        objEphem,
                                                        donefld)
        moonFlag = calcMoonFlag(fld, startTime, endTime, mmt)

        obsWeight['fitWeight'] = fitWeight
        obsWeight['obsTime'] = (endTime-startTime).total_seconds()
        obsWeight['moonFlag'] = moonFlag
        obsWeight['nVisits'] = fitVisits

        # Priority Weight
        priorFlag = 1.0 / float(fld['priority'])

        # Now combine the weights
        weightTAC = 1.0 - 1.0 * donefld['doneTime'].values[0] \
            / donefld['totalTime'].values[0]
        if weightTAC < 0:
            weightTAC = 0.001
        totalWeight = fitWeight * moonFlag * (weightTAC)**2*priorFlag
        # We need to account for a few extra things.
        # 1. I want fields that have had previously observed
        #    fields to be the top priority if they are observable.
        # This modification will weight fields with no observations
        # at one (1+0) and those partially observed at 10
        partCompWeight = int(donefld['doneVisit'].values[0] > 0)
        totalWeight = totalWeight * (1+9.0*partCompWeight)

        # Now account for weighting from preivous iteration of the
        # code. This should smooth things out
        totalWeight = totalWeight * donefld['prevWeight'].values[0]
        obsWeight['totalWeight'] = totalWeight
        # Append to weight listing
        weightList.append(obsWeight)

    obsWeights = pd.DataFrame(weightList)

    # Doing the O(n) problem here
    hasMax = []

    # Check to see there are actually targets!
    if len(obsWeights) == 0:
        return None

    maxWeight = max(obsWeights['totalWeight'])
    for ii in range(len(obsWeights)):
        if obsWeights['totalWeight'].values[ii] == maxWeight:
            hasMax.append(ii)

    # Now, check there are some with non-zero weight and
    # Observe a random one
    if maxWeight == 0:
        # No targets were done.
        return None
    else:
        randIndex = hasMax[randint(0, len(hasMax)-1)]

        diffTime = obsWeights['obsTime'].values[randIndex]
        # Increment the schedule
        schedule = []
        schedule.append(startTime)
        schedule.append(diffTime)
        schedule.append(obsWeights['objid'].values[randIndex])

        # Update the done masks
        index = [i for i, x in enumerate(donePar['objid'] ==
                                         schedule[2]) if x]

        # Donetime tracks total time for this PI.  We weight by
        # this, so let's increment all entries in donefld

        PI = donePar.loc[index, 'PI'].values[0]
        for ii in range(len(fldPar)):
            if donePar["PI"].values[ii] == PI:
                donePar.loc[ii, 'doneTime'] += diffTime/3600.
                donePar.loc[ii, 'currentWeight'] = \
                    donePar.loc[ii, 'doneTime'] / \
                    donePar.loc[ii, 'totalTime']

        donePar.loc[index, 'doneVisit'] += \
            obsWeights['nVisits'].values[randIndex]

        # Check to see if the target is now done
        requested = int(fldPar[fldPar['objid'] ==
                               schedule[2]]["repeats"].values[0])
        if int(donePar.loc[index, 'doneVisit']) >= requested:
            donePar.loc[index, 'complete'] = 1
        return schedule


def obsOneNight(fldPar, donePar, date):
    """Fully Schedules one night."""
    # Get the ephemeris for the night
    mmt = MMTEphem.ephem(date)
    startTime = mmt.eveningTwilight

    # Now, start at twilight and add observervations. Each
    # time increment currentTime by the observation time
    # If no observations are found, add 5 minutes and check again
    currentTime = startTime
    schedule = []  # Will contain scheduled fields
    allDone = False
    while (currentTime < mmt.morningTwilight) & (allDone is False):
        newSched = obsUpdateRow(fldPar, donePar, currentTime, mmt)
        # Check to see if something was observed
        if newSched is None:
            # Increment time and continue
            currentTime += datetime.timedelta(minutes=20)
        else:
            # Append new entry to schedule, increment currentTime
            # by exposure time.
            schedule.append(newSched)
            currentTime += datetime.timedelta(seconds=newSched[1])
            completed = donePar['complete'].values
            if min(completed) == 1:
                allDone = True

    return schedule


def main(args):
    """Main module where the bulk of the work is completed."""
    # Get all of the observations for this dataset
    if len(args) > 1:
        trimester = args[1]
    else:
        trimester = None  # The default is handled in queueTools
    obsPars = queueTools.readAllFLDfiles(trimester)


    # Create the blank done file if it doesn't exist
    donePar = queueTools.createBlankDoneMask(obsPars)

    donefile = 'donefile.dat'
    if os.path.isfile(donefile):
        # Read in the existing file
        f = open(donefile)
        for line in f.readlines():
            if line[0] != "#":
                split = line.split()
                id = split[0]
                pi = split[1]
                visits = float(split[2])
                donetime = float(split[3])


                donePar.loc[donePar['objid']==id, 'doneVisit'] = visits
                donePar.loc[donePar['PI']==pi, 'doneTime'] += donetime


    # Run one call of obsUpdateRow as a test
    allDates = ["2016/03/18",
                "2016/03/19",
                "2016/03/20",
                "2016/03/21",
                "2016/03/24",
                "2016/03/25",
                "2016/03/26",
                "2016/03/27",
                "2016/03/28"]

    finishedFlag = False
    for iter in range(5):

        # Check to see if I reached a success mark
        if finishedFlag is True:
            continue

        print("")
        print("**** ITERATION # %i ********" % (iter+1))
        fullSched = []
        for date in allDates:
            sys.stdout.write("\r Working on Date %s   " % date)
            schedule = obsOneNight(obsPars, donePar, date)
            for line in schedule:
                fullSched.append(line)

        # Now, fill in previous weights and zero out entries in a
        # new donePar
        newDone = donePar.copy()
        newDone.loc[:, 'complete'] = 0
        newDone.loc[:, 'doneTime'] = 0.0
        newDone.loc[:, 'doneVisit'] = 0.0

        # Make a list of PIS and see if they got all the fields they
        # asked for.
        allDone = {}
        for ii in range(len(newDone)):
            PI = newDone["PI"].values[ii]
            if PI not in allDone:
                completed = donePar[donePar['PI'] == PI]['complete'].values
                if min(completed) == 1:
                    allDone[PI] = True
                else:
                    allDone[PI] = False
            weightMod = int(allDone[PI])
            newDone.loc[ii, 'prevWeight'] = \
                donePar.loc[ii, 'currentWeight']*(1-weightMod*0.9)

        if min(allDone.values()) == 1:
            # We've scheduled everything, stop iterating
            finishedFlag = True
        newDone.loc[:, 'currentWeight'] = 0.0
        pprint(donePar)
        if iter != 5:
            donePar = newDone

    # Format the schedule for output
    f = open('schedule.dat', 'w')
    for sched in fullSched:
        startTime = sched[0]
        duration = sched[1]
        field = sched[2]

        FORMAT = "%Y/%m/%d %H:%M:%S"
        outStart = startTime.strftime(FORMAT)
        endTime = startTime + datetime.timedelta(seconds=duration)
        outEnd = endTime.strftime(FORMAT)
        f.write("%s %s %s\n" % (outStart, outEnd, field))

    f.close()
    plotSchedule.main([])

if __name__ == "__main__":
    main(sys.argv)
