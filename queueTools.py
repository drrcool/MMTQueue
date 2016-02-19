"""A collection of tools to augment the MMT queue scheduler.

These tools are ones that aren't explicitly handling the ephemeris for
targets (found in MMTEphem) or the queue calculations themselves.

These include wrappers on common tasks, routines to read input
and format output.
"""

# Imports
from pandas import DataFrame
from os import walk
import MMTEphem
import ephem as pyEphem

def readAllocatedTime(startDay="1900/1/1", endDay="3000/1/1"):
    """Read a log file to determine how much each PI was allocated."""
    filename = "AllocatedTime.dat"
    f = open(filename, 'r')

    if type(startDay) == str:
        startDay = pyEphem.date(startDay).datetime()
    if type(endDay) == str:
        endDay = pyEphem.date(endDay).datetime()


    allocatedTime = {}
    for line in f.readlines():
        if line[0] == '#':
            # Comment string, skip
            continue

        date, PI = line.strip().split()
        date = pyEphem.date(date).datetime()
        mmt = MMTEphem.ephem(date)
        if (abs((date-startDay).total_seconds()) < 24*3600.) | \
           (abs((date-endDay).total_seconds()) < 24*3600.):
            # This night is not in the current run
            continue

        nightLength = (mmt.morningTwilight - mmt.eveningTwilight)
        nightLength = nightLength.total_seconds() / 3600.0

        if PI in allocatedTime:
            allocatedTime[PI] += nightLength
        else:
            allocatedTime[PI] = nightLength

    return allocatedTime


def readFLDfile(fileName):
    """Read FLD file and return a dictionary with relavent information."""
    # Intiialize output dictionary
    obspars = {}

    f = open(fileName, 'r')

    # Read the PI
    junk, piName = f.readline().strip().split()
    obspars['PI'] = piName

    # Program ID
    junk, progID = f.readline().strip().split()
    obspars['progID'] = progID

    # Add filename
    obspars['fldfile'] = fileName

    # Parse the column keywords
    keywords = f.readline().strip().split()
    f.readline()  # line of ---
    values = f.readline().strip().split()

    for key, val in zip(keywords, values):
        obspars[key] = val

    return obspars


def readAllFLDfiles(path=None):
    """Read all FLD files in path (or a walk through path)."""

    # Path
    if path is None:
        path = '/Users/rcool/MMTQueue/mmirs_catalogs/'
    f = []
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend([dirpath+'/'+f for f in filenames])

    outfiles = [xx for xx in f if xx[-4:] == '.fld']

    dictList = []
    for file in outfiles:
        df = readFLDfile(file)
        dictList.append(df)

    outFrame = DataFrame(dictList)

    return outFrame


def createBlankDoneMask(obsFrame, startDay="1900/1/1", endDay="3000/1/1"):
    """Create a blank structure to hold information about completed data.

    Input:
        obsFrame : DF of all observations being considered
    Output:
        doneFrame : DF with one line per observation with all "done"
                    information set to not done.
    """
    nObs = len(obsFrame)
    allocatedTime = readAllocatedTime(startDay=startDay, endDay=endDay)

    # Create the blank format
    blankDict = {}
    blankDict['complete'] = 0  # Is this observation completely done?
    blankDict['doneVisit'] = 0  # How many visits have been completed?
    blankDict['totalTime'] = 0.0
    blankDict['doneTime'] = 0.0  # Total (for PI) completed time
    blankDict['prevWeight'] = 1.0  # For tracking weights from previous pass
    blankDict['currentWeight'] = 0.0
    blankDict['PI'] = ''
    blankDict['objid'] = ''

    dictList = []
    for ii in range(0, nObs):
        copyDict = blankDict.copy()
        PI = obsFrame['PI'].values[ii]
        copyDict['objid'] = obsFrame['objid'].values[ii]
        copyDict['PI'] = obsFrame["PI"].values[ii]
        copyDict['totalTime'] = allocatedTime[PI]
        dictList.append(copyDict)

    outFrame = DataFrame(dictList)
    return outFrame


if __name__ == "__main__":
    pass
