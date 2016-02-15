"""A collection of tools to augment the MMT queue scheduler.

These tools are ones that aren't explicitly handling the ephemeris for
targets (found in MMTEphem) or the queue calculations themselves.

These include wrappers on common tasks, routines to read input
and format output.
"""

# Imports
from pandas import DataFrame
from os import walk


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


def createBlankDoneMask(obsFrame):
    """Create a blank structure to hold information about completed data.

    Input:
        obsFrame : DF of all observations being considered
    Output:
        doneFrame : DF with one line per observation with all "done"
                    information set to not done.
    """
    nObs = len(obsFrame)

    # Create the blank format
    blankDict = {}
    blankDict['complete'] = 0  # Is this observation completely done?
    blankDict['doneVisit'] = 0  # How many visits have been completed?
    blankDict['visitsScheduled'] = 0  # How many visits have been scheduled?
    blankDict['weightBK'] = 1e6  # This is a bookkeeping weight.
    blankDict['weightPartial'] = 0  # This says some visits were taken.

    dictList = []
    for ii in range(0, nObs):
        dictList.append(blankDict)

    outFrame = DataFrame(dictList)
    return outFrame

if __name__ == "__main__":
    pass
