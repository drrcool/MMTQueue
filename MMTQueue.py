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
from pprint import pprint
import queueTools
import MMTEphem

def mmirsOberhead(fld):
    """Return the expected overhead time (in seconds) for MMIRS observation."""

    if fld['obstype'] == 'mask':
        return 600.0
    else if fld['obstype'] == ''

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
    print(len(fld))
    # First, get the end of the night time
    nightEnd = mmt.morningTwilight
    timeRemaining = (nightEnd - startTime).total_seconds()  # In seconds

    # This loop may be specific for MMIRS as BINOSPEC
    # may not have the nvisit/nexposure definitions.
    exptime = fld['exptime'].values * 60  # Exposure time is stored in minutes
    repeats = fld['repeats'].values
    nexp = fld['nexp'].values

    # Calculate the number of visits that fit (keep it whole)
    # We don't want to worry about partial visits
    expPerVisit = exptime * nexp
    possibleFinishedVisits = (timeRemaining - mmirsOverhead(fld)) / expPerVisit


    return return possibleFinishedVisits

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
        fitWeight = willItFitWeight(fld, startTime, mmt)
        print(objID, fitWeight)


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
    date = "2016/02/15"
    mmt = MMTEphem.ephem(date)
    startTime = mmt.eveningTwilight

    obsUpdateRow(obsPars, donePar, startTime, mmt)


if __name__ == "__main__":
    main(sys.argv)
