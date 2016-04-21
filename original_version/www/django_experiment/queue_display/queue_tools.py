import pandas as pd
from os import walk


def parse_mask_position_angle(mask, runname):

    """Parse the .msk file for a mask to get it's position angle.

    For the March 2016 run, the position angle was not written to the FLD
    files for any mask observations. This makes checking the rotator limits
    impossible. This code parses the .msk file to get the needed position
    angle to add to the fldPar.
    """
    # Read the mask file
    maskfile = '/Users/rcool/MMTQueue/experiment/mmirs_catalogs/' + \
        runname + '/' + mask + '.msk'
    f = open(maskfile, 'r')

    # Check each line in the maskfile and find the line that starts with 'pa'
    mask_position_angle = False
    for line in f.readlines():
        sline = line.strip().split()
        if len(sline) > 1 and sline[0] == 'pa':
            mask_position_angle = float(sline[1])
    f.close()

    return mask_position_angle


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

    return obspars


def read_all_fld_files(runname):
    """Read all of the fld files for a run and output a dataframe."""

    # Set the path
    path = '/Users/rcool/MMTQueue/experiment/mmirs_catalogs/' + \
        runname

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
