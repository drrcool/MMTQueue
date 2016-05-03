import ipdb
"""Verification tools for MMIRS catalogs before they are scheduled.

We focus here on the "low lying fruit" of errors we found were pretty common
during the runs early in 2016. There are likely many more we could go after,
and those should be added as separate modules.

The design philosophy is that when run with a trimester name, ever FLD file
will be verified.
"""
import sys
import os
from MMTQueue import read_all_fld_files


def get_fld_file_list(trimester):
    """Return a list of all FLD files for the given trimester."""
    # Set the path
    path = 'catalogs/' + trimester
    filelist = []
    for (dirpath, dirnames, filenames) in os.walk(path):
        filelist.extend([dirpath+'/' + f
                         for f in filenames if f[-4:] == '.fld'])

    return filelist


def do_maskfiles_exist(trimester):
    """Identify any missing msk files for this run."""
    fldfiles = get_fld_file_list(trimester)

    needed_masks = []
    # Loop through the fldfiles and check for obstype=mask
    for ifile in fldfiles:
        f = open(ifile, 'r')
        # Read off the header
        _ = f.readline()
        _ = f.readline()
        param_names = f.readline().strip().split()
        _ = f.readline()
        param_values = f.readline().strip().split()
        f.close()

        fldpar = {}
        for x, y in zip(param_names, param_values):
            fldpar[x] = y
        if fldpar['obstype'] == 'mask':
            needed_masks.append(fldpar['mask'])

    missing_msk_files = []
    for mask in needed_masks:
        mask_file = 'catalogs/' + trimester + '/masks/' + mask + '.msk'
        if os.path.isfile(mask_file) == 0:
            missing_msk_files.append(mask)

    if len(missing_msk_files) > 0:
        print("MISSING .msk FILES IDENTIFIED. PLEASE FIND .msk FILES " +
              "FOR THE FOLLOWING MASKS BEFORE CONTINUING")
        print("".join(["*"]*80))
        for mask in missing_msk_files:
            print(mask)
        return 0

    return 1


def add_error(error_log, errorstring, PI, field):
    """Add error to the error_log[PI][field] element."""
    # If no errors yet, just create the error array,
    # otherwise, append
    if error_log[PI][field] is not None:
        error_log[PI][field].append(errorstring)
    else:
        error_log[PI][field] = [errorstring]

    return error_log


def verify_exposure_times(fldpar, error_log):
    """Check each fld file for exposure times that make sense.

    Primarily, we check that all exposure times are < 300s and > 1s.
    If any of these fail, a dictionary containing errors will be generated
    that will be used later in the code to generate email responses to the
    PI
    """
    error_count = 0
    # We need to loop through all the fld files
    for ii in range(len(fldpar)):
        objpar = fldpar.iloc[ii]
        PI = objpar['PI']
        field = objpar['objid']
        exptime = float(objpar['exptime'])*60  # Convert to seconds

        if (exptime < 1) or (exptime > 300):
            error_count += 1

            errorstring = \
                ("* Has a specified exposure time "
                 "of {1} seconds per image."
                 .format(field, round(exptime, 1)))

            error_log = add_error(error_log, errorstring, PI, field)

    return error_count


def create_error_log(fldpar):
    """Create a blank dictionary to hold any errors for the file.

    The organization is by PI, then list of fields, and a list of
    errors. We'll parse this out later.
    """
    error_dict = {}
    for ii in range(len(fldpar)):
        PI = fldpar.iloc[ii]['PI']
        if PI not in error_dict:
            error_dict[PI] = {fldpar.iloc[ii]['objid']: None}
        else:
            error_dict[PI][fldpar.iloc[ii]['objid']] = None

    return error_dict


def verify_nvisits(fldpar, error_log):
    """Check visits parameter to ensure that visits > 0"""
    error_count = 0
    for ii in range(len(fldpar)):
        PI = fldpar.iloc[ii]['PI']
        nvisits = fldpar.iloc[ii]['repeats']
        field = fldpar.iloc[ii]['objid']

        if nvisits == 0:
            error_count += 1
            errorstring = \
                ("* Has a specified number of repeats set \n"
                 "set to zero.  Remember that the total number of  \n"
                 "exposures is defined to be repeats*nexp where repeats \n"
                 "is the number of time a dither pattern with nexp \n"
                 "positions should be repeated.".format(field))
            error_log = add_error(error_log, errorstring, PI, field)
    return error_count


def verify_nexp(fldpar, error_log):
    """Check nexp to ensure it's valid.

    Nexp *must* be even for spectroscopy.

    If longslit and not 2, then note that additional info will be needed.

    If mask and not 4, then note that additional info will be needed.
    """
    error_count = 0
    for ii in range(len(fldpar)):
        PI = fldpar.iloc[ii]['PI']
        field = fldpar.iloc[ii]['objid']
        nexp = int(fldpar.iloc[ii]['nexp'])
        obstype = fldpar.iloc[ii]['obstype']

        if obstype == 'mask':
            if nexp != 4:
                error_count += 1
                errorstring = \
                    ("* Has nexp set to {1}. Typically, for masks "
                     "a setting 4-position \n "
                     "dither pattern is used.".format(field, nexp))
                error_log = add_error(error_log, errorstring, PI, field)
        elif obstype == 'longslit':
            if (nexp % 2) != 0:
                error_count += 1
                errorstring = \
                    ("* Has a nexp set to {1}.  For proper sky "
                     "subtraction, an *even* number\n"
                     " of exposures is "
                     "required.".format(field, nexp))
                error_log = add_error(error_log, errorstring, PI, field)

    return error_count


def verify_readtab(fldpar, error_log):
    """Ensure that the proper readtab is used for long and short exposures.

    Typicall, we want to use ramp_4.426 for exposure times > 30s and
    ramp_1.475 for < 30s.
    """
    error_count = 0
    for ii in range(len(fldpar)):
        PI = fldpar.iloc[ii]['PI']
        field = fldpar.iloc[ii]['objid']
        exptime = float(fldpar.iloc[ii]['exptime'])*60
        readtab = fldpar.iloc[ii]['readtab']
        error_string1 = \
            ("* Has a readtab setting of {1}. Typically "
             "when exposure times \n"
             "are longer than 30s, we use "
             "ramp_4.426.".format(field, readtab))
        error_string2 = \
            ("* Has a readtab setting of {1}. Typically "
             "when exposure times \n"
             "are shorter than 30s, we use "
             "ramp_1.475.".format(field, readtab))
        if (exptime > 30) and (readtab == 'ramp_1.475'):
            error_count += 1
            error_log = add_error(error_log, error_string1, PI, field)
        if (exptime < 30) and (readtab == 'ramp_4.426'):
            error_count += 1
            error_log = add_error(error_log, error_string2, PI, field)

    return error_count


def verify_nondark(fldpar, error_log):
    """Check that mmirs catalogs only are asking for bright time."""
    error_count = 0
    for ii in range(len(fldpar)):
        PI = fldpar.iloc[ii]['PI']
        field = fldpar.iloc[ii]['objid']
        moon = fldpar.iloc[ii]['moon']

        if moon != 'bright':
            error_count += 1
            error_string = \
                ("* Specified to have a {1} lunar phase or "
                 "darker.  Please verify if \n"
                 "this is a requirement as it "
                 "will dramatically limit the flexibility of the MMIRS \n"
                 "queue observations.".format(field, moon))
            error_log = add_error(error_log, error_string, PI, field)
    return error_count


def count_errors_for_pi(error_log, pi):
    """Count the number of errors listed for a single PI."""
    error_count = 0

    # Check to be sure the PI is listed
    if pi not in error_log.keys():
        return 0

    # Now, count for each field
    for field, error in error_log[pi].items():
        if error is not None:
            error_count += len(error)

    return error_count


def create_email_templates(error_log, outroot):
    """Based on the error_log create an email to be sent to PIs.

    Email files are written to $outroot_PINAME_queuemail
    """
    # Loop through each of the PIs
    for pi in error_log.keys():
        # Only do this step if the PI actually has errors
        if count_errors_for_pi(error_log, pi) > 0:
            f = open('{0}_{1}_catalogerrors'.format(outroot, pi), 'w')
            email = 'Dear {0},'.format(pi)
            email = email + "\n\n" + \
                "Below, you will find a list of errors found in your \n" \
                "MMIRS catalog submission.  Note, some may be \n" \
                "intentional, in which case, please verify with us.\n\n"
            for field in error_log[pi].keys():
                if error_log[pi][field] is not None:
                    email = email + '\n' + '{0}:'.format(field) + '\n'
                    for error in error_log[pi][field]:
                        email = email + error + '\n'
            email = email + '\n \n' + ('').join(['-']*70) + '\n\n'
            f.write(email)
            f.close()


def main(args):
    """Pipeline the user experience."""
    # Check to ensure that a trimester was specified
    try:
        trimester = args[1]
        outroot = args[2]
    except IndexError:
        print("Syntax: ./verify_mmirs_catalogs.py trimester_name outroot")
        sys.exit(0)



    # First, check for missing files
    found_all_masks = do_maskfiles_exist(trimester)
    if found_all_masks == 0:
        return

    # Now, read all of the FLD files
    error_count = 0  # Initialize
    fldpar = read_all_fld_files(trimester)

    # Create the error log
    error_log = create_error_log(fldpar)

    # Verify the exposure times
    error_count += \
        verify_exposure_times(fldpar, error_log)

    # Check visits to ensure that it was actually set
    error_count += \
        verify_nvisits(fldpar, error_log)

    # Check to be sure the nexp set was feasible
    error_count += \
        verify_nexp(fldpar, error_log)

    # Check to be sure the proper readtabs are being used
    error_count += \
        verify_readtab(fldpar, error_log)

    # Check that if grey or bright is asked for, that this is right
    error_count += \
        verify_nondark(fldpar, error_log)


    create_email_templates(error_log, outroot)

if __name__ == "__main__":
    main(sys.argv)
