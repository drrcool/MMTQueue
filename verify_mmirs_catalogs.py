"""Verification tools for MMIRS catalogs before they are scheduled.

We focus here on the "low lying fruit" of errors we found were pretty common
during the runs early in 2016. There are likely many more we could go after,
and those should be added as separate modules.

The design philosophy is that when run with a trimester name, ever FLD file
will be verified.
"""
import sys
import os


def main(args):
    """Pipeline the user experience."""
    # Check to ensure that a trimester was specified
    try:
        trimester = args[1]
    except IndexError:
        print("Syntax: ./verify_mmirs_catalogs.py trimester_name")
        print("No trimester name was given")
    

if __name__ == "__main__":
    main(sys.argv)
