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


def main(args):
    """Main module where the bulk of the work is completed."""
    pass

if __name__ == "__main__":
    main(sys.argv)
