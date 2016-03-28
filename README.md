# MMTQueue

### Purpose:
This software is intended to run the MMT queue operations.
For now, this specifically is for MMIRS, but that will
eventually change to include Binospec.

Please contact Richard Cool for questions or concerns (rcool@mmto.org).

This code is **not** production level at this point.  Work is ongoing.

## Current Needs

1. Better input from the astronomers
  * Getting faulty or incomplete details was the biggest issue
  during the March run.  We need to overhaul how astronomers
  submit catalogs (and maybe masks).
  * This interface needs a number of error checks for common errors:
    * for imaging are they requesting a spectroscopic filter?
    * for longslit and masks did they provide an even number
    of positions per dither
    * Are the exposure times reasonable? For most things this
    should be less than 5 minutes.
    * Are the ramps set properly for the exposure time?
  * We need a way to get the dither information for imaging
  * Take some freedom away from the form -- set default dither
  patterns and make them override this if they really want to
  * Probably want to feed a database. This makes querrying so much easier
  * Better ways to communicate conditional priorities
2. Implement priorities for a PI. Right now, priority 2 is set to
   globally lower priority than priority 1. This should change.
3. Reevaluate how we do fractional fields. This may be biasing us?
4. Look at re-implementing partial completeness as a flag for ranking within
   a single PI. 
