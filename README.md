# MMTQueue

### Purpose:
This software is intended to run the MMT queue operations. For now, this specifically is for MMIRS, but that will eventually change to include Binospec.

Please contact Richard Cool for questions or concerns (rcool@mmto.org).

This code is **not** production level at this point.  Work is ongoing.

## Todos:
1. <s> Create functions for sunrise and sunset </s>
2. <s> Create function to calculate to compute moon position and brightness at given time</s>
3. <s> Create Data Structure / Model for storing "DONE" observations</s>
4. <s> Create Data Structure / Model for observing parameters</s>
5. Create function to check for time until transit
6. <s> Create function that edits TAC weight for each target given observed status</s>
7. <s> Create visualization tools to see what's been scheduled on a given night</s>
8. <s> Create schedule output for observer</s>
9. Create "poor weather" override mode



## Lunar age notes:

Typical definitions:

Dark : 0.0 < FLI < 0.25
Grey : 0.25 < FLI < 0.65
Bright : 0.65 < FLI < 1.00

Typically FLI moves about 0.1 / day

Trying to split fairly evenly, I make each window about 9 days

Dark :  abs(age) < 4.5
Grey : 4.5 < abs(age) < 9.0
Bright : 9.0 < abs(age) < 14.0

Note that we consider a night dark even if the night is "grey" but the observation can be fit in before moonrise.


# Big things to consider:
* How do we account for priority within a PI's program?  If I give you field 1 and field 2, and field 2 has lower priority, do I only do that if field 1 is done?
* What about about parity? Like I only want J band if you already have H
  
