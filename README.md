# MMTQueue

### Purpose:
This software is intended to run the MMT queue operations. For now, this specifically is for MMIRS, but that will eventually change to include Binospec.

Please contact Richard Cool for questions or concerns (rcool@mmto.org).

This code is **not** production level at this point.  Work is ongoing.

## Todos:
1. <s> Create functions for sunrise and sunset </s>
2. Create function to calculate to compute moon position and brightness at given time
3. Create Data Structure / Model for storing "DONE" observations
4. Create Data Structure / Model for observing parameters
5. Create function to check for time until transit
6. Create function that edits TAC weight for each target given observed status
7. Create visualization tools to see what's been scheduled on a given night
8. Create schedule output for observer
9. Create "poor weather" override mode


####Current progress:

I'll try and keep this updated, but no promises.  

Feb 9 -- Repo created.


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

Need to write a routine that checks if the moon is up at start and at end.  It's not really possible that the moon is down for both but was up in the middle.  