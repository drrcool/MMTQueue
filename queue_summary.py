"""Creates human-readable "todo" list for the night.

Takes the queue schedule and the parameter files and creates a file for humans
to read to use as a "checklist" for the night.

"""
import queueTools


def schedule_to_text(trimester):
    """Parse the schedule.dat file to text"""

    print(trimester)
    obsPars = queueTools.readAllFLDfiles(trimester)

    schedule_file = "schedule.dat"
    f = open(schedule_file, 'r')

    currentNight = 0

    for line in f.readlines():
        startD, startT, endD, endT, afield, nvisit = line.strip().split()

        if startD != currentNight:
            print(''.join(['*']*80))
            print(''.join(['*']*80))
            print("Starting UT Night of %s" % startD)
            print(''.join(['*']*80))
            print(''.join(['*']*80))
            print("")
            currentNight = startD


        # Get this data
        par = obsPars[obsPars["objid"] == afield]

        print("".join(['*']*80))
        print("")
        print("Field: %s" % afield)
        print("Obs Type: %s" % par['obstype'].values[0])
        print("PI: %s" % par['PI'].values[0])
        print("Start Time (UT): %s" % startT)

        print("")
        print("RA: %s" % par['ra'].values[0])
        print("DEC: %s" % par['dec'].values[0])
        print("Estimated Mag: %s" % par['mag'].values[0])
        print("")
        print("Mask: %s" % par['mask'].values[0])
        print("")
        exptime = [float(x)*60.0 for x in par['exptime'].values][0]
        print("Exposure Time: %s" % exptime)
        print("Number of points per dither pattern: %s" %
              par['nexp'].values[0])
        print("Number of Dither Patterns to Run Tonight: %s" % nvisit)
        print("")
        print("Dithersize: %s" % par['dithersize'].values[0])
        print("Filter: %s" % par['filter'].values[0])
        print("Grism: %s" % par['grism'].values[0])
        print("Gain: %s" % par['gain'].values[0])
        print("ReadTab: %s" % par['readtab'].values[0])
        print("")
        print("Requested Seeing: %s" % par['seeing'].values[0])
        print("Photometric Required?: %s" % par['photometric'].values[0])
        for _ in range(3):
            print("")


def main(trimester):
    """Parse the inputs and create a text schedule for each date.

    Inputs : string of whitespace separated dates
    """

    schedule_to_text(trimester)




if __name__ == "__main__":
    main('2016a')
