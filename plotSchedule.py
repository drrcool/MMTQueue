""" Testbed for schedule plotting code."""

import matplotlib.pyplot as plt

import datetime
import numpy as np
import sys
import re


def strings2datetime(date, time):
    """Convert a date and time string to a datetime."""
    y, m, d = map(int, date.split('/'))
    H, M, S = map(int, time.split(':'))
    dt = datetime.datetime(y, m, d, H, M, S)
    return dt


def string2decTime(time):
    """Convert a time to decimal."""
    h, m, s = map(float, time.split(':'))
    seconds = h*3600.0 + m*60.0 + s
    return seconds / 3600.


def main(argv):
    """Test code."""

    fileName = 'schedule.dat'
    f = open(fileName, 'r')
    startTime = []
    endTime = []
    field = []
    dates = []

    for line in f.readlines():
        startD, startT, endD, endT, afield, nvisit = line.strip().split()

        dates.append(startD)
        startTime.append(string2decTime(startT))
        endTime.append(string2decTime(endT))

        reString = "^[a-zA-Z]+-[A-Za-z0-9]+_(.*)$"
        m = re.search(reString, afield)
        field.append(m.group(1))

    print(min(startTime), max(endTime))
    udates = sorted(set(dates))
    mindate = strings2datetime(min(udates), "00:00:00")
    maxdate = strings2datetime(max(udates), '00:00:00')

    nDays = (maxdate-mindate).days

    plt.figure(figsize=(12, 12))
    plt.plot([0, 0], [0, 0], visible=False)
    plt.rcParams['font.size'] = 8

    xpos = []
    for ii in range(len(dates)):
        idate = strings2datetime(dates[ii], "00:00:00")
        index = (idate-mindate).days
        xpos.append(index)
        # Package the rectangle
        xval = [index - 0.45, index + 0.45]
        y1 = [startTime[ii], startTime[ii]]
        y2 = [endTime[ii], endTime[ii]]

        plt.fill_between(xval, y1, y2, alpha=0.15)
        plt.text(index, startTime[ii]+0.5*(endTime[ii]-startTime[ii]),
                 field[ii], rotation=90, va='center')

    plt.title("MMIRS Schedule March 2016")
    plt.ylabel("Time (UT hour)")

    plt.xticks(xpos, dates, rotation=90)
    plt.xlabel("UT Date (local day + 1)")

    plt.xlim(-0.8, nDays+0.8)
    plt.ylim(max(endTime)+1, min(startTime)-1)
    plt.grid(0)

    #plt.show()
    plt.savefig('schedule.pdf')

if __name__ == "__main__":
    main(sys.argv)
