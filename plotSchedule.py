""" Testbed for schedule plotting code."""

import matplotlib.pyplot as plt

import datetime
import os
import sys
import MMTEphem


def strings2datetime(date, time):
    y, m, d = map(int, date.split('/'))
    H, M, S = map(int, time.split(':'))
    dt = datetime.datetime(y, m, d, H, M, S)
    return dt

def main(argv):
    """Test code."""

    fileName = 'schedule.dat'
    f = open(fileName, 'r')
    startTime = []
    endTime = []
    field = []
    date = []
    initstartD = None
    ii= 0
    index = []
    for line in f.readlines():
        startD, startT, endD, endT, iiField = line.strip().split()

        date.append(startD)
        startTime.append(strings2datetime(startD, startT))
        endTime.append(strings2datetime(endD, endT))
        field.append(iiField)
        index.append(ii+0.5)
        ii += 1

    udates = sorted(set(date))

    plt.rcParams['font.size'] = 12
    plt.plot(startTime, index, '>', color='darkblue')
    plt.plot(endTime, index, '<', color='darkblue')

    plt.ylim(max(index)+0.5, 0)
    plt.xlim(min(startTime)-datetime.timedelta(hours=12),
             max(endTime)+datetime.timedelta(hours=12))
    plt.gcf().autofmt_xdate()

    mindate = strings2datetime(min(udates), '00:00:00')
    maxdate = strings2datetime(max(udates), '00:00:00')
    idate = mindate - datetime.timedelta(days=2)
    while idate < maxdate+datetime.timedelta(days=2):
        mmt = MMTEphem.ephem(idate)
        plt.fill_between([mmt.eveningTwilight, mmt.morningTwilight],
                        [-1, 100], alpha=0.25, color='darkblue')

        idate += datetime.timedelta(days=1)

    for ii in range(len(startTime)):
        diff = (endTime[ii]-startTime[ii]).total_seconds()/2.0
        xlabel = startTime[ii]-datetime.timedelta(seconds=diff*2)
        ylabel = index[ii]-0.1
        plt.annotate(field[ii], xy=(xlabel, ylabel))

    plt.title("MMIRS Schedule March 2016")
    plt.xlabel("UT Date")
    plt.ylabel("Obs Index")
    plt.show()
    plt.savefig('schedule.png')

if __name__ == "__main__":
    main(sys.argv)
