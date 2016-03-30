"""
Generate an interactive BarChart for the MMT Queue.

Based on the Bokeh library for python.

Assumes the schedule resides in schedule.dat
"""

from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import figure, show, output_file
import queueTools
import pandas as pd
import datetime
import re
from tableau_colormap import tableau20


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


def tuple_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


fldPar = queueTools.readAllFLDfiles()
schedfile = 'schedule.dat'
f = open(schedfile, 'r')

# We will use these colors to color each group by a different color
colormap = tableau20(raw=True)

startDate = []
startTime = []
startTimeOrig = []
endTime = []
field = []
repeats = []
zeroDateOrig = None
zeroDate = None
duration = []
midtime = []
objid = []
color = [tuple_to_hex(colormap[x % len(colormap)]) for x in range(len(fldPar))]
fldPar['color'] = color
idx = 0

for line in f.readlines():
    startD, startT, endD, endT, afield, nvisit = line.strip().split()

    zeroDateOrig = zeroDateOrig or startD
    zeroDate = zeroDate or strings2datetime(startD, "00:00:00")

    startTimeOrig.append(startT)
    startDate.append((strings2datetime(startD, '00:00:00')-zeroDate).days)
    startTime.append(string2decTime(startT))
    endTime.append(string2decTime(endT))
    dur = string2decTime(endT)-string2decTime(startT)
    duration.append(dur)
    midtime.append(string2decTime(startT)+dur/2.)

    reString = "^[a-zA-Z]+-[A-Za-z0-9]+_(.*)$"
    m = re.search(reString, afield)
    objid.append(afield)
    field.append(m.group(1))
    repeats.append(nvisit)
    color.append(colormap[idx % len(colormap)])
    idx += 1

fontsize = [str(12 - 2*int(len(x)/4))+'pt' for x in field]

schedule_df = pd.DataFrame({
    'startDate': startDate,
    'startTime': startTime,
    'endTime': endTime,
    'duration': duration,
    'objid': objid,
    'field': field,
    'midTime': midtime,
    'repeatsThisPass': repeats,
    'fontsize': fontsize,
    'startTimeOrig': startTimeOrig})

merged = pd.merge(schedule_df, fldPar, on='objid')
dateRange = [schedule_df['startDate'].min()-0.75,
             schedule_df['startDate'].max()+0.75]
timeRange = [schedule_df['endTime'].max()+0.5,
             schedule_df['startTime'].min()-5.5]

data = ColumnDataSource(merged)


p = figure(title='MMT Queue', tools='hover',
           x_range=dateRange, y_range=timeRange)
p.plot_width = 1000
p.grid.grid_line_color = None
p.xaxis.axis_label = "UT Date. Starting Night %s" % zeroDateOrig
p.yaxis.axis_label = "UT Time"

text_props = {
    "source": data,
    "angle": 0,
    "color": "black",
    "text_align": "center",
    "text_baseline": "middle"
}

p.text(x='startDate', y='midTime', text='field', **text_props,
       text_font_size='fontsize')

p.rect('startDate', 'midTime', 1.0, 'duration', source=data, fill_alpha=0.4,
       color="color")


p.select_one(HoverTool).tooltips = [
    ("Field", "@field"),
    ("Obstype", "@obstype"),
    ("PI", "@PI"),
    ("", ""),
    ("Start Time", "@startTimeOrig"),
    ("", ""),
    ("RA", "@ra"),
    ("DEC", "@dec"),
    ("Estimated Mag", "@mag"),
    ("", ""),
    ("Mask", "@mask"),
    ("", ""),
    ("Exposure Time", "@exptime"),
    ("Points per Dither", "@nexp"),
    ("Dither Passes", "@repeatsThisPass"),
    ("", ""),
    ("Dithersize", "@dithersize"),
    ("Filter", "@filter"),
    ("Grisms", "@grism"),
    ("Readtab", "@readtab"),
    ("", ""),
    ("Seeing", "@seeing"),
    ("Photometric", "@photometric")
]


output_file("test.html")
show(p)
