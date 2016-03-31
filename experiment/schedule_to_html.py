import ipdb
"""
Generate an interactive BarChart for the MMT Queue.

Based on the Bokeh library for python.

Assumes the schedule resides in schedule.dat
"""

from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import figure, show, output_file
import pandas as pd
import datetime
from tableau_colormap import tableau20
import MMTQueue
import math
import ipdb
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


def tuple_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


fldPar = MMTQueue.read_all_fld_files('2016a')
# The ephemeris kills things later, drop it
fldPar.drop('ephem', 1, inplace=True)

schedule_df = pd.read_csv('schedule.csv', parse_dates=['start_time'])
# We will use these colors to color each group by a different color
colormap = tableau20(raw=True)
color = [tuple_to_hex(colormap[x % len(colormap)]) for x in range(len(fldPar))]
fldPar['color'] = color

schedule_df['end_time'] = [x + datetime.timedelta(seconds=float(y))
                           for x, y in zip(schedule_df['start_time'],
                                           schedule_df['duration'])]
schedule_df['mid_time'] = [x + datetime.timedelta(seconds=float(0.5*y))
                           for x, y in zip(schedule_df['start_time'],
                           schedule_df['duration'])]
schedule_df['mid_time'] = [x.hour + x.minute/60.
                           for x in schedule_df['mid_time']]
print(schedule_df['mid_time'])
mindate = schedule_df['start_time'].min()
schedule_df['start_date'] = [
    math.floor((x - mindate).total_seconds()/24.0/3600.0)
    for x in schedule_df['start_time']]


dateRange = [schedule_df['start_date'].min()-1,
             schedule_df['start_date'].max()+1]
timeRange = [max(schedule_df['end_time'].dt.hour+2),
             min(schedule_df['start_time'].dt.hour-2)]
schedule_df['duration'] /= 3600.
reString = "^[a-zA-Z]+-[A-Za-z0-9]+_(.*)$"
schedule_df['field'] = [re.search(reString, x).group(1)
                        for x in schedule_df['objid']]
schedule_df['fontsize'] = [str(12 - 2*int(len(x)/4))+'pt' for x in schedule_df['field']]

merged = pd.merge(schedule_df, fldPar, on='objid')
data = ColumnDataSource(merged)
p = figure(title='MMT Queue', tools='hover',
           x_range=dateRange, y_range=timeRange)
p.plot_width = 1000
p.grid.grid_line_color = None
p.xaxis.axis_label = "UT Date."
p.yaxis.axis_label = "UT Time"

text_props = {
    "source": data,
    "angle": 0,
    "color": "black",
    "text_align": "center",
    "text_baseline": "middle"
}

p.text(x='start_date', y='mid_time', text='field', **text_props,
       text_font_size='fontsize')


p.rect('start_date', 'mid_time', 1.0, 'duration', source=data, fill_alpha=0.4,
       color="color")


# p.select_one(HoverTool).tooltips = [
#     ("Field", "@field"),
#     ("Obstype", "@obstype"),
#     ("PI", "@PI"),
#     ("", ""),
#     ("Start Time", "@startTimeOrig"),
#     ("", ""),
#     ("RA", "@ra"),
#     ("DEC", "@dec"),
#     ("Estimated Mag", "@mag"),
#     ("", ""),
#     ("Mask", "@mask"),
#     ("", ""),
#     ("Exposure Time", "@exptime"),
#     ("Points per Dither", "@nexp"),
#     ("Dither Passes", "@repeatsThisPass"),
#     ("", ""),
#     ("Dithersize", "@dithersize"),
#     ("Filter", "@filter"),
#     ("Grisms", "@grism"),
#     ("Readtab", "@readtab"),
#     ("", ""),
#     ("Seeing", "@seeing"),
#     ("Photometric", "@photometric")
# ]


output_file("test.html")
show(p)
