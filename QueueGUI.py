import matplotlib
import re
import queueTools
import datetime
import tkinter
import matplotlib.backends.backend_tkagg as tkagg
import matplotlib.patches as mpatches
import matplotlib.figure as pltfig
import pandas as pd
import tableau_colormap
matplotlib.use("TkAgg")

"""Create a user interface for the MMT Queue software.

This gui provides an interface for the queue system to hopefully automate
much of the more annoying aspects of the queue and to great an all in one
interface to work with the queue.

Some goals :
    * Not easy to get into a fubar state
    * Ease of use
    * Ease of updating

The interface will be created in tkinter, so should feel TCL-ish.
"""


def placeholder():
    pass


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

class QueueGUI():
    """Create the main GUI frame."""
    def __init__(self, trimester=None):
        """Initialize the QueueGUI."""
        self.root = tkinter.Tk()

        # Create the Drop Down Menu
        self.build_drop_down_menu()

        # Create the container Frames
        self.build_container_frames()

        # Create the Toolbar
        self.build_left_toolbar(self.leftFrame)

        # On Initialize, set the plot limits to None
        self.plot_limits = [None, None, None, None]

        # Read in the schedule and all the FLD files
        self.read_field_information(trimester=trimester)
        self.read_schedulefile()

        # Initial Plot of Schedule
        self.plot_schedule()

    def set_plot_limits(self, x1, x2, y1, y2):
        """Given an input from the user, change the plot range."""
        self.plot_limits = [x1, x2, y1, y2]

    def build_container_frames(self):
        """Build containters for the components in the frame."""
        # Top frame
        self.topFrame = tkinter.Frame(self.root)
        self.topFrame.grid(row=0)

        # Buttom Frame
        self.botFrame = tkinter.Frame(self.root)
        self.botFrame.grid(row=100, column=0)

        # Left side buttons
        self.leftFrame = tkinter.Frame(self.root)
        self.leftFrame.grid(row=1, column=0)

        # Middle frame
        self.midFrame = tkinter.Frame(self.root)
        self.midFrame.grid(row=1, column=1, columnspan=99, rowspan=99)
        self.fig = None


    def build_drop_down_menu(self):
        """Create Drop Down Menu and submenus."""
        # Create the scope and fill the menu

        self.drop_down_menu = tkinter.Menu(self.root)
        self.root.config(menu=self.drop_down_menu)

        # Create the File menu in the drop down
        fileMenu = tkinter.Menu(self.drop_down_menu)
        self.drop_down_menu.add_cascade(label="File", menu=fileMenu)
        # TODO : Implement Calculate Queue command
        fileMenu.add_command(label='Calculate Queue', command=placeholder)
        fileMenu.add_separator()
        # TODO : Implement Quit Queue Command
        fileMenu.add_command(label='Quit...', command=placeholder)

        # Create a Queue Settings Menu
        qMenu = tkinter.Menu(self.drop_down_menu)
        self.drop_down_menu.add_cascade(label="Queue Settings", menu=qMenu)
        # Implement subcommands
        # TODO Create function for fit days
        qMenu.add_command(label="Edit Days to Fit", command=placeholder)
        # TODO create function for done fields
        qMenu.add_command(label='Edit Previously Obstained Data',
                          command=placeholder)
        # TODO create function for modifying allocation
        qMenu.add_command(label='Edit Allocated Time', command=placeholder)

        # Create a Utilities Menu
        utilMenu = tkinter.Menu(self.drop_down_menu)
        self.drop_down_menu.add_cascade(label='Utilities', menu=utilMenu)
        # Implement subcommands
        # TODO create function for find alternative
        utilMenu.add_command(label="Find Alternate Field", command=placeholder)

        # Create an About Menu
        aboutMenu = tkinter.Menu(self.drop_down_menu)
        self.drop_down_menu.add_cascade(label='About', menu=aboutMenu)
        # TODO Create about function
        aboutMenu.add_command(label='About MMTQueue', command=placeholder)

    def build_left_toolbar(self, frame):
        """Create a toolbar at the left of the window for common functions."""
        # Create the holder
        topToolbar = tkinter.Frame(frame, bg="Blue")

        # TODO function for running the queue from button
        runQueueButton = tkinter.Button(topToolbar, text="Run Queue",
                                        command=placeholder)
        runQueueButton.pack(side='left', padx=2, pady=2)

        topToolbar.pack(side='top', fill='x')

    def read_schedulefile(self):


        schedulefile = 'schedule.dat'
        f = open(schedulefile, 'r')

        startTime = []
        endTime = []
        field = []
        dates = []
        fullfield = []
        nvisit_scheduled = []
        duration = []
        rects = []

        # Read the schedule file
        for line in f.readlines():
            startD, startT, _, \
                endT, afield, nvisit = line.strip().split()

            dates.append(startD)
            startTime.append(string2decTime(startT))
            endTime.append(string2decTime(endT))
            duration.append(endTime[-1]-startTime[-1])
        # Parse the field name down to something readable
            reString = "^[a-zA-Z]+-[A-Za-z0-9]+_(.*)$"
            m = re.search(reString, afield)
            field.append(m.group(1))

            # Last two bookkeepers
            fullfield.append(afield)
            nvisit_scheduled.append(nvisit)
        f.close()

        # find the minimum date
        mindate = strings2datetime(min(dates), '00:00:00')
        dates = [(strings2datetime(x, '00:00:00')-mindate).days
                 for x in dates]

        artist = []

        for ii, date in enumerate(dates):
            rectval = [dates[ii]-0.5, startTime[ii], 1.0, duration[ii]]
            rects.append(rectval)
            artist.append(
                (rectval[0], rectval[1])
            )

        # Create a dataframe
        self.sched_frame = pd.DataFrame({'fieldID': field,
                                         'objid': fullfield,
                                         'startTime': startTime,
                                         'endTime': endTime,
                                         'date': dates,
                                         'rect': rects,
                                         'artist': artist,
                                         'duration': duration,
                                         'nvisit_sched': nvisit_scheduled})

    def click_schedule_block(self, event):
        # Only do this on down click
        if event.mouseevent.button == 1:
            artist = event.artist
            artist.set_visible(False)
            self.plot_schedule(selected_artist=artist.xy)

    def read_field_information(self, trimester=None):
        """Read in all the FLD files."""
        self.fldPar = queueTools.readAllFLDfiles(trimester)
        # Read in the colormap
        colormap = tableau_colormap.tableau20(raw=True)
        color = [tuple_to_hex(colormap[x % len(colormap)])
                 for x in range(len(self.fldPar))]
        self.fldPar['color'] = color

    def plot_schedule(self, selected_artist=None):
        """Update the plotted schedule.

        Selected artist (defaults to none), can be
        set if an item is clicked which will mark it
        (and more later).
        """
        if self.fig is None:
            self.fig = pltfig.Figure(figsize=(5, 4))
            self.ax = self.fig.add_subplot(111)
            # Create the canvas
            self.canvas = tkagg.FigureCanvasTkAgg(self.fig,
                                                  master=self.midFrame)
            self.canvas.show()
            self.canvas.get_tk_widget().pack()
            self.fig.canvas.mpl_connect('pick_event',
                                        self.click_schedule_block)
        # Set the font for labels
        matplotlib.rcParams.update({'font.size': 6})

        # Plot each of the rects
        for ii in range(len(self.sched_frame)):

            rectval = self.sched_frame.loc[ii, 'rect']
            edgecolor = 'none'
            fldMatch = self.fldPar[self.fldPar['objid'] ==
                                   self.sched_frame.loc[ii, 'objid']]
            colorval = fldMatch['color'].values[0]
            if (selected_artist is not None):
                artist = self.sched_frame.loc[ii, 'artist']
                if ((selected_artist[0] == artist[0]) and
                    (abs(selected_artist[1]-artist[1]) == 0)):
                    edgecolor = 'black'

            rect = mpatches.Rectangle(
                (rectval[0], rectval[1]), rectval[2], rectval[3],
                edgecolor=edgecolor,
                facecolor=colorval,
                picker=True,
                linewidth=10,
                alpha=0.1
            )


            self.ax.add_artist(rect)
            self.ax.draw_artist(rect)
            self.ax.text(self.sched_frame.loc[ii, 'date'],
                         self.sched_frame.loc[ii, 'startTime'] +
                         0.5*self.sched_frame.loc[ii, 'duration'],
                         self.sched_frame.loc[ii, 'fieldID'],
                         va='center', ha='center')
        plot_limits = self.plot_limits.copy()

        # Set y limits
        if plot_limits[2] is None:
            plot_limits[2] = self.sched_frame['endTime'].max()+1.0
        if plot_limits[3] is None:
            plot_limits[3] = self.sched_frame['startTime'].min() - 1.0

        # Set x limits
        if plot_limits[0] is None:
            plot_limits[0] = self.sched_frame['date'].min()-0.6
        if plot_limits[1] is None:
            plot_limits[1] = self.sched_frame['date'].max()+0.6

        self.ax.set_xlim(plot_limits[0:2])
        self.ax.set_ylim(plot_limits[2:4])

def main():
    """Code that is run on startup. Create the GUI and go."""
    Queue = QueueGUI()
    Queue.root.mainloop()


if __name__ == "__main__":
    main()
