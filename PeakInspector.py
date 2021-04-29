#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Rafael Barrero Rodriguez"
__credits__ = ["Rafael Barrero Rodriguez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "rbarreror@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# Import modules
import os
import sys
from pathlib import Path
import argparse
import configparser
import logging
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.ticker
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import SingleIntervalTicker, LinearAxis
from bokeh.layouts import gridplot
from bokeh.palettes import d3
import tkinter as tk
from tkinter import font as tkfont
from tkinter import filedialog as fd

import pdb


########################
########################
## Plotting functions ##
########################
########################

class PlotObject():
    """
    PlotObject stores the information used to plot
    """

    def __init__(self):

        self.path = ""  # path to histogram
        self.pathPeakList = "" # path to Peak List
        self.data = None    # pandas dataframe
        self.read = False   # True when dataframe is charged
        self.nPlots = 0     # number of graphs
        self.presentPlot = 0    # number of plot being customized in GUI
        self.plots = {}     # key = Plot section
                            # value = dictionary with 4 keys: columns (list), plotType (list), threshold (list), x_axis (string)
        self.plotSections = []  # List with plot sections
        self.color = 0  
        self.peaks = []     # List of pairs (float, string)
    
    def readData(self, path):
        """
        """
        logging.info("Reading histogram")

        self.path = path
        
        try:
            self.data = pd.read_csv(self.path, sep="\t", float_precision='high')
            self.read = True
        
        except:
            logging.info(f"Unexpected error: {sys.exc_info()[0]}")
            sys.exit()
    

    def getPathFromGUI(self, firstTab):
        """
        """
        self.path = fd.askopenfilename()
        self.readData(self.path)
        
        msg = tk.Label(firstTab, text="File uploaded!", font=('Helvetica', 10))
        msg.place(x=208, y=210)
    

    def readPlotsFromConfig(self, config):
        """
        Search sections containing [Plot x], where x is a number.
        From each section get:
            - columns: List of columns in Y axis
            - plotType: List with type of plot (line, scatter)
            - threshold: List with thresholds 
            - x_axis: String with column in X axis
        
        This dictionary is stored in plots with section name as key
        """

        logging.info("Reading config file")

        self.plotSections = [section for section in config.sections() if re.search(r"^Plot\s?\w*$", section, re.IGNORECASE)]
        
        self.nPlots = len(self.plotSections)

        for section in self.plotSections:
            self.plots[section] = {
                'columns': [i.strip() for i in re.split(r",\s?|;\s?|\s", config.get(section, 'columnName'))],
                'plotType': [i.strip() for i in re.split(r",\s?|;\s?|\s", config.get(section, 'plotType'))],
                'threshold': [float(i.strip()) for i in re.split(r",\s?|;\s?|\s", config.get(section, 'thresholds'))\
                                if re.search(r"^-?\d+(\.\d*)?([eE]-?\d+)?$", i.strip())],
                'x_axis': config.get(section, 'x_axis')
            }


    def guiSelection(self, user_selection, thr1, thr2):
        """
        Parameters selected by user through GUI are stored as dictionary in
        plots. The same structure as readPlotsFromConfig
        """

        # Create name of section, which will be the key
        self.plotSections.append(f"Plot {self.presentPlot}")


        self.plots[self.plotSections[-1]] = {
            "columns": [col for col in user_selection['y_axis'] if user_selection['y_axis'][col]["selected"].get()],
            "plotType": [user_selection['y_axis'][col]['type'].get() \
                for col in user_selection['y_axis'] if user_selection['y_axis'][col]["selected"].get()],
            "threshold": [float(thr) for thr in [thr1, thr2] if re.search(r"^-?\d+(\.\d*)?([eE]-?\d+)?$", thr)],
            "x_axis": user_selection['x_axis'].get()
        }


    def getPeaks(self, path):
        '''
        Input:
            - path: String containing path to peaks list
        Effect:
            - self.peaks: List of pairs. The first element of each pair is a float with the theoretical DM and the
            second element is its associated name
        '''

        logging.info("Reading peaks list")

        try:
            df_peaks = pd.read_csv(path, sep="\t", float_precision="high")
        
        except:
            logging.info(f"Error reading peaks list: {sys.exc_info()[0]}")
            sys.exit()

        self.peaks = [[float(dm), str(name)] for dm, name in zip(df_peaks.loc[:, 'DM'].to_list(), df_peaks.loc[:, 'Name'].to_list())]

    
    def getPeakListPathFromGUI(self, firstTab):
        """
        """
        self.pathPeakList = fd.askopenfilename()
        self.getPeaks(self.pathPeakList)
        
        msg = tk.Label(firstTab, text="File uploaded!", font=('Helvetica', 10))
        msg.place(x=208, y=325)


    def reset(self):
        """
        Reset object, when graph was plotted
        """
        self.__init__()
        


def plot_bottom_graph(main_plot, section):
    '''
    Represent bottom graph using main_plot as reference. It will use its
    x axis, so they are coupled. The function will return the bottom plot figure
    object.
    '''

    # Build bottom plot
    bottom_plot = figure(title=f"{section}",\
         x_axis_label=plotObject.plots[section]['x_axis'], \
         width=1300, height=400, x_range=main_plot.x_range,\
         tools = "pan,xzoom_in,xzoom_out,ywheel_zoom,box_zoom,reset,save,undo,hover", tooltips=[("Name", "$name")])
		
    bottom_plot = addPlotsToFigure(bottom_plot, section)

    return bottom_plot


def plot_pleak(figure, section):
    '''
    Plot list of peaks given by the user
    '''

    # Get all values from columns plotted by user to find minimum value
    all_values = [plotObject.data.loc[:, col] for col in plotObject.plots[section]['columns']]
    
    all_values = [j for i in all_values for j in i if not pd.isna(j)]

    min_value = np.min(all_values)
    min_value = min_value - abs(0.1*min_value)
    max_value = np.max(all_values)*1.1

    y_axis = (min_value, 0, max_value)
    for peak, peak_name in plotObject.peaks:
        x_axis = np.ones_like(y_axis)*peak
        figure.line(x_axis, y_axis, line_color='green', line_width=2, name=peak_name)
    
    return figure


def plot_threshold(figure, section):
    '''
    Plot thresholds
    '''

    min_mz = np.min(plotObject.data.loc[:, plotObject.plots[section]['x_axis']])
    max_mz = np.max(plotObject.data.loc[:, plotObject.plots[section]['x_axis']])

    x_axis = (min_mz, max_mz)

    for threshold_i in plotObject.plots[section]['threshold']:
        y_axis = np.ones_like(x_axis)*threshold_i
        figure.line(x_axis, y_axis, line_color='black', line_dash="4 4", name="Threshold")

    return figure


def addPlotsToFigure(figure, section):
    """
    Add lines or scatter to the figure
    """

    try:
        figure.xaxis.ticker.desired_num_ticks = 30
	
    except AttributeError as err:
        logging.info(f"AttributeError: {err}")
        print("bokeh package needs to be updated (pip install bokeh -U)")
        sys.exit()


    for i, column in enumerate(plotObject.plots[section]["columns"]):

        # assert that it is present in dataframe
        if column not in plotObject.data.columns:
            logging.info(f"Column not found: {column}")
            continue
        

        # get boolean with rows without na. Only those will be plotted
        bool_not_na = (~pd.isna(plotObject.data.loc[:, column])).to_list()

        if re.search(r"^line$", plotObject.plots[section]['plotType'][i], re.IGNORECASE):
            # If line...
            figure.line(plotObject.data.loc[bool_not_na, plotObject.plots[section]['x_axis']],
                        plotObject.data.loc[bool_not_na, column],
                        line_width=2,
                        color=d3['Category20'][20][plotObject.color],
                        legend_label=column,
                        name = column)

        elif re.search(r"^scatter$", plotObject.plots[section]['plotType'][i], re.IGNORECASE):
            # If scatter...
            figure.circle(plotObject.data.loc[bool_not_na, plotObject.plots[section]['x_axis']],
                        plotObject.data.loc[bool_not_na, column],
                        size=1,
                        color=d3['Category20'][20][plotObject.color],
                        legend_label=column,
                        name = column)
        
        # Change color for the next plot
        plotObject.color += 2 if plotObject.color < 18 else 17
    
    # plot threshold
    figure = plot_threshold(figure, section)

    # plot peaks of interest
    figure = plot_pleak(figure, section)

    return figure


def plot_graphs():
    '''
    Genera function for plotting graphs
    '''
    
    logging.info('Plotting graphs')

    # Save graphs in html file
    output_file(os.path.splitext(plotObject.path)[0] + '_plot.html')
    
    # Build the first plot
    p1 = figure(title="Plot 1",\
         x_axis_label=plotObject.plots[plotObject.plotSections[0]]['x_axis'],\
         width=1300, height=400, tools = "pan,xzoom_in,xzoom_out,ywheel_zoom,box_zoom,reset,save,undo,hover", tooltips=[('Name', '$name')])
    
    p1 = addPlotsToFigure(p1, plotObject.plotSections[0])

    # If there are more plots, these are represented below using plot_bottom_graph function
    if plotObject.nPlots > 1:

        bottom_graphs_list = [[plot_bottom_graph(p1, section)] for section in plotObject.plotSections[1:]]
        bottom_graphs_list = [plot for plot in bottom_graphs_list if plot != [""]]

        all_graphs_list = [[p1]] + bottom_graphs_list

        # Show the plot
        plot = gridplot(all_graphs_list)
        show(plot)
        # savePlot(plot)

    else: 
        show(p1)

    
    
    logging.info('Graphs plotted')

    return 0
    

####################
####################
## GUI functions ##
####################
####################

def showTab(root, container):
    """
    GUI Tab to customize plot
    Columns showed are taken from dataframe
    """
    # Trial with one extra plottab
    trialTab = tk.Frame(container)
    trialTab.grid(row=0, column=0, sticky="nsew")

    label_TT_1 = tk.Label(trialTab, text=f"Plot {plotObject.presentPlot}", font=("comicsans", 16))
    label_TT_1.pack(side="top", fill="x", pady=10) 

    # Get column names with float or integers
    column_names = [col for col in plotObject.data.columns if pd.api.types.is_numeric_dtype(plotObject.data.loc[:, col])]

    ########################
    # USER PLOTS SELECTION #
    ########################

    font_selections = ("comicsans", 10)
    checkbox_x = 100
    checkbox_y  = 130

    label_TT_12 = tk.Label(trialTab, text="Select columns to be plotted", font=("comicsans", 12))
    label_TT_12.place(relx=0.30, y=checkbox_y-70)

    label_TT_13 = tk.Label(trialTab, text="X", font=("comicsans", 11))
    label_TT_13.place(x=checkbox_x-6, y=checkbox_y-30)

    label_TT_14 = tk.Label(trialTab, text="Y", font=("comicsans", 11))
    label_TT_14.place(x=checkbox_x+34, y=checkbox_y-30)

    label_TT_15 = tk.Label(trialTab, text="Column", font=("comicsans", 11))
    label_TT_15.place(x=checkbox_x+68, y=checkbox_y-30)

    # Loop over each possible column, to show buttons...
    # store user selection in a dictionary of dictionaries. Each one stores column info
    user_selection = {'x_axis': tk.StringVar(value='midpoint'), 'y_axis': {}}

    for i, col in enumerate(column_names):
        user_selection['y_axis'][col] = {
            'selected': tk.BooleanVar(), 
            'type':tk.StringVar(value='line')
            }

        # Radio Button X axis #
        radioButton = tk.Radiobutton(trialTab, variable=user_selection['x_axis'], value=col, font=font_selections)
        radioButton.pack()
        radioButton.place(x=checkbox_x-10, y=checkbox_y+30*i)

        # CheckBox Y axis #
        col_check = tk.Checkbutton(trialTab, text=f"    {col}", var=user_selection['y_axis'][col]['selected'], font=font_selections)
        col_check.place(x=checkbox_x+30, y=checkbox_y+30*i)

        # Radio Button Line #
        radioButton_1 = tk.Radiobutton(trialTab, text="Line", variable=user_selection['y_axis'][col]['type'], value='line', font=font_selections)
        radioButton_1.pack()
        radioButton_1.place(x=checkbox_x+160, y=checkbox_y+30*i)

        # Radio Button Scatter #
        radioButton_2 = tk.Radiobutton(trialTab, text="Scatter", variable=user_selection['y_axis'][col]['type'], value='scatter', \
            font=font_selections)
        radioButton_2.pack()
        radioButton_2.place(x=checkbox_x+220, y=checkbox_y+30*i)


    #############################
    # USER THRESHOLDS SELECTION #
    #############################
    y_label = checkbox_y+30*len(column_names) + 20
    y_entry = y_label + 30

    # Create second section with entry spaces to introduce thresholds
    label_TT_2 = tk.Label(trialTab, text="Define Thresholds (e.g. 5e3)", font=("comicsans", 11))
    label_TT_2.place(relx=0.3, y=y_label)
  

    # Threshold 1 #
    ###############

    label_threshold_1 = tk.Label(trialTab, text="Threshold 1", font=font_selections)
    label_threshold_1.place(relx=0.3, y=y_label+30)

    threshold_entry1 = tk.Entry(trialTab, width=10, justify="center")
    threshold_entry1.place(relx=0.5, y=y_label+30)

    
    # Threshold 2 #
    ###############

    label_threshold_2 = tk.Label(trialTab, text="Threshold 2", font=font_selections)
    label_threshold_2.place(relx=0.3, y=y_label+60)

    threshold_entry2 = tk.Entry(trialTab, width=10, justify="center")
    threshold_entry2.place(relx=0.5, y=y_label+60)

    ################
    # PLOT OR NEXT #
    ################

    if plotObject.nPlots == plotObject.presentPlot:
        # if it is the last plot, show "Plot" option
        plot_button = tk.Button(trialTab, text="Plot", font=font_selections, pady=5, padx=10, bd=2,
                                command=lambda: moveToPlot(root, container, user_selection, threshold_entry1.get(), threshold_entry2.get()))
        plot_button.pack()
        plot_button.place(relx=0.40, y=y_label+100, height=40, width=120)
    
    else:
        # else, show "Next" option
        plot_button = tk.Button(trialTab, text="Next", font=font_selections, pady=5, padx=10, bd=2,
                                command=lambda: moveToNext(root, container, user_selection, threshold_entry1.get(), threshold_entry2.get()))
        plot_button.pack()
        plot_button.place(relx=0.40, y=y_label+100, height=40, width=120)

    trialTab.tkraise()


def moveToPlot(root, container, user_selection, threshold_entry1, threshold_entry2):
    """
    Plot button was pressed...
    """
    # Store in plotObject user parameters
    plotObject.guiSelection(user_selection, threshold_entry1, threshold_entry2)
   
    # plot graphs
    plot_graphs()

    # show main
    showMain(root, container)


def moveToNext(root, container, user_selection, threshold_entry1, threshold_entry2):
    """
    Store user selected parameters and show
    next tab
    """
    # Store in plotObject user parameters
    plotObject.guiSelection(user_selection, threshold_entry1, threshold_entry2)
    
    # Show next tab
    plotObject.presentPlot += 1
    showTab(root, container)


def showMain(root, container):
    """
    Home tab
    """
    # reset plotObject
    plotObject.reset()

    firstTab = tk.Frame(container)
    firstTab.grid(row=0, column=0, sticky="nsew")
    
    label_FT_1 = tk.Label(firstTab, text="PeakInspector", font=('Helvetica', 18))
    label_FT_1.pack(side="top", fill="x", pady=15)

    # Number of plots #
    label_FT_1 = tk.Label(firstTab, text="Enter number of plots", font=root.title_font)
    label_FT_1.pack(side="top", fill="x", pady=10)
    
    entry_FT_1 = tk.Entry(firstTab, width=3, justify="center")
    entry_FT_1.pack()

    # Input File #
    label_FT_1 = tk.Label(firstTab, text="Select Input Table", font=root.title_font)
    label_FT_1.pack(side="top", fill="x", pady=20)

    button_FT_1 = tk.Button(firstTab, text="Click to select file", command=lambda: plotObject.getPathFromGUI(firstTab))
    button_FT_1.pack()
    button_FT_1.place(x=200, y=180)

    # File with peaks #
    label_FT_1 = tk.Label(firstTab, text="Select Peak File", font=root.title_font)
    label_FT_1.pack(side="top", fill="x", pady=60)

    button_FT_1 = tk.Button(firstTab, text="Click to select file", command=lambda: plotObject.getPeakListPathFromGUI(firstTab))
    button_FT_1.pack()
    button_FT_1.place(x=200, y=290)

    # Go to next operation
    button_FT_2 = tk.Button(firstTab, text="Customize Plots", pady=5, width=20,
                            command=lambda: startCustom(entry_FT_1.get(), root, container))
    button_FT_2.pack()
    button_FT_2.place(x=175, y=390)

    firstTab.tkraise()


def isNumber(entry):
    """
    If user entered a number (number of plots), it is accepted
    """
    if re.search("^[1-5]$", entry):
        return True
    else:
        return False


def startCustom(nPlots, root, container):
    """
    If user entered a number, go to next plot
    """
    if isNumber(nPlots) and plotObject.read:
        plotObject.nPlots = int(nPlots)
        plotObject.presentPlot += 1
        showTab(root, container)


def startGUI():
    '''
    Execution using GUI
    '''

    ###############
    # CREATE ROOT #
    ###############

    root = tk.Tk(className="PeakInspector")
    root.title_font = tkfont.Font(family='Helvetica', size=12)
    root.geometry("500x550")
    root.frames = {}

    ##################
    # MAIN CONTAINER #
    ##################

    container = tk.Frame(root)
    container.pack(side="top", fill="both", expand=True)
    container.grid_rowconfigure(0, weight=1)
    container.grid_columnconfigure(0, weight=1)

    # Main Page #
    #############
    showMain(root, container)

    root.mainloop()


####################
####################
## Main functions ##
####################
####################

def main(args):
    """
    main function
    """

    if args.gui:
        # if user selected gui, execute it...
        startGUI()
    
    else:
        # otherwise, get values from config.ini
        
        # read dataframe
        plotObject.readData(args.infile)
        
        # get peaks of interest
        plotObject.getPeaks(args.peaks)

        plotObject.readPlotsFromConfig(config)
        plot_graphs()



if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PeakInspector',
        epilog='''
        Example:
            python PeakInspector.py
        
        '''
    )

    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakInspector.ini")

    parser.add_argument('-i', '--infile', help='Path to input file', type=str)
    parser.add_argument('-p', '--peaks', help='Path to peaks list', type=str)
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file', type=str)
    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-gui', action='store_true', help='Read data from GUI', default=False)

    args = parser.parse_args()

    # If user use config file, parse it
    if args.config and not args.gui:
        config  = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(args.config)

    if not args.infile and not args.gui:
        args.infile = config['Parameters']['infile']
    
    if not args.peaks and not args.gui:
        args.peaks = config['Parameters']['peaksList']

    # logging debug level. By default, info level
    if args.infile:
        log_file = outfile = args.infile[:-4] + '_log.txt'
        log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
    
    else:
        log_file = outfile = 'log.txt'
        log_file_debug = outfile = 'log_debug.txt'


    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])
    

    # global plot object
    plotObject = PlotObject()

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')
