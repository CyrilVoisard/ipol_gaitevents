#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from scipy.signal import butter, filtfilt
from scipy import interpolate
import sys

from package import import_data, dtw_detection, seg_detection, plot_stepdetection

# if you need to access a file next to the source code, use the variable ROOT
ROOT = os.path.dirname(os.path.realpath(__file__))

# Save the current CWD
data_WD = os.getcwd()

# Change the CWD to ROOT
os.chdir(ROOT)


def print_quality_index(steps_lim_full, seg_lim):
    """Compute the quality index of the trial gait events detection (between 0 and 100) and produce a picture of the number surrounded by an appropriately colored circle. 
    Add quality index formula ? 

    Parameters
    ----------
        steps_lim_full {dataframe} -- pandas dataframe with all the detected gait events
        seg_lim {dataframe} -- pandas dataframe with phases events 

    Returns
    -------
        qi {int} -- quality index 
        steps_lim_corrected {dataframe} -- pandas dataframe with gait events after elimination of the extra trial steps
    """
    
    steps_lim_corrected = steps_lim_full
    qi = 50
    max_qi=100
    xval = np.arange(0, 2*np.pi*(.05+0.90*(qi/max_qi)), 0.01)
    colormap = plt.get_cmap("RdYlGn")
    norm = mpl.colors.Normalize(0.0, 2*np.pi)
    f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4),subplot_kw=dict(projection='polar'))
    #Scatter version
    yval = np.ones_like(xval)
    ax.scatter(xval, yval, c=xval, s=300, cmap=colormap, norm=norm, linewidths=0)

    ax.set_axis_off()
    ax.set_ylim(0,1.5)
    if qi<10:
        ax.annotate(qi,xy=( 1.25*np.pi, .3),color=colormap(.05+0.90*(qi/max_qi)),fontsize=50)
    else :
        ax.annotate(qi,xy=( 1.18*np.pi, .5),color=colormap(.05+0.90*(qi/max_qi)),fontsize=50)

    path_out = os.path.join(data_WD, "quality_index.svg")
    plt.savefig(path_out, dpi=80,
                    transparent=True, bbox_inches="tight")
    
    return qi, steps_lim_corrected
            

def print_seg_detection(seg_lim, freq):
    """Dump the phase segmentation computed from the trial. 

    Parameters
    ----------
        seg_lim {dataframe} -- pandas dataframe with phases events 
    """

    seg_lim_dict = {'Start': seg_lim[0],
                    'U-Turn start': seg_lim[1],
                    'U-Turn end': seg_lim[2],
                    'End': seg_lim[3]}

    display_dict = {'Start_title': "Trial start",
                    'Start': "{Start}".format(**seg_lim_dict),
                    'Start_sec': "{}".format(round(seg_lim_dict['Start']/100, 2)),
                    'U-Turn start_title': "U-turn start",
                    'U-Turn start': "{U-Turn start}".format(**seg_lim_dict),
                    'U-Turn start_sec': "{}".format(round(seg_lim_dict['U-Turn start']/100, 2)),
                    'U-Turn end_title': "U-turn end",
                    'U-Turn end': "{U-Turn end}".format(**seg_lim_dict),
                    'U-Turn end_sec': "{}".format(round(seg_lim_dict['U-Turn end']/100, 2)),
                    'End_title': "Trial end",
                    'End': "{End}".format(**seg_lim_dict), 
                    'End_sec': "{}".format(round(seg_lim_dict['End']/100, 2))}
        
    info_msg = """
    {Start_title:<15}| {Start:<10}| {Start_sec:<10}
    {U-Turn start_title:<15}| {U-Turn start:<10}| {U-Turn start_sec:<10}
    {U-Turn end_title:<15}| {U-Turn end:<10}| {U-Turn end_sec:<10}
    {End_title:<15}| {End:<10}| {End_sec:<10}
    """

    # Dump information
    os.chdir(data_WD) # Get back to the normal WD

    with open("seg_lim.txt", "wt") as f:
        print(info_msg.format(**display_dict), file=f)
        

def print_steps_detection(steps_lim_corrected):
    """Dump the trial parameters computed from the gait events detection.  

    Parameters
    ----------
        steps_lim_corrected {dataframe} -- pandas dataframe with gait events after elimination of the extra trial steps
    """

    steps_dict = {"TrialDuration": 1000, 
                  "LeftGaitCycles": 1000, 
                  "RightGaitCycles": 1000, 
                  "WalkingSpeed": 1000, 
                  "LeftGaitCyclesOk": 1000, 
                  "RightGaitCyclesOk": 1000
                 }

    display_dict = {'Raw': "Raw data",
                    'Corrected': "Corrected data",
                    'Number': "Number of footsteps:",
                    'TrialDuration': "Trial duration (s): {TrialDuration}".format(**steps_dict),
                    'LeftGaitCycles': "    - Left foot: {LeftGaitCycles}".format(**steps_dict),
                    'RightGaitCycles': "    - Right foot: {RightGaitCycles}".format(**steps_dict),
                    'WalkingSpeed': "WalkingSpeed (m/s): {WalkingSpeed}".format(**steps_dict),
                    'LeftGaitCyclesOk': '    - Left foot: {LeftGaitCyclesOk}'.format(**steps_dict),
                    'RightGaitCyclesOk': '    - Right foot: {RightGaitCyclesOk}'.format(**steps_dict)
                    }
    info_msg = """
    {Raw:^30}|{Corrected:^30}
    ------------------------------+------------------------------
    {TrialDuration:<30}| {WalkingSpeed:<30}
    {Number:<30}| Number of validated footsteps:
    {LeftGaitCycles:<30}| {LeftGaitCyclesOk:<30}
    {RightGaitCycles:<30}| {RightGaitCyclesOk:<30}
    """

    # dump information
    os.chdir(data_WD) # Get back to the normal WD

    with open("gait_events.txt", "wt") as f:
        print(info_msg.format(**display_dict), file=f)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Return a semiogram for a given trial.')
    parser.add_argument('-i0', metavar='data_lb', help='Time series for the lower back sensor.')
    parser.add_argument('-i1', metavar='data_rf', help='Time series for the right foot sensor.')
    parser.add_argument('-i2', metavar='data_lf', help='Time series for the left foot sensor.')

    
    parser.add_argument('-freq', metavar='freq',
                        help='Acquistion frequency.')
    args = parser.parse_args()

    freq = int(args.freq)
    
    # load data
    data_lb = import_data.import_XSens(os.path.join(data_WD, args.i0))
    data_rf = import_data.import_XSens(os.path.join(data_WD, args.i1))
    data_lf = import_data.import_XSens(os.path.join(data_WD, args.i2))
    
    # gait events and steps detection
    steps_rf, steps_lf, steps_lim_full = dtw_detection.steps_detection_full(data_rf, data_lf, freq)
    
    # phase boundaries detection and figure
    seg_lim = seg_detection.seg_detection(data_lb, steps_lim_full, freq)

    # quality index and 
    qi, steps_lim_corrected = print_quality_index(steps_lim_full, seg_lim)

    # print phases and figure
    print_seg_detection(seg_lim, freq)
    seg_detection.plot_seg_detection(seg_lim, steps_lim_full, data_lb, freq, output=data_WD)

    # print validated gait events and figure 
    print_steps_detection(steps_lim_corrected)
    #plot_stepdetection.plot_stepdetection(steps_lim_corrected, data_rf, data_lf, freq, output=data_WD, corrected=True)
    plot_stepdetection.plot_stepdetection(steps_rf, steps_lf, data_rf, data_lf, freq, output=data_WD, corrected=True)

    print("ok charge")
    sys.exit(0)
