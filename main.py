#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import json
import sys
from package import import_data, dtw_detection, seg_detection, plot_stepdetection, quality, download

# if you need to access a file next to the source code, use the variable ROOT
ROOT = os.path.dirname(os.path.realpath(__file__))

# Save the current CWD
data_WD = os.getcwd()

# Change the CWD to ROOT
os.chdir(ROOT)
            

def print_seg_detection(seg_lim, freq):
    """Dump the phase segmentation computed from the trial. 

    Parameters
    ----------
        seg_lim {dataframe} -- pandas dataframe with phases events 
    """

    seg_lim_dict = {'Start': int(seg_lim[0]),
                    'U-Turn start': int(seg_lim[1]),
                    'U-Turn end': int(seg_lim[2]),
                    'End': int(seg_lim[3])}

    display_dict = {'Start_title': "Trial start",
                    'Start': "{Start}".format(**seg_lim_dict),
                    'Start_sec': "{}".format(round(seg_lim_dict['Start']/freq, 2)),
                    'U-Turn start_title': "U-turn start",
                    'U-Turn start': "{U-Turn start}".format(**seg_lim_dict),
                    'U-Turn start_sec': "{}".format(round(seg_lim_dict['U-Turn start']/freq, 2)),
                    'U-Turn end_title': "U-turn end",
                    'U-Turn end': "{U-Turn end}".format(**seg_lim_dict),
                    'U-Turn end_sec': "{}".format(round(seg_lim_dict['U-Turn end']/freq, 2)),
                    'End_title': "Trial end",
                    'End': "{End}".format(**seg_lim_dict), 
                    'End_sec': "{}".format(round(seg_lim_dict['End']/freq, 2))}
        
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
        

def print_steps_detection(seg_lim_full, seg_lim_corrected, steps_lim, steps_lim_corrected, freq):
    """Dump the trial parameters computed from the gait events detection.  

    Parameters
    ----------
        steps_lim_corrected {dataframe} -- pandas dataframe with gait events after elimination of the extra trial steps
    """

    steps_dict = {"TrialDuration": (seg_lim_full[3] - seg_lim_full[0])/freq, 
                  "LeftGaitCycles": len(steps_lim[steps_lim["Foot"]==0]), 
                  "RightGaitCycles": len(steps_lim[steps_lim["Foot"]==1]), 
                  "WalkingSpeed": (seg_lim_corrected[3] - seg_lim_corrected[0] - seg_lim_corrected[2] + seg_lim_corrected[1])/freq, 
                  "LeftGaitCyclesOk": len(steps_lim_corrected[(steps_lim_corrected["Foot"]==0) & (steps_lim_corrected["Correct"]==1)]), 
                  "RightGaitCyclesOk": len(steps_lim_corrected[(steps_lim_corrected["Foot"]==1) & (steps_lim_corrected["Correct"]==1)])
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
    parser.add_argument('-i0', metavar='data_rf', help='Time series for the right foot sensor.')
    parser.add_argument('-i1', metavar='data_lf', help='Time series for the left foot sensor.')

    
    parser.add_argument('-freq', metavar='freq',
                        help='Acquistion frequency.')
    args = parser.parse_args()

    freq = int(args.freq)

    # quality index tab
    q1 = [0, 0]  # intrinsic detection quality [autocorrelation coefficient, DTW coefficient]
    q2 = [0]  # extrinsic detection quality [right-left step alternation]
    
    # load data
    data_rf = import_data.import_XSens(os.path.join(data_WD, args.i0), freq)
    data_lf = import_data.import_XSens(os.path.join(data_WD, args.i1), freq)
    
    # gait events and steps detection
    steps_lim, q1 = dtw_detection.steps_detection_full(data_rf, data_lf, freq, output=data_WD)

    # quality index
    q2 = quality.correct_detection(steps_lim)
    quality.print_all_quality_index(q1, q2, output=data_WD)

    # print validated gait events and figure 
    print_steps_detection(steps_lim, freq)
    plot_stepdetection.plot_stepdetection(steps_lim, data_rf, data_lf, freq, output=data_WD)
    plot_stepdetection.plot_stepdetection_construction(steps_lim, data_rf, data_lf, freq, output=data_WD, corrected=True)

    # load file to be download
    download.json_report(steps_lim, freq, output=data_WD)

    print("ok charge")
    sys.exit(0)
