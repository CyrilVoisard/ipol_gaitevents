#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
from package import import_data, dtw_detection, plot_stepdetection, quality, download

# if you need to access a file next to the source code, use the variable ROOT
ROOT = os.path.dirname(os.path.realpath(__file__))

# Save the current CWD
data_WD = os.getcwd()

# Change the CWD to ROOT
os.chdir(ROOT)
            

def print_steps_detection(steps_lim, t_full, freq):
    """Dump the trial parameters computed from the gait events detection.  

    Parameters
    ----------
        steps_lim {dataframe} -- pandas dataframe with gait events
        t_full {int} -- total number of samples in the trial
        freq {int} -- acquisition frequency in Herz
    """

    steps_dict = {"TrialDuration": t_full/freq, 
                  "LeftGaitCycles": len(steps_lim[steps_lim["Foot"]==0]), 
                  "RightGaitCycles": len(steps_lim[steps_lim["Foot"]==1])
                 }

    display_dict = {'Raw': "Raw data",
                    'Number': "Number of footsteps:",
                    'TrialDuration': "Trial duration (s): {TrialDuration}".format(**steps_dict),
                    'LeftGaitCycles': "    - Left foot: {LeftGaitCycles}".format(**steps_dict),
                    'RightGaitCycles': "    - Right foot: {RightGaitCycles}".format(**steps_dict)
                    }
    info_msg = """
    {Raw:^30}
    ------------------------------
    {TrialDuration:<30}
    {Number:<30}
    {LeftGaitCycles:<30}
    {RightGaitCycles:<30}
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
    q2 = quality.compute_extrinsic_quality(steps_lim)
    quality.print_all_quality_index(q1, q2, output=data_WD)

    # print validated gait events and figure 
    print_steps_detection(steps_lim, len(data_rf), freq)
    plot_stepdetection.plot_stepdetection(steps_lim, data_rf, data_lf, freq, output=data_WD)
    plot_stepdetection.plot_stepdetection_construction(steps_lim, data_rf, data_lf, freq, output=data_WD)

    # load file to be download
    download.json_report(steps_lim, freq, output=data_WD)

    print("ok charge")
    sys.exit(0)
