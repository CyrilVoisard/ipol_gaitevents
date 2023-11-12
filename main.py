#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, filtfilt
from scipy import interpolate
import sys

from package import import_data, compute_semio_val, radar_design

# if you need to access a file next to the source code, use the variable ROOT
ROOT = os.path.dirname(os.path.realpath(__file__))

# Save the current CWD
data_WD = os.getcwd()

# Change the CWD to ROOT
os.chdir(ROOT)
        

def print_semio_criteria(criteria_dict):
    """Dump the parameters computed from the trial in a text file (trial_info.txt)

    Parameters
    ----------
    parameters_dict : dict
        Parameters of the trial.
    """

    display_dict = {'Average Speed': "Average Speed: {Average Speed}".format(**criteria_dict),
                    'Springiness': "Springiness: {Springiness}".format(**criteria_dict),
                    'Sturdiness': "Sturdiness: {Sturdiness}".format(**criteria_dict),
                    'Smoothness': "Smoothness: {Smoothness}".format(**criteria_dict),
                    'Steadiness': "Steadiness: {Steadiness}".format(**criteria_dict),
                    'Stability': "Stability: {Stability}".format(**criteria_dict),
                    'Symmetry': "Symmetry: {Symmetry}".format(**criteria_dict),
                    'Synchronisation': "Synchronisation: {Synchronisation}".format(**criteria_dict)
                    }
    info_msg = """
    Z-Scores
    --------------------------------------------------+--------------------------------------------------
    {Average Speed:<50}| {Steadiness:<50}
    {Springiness:<50}| {Stability:<50}
    {Sturdiness:<50}| {Symmetry:<50}
    {Smoothness:<50}| {Synchronisation:<50}
    """

    # Dump information
    os.chdir(data_WD) # Get back to the normal WD

    with open("trial_criteria.txt", "wt") as f:
        print(info_msg.format(**display_dict), file=f)


def print_quality_index(steps_lim_full, seg_lim):
    """Dump the parameters computed from the trial in a text file (trial_info.txt)

    Parameters
    ----------
    parameters_dict : dict
        Parameters of the trial.
    """

    display_dict = {'Average Speed': "Average Speed: {Average Speed}".format(**criteria_dict),
                    'Springiness': "Springiness: {Springiness}".format(**criteria_dict),
                    'Sturdiness': "Sturdiness: {Sturdiness}".format(**criteria_dict),
                    'Smoothness': "Smoothness: {Smoothness}".format(**criteria_dict),
                    'Steadiness': "Steadiness: {Steadiness}".format(**criteria_dict),
                    'Stability': "Stability: {Stability}".format(**criteria_dict),
                    'Symmetry': "Symmetry: {Symmetry}".format(**criteria_dict),
                    'Synchronisation': "Synchronisation: {Synchronisation}".format(**criteria_dict)
                    }
    info_msg = """
    Z-Scores
    --------------------------------------------------+--------------------------------------------------
    {Average Speed:<50}| {Steadiness:<50}
    {Springiness:<50}| {Stability:<50}
    {Sturdiness:<50}| {Symmetry:<50}
    {Smoothness:<50}| {Synchronisation:<50}
    """

    # Dump information
    os.chdir(data_WD) # Get back to the normal WD

    with open("trial_criteria.txt", "wt") as f:
        print(info_msg.format(**display_dict), file=f)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description='Return a semiogram for a given trial.')
    parser.add_argument('-i0', metavar='data_lb', help='Time series for the lower back sensor.')
    parser.add_argument('-i1', metavar='data_rf', help='Time series for the right foot sensor.')
    parser.add_argument('-i2', metavar='data_lf', help='Time series for the left foot sensor.')
    
    parser.add_argument('-freq', metavar='freq',
                        help='Acquistion frequency.')
    args = parser.parse_args()

    freq = int(args.freq)
    
    # load data (only lower back in this demo)
    data_lb = import_data.import_XSens(os.path.join(data_WD, args.i0))
    data_rf = import_data.import_XSens(os.path.join(data_WD, args.i1))
    data_lf = import_data.import_XSens(os.path.join(data_WD, args.i2))
    
    # gait events and steps detection
    steps_lim_full = dtw_detection.steps_detection(data_rf, data_lb, freq)
    
    # phase boundaries detection and figure
    seg_lim = seg_detection.seg_detection(data_lb, steps_lim_full, freq)

    # quality index and 
    qi, steps_lim_corrected = print_quality_index(steps_lim_full, seg_lim)

    # print phases and figure
    print_seg_detection(seg_lim)
    seg_detection.plot_seg_detection(seg_lim, data_lb, freq)

    # print validated gait events and figure 
    print_steps_detection(steps_lim_corrected)
    dtw_detection.plot_steps_detection(steps_lim_corrected, data_rf, data_lf, freq, corrected=True)
    



    
    
    parameters_dict = dict(zip(parameters_names, parameters))
    print_semio_parameters(parameters_dict)

    criteria_dict = dict(zip(criteria_names, criteria))
    print_semio_criteria(criteria_dict)

    # semiogram design
    radar_design.new_radar_superpose({"unique": criteria}, min_r=int(args.min_z), max_r=int(args.max_z), output=data_WD, name="semio")
    if compare : 
        radar_design.new_radar_superpose({"unique": ref_criteria}, min_r=int(args.min_z), max_r=int(args.max_z), output=data_WD, name="semio_ref")
        radar_design.new_radar_superpose({"ref": ref_criteria, "new": criteria}, min_r=int(args.min_z), max_r=int(args.max_z), output=data_WD, name="semio_sup")
    print("ok charge")
    sys.exit(0)
