import os
import numpy as np
import pandas as pd
from scipy import stats
from tslearn import metrics

from package import find_stride, deal_stride, plot_stepdetection


def steps_detection_full(data_rf, data_lf, freq, output):
    """Detection of all gait events.  

    Parameters
    ----------
        data_rf {dataframe} -- pandas dataframe with data from the right foot sensor
        data_lf {dataframe} -- pandas dataframe with data from the left foot sensor
        freq {int} -- acquisition frequency 
        output {str} -- folder path for output fig

    Returns
    -------
       steps_lim {dataframe} -- pandas dataframe with the detected gait events
       q {int} -- extrinsic quality index 
    """
    
    steps_rf, q_rf = steps_detection(data_rf, data_lf, 1, freq, output)
    steps_lf, q_lf = steps_detection(data_lf, data_rf, 0, freq, output)
    
    full = np.concatenate((steps_rf, steps_lf))
    steps_lim = pd.DataFrame(full, columns=["Foot", "Phase", "HO", "TO", "HS", "FF", "Score"])
    q = [min(q_rf[0], q_lf[0]), min(q_rf[1], q_lf[1])]

    return steps_lim, q


def steps_detection(data_1, data_2, foot, freq, output):
    """Detection of gait events for one foot.  

    Parameters
    ----------
        data_1 {dataframe} -- pandas dataframe with data from the considered foot sensor
        data_2 {dataframe} -- pandas dataframe with data from the other foot sensor
        foot {int} -- considered foot
        freq {int} -- acquisition frequency 
        output {str} -- folder path for output fig

    Returns
    -------
       steps_side {array} -- numpy array with the detected gait events
       q_side {int} -- extrinsic quality index for
    """

    # data 
    x = data_1["Gyr_Y"]
    z = deal_stride.calculate_jerk_tot(data_1, freq)

    # template for the template based detection
    gyr_ref, jerk_ref, stride_annotations_ref, q_side = find_stride.annotate_ref_stride(data_1, data_2, foot, freq, output=output)

    # matrix cost and gait event intuition with basic correlation
    cost = matrix_cost(x, z, gyr_ref, jerk_ref)
    
    pic_correl_start = find_stride.indexes(cost, 0.35, min_dist=len(gyr_ref) // 2, thres_abs=True)
    pic_correl_start = pic_correl_start[np.argsort(-cost[pic_correl_start])]
    F = [0] * len(x)  # same size as the signal, allows for counting if the steps are identified.

    # gait events exact determination with DTW
    starts = []
    ends = []
    sims = []
    annotations = []
    steps_list = []
    
    for i in range(len(pic_correl_start)):
        step = []
        start_min, end_min, path_min, sim_min, annotations_min = affine_annotate_dtw(x, z, pic_correl_start[i],
                                                                                     gyr_ref, jerk_ref,
                                                                                     stride_annotations_ref)

        add = False
        ho = start_min + annotations_min["HO"]
        to = start_min + annotations_min["TO"]
        hs = start_min + annotations_min["HS"]
        ff = start_min + annotations_min["FF"]

        if (to < hs) & (np.sum(F[to:hs]) == 0):
            for k in range(min(ff, ho), max(ff, ho) + 1):
                F[k] = 1
            add = True
        if (to > hs) & (np.sum(F[hs:to]) == 0):
            for k in range(min(ff, ho), max(ff, ho) + 1):
                F[k] = 1
            add = True

        if add:
            starts.append(start_min)
            ends.append(end_min)
            sims.append(sim_min)
            annotations.append(annotations_min)

            step.append(foot)
            step.append(100)
            step.append(ho)
            step.append(to)
            step.append(hs)
            step.append(ff)
            step.append(sim_min)
            steps_list.append(step)

    steps_side = np.array(steps_list)
    steps_side = steps_side[steps_side[:, 3].argsort()]

    return steps_side, q_side


def affine_annotate_dtw(x, y, start, gyr, jerk, stride_annotations):
    """Matrix cost for the template based detection.

    Parameters
    ----------
        x {vector} -- time series with gyration signal of the entire signal
        y {vector} -- time series with jerk signal of the entire signal
        start {int} -- time series with jerk signal of the entire signal
        gyr {vector} -- time series with gyration signal of the reference stride
        jerk {vector} -- time series with jerk signal of the reference stride
        stride_annotations {dict} -- dictionary with the index for the 4 gait events of the model sride: HS, FF, HO, TO. 

    Returns
    -------
       numpy array -- numpy array with the matrix cost
    """
    
    L = len(gyr)
    end = start + L
    
    # parameters
    s_y1 = np.array([y[start:end] / np.max(y[start:end]),
                     x[start:end] / np.max(abs(x[start:end]))])
    s_y1 = s_y1.transpose()

    s_y2 = np.array([jerk / np.max(jerk),
                     gyr / np.max(abs(gyr))])
    s_y2 = s_y2.transpose()

    r = 2
    path_min, sim_min = metrics.dtw_path(s_y1, s_y2, global_constraint="itakura", itakura_max_slope=r)
    start_min = start
    end_min = start + L

    annotations_min = deal_stride.annotate(path_min, stride_annotations)

    return start_min, end_min, path_min, sim_min, annotations_min


def matrix_cost(x, z, gyr_ref, jerk_ref, mu=0.1):
    """Matrix cost for the template based detection.

    Parameters
    ----------
        x {vector} -- time series with gyration signal of the entire signal
        z {vector} -- time series with jerk signal of the entire signal
        gyr_ref {vector} -- time series with gyration signal of the reference stride
        jerk_ref {vector} -- time series with jerk signal of the reference stride
        mu {float} -- coefficient of expansion and compression above which 2 signals are considered uncorrelated

    Returns
    -------
       numpy array -- numpy array with the matrix cost
    """
    
    Nx = len(x)

    matrix = [0 for i in range(Nx)]
    matrix_2 = [0 for i in range(Nx)]

    u = gyr_ref
    v = jerk_ref
    Nd = len(u)

    for j in range(0, Nx - Nd + 1):  

        # filter with amplitude
        cx = np.std(x[j:j + Nd])
        cu = np.std(u)
        cz = np.std(z[j:j + Nd])
        cv = np.std(v)

        # correlation estimation
        w = stats.pearsonr(x[j:j + Nd], u)[0] / 2 + stats.pearsonr(z[j:j + Nd], v)[0] / 2
        w_2 = max(stats.pearsonr(x[j:j + Nd], u)[0], stats.pearsonr(z[j:j + Nd], v)[0])
        if (cx > mu * cu) & (cz > mu * cv):
            matrix[j] = w
            matrix_2[j] = w_2
        else:
            matrix[j] = -abs(w) / 10  
            matrix_2[j] = -abs(w_2) / 10

    return np.array(matrix, dtype=float)
