import os
import numpy as np
import pandas as pd
from scipy import stats
from tslearn import metrics

from package import find_stride, deal_stride, plot_stepdetection


def steps_detection_full(data_rf, data_lf, freq):
    steps_rf = steps_detection(data_rf, data_lf, 1, freq)
    steps_lf = steps_detection(data_lf, data_rf, 0, freq)
    
    full = np.concatenate((steps_rf, steps_lf))

    return pd.DataFrame(full, columns=["Foot", "Phase", "HO", "TO", "HS", "FF", "Score"])


def steps_detection(data, data_other_side, foot, freq):
    x = data["Gyr_Y"]
    z = deal_stride.calculate_jerk_tot(data, freq)

    gyr_ok, acc_ok, stride_annotations_ok, comp = find_stride.annotate_stride_estimation(data, data_other_side, foot, freq)

    cost = matrix_cost(x, z, gyr_ok, acc_ok)

    pic_correl_start = find_stride.indexes(cost, 0.35, min_dist=len(gyr_ok) // 2, thres_abs=True)
    pic_correl_start = pic_correl_start[np.argsort(-cost[pic_correl_start])]
    F = [0] * len(x)  # same size as the signal, allows for counting if the steps are identified.

    starts = []
    ends = []
    sims = []
    annotations = []
    steps_list = []

    for i in range(len(pic_correl_start)):
        step = []
        start_min, end_min, path_min, sim_min, annotations_min = affine_annotate_dtw(x, z, pic_correl_start[i],
                                                                                     gyr_ok, acc_ok,
                                                                                     stride_annotations_ok)

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
        #else:
         #   print("Step déjà rentré : ", ho, to, hs, ff)
          #  print(F[to:hs], F[hs:to])

    steps_list = np.array(steps_list)
    steps_list = steps_list[steps_list[:, 3].argsort()]

    return steps_list


def affine_annotate_dtw(x, y, start, gyr, acc, stride_annotations, disp=False):
    L = len(gyr)
    end = start + L

    dx = np.array(np.diff(x).tolist() + [0])
    dgyr = np.array(np.diff(gyr).tolist() + [0])

    # parameters
    s_y1 = np.array([y[start:end] / np.max(y[start:end]),
                     x[start:end] / np.max(abs(x[start:end]))])
    s_y1 = s_y1.transpose()

    s_y2 = np.array([acc / np.max(acc),
                     gyr / np.max(abs(gyr))])
    s_y2 = s_y2.transpose()

    r = 2
    path_min, sim_min = metrics.dtw_path(s_y1, s_y2, global_constraint="itakura", itakura_max_slope=r)
    start_min = start
    end_min = start + L

    annotations_min = deal_stride.annotate(path_min, stride_annotations)

    return start_min, end_min, path_min, sim_min, annotations_min


def matrix_cost(x, z, gyr_ok, jerk_ok, mu=0.1):
    Nx = len(x)

    matrix = [0 for i in range(Nx)]
    matrix_2 = [0 for i in range(Nx)]

    u = gyr_ok
    v = jerk_ok
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
