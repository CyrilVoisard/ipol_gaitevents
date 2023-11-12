import os
import numpy as np
import pandas as pd
from scipy import stats
from tslearn import metrics

from package import find_stride, deal_stride, plot_stepdetection


def steps_detection_full(data_rf, data_lf, id_exp, freq=100, gr=False, exo=[], download=False, output_files=0):
    steps_rf = steps_detection(data_rf, data_lf, 1, id_exp, freq=freq, gr=gr, exo=exo, download=download,
                               output_files=output_files)
    steps_lf = steps_detection(data_lf, data_rf, 0, id_exp, freq=freq, gr=gr, exo=exo, download=download,
                               output_files=output_files)

    full = np.concatenate((steps_rf, steps_lf))

    if download:
        os.chdir(output_files)
        name_file = id_exp + '_dtw_steps.txt'
        print(name_file)
        np.savetxt(name_file, full, delimiter=';', fmt='%d')

    plot_stepdetection.plot_stepdetection_dtw(steps_rf, steps_lf, data_rf, data_lf, id_exp=id_exp, download=True,
                                              output_files=output_files)
    print('RAS image steps')

    return pd.DataFrame(full, columns=["Foot", "Phase", "HO", "TO", "HS", "FF", "Score"])


def steps_detection(data, data_autre, pied, id_exp, freq=100, gr=False, exo=None, download=False, output_files=0):
    x = data["Gyr_Y"]
    z = deal_stride.calculate_jerk_tot(data, freq)

    gyr_ok, acc_ok, stride_annotations_ok, comp = find_stride.annotate_stride_estimation(data, data_autre, pied,
                                                                                         id_exp=id_exp, freq=freq,
                                                                                         gr=gr, exo=exo,
                                                                                         download=download,
                                                                                         output=output_files)

    cost = matrix_cost(x, z, gyr_ok, acc_ok)

    pic_correl_debut = indexes(cost, 0.35, min_dist=len(gyr_ok) // 2, thres_abs=True)
    pic_correl_debut = pic_correl_debut[np.argsort(-cost[pic_correl_debut])]
    F = [0] * len(x)  # Same size as the signal, allows for counting if the steps are identified.

    starts = []
    ends = []
    sims = []
    annotations = []
    steps_list = []

    for i in range(len(pic_correl_debut)):
        step = []
        start_min, end_min, path_min, sim_min, annotations_min = affine_annotate_dtw(x, z, pic_correl_debut[i],
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

            step.append(pied)
            step.append(100)
            step.append(ho)
            step.append(to)
            step.append(hs)
            step.append(ff)
            step.append(sim_min)
            steps_list.append(step)
        else:
            print("Step déjà rentré : ", ho, to, hs, ff)
            print(F[to:hs], F[hs:to])

    steps_list = np.array(steps_list)
    steps_list = steps_list[steps_list[:, 3].argsort()]

    return steps_list


def affine_annotate_dtw(x, y, start, gyr, acc, stride_annotations, disp=False):
    L = len(gyr)
    end = start + L

    dx = np.array(np.diff(x).tolist() + [0])
    dgyr = np.array(np.diff(gyr).tolist() + [0])

    # Paramètres de départ
    s_y1 = np.array([1 * y[start:end] / np.max(y[start:end]),
                     1 * x[start:end] / np.max(abs(x[start:end])),
                     0 * abs(dx[start:end]) / (np.max(abs(dx[start:end])))])
    s_y1 = s_y1.transpose()

    s_y2 = np.array([1 * acc / np.max(acc),
                     1 * gyr / np.max(abs(gyr)),
                     0 * abs(dgyr) / (np.max(abs(dgyr)))])
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

    for j in range(0, Nx - Nd + 1):  # À quelle point le signal qui suit correspond au pas du dictionnaire ?

        # Sélection avec l'amplitude mu
        cx = np.std(x[j:j + Nd])
        cu = np.std(u)
        cz = np.std(z[j:j + Nd])
        cv = np.std(v)

        # Calcul de la corrélation
        w = stats.pearsonr(x[j:j + Nd], u)[0] / 2 + stats.pearsonr(z[j:j + Nd], v)[0] / 2
        w_2 = max(stats.pearsonr(x[j:j + Nd], u)[0], stats.pearsonr(z[j:j + Nd], v)[0])
        if (cx > mu * cu) & (cz > mu * cv):
            matrix[j] = w
            matrix_2[j] = w_2
        else:
            matrix[j] = -abs(w) / 10  # si on met des 0, c'est sélectionné dans les >=...
            matrix_2[j] = -abs(w_2) / 10  # si on met des 0, c'est sélectionné dans les >=...

    return np.array(matrix, dtype=float)


def indexes(y, thres=0.3, min_dist=1, thres_abs=False):
    """Peak detection routine.

    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.

    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    thres_abs: boolean
        If True, the thres value will be interpreted as an absolute value, instead of
        a normalized threshold.

    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected.
        When using with Pandas DataFrames, iloc should be used to access the values at the returned positions.
    """
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    if not thres_abs:
        thres = thres * (np.max(y) - np.min(y)) + np.min(y)

    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)
    # print("dy", dy)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    # print("where", np.where(dy == 0))
    zeros, = np.where(dy == 0)

    # check if the signal is totally flat
    if len(zeros) == len(y) - 1:
        print("Signal seems to be flat !")
        return np.array([])

    if len(zeros):
        # compute first order difference of zero indexes
        zeros_diff = np.diff(zeros)
        # check when zeros are not chained together
        zeros_diff_not_one, = np.add(np.where(zeros_diff != 1), 1)
        # make an array of the chained zero indexes
        zero_plateaus = np.split(zeros, zeros_diff_not_one)

        # fix if leftmost value in dy is zero
        if zero_plateaus[0][0] == 0:
            dy[zero_plateaus[0]] = dy[zero_plateaus[0][-1] + 1]
            zero_plateaus.pop(0)

        # fix if rightmost value of dy is zero
        if len(zero_plateaus) and zero_plateaus[-1][-1] == len(dy) - 1:
            dy[zero_plateaus[-1]] = dy[zero_plateaus[-1][0] - 1]
            zero_plateaus.pop(-1)

        # for each chain of zero indexes
        for plateau in zero_plateaus:
            median = np.median(plateau)
            # set leftmost values to leftmost non zero values
            dy[plateau[plateau < median]] = dy[plateau[0] - 1]
            # set rightmost and middle values to rightmost non zero values
            dy[plateau[plateau >= median]] = dy[plateau[-1] + 1]

    # find the peaks by using the first order difference
    peaks = np.where(
        (np.hstack([dy, 0.0]) < 0.0)
        & (np.hstack([0.0, dy]) > 0.0)
        & (np.greater(y, thres))
    )[0]

    # handle multiple peaks, respecting the minimum distance
    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks