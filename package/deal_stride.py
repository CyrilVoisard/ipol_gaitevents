import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import interpolate


def stride_sain_decal(len_estimation, freq=100):
    gyr_ref_decal = dict()
    acc_ref_decal = dict()
    stride_ref_decal_annotations = dict()
    gyr_ref, acc_ref, stride_ref_annotations = stride_sain(len_estimation, freq)
    gyr_ref_decal["0"] = gyr_ref
    acc_ref_decal["0"] = acc_ref
    stride_ref_decal_annotations["0"] = stride_ref_annotations

    for i in range(1, len(gyr_ref)):
        gyr_ref_decal[str(i)] = np.concatenate((gyr_ref[i:], gyr_ref[:i]), axis=0)
        acc_ref_decal[str(i)] = np.concatenate((acc_ref[i:], acc_ref[:i]), axis=0)
        annotations = dict()
        annotations["HS"] = (stride_ref_annotations["HS"] - i) % len(gyr_ref)
        annotations["FF"] = (stride_ref_annotations["FF"] - i) % len(gyr_ref)
        annotations["HO"] = (stride_ref_annotations["HO"] - i) % len(gyr_ref)
        annotations["TO"] = (stride_ref_annotations["TO"] - i) % len(gyr_ref)
        stride_ref_decal_annotations[str(i)] = annotations

    return gyr_ref_decal, acc_ref_decal, stride_ref_decal_annotations


def stride_sain(len_estimation, freq=100):
    
    gyr_ref_100 = np.array([1.26760644e-01, 1.11668357e-01, 8.32307508e-02, 5.23749713e-02,
                            2.49340487e-02, 8.49469197e-04, -2.21600981e-02, -4.43579484e-02,
                            -6.35252125e-02, -7.98064015e-02, -1.02934550e-01, -1.55891377e-01,
                            -2.70143688e-01, -4.72886234e-01, -7.73106030e-01, -1.15611776e+00,
                            -1.59250945e+00, -2.05840201e+00, -2.55496150e+00, -3.11276623e+00,
                            -3.77422420e+00, -4.56166497e+00, -5.45103567e+00, -6.37140923e+00,
                            -7.23515596e+00, -7.98030738e+00, -8.59202565e+00, -9.07916982e+00,
                            -9.41461364e+00, -9.48374908e+00, -9.09457055e+00, -8.06869603e+00,
                            -6.37173223e+00, -4.19714027e+00, -1.93070554e+00, 5.27506218e-03,
                            1.34410524e+00, 2.08519285e+00, 2.45319284e+00, 2.73593271e+00,
                            3.11093218e+00, 3.47143548e+00, 3.76652520e+00, 4.10501159e+00,
                            4.38901550e+00, 4.62756335e+00, 4.84052117e+00, 5.05616942e+00,
                            5.24828300e+00, 5.34637901e+00, 5.40054707e+00, 5.40620508e+00,
                            5.40620508e+00, 5.40620508e+00, 5.40620508e+00, 5.40620508e+00,
                            5.40620508e+00, 5.40620508e+00, 5.40620508e+00, 5.40620508e+00,
                            5.40620508e+00, 5.36620508e+00, 5.29285011e+00, 5.16633289e+00,
                            5.01176586e+00, 4.75517299e+00, 4.29232979e+00, 3.76639946e+00,
                            3.15159076e+00, 2.43453909e+00, 1.59675984e+00, 6.60706177e-01,
                            -3.83399033e-01, -1.57999761e+00, -2.75348364e+00, -3.88444378e+00,
                            -4.62385631e+00, -5.22448655e+00, -5.52268047e+00, -5.49605324e+00,
                            -5.16153231e+00, -4.72743923e+00, -4.18226282e+00, -3.55346652e+00,
                            -2.87584621e+00, -2.24011482e+00, -1.63418599e+00, -1.08828296e+00,
                            -6.06997277e-01, -1.61468643e-01, 1.18824624e-01, 3.40925762e-01,
                            3.90451580e-01, 3.09668516e-01, 1.70005975e-01, 3.26757531e-02,
                            -7.51470858e-02, -1.57077914e-01, -2.26225343e-01, -2.84925791e-01,
                            -3.19621130e-01, -3.11547167e-01, -2.52768034e-01, -1.55177074e-01,
                            -4.67644606e-02, 4.13623151e-02, 8.84683220e-02])

    t_ref_100 = np.linspace(0, len(gyr_ref_100) - 1, len(gyr_ref_100))
    t_ref_freq = np.linspace(0, len(gyr_ref_100) - 1, round(len(gyr_ref_100) * freq / 100))
    interp_gyr = interpolate.interp1d(t_ref_100, gyr_ref_100, kind='cubic')
    gyr_ref_freq = interp_gyr(t_ref_freq).tolist()

    n_ajout = max(0, len_estimation - len(t_ref_freq))
    ajout = [1.72042875e-01]

    gyr_ref = np.array(ajout * n_ajout + gyr_ref_freq)

    jerk_ref_100 = np.array([7.33013844e-01, 5.08680015e-01, 2.83958884e-01, 2.13838951e-01,
                                           2.62533483e-01, 3.39655179e-01, 3.95763084e-01, 4.42510601e-01,
                                           5.88738895e-01, 9.28792999e-01, 1.35628149e+00, 1.67269301e+00,
                                           1.84836619e+00, 1.91495327e+00, 1.73255119e+00, 1.30569317e+00,
                                           1.14655260e+00, 1.61658351e+00, 2.14763725e+00, 2.34222221e+00,
                                           3.65667316e+00, 8.04565229e+00, 1.53677445e+01, 2.60056893e+01,
                                           4.61409381e+01, 8.31398263e+01, 1.31398301e+02, 1.69812792e+02,
                                           1.79394619e+02, 1.94540100e+02, 1.98540100e+02, 1.98540100e+02,
                                           1.98540100e+02, 1.94540100e+02, 1.79394619e+02, 1.69812792e+02,
                                           1.51398301e+02, 1.23737150e+02, 6.45552726e+01, 3.60850688e+01,
                                           2.45618215e+01, 1.71648442e+01, 1.49930182e+01, 1.87042949e+01,
                                           2.23649762e+01, 2.31467738e+01, 2.23968317e+01, 1.92686082e+01,
                                           1.30519897e+01, 7.02659614e+00, 4.14744285e+00, 3.60903171e+00,
                                           3.66007094e+00, 3.80497913e+00, 4.09197081e+00, 4.48508542e+00,
                                           4.70672553e+00, 4.64844313e+00, 5.11775977e+00, 7.24105098e+00,
                                           1.15424821e+01, 1.77197928e+01, 2.37151860e+01, 2.68721413e+01,
                                           2.79508998e+01, 2.97135396e+01, 3.34259338e+01, 4.43466060e+01,
                                           6.69896117e+01, 8.36166879e+01, 7.83897281e+01, 1.03591705e+02,
                                           2.38669992e+02, 4.44851476e+02, 5.82533481e+02, 6.09298880e+02,
                                           5.95001660e+02, 5.39271190e+02, 4.00169700e+02, 2.67868961e+02,
                                           2.51055605e+02, 2.70900752e+02, 2.09216231e+02, 1.20965555e+02,
                                           9.74239139e+01, 9.72068879e+01, 6.32176929e+01, 3.12584733e+01,
                                           3.36207832e+01, 4.03751047e+01, 3.18245086e+01, 2.41630004e+01,
                                           2.24790693e+01, 1.62291136e+01, 7.34778457e+00, 4.39378791e+00,
                                           5.41511294e+00, 4.97150176e+00, 3.74616344e+00, 3.67267255e+00,
                                           3.69761117e+00, 2.90919078e+00, 2.04212959e+00, 1.61950796e+00,
                                           1.41944068e+00, 1.22664442e+00, 9.52457080e-01])

    t_ref_100 = np.linspace(0, len(jerk_ref_100) - 1, len(jerk_ref_100))
    t_ref_freq = np.linspace(0, len(jerk_ref_100) - 1, round(len(jerk_ref_100) * freq / 100))
    interp_jerk = interpolate.interp1d(t_ref_100, jerk_ref_100, kind='cubic')
    jerk_ref_freq = interp_jerk(t_ref_freq).tolist()

    n_ajout = max(0, len_estimation - len(t_ref_freq))
    ajout = [1.72042875e-01]

    jerk_ref = np.array(ajout * n_ajout + jerk_ref_freq)

    stride_ref_annotations = {'HS': round(74*freq/100) + len(ajout * n_ajout),
                              'FF': round(89*freq/100) + len(ajout * n_ajout),
                              'HO': round(10*freq/100) + len(ajout * n_ajout),
                              'TO': round(32*freq/100) + len(ajout * n_ajout)}

    return gyr_ref, np.sqrt(jerk_ref), stride_ref_annotations


def plot_annotate_stride(gyr, acc, stride_annotations, output):

    fig, ax = plt.subplots(2, figsize=(5, 8))
    ax[0].plot(gyr / np.max(abs(gyr)))

    mi, ma = -1, 1
    ax[0].vlines(stride_annotations["HS"], mi, ma, 'black', label="Heel Strike")
    ax[0].vlines(stride_annotations["FF"], mi, ma, 'green', label="Foot Flat")
    ax[0].vlines(stride_annotations["HO"], mi, ma, 'green', label="Heel Off")
    ax[0].vlines(stride_annotations["TO"], mi, ma, 'red', label="Toe Off")
    ax[0].legend()
    ax[0].grid()
    ax[0].set_ylabel("Gyr_Y")
    ax[1].plot(acc / np.max(abs(acc)))
    mi, ma = 0, 1
    ax[1].vlines(stride_annotations["HS"], mi, ma, 'black', label="Heel Strike")
    ax[1].vlines(stride_annotations["FF"], mi, ma, 'green', label="Foot Flat")
    ax[1].vlines(stride_annotations["HO"], mi, ma, 'green', label="Heel Off")
    ax[1].vlines(stride_annotations["TO"], mi, ma, 'red', label="Toe Off")
    ax[1].legend()
    ax[1].set_ylabel("Jerk_tot")
    ax[1].grid()

    titre = "stride_reference.png"
    os.chdir(output)
    plt.savefig(titre, bbox_inches="tight")
    plt.close()

    return None


def annotate(path, stride_annotations):
    new_stride_annotations = {'HS': float("inf"), 'FF': float("inf"), 'HO': float("inf"), 'TO': float("inf")}

    for i in range(len(path)):
        if path[i][1] == stride_annotations["HS"]:
            new_stride_annotations["HS"] = path[i][0]
        if path[len(path) - i - 1][1] == stride_annotations["FF"]:
            new_stride_annotations["FF"] = path[len(path) - i - 1][0]
        if path[i][1] == stride_annotations["HO"]:
            new_stride_annotations["HO"] = path[i][0]
        if path[len(path) - i - 1][1] == stride_annotations["TO"]:
            new_stride_annotations["TO"] = path[len(path) - i - 1][0]

    return new_stride_annotations


def calculate_jerk_tot(data, freq=100):
    
    jerk_tot = np.sqrt(
        np.diff(data["FreeAcc_X"]) ** 2 + np.diff(data["FreeAcc_Z"]) ** 2 + np.diff(data["FreeAcc_Y"]) ** 2)
    jerk_tot = np.array(jerk_tot.tolist() + [0])
    y = pd.DataFrame(jerk_tot)
    y_mean = y.rolling(9, center=True, win_type='boxcar').sum()
    y_mean = y_mean.fillna(0)
    z = y_mean.to_numpy().transpose()[0]
    return z
