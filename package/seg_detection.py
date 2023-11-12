# Objective: Detection of gait phase segmentation based on angular position in the plane orthogonal to X for U-turn
# and movement onset detection. The XSens experiment is divided into 3 phases: go, U-turn, back.

# Output file includes 2 timestamps delineating the U-turn.

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from semio_package import features as ft


def disp_seg(seg):
    seg.columns = ['Index de temps (en 0.01s)']
    seg.index = ['Début marche', 'Début demi-tour', 'Fin demi-tour', 'Fin marche']
    print(seg)


def plot_segdetection_xs(data_tr, steps_lim, id_exp=0, freq=100, download=False, output_files=0):
    # segmentation
    seg = seg_detection(data_tr, steps_lim, id_exp=id_exp, freq=freq, download=download, output_files=output_files)

    # Graphic signals
    t_full, angle_x_full = signals_for_seg(data_tr, steps_lim)

    # Figure initialisation et signal brut
    plt.rcParams["figure.figsize"] = (20, 10)
    fig, ax = plt.subplots()
    ax.plot(t_full, angle_x_full, linewidth=3, label='angular position')
    ax.grid()
    ax.set_yticks([0, 90, 180])
    ax.yaxis.set_tick_params(labelsize=12)
    ax.set_ylabel('Angular position (°)', fontsize=15)
    ax.set_xlabel('Time (s)', fontsize=15)
    ax.xaxis.set_tick_params(labelsize=12)
    ax.set(title="Position angulaire selon X : " + id_exp)

    # Marqueurs de la segmentation de la marche en rouge pour le demi-tour
    ax.vlines(seg[1] / freq, -50, 230, 'red', '-', linewidth=2, label="$u_{go}$ and $u_{back}$")
    ax.vlines(seg[2] / freq, -50, 230, 'red', '-', linewidth=2)
    fig.legend(fontsize=15)

    if download:
        # On enregistre la nouvelle figure
        titre = id_exp + "_auto_seg_image.png"
        os.chdir(output_files)
        plt.savefig(titre, bbox_inches="tight")
        # plt.show()
        plt.close()

    return seg


def seg_detection(data_tr, steps_lim, id_exp=0, freq=100, download=False, output_files=0):
    start = int(np.min(steps_lim["TO"]))
    end = int(np.max(steps_lim["HS"]))
    # Useful signals
    t_full, angle_x_full = signals_for_seg(data_tr, steps_lim)

    # Middle argument
    mid_index = find_nearest(angle_x_full, 90)
    l_search = min(abs(end - mid_index), abs(start - mid_index))
    start = int(mid_index - l_search)
    end = int(mid_index + l_search)

    t = data_tr["PacketCounter"][start:end]
    angle_x = angle_x_full[start:end]

    L = len(t)
    # Affine approximation of go phase
    a_go, b_go, r_go, p_value_go, std_err_go = linregress(t[0:2*L//5], angle_x[0:2*L//5])
    # Affine approximation of back phase
    a_back, b_back, r_back, p_value_back, std_err_back = linregress(t[3*L//5:], angle_x[3*L//5:])
    # Affine approximation of U-turn phase
    a_u, b_u, r_u, p_value_u, std_err_u = linregress(t_full[mid_index - freq:mid_index + freq],
                                                     angle_x_full[mid_index - freq:mid_index + freq])

    # Intersection points
    x_inter_go = (b_go - b_u) / (a_u - a_go)
    x_inter_back = (b_back - b_u) / (a_u - a_back)
    approx_seg_lim = [start, freq * x_inter_go, freq * x_inter_back, end]

    # U-Turn boundaries
    # strT = ft.stride_time(data_tr, approx_seg_lim, steps_lim, freq=freq)
    # print(x_inter_go, strT, start, approx_seg_lim)
    # start_uturn = 100 * (x_inter_go - strT) + np.argmin(
    # angle_x[int(100 * x_inter_go - 100 * strT - start):int(100 * x_inter_go - start)])
    # end_uturn = int(100 * x_inter_back) + np.argmax(
    # angle_x[int(100 * x_inter_back - start):int(100 * x_inter_back + 100 * strT - start)])

    # seg = [0, start_uturn, end_uturn, len(data_tr["Gyr_X"])]
    seg = [0, int(100*x_inter_go), int(100*x_inter_back), len(data_tr["Gyr_X"])]

    if download:
        os.chdir(output_files)
        name_file = id_exp + '_auto_seg.txt'
        np.savetxt(name_file, seg, delimiter=';', fmt='%d')

    return seg


def signals_for_seg(data_tr, steps_lim):

    gyr_x = data_tr['Gyr_X']
    angle_x_full = np.cumsum(gyr_x)
    a = np.median(angle_x_full[0:len(angle_x_full) // 2])  # Tout début du signal
    z = np.median(angle_x_full[len(angle_x_full) // 2:len(angle_x_full)])  # Fin du signal

    angle_x_full = np.sign(z) * (angle_x_full - a) * 180 / abs(z)
    t_full = data_tr["PacketCounter"]

    return t_full, angle_x_full


def find_nearest(array, value):
    i = 0
    while (value - array[i]) * value > 0:  # Tant qu'on est du même côté, c'est à dire qu'ils ont le même signe
        i += 1
    if abs(value - array[i]) < abs(value - array[i-1]):
        return i
    else:
        return i-1
