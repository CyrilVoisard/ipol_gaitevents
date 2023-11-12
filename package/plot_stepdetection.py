# Objectif :

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import patches
import os

from mp_dtw_package import deal_stride


def plot_stepdetection_dtw(steps_rf, steps_lf, data_pied_droit, data_pied_gauche, id_exp=0,
                           download=False, output_files=0, *args, **kwargs):
    nom = id_exp + "_Détection de pas "

    fig, ax = plt.subplots(4, figsize=(20, 12), sharex=True, sharey=False)

    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[3].grid()

    ax[0].set(ylabel='Jerk total', title=nom + "Left Foot")
    ax[0].yaxis.set_tick_params(labelsize=12)
    ax[1].set(xlabel='Time (s)', ylabel='Gyr ML')
    ax[1].xaxis.set_tick_params(labelsize=12)
    ax[2].set(ylabel='Jerk total', title=nom + "Right Foot")
    ax[2].yaxis.set_tick_params(labelsize=12)
    ax[3].set(xlabel='Time (s)', ylabel='Gyr ML')
    ax[3].xaxis.set_tick_params(labelsize=12)

    # ---------------------------- Données pour le pied gauche ---------------------------------------------
    t_lf = data_pied_gauche["PacketCounter"]
    gyr_lf = data_pied_gauche["Gyr_Y"]
    jerk_lf = deal_stride.calculate_jerk_tot(data_pied_gauche)
    ma_lf = max(gyr_lf)
    mi_lf = min(gyr_lf)
    ax[0].plot(t_lf, jerk_lf)
    ax[1].plot(t_lf, gyr_lf)

    # ---------------------------- Données pour le pied droit ---------------------------------------------
    t_rf = data_pied_droit["PacketCounter"]
    gyr_rf = data_pied_droit["Gyr_Y"]
    jerk_rf = deal_stride.calculate_jerk_tot(data_pied_droit)
    ma_rf = max(gyr_rf)
    mi_rf = min(gyr_rf)
    ax[2].plot(t_rf, jerk_rf)
    ax[3].plot(t_rf, gyr_rf)

    for i in range(len(steps_lf)):
        # ----------------------- Remplissage pour le pied gauche -------------------------------------------
        ax[0].plot(t_lf[int(steps_lf[i][2])], jerk_lf[int(steps_lf[i][2])], 'g', marker='x')
        ax[0].plot(t_lf[int(steps_lf[i][3])], jerk_lf[int(steps_lf[i][3])], 'r', marker='x')
        ax[0].plot(t_lf[int(steps_lf[i][4])], jerk_lf[int(steps_lf[i][4])], 'k', marker='x')
        ax[0].plot(t_lf[int(steps_lf[i][5])], jerk_lf[int(steps_lf[i][5])], 'g', marker='x')

        to = int(steps_lf[i][3])
        ax[1].vlines(t_lf[to], mi_lf, ma_lf, 'k', '--')
        hs = int(steps_lf[i][4])
        ax[1].vlines(t_lf[hs], mi_lf, ma_lf, 'k', '--')

        if to < hs:
            ax[1].add_patch(patches.Rectangle((t_lf[to], mi_lf),  # (x,y)
                                              t_lf[hs] - t_lf[to],  # width
                                              ma_lf - mi_lf,  # height
                                              alpha=0.1,
                                              facecolor='red', linestyle='dotted'))
            if i < len(steps_lf) - 1:
                to_ap = int(steps_lf[i + 1][3])
                ax[1].add_patch(patches.Rectangle((t_lf[hs], min(gyr_lf)),  # (x,y)
                                                  t_lf[to_ap] - t_lf[hs],  # width
                                                  max(gyr_lf) - min(gyr_lf),  # height
                                                  alpha=0.1,
                                                  facecolor='green', linestyle='dotted'))
        else:
            ax[1].add_patch(patches.Rectangle((t_lf[hs], min(gyr_lf)),  # (x,y)
                                              t_lf[to] - t_lf[hs],  # width
                                              max(gyr_lf) - min(gyr_lf),  # height
                                              alpha=0.1,
                                              facecolor='green', linestyle='dotted'))
            if i < len(steps_lf) - 1:
                hs_ap = int(steps_lf[i + 1][4])
                ax[1].add_patch(patches.Rectangle((t_lf[to], min(gyr_lf)),  # (x,y)
                                                  t_lf[hs_ap] - t_lf[to],  # width
                                                  max(gyr_lf) - min(gyr_lf),  # height
                                                  alpha=0.1,
                                                  facecolor='red', linestyle='dotted'))

        # ----------------------- Remplissage pour le pied droit -------------------------------------------
    for i in range(len(steps_rf)):
        ax[2].plot(t_rf[int(steps_rf[i][2])], jerk_rf[int(steps_rf[i][2])], 'g', marker='x')
        ax[2].plot(t_rf[int(steps_rf[i][3])], jerk_rf[int(steps_rf[i][3])], 'r', marker='x')
        ax[2].plot(t_rf[int(steps_rf[i][4])], jerk_rf[int(steps_rf[i][4])], 'k', marker='x')
        ax[2].plot(t_rf[int(steps_rf[i][5])], jerk_rf[int(steps_rf[i][5])], 'g', marker='x')

        to = int(steps_rf[i][3])
        ax[3].vlines(t_rf[to], mi_rf, ma_rf, 'k', '--')
        hs = int(steps_rf[i][4])
        ax[3].vlines(t_rf[hs], mi_rf, ma_rf, 'k', '--')

        if to < hs:
            ax[3].add_patch(patches.Rectangle((t_rf[to], mi_rf),  # (x,y)
                                              t_rf[hs] - t_rf[to],  # width
                                              ma_rf - mi_rf,  # height
                                              alpha=0.1,
                                              facecolor='red', linestyle='dotted'))
            if i < len(steps_rf) - 1:
                to_ap = int(steps_rf[i + 1][3])
                ax[3].add_patch(patches.Rectangle((t_rf[hs], mi_rf),  # (x,y)
                                                  t_rf[to_ap] - t_rf[hs],  # width
                                                  ma_rf - mi_rf,  # height
                                                  alpha=0.1,
                                                  facecolor='green', linestyle='dotted'))
        else:
            ax[3].add_patch(patches.Rectangle((t_rf[hs], min(gyr_rf)),  # (x,y)
                                              t_rf[to] - t_rf[hs],  # width
                                              max(gyr_rf) - min(gyr_rf),  # height
                                              alpha=0.1,
                                              facecolor='green', linestyle='dotted'))
            if i < len(steps_rf) - 1:
                hs_ap = int(steps_rf[i + 1][4])
                ax[3].add_patch(patches.Rectangle((t_rf[to], min(gyr_rf)),  # (x,y)
                                                  t_rf[hs_ap] - t_rf[to],  # width
                                                  max(gyr_rf) - min(gyr_rf),  # height
                                                  alpha=0.1,
                                                  facecolor='red', linestyle='dotted'))

    # Légende pour les couleurs des pas
    red_patch = mpatches.Patch(color='red', alpha=0.1, label='swing')
    green_patch = mpatches.Patch(color='green', alpha=0.1, label='stance')
    red_circle = mpatches.Patch(color='red', label='Toe Off')
    black_circle = mpatches.Patch(color='black', label='Heel Strike')
    green_circle = mpatches.Patch(color='green', label='Flat Foot & Heel Off')

    ax[0].legend(handles=[black_circle, green_circle, red_circle], loc="upper left")
    ax[1].legend(handles=[red_patch, green_patch], loc="upper left")
    ax[2].legend(handles=[black_circle, green_circle, red_circle], loc="upper left")
    ax[3].legend(handles=[red_patch, green_patch], loc="upper left")

    if download:
        # On enregistre la nouvelle figure
        titre = id_exp + "_dtw_steps_image.png"
        os.chdir(output_files)
        plt.savefig(titre)
        # plt.show()
        plt.close()
