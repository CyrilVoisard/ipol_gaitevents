import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from load_package import seg_steps


def compute_precision(precision, data, steps_full, exo, id_exp):
    precision_full = precision
    precision_lf = precision_one_foot(data[0], exo, steps_full, 0, id_exp=id_exp)
    precision_rf = precision_one_foot(data[1], exo, steps_full, 1, id_exp=id_exp)

    # Pour les pas gauches
    # print(len(precision_lf[0]), len(precision_lf[1]), len(precision_lf[2]), len(precision_lf[3]))
    for i in range(len(precision_lf[0])):
        # print(i)
        precision_one = [id_exp[:min(10, len(id_exp))], id_exp[:min(19, len(id_exp))], id_exp, 0, precision_lf[4][i],
                         precision_lf[0][i], precision_lf[1][i], precision_lf[2][i], precision_lf[3][i]]
        precision_full.append(precision_one)

    # Pour les pas gauches
    for i in range(len(precision_rf[0])):
        precision_one = [id_exp[:min(10, len(id_exp))], id_exp[:min(19, len(id_exp))], id_exp, 1, precision_rf[4][i],
                         precision_rf[0][i], precision_rf[1][i], precision_rf[2][i], precision_rf[3][i]]
        precision_full.append(precision_one)

    return precision_full


def precision_one_foot(data, exo, steps_full, pied, id_exp=0, download=False, output_files=0):
    exo_float = np.array(exo[:, (0, 1, 2, 3, 4)], dtype=float)
    exo_p = exo[exo_float[:, 4] == pied]
    exo_p_float = np.array(exo_p[:, (0, 1, 2, 3, 4)], dtype=float)

    # Extraction du demi tour
    u_turn_start = 0
    u_turn_end = 0
    if exo[0][6] == "U_turn":
        u_turn_start = exo_float[0][0]
        u_turn_end = exo_float[0][1]

    # print("exo_p", exo_p)
    # print("exo_p_float", exo_p_float)
    len_str_estimation, len_swg_estimation = extract_len(exo_p_float, u_turn_start, u_turn_end)
    # print("len_str_swg_estimation", len_str_estimation, len_swg_estimation)

    to_tapis = []
    hs_tapis = []
    #print("to_tapis_avant", exo_p_float[exo_p[:, 6] == "Validated", 1])
    #print("hs_tapis_avant", exo_p_float[exo_p[:, 5] == "Validated", 0])
    for i in range(len(exo_p)):
        if (exo_p[i, 6] == "Validated") & ((u_turn_start > exo_p_float[i, 1]) | (exo_p_float[i, 1] > u_turn_end)):
            to_tapis.append(100*exo_p_float[i, 1])
        if (exo_p[i, 5] == "Validated") & ((u_turn_start > exo_p_float[i, 0]) | (exo_p_float[i, 0] > u_turn_end)):
            hs_tapis.append(100*exo_p_float[i, 0])
    # print("to_tapis_après", u_turn_start, u_turn_end, to_tapis)
    # print("hs_tapis_après", u_turn_start, u_turn_end, hs_tapis)

    miss = exo_p[exo_p[:, 5] == "Missing"]  # U_turn est pas toujours compris dedans
    # print("miss", miss)

    steps = steps_full[steps_full[:, 0] == pied]
    to_detect = []
    hs_detect = []
    for i in range(len(steps)):
        to_moment = steps[i][3]
        hs_moment = steps[i][4]
        j_to = 0
        j_hs = 0
        if 100 * u_turn_start < to_moment < 100 * u_turn_end:
            j_to = 1
        if 100 * u_turn_start < hs_moment < 100 * u_turn_end:
            j_hs = 1
        for j in range(len(miss)):
            min_miss = min(float(miss[j][0]), float(miss[j][1])) + 0.05
            max_miss = max(float(miss[j][0]), float(miss[j][1])) - 0.05
            if 100 * min_miss < int(to_moment) < 100 * max_miss:
                # print("Missing phase :", miss[j])
                # print(100 * min_miss, int(to_moment), 100 * max_miss)
                j_to = 1
            if 100 * min_miss < int(hs_moment) < 100 * max_miss:
                # print("Missing phase :", miss[j])
                # print(100*min_miss, int(hs_moment), 100*max_miss)
                j_hs = 1
        if j_to == 0:
            to_detect.append(to_moment)
        if j_hs == 0:
            hs_detect.append(hs_moment)

    hs_tapis_copie = hs_tapis.copy()
    to_tapis_copie = to_tapis.copy()

    start = min(min(to_tapis), min(hs_tapis)) - 10
    end = max(max(to_tapis), max(hs_tapis)) + 10

    # Reconstituer les couples en prenant en base le GR
    event_detect_type = []
    event_detect_value = []
    event_tapis_corresp_type = []
    event_tapis_corresp_value = []
    validate = []

    # print("HS detect : ", hs_detect)
    for i in range(len(hs_detect)):
        if (hs_detect[i] > start) & (hs_detect[i] < end):
            event_detect_type.append("HS")
            event_detect_value.append(hs_detect[i])
            to_nearest = find_nearest(to_tapis, hs_detect[i])
            if len(hs_tapis) == 0:
                print("len_0 HS", hs_tapis_copie, hs_tapis, hs_detect[i], len(hs_detect), hs_detect[i:])
                hs_nearest = -10000
            else:
                hs_nearest = find_nearest(hs_tapis, hs_detect[i])
            event_tapis_corresp_type.append("HS")
            event_tapis_corresp_value.append(hs_nearest)
            if abs(hs_nearest - hs_detect[i]) < 0.2 * len_str_estimation:
                if abs(hs_nearest - hs_detect[i]) <= abs(to_nearest - hs_detect[i]):
                    # print("Ok precision TO : ", to_detect[i], " correspond à ", to_nearest)
                    hs_tapis.remove(hs_nearest)
                    validate.append(1)
                else:
                    print("Problème de précision HS :", hs_detect[i], " plus proche de TO ", hs_nearest, " que de HS ",
                          to_nearest)
                    hs_tapis.remove(hs_nearest)
                    validate.append(1)
            else:
                print("Problème, trop loin precision HS :", hs_detect[i], " pour ", hs_nearest, " écart max ",
                      0.2 * len_str_estimation)
                hs_nearest_copie = find_nearest(hs_tapis_copie, hs_detect[i])
                if abs(hs_nearest_copie - hs_detect[i]) < 0.2 * len_str_estimation :
                    print("Pourtant dans hs_copie on trouve un pas :", hs_nearest_copie, hs_tapis_copie)
                validate.append(0)
        else:
            l = 1
            # print("HS en dehors du tapis ", hs_detect[i], start, end)

    # print("TO detect : ", to_detect)
    for i in range(len(to_detect)):
        if (to_detect[i] > start) & (to_detect[i] < end):
            event_detect_type.append("TO")
            event_detect_value.append(to_detect[i])
            if len(to_tapis) == 0:
                print("len_0 TO", to_tapis_copie, to_tapis, to_detect[i], len(to_detect), to_detect[i:])
                to_nearest = -10000
            else:
                to_nearest = find_nearest(to_tapis, to_detect[i])
            hs_nearest = find_nearest(hs_tapis_copie, to_detect[i])
            event_tapis_corresp_type.append("TO")
            event_tapis_corresp_value.append(to_nearest)
            if abs(to_nearest - to_detect[i]) < 0.2 * len_str_estimation:
                if abs(to_nearest - to_detect[i]) <= abs(hs_nearest - to_detect[i]):
                    # print("Ok precision TO : ", to_detect[i], " correspond à ", to_nearest)
                    #print("TO avant", to_nearest, to_tapis)
                    to_tapis.remove(to_nearest)
                    validate.append(1)

                else:
                    print("Problème de précision TO :", to_detect[i], " plus proche de HS ", hs_nearest, " que de TO ",
                          to_nearest)
                    to_tapis.remove(to_nearest)
                    validate.append(1)
            else:
                print("Problème, trop loin precision TO :", to_detect[i], " pour ", to_nearest, " écart max ",
                      0.2 * len_str_estimation)
                to_nearest_copie = find_nearest(to_tapis_copie, to_detect[i])
                if abs(to_nearest_copie - to_detect[i]) < 0.2 * len_str_estimation:
                    print("Pourtant dans hs_copie on trouve un pas :", to_nearest_copie, to_tapis_copie)
                validate.append(0)
        else:
            l = 1
            # print("TO en dehors du tapis ", to_detect[i], start, end)
    if len(hs_tapis) > 0:
        print("HS après", hs_tapis)
    if len(hs_tapis) > 1:
        print("TO après", to_tapis)

    return np.array([event_detect_type, event_detect_value, event_tapis_corresp_type, event_tapis_corresp_value, validate])


def compute_recall(recall, data, steps_full, exo, id_exp):
    recall_full = recall

    recall_lf = recall_one_foot(data[0], exo, steps_full, 0, id_exp=id_exp)
    recall_rf = recall_one_foot(data[1], exo, steps_full, 1, id_exp=id_exp)
    # print("Len recall", len(recall_rf[0]), len(recall_rf[1]), len(recall_rf[2]), len(recall_rf[3]), len(recall_rf[4]))

    # Pour les pas gauches
    for i in range(len(recall_lf[0])):
        recall_one = [id_exp[:min(10, len(id_exp))], id_exp[:min(19, len(id_exp))], id_exp, 0, recall_lf[4][i],
                      recall_lf[0][i], recall_lf[1][i], recall_lf[2][i], recall_lf[3][i]]
        recall_full.append(recall_one)

    # Pour les pas gauches
    for i in range(len(recall_rf[0])):
        recall_one = [id_exp[:min(10, len(id_exp))], id_exp[:min(19, len(id_exp))], id_exp, 1, recall_rf[4][i],
                      recall_rf[0][i], recall_rf[1][i], recall_rf[2][i], recall_rf[3][i]]
        recall_full.append(recall_one)

    return recall_full


def recall_one_foot(data, exo, steps_full, pied, id_exp=0, download=False, output_files=0):
    exo_float = np.array(exo[:, (0, 1, 2, 3, 4)], dtype=float)
    exo_p = exo[exo_float[:, 4] == pied]
    exo_p_float = np.array(exo_p[:, (0, 1, 2, 3, 4)], dtype=float)
    # Extraction du demi tour
    u_turn_start = 0
    u_turn_end = 0
    if exo[0][6] == "U_turn":
        u_turn_start = exo_float[0][0]
        u_turn_end = exo_float[0][1]
        # print("U_turn find", u_turn_start, exo_p[0][6], u_turn_end)
    # print("U_turn tot", u_turn_start, exo_p[0][6], u_turn_end)

    # print("exo_p", exo_p)
    # print("exo_p_float", exo_p_float)
    len_str_estimation, len_swg_estimation = extract_len(exo_p_float, u_turn_start, u_turn_end)
    # print("len_str_swg_estimation", len_str_estimation, len_swg_estimation)

    steps = steps_full[steps_full[:, 0] == pied]
    to_detect = steps[:, 3]
    hs_detect = steps[:, 4]

    # Reconstituer les couples en prenant en base le GR
    event_tapis_type = []
    event_tapis_value = []
    event_detect_corresp_type = []
    event_detect_corresp_value = []
    validate = []

    for i in range(len(exo_p)):
        if exo_p[i][5] != "Missing":

            # print("HS tapis : ", 100*exo_p_float[i][0], exo_p[i])
            event_tapis_type.append("HS")
            event_tapis_value.append(100 * exo_p_float[i][0])
            to_nearest = find_nearest(to_detect, 100 * exo_p_float[i][0])
            hs_nearest = find_nearest(hs_detect, 100 * exo_p_float[i][0])
            event_detect_corresp_type.append("HS")
            event_detect_corresp_value.append(hs_nearest)
            # print("U_turn HS", u_turn_start, exo_p_float[i][0], u_turn_end)
            # print(exo_p_float[i][0], u_turn_start, u_turn_end)
            if (exo_p_float[i][0] < u_turn_start) | (exo_p_float[i][0] > u_turn_end):
                if (100 * exo_p_float[i][0] > len_str_estimation) & (
                        100 * exo_p_float[i][0] < len(data) - len_str_estimation):
                    if (exo_p[i][5] == "Validated") | (exo_p[i][5] == "Not_validated"):
                        # print("Ok recall HS : ", 100*exo_p_float[i][0], " correspond à ", to_nearest)
                        if abs(hs_nearest - 100 * exo_p_float[i][0]) < 0.2 * len_str_estimation:
                            if exo_p[i][5] == "Validated":
                                validate.append(1)
                                # print("validate", len(validate), exo_p[i][5], validate)
                            else:
                                print("Intermédiaire :", exo_p[i][5], 0.5)
                                validate.append(0.5)
                                # print("validate", len(validate), exo_p[i][5], validate)
                        else:
                            print("Problème, trop loin recall HS :", hs_nearest, 100 * exo_p_float[i][0], 0.2 * len_str_estimation)
                            validate.append(0)
                            # print("validate", len(validate), exo_p[i][5], validate)

                        if not (abs(hs_nearest - 100 * exo_p_float[i][0]) <= abs(
                                to_nearest - 100 * exo_p_float[i][0])) | (
                                       (hs_nearest - 100 * exo_p_float[i][0]) * (
                                       to_nearest - 100 * exo_p_float[i][0]) < 0):
                            print("Recall HS :", 100 * exo_p_float[i][0], " plus proche de TO ", to_nearest,
                                  " que de HS ", hs_nearest)
                            print("Mais dans les limites :", len_str_estimation, 100 * exo_p_float[i][0],
                                  len(data) - len_str_estimation)
                            # event_detect_corresp_type.append("TO")
                            # event_detect_corresp_value.append(to_nearest)
                            # validate.append(0)
                    else:
                        if exo_p[i][5] == "Deleted":
                            #print("Suppression du pas de mauvaise qualité HS : ", exo_p_float[i][0])
                            validate.append(-1)
                            # print("validate", len(validate), exo_p[i][5], validate)
                else:
                    # print("Attention trop proche des limites HS", 100*exo_p_float[i][0], len(data), len_str_estimation)
                    validate.append(-2)
                    # print("validate", len(validate), exo_p[i][5], validate)
            else:
                # print("U_turn HS", u_turn_start, exo_p_float[i][0], u_turn_end)
                validate.append(-3)
                # print("validate", len(validate), exo_p[i][5], validate)

            # print("TO tapis : ", 100 * exo_p_float[i][1])
            event_tapis_type.append("TO")
            event_tapis_value.append(100 * exo_p_float[i][1])
            to_nearest = find_nearest(to_detect, 100 * exo_p_float[i][1])
            hs_nearest = find_nearest(hs_detect, 100 * exo_p_float[i][1])
            event_detect_corresp_type.append("TO")
            event_detect_corresp_value.append(to_nearest)
            if (exo_p_float[i][1] < u_turn_start) | (exo_p_float[i][1] > u_turn_end):
                if (100 * exo_p_float[i][1] > len_str_estimation) & (
                        100 * exo_p_float[i][1] < len(data) - len_str_estimation):
                    if (exo_p[i][6] == "Validated") | (exo_p[i][6] == "Not_validated"):
                        # print("Ok recall TO : ", 100*exo_p_float[i][1], " correspond à ", to_nearest)

                        if abs(to_nearest - 100 * exo_p_float[i][1]) < 0.2 * len_str_estimation:
                            if exo_p[i][6] == "Validated":
                                validate.append(1)
                                # print("validate", len(validate), exo_p[i][6], validate)
                            else:
                                print("Intermédiaire :", exo_p[i][6], 0.5)
                                validate.append(0.5)
                                # print("validate", len(validate), exo_p[i][6], validate)
                        else:
                            print("Problème, trop loin recall TO :", to_nearest, 100 * exo_p_float[i][1], 0.2 * len_str_estimation)
                            validate.append(0)
                            # print("validate", len(validate), exo_p[i][6], validate)

                        if not (abs(to_nearest - 100 * exo_p_float[i][1]) <= abs(
                                hs_nearest - 100 * exo_p_float[i][1])) | (
                                       (hs_nearest - 100 * exo_p_float[i][1]) * (
                                       to_nearest - 100 * exo_p_float[i][1]) < 0):
                            print("Recall TO :", 100 * exo_p_float[i][1], " plus proche de HS ", hs_nearest,
                                  " que de TO ",
                                  to_nearest)
                            print(abs(to_nearest - 100 * exo_p_float[i][1]), " plus grand que ",
                                  abs(hs_nearest - 100 * exo_p_float[i][1]),
                                  " et ", hs_nearest - 100 * exo_p_float[i][1], " de même signe que ",
                                  to_nearest - 100 * exo_p_float[i][1])
                            print("Mais dans les limites :", len_str_estimation, 100 * exo_p_float[i][1],
                                  len(data) - len_str_estimation)
                            # event_detect_corresp_type.append("HS")
                            # event_detect_corresp_value.append(hs_nearest)
                            # validate.append(0)
                    else:
                        if exo_p[i][6] == "Deleted":
                            #print("Suppression du pas de mauvaise qualité TO : ", exo_p_float[i][1])
                            validate.append(-1)
                            # print("validate", len(validate), exo_p[i][6], validate)

                else:
                    # print("Attention trop proche des limites TO", 100 * exo_p_float[i][1], len(data), len_str_estimation)
                    validate.append(-2)
                    # print("validate", len(validate), exo_p[i][6], validate)
            else:
                # print("U_turn TO", u_turn_start, exo_p_float[i][1], u_turn_end)
                validate.append(-3)
                # print("validate", len(validate), exo_p[i][6], validate)

    # print("len", len(event_tapis_type), len(event_tapis_value), len(event_detect_corresp_type), len(event_detect_corresp_value), len(validate))

    return np.array(
        [event_tapis_type, event_tapis_value, event_detect_corresp_type, event_detect_corresp_value, validate])


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def extract_len(exo, u_turn_start=0, u_turn_end=0):
    len_str = []
    len_swg = []
    for j in range(len(exo) - 1):
        if (exo[j + 1][0] < u_turn_start) | (exo[j + 1][0] > u_turn_end):
            if (exo[j][0] < u_turn_start) | (exo[j][0] > u_turn_end):
                len_str.append(exo[j + 1][0] - exo[j][0])
            if (exo[j][1] < u_turn_start) | (exo[j][1] > u_turn_end):
                len_swg.append(exo[j + 1][0] - exo[j][1])
        if (exo[j + 1][1] < u_turn_start) | (exo[j + 1][1] > u_turn_end):
            if (exo[j][1] < u_turn_start) | (exo[j][1] > u_turn_end):
                len_str.append(exo[j + 1][1] - exo[j][1])

    len_str_estimation = round(np.percentile(len_str, 50) * 100, 0)
    len_swg_estimation = max(round(np.percentile(len_swg, 50) * 100, 0), 40)

    return len_str_estimation, len_swg_estimation
