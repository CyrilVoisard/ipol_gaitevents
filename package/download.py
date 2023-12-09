from datetime import datetime


def json_report(seg_lim, steps_lim, freq, output):
    dict_events = dict()
    dict_events["Detection date"] = datetime.now()
    dict_events["Freq"] = freq
    
    dict_events["TrialBoundaries"] = [int(seg_lim[0]), int(seg_lim[3])]
    dict_events["UTurnBoundaries"] = [int(seg_lim[1]), int(seg_lim[2])]

    steps_lim_lf = []
    steps_lim_rf = []
    for i in range(len(steps_lim)):
        if steps_lim["Corrected"][i]==1:
            if steps_lim["Foot"][i]==1:
                steps_lim_rf.append([int(steps_lim["TO"][i]), int(steps_lim["HS"][i])])
            else:
                steps_lim_lf.append([int(steps_lim["TO"][i]), int(steps_lim["HS"][i])])
    dict_events["LeftFootEvents"] = steps_lim_lf
    dict_events["RightFootEvents"] = steps_lim_rf
    
    with open(os.path.join(output, "gait_events.json"), "w") as f:
        json.dump(dict_events, f)
