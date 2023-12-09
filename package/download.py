from datetime import datetime


def json_report(seg_lim, steps_lim, freq, output):
    dict_events = dict()
    dict_events["Detection date"] = datetime.now()
    dict_events["Freq"] = freq
    
    dict_events["TrialBoundaries"] = [int(seg_lim[0]), int(seg_lim[3])]
    dict_events["UTurnBoundaries"] = [int(seg_lim[1]), int(seg_lim[2])]

    strides = steps_lim[steps_lim["Correct"] == 1]
    print(len(strides)) #, strides["Foot"][9])
    strides_lf = []
    strides_rf = []
    for i in range(len(strides)):
        print("i", i)
        if strides["Foot"][i]==1:
            strides_rf.append([int(strides["TO"][i]), int(strides["HS"][i])])
        else:
            strides_lf.append([int(strides["TO"][i]), int(strides["HS"][i])])
    dict_events["LeftFootEvents"] = strides_lf
    dict_events["RightFootEvents"] = strides_rf
    
    with open(os.path.join(output, "gait_events.json"), "w") as f:
        json.dump(dict_events, f)
