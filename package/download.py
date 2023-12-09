from datetime import datetime


def json_report(seg_lim_corrected, steps_lim_corrected, freq, output):
    dict_events = dict()
    dict_events["Detection date"] = datetime.now()
    dict_events["Freq"] = freq
    
    dict_events["TrialBoundaries"] = [int(seg_lim_corrected[0]), int(seg_lim_corrected[3])]
    dict_events["UTurnBoundaries"] = [seg_lim_corrected[1]), seg_lim_corrected[2])]

    strides = steps_lim_corrected[steps_lim_corrected["Correct"] == True]
    strides_lf = []
    strides_rf = []
    for i in range(len(strides)):
        if strides["Foot"][i]==1:
            strides_rf.append([int(strides["TO"][i]), int(strides["HS"][i])])
        else:
            strides_lf.append([int(strides["TO"][i]), int(strides["HS"][i])])
    dict_events["LeftFootEvents"] = strides_lf
    dict_events["RightFootEvents"] = strides_rf
    
    with open(os.path.join(output, "gait_events.json"), "w") as f:
        json.dump(dict_events, f)
