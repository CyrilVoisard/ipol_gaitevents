import os 
import json
from datetime import datetime


def json_report(steps_lim, freq, output):
    dict_events = dict()
    dict_events["Detection date"] = str(datetime.now())
    dict_events["Freq"] = freq
    
    steps_lim_lf = []
    steps_lim_rf = []
    for i in range(len(steps_lim)):
        if steps_lim["Foot"][i]==1:
            steps_lim_rf.append([int(steps_lim["TO"][i]), int(steps_lim["HS"][i])])
        else:
            steps_lim_lf.append([int(steps_lim["TO"][i]), int(steps_lim["HS"][i])])
    dict_events["LeftFootEvents"] = steps_lim_lf
    dict_events["RightFootEvents"] = steps_lim_rf
    
    with open(os.path.join(output, "gait_events.json"), "w") as f:
        json.dump(dict_events, f)
