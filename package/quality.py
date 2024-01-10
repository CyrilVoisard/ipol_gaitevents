import os 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats


def print_all_quality_index(q1, q2, output):
    """Compute the quality index of the trial gait events detection (between 0 and 100) and produce a picture of the number surrounded by an appropriately colored circle. 
    Add quality index formula ? 

    Parameters
    ----------
        steps_lim_full {dataframe} -- pandas dataframe with all the detected gait events
        seg_lim {dataframe} -- pandas dataframe with phases events 
        output {str} -- folder path for output fig

    Returns
    -------
        qi {int} -- corrected quality index 
        steps_lim_corrected {dataframe} -- pandas dataframe with gait events after elimination of the extra trial steps
    """

    q_mean = round(1/2 * (np.mean(q1) + np.mean(q2)))

    fig = plt.figure(figsize=(6, 5))
    gs = GridSpec(nrows=2, ncols=2, width_ratios = [2, 1])
    ax0 = fig.add_subplot(gs[:, 0], projection='polar')
    ax0 = plot_quality_index(q_mean, ax0, scale = 1)
    ax0.text(0.22, 0.79, 'Global quality score', fontsize = 14, fontweight='bold', transform=plt.gcf().transFigure)
    ax1 = fig.add_subplot(gs[0, 1], projection='polar')
    ax1 = plot_quality_index(round(np.mean(q1)), ax1, scale = 2)
    ax1.text(0.7, 0.85, 'Protocol quality score', fontsize = 11, fontweight='bold', transform=plt.gcf().transFigure)
    ax2 = fig.add_subplot(gs[1, 1], projection='polar')
    ax2 = plot_quality_index(round(np.mean(q2)), ax2, scale = 2)
    ax2.text(0.7, 0.425, 'Intrinsec quality score', fontsize = 11, fontweight='bold', transform=plt.gcf().transFigure)

    path = os.path.join(output, "quality_index.svg")

    plt.savefig(path, dpi=80, transparent=True, bbox_inches="tight")
  

def plot_quality_index(q, ax, scale):
    """Compute the quality index of the trial gait events detection (between 0 and 100) and produce a picture of the number surrounded by an appropriately colored circle. 

    Parameters
    ----------
        q {int} -- quality index 
        ax {} -- ax in output fig
    """
  
    # plot qi
    max_q=100
    xval = np.arange(0, 2*np.pi*(.05+0.90*(q/max_q)), 0.01)
    colormap = plt.get_cmap("Greens")
    norm = mpl.colors.Normalize(0.0, 2*np.pi)
    # f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4),subplot_kw=dict(projection='polar'))
    #Scatter version
    yval = np.ones_like(xval)
    ax.scatter(xval, yval, c=xval, s=300/(2*scale), cmap=colormap, norm=norm, linewidths=1)
  
    ax.set_axis_off()
    ax.set_ylim(0,1.5)
    if q<10:
        ax.annotate(q, xy=( 1.25*np.pi, .3), color=colormap(.05+0.90*(q/max_q)), fontsize=50/scale)
    else :
        if  q == 100: 
            ax.annotate(q, xy=(1.11*np.pi, .7), color=colormap(.05+0.90*(q/max_q)), fontsize=50/scale)
        else:
            ax.annotate(q, xy=(1.18*np.pi, .5), color=colormap(.05+0.90*(q/max_q)), fontsize=50/scale)
  
    return ax
    

def compute_extrinsic_quality(steps_lim):
  """Compute the quality index of the trial gait events detection (between 0 and 100) referring to extrinsic quality: right-left step alternation

    Parameters
    ----------
        steps_lim {dataframe} -- pandas dataframe with gait events

    Returns
    -------
        qi {int} -- quality index 
    """
  
  # estimation of stride alternation
  steps_lim_sort = steps_lim.sort_values(by = ['HS', 'TO'])
  foot_list = steps_lim_sort['Foot'].tolist()
  i = 0
  for k in range(len(foot_list)-1):
      i = i + abs(foot_list[k+1]-foot_list[k])
  q3 = round(100*i/(len(foot_list) -1))

  return q3
