import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def compute_quality_index(steps_lim, seg_lim):
  """Compute the quality index of the trial gait events detection (between 0 and 100)
  Add quality index formula ? 

    Parameters
    ----------
        steps_lim {dataframe} -- pandas dataframe with gait events
        seg_lim {dataframe} -- pandas dataframe with phases events 

    Returns
    -------
        qi {int} -- quality index 
    """
  
  # estimation of stride alternation
  steps_lim_sort = steps_lim.sort_values(by = ['HS', 'TO'])
  alt_go = steps_lim_sort[steps_lim_sort['HS'] < seg_lim[1]]['Foot'].tolist()
  alt_back = steps_lim_sort[steps_lim_sort['HS'] > seg_lim[2]]['Foot'].tolist()
  i = 0
  for k in range(len(alt_go)-1):
      i = i + abs(alt_go[k+1]-alt_go[k])
  for k in range(len(alt_back)-1):
      i = i + abs(alt_back[k+1]-alt_back[k])
  qi = round(100*i/(len(alt_go) + len(alt_back)-2))

  return qi
    

def print_quality_index(steps_lim_full, seg_lim, output):
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

    qi_raw = compute_quality_index(steps_lim_full, seg_lim)
    path_raw = os.path.join(output, "quality_index_raw.svg")
    plot_quality_index(qi_raw, path_raw)
  
    steps_lim_corrected = steps_lim_full
  
    qi_corrected = compute_quality_index(steps_lim_corrected, seg_lim)
    path_corrected = os.path.join(output, "quality_index_raw.svg")
    plot_quality_index(qi_corrected, path_corrected)
  
    return qi_corrected, steps_lim_corrected


def plot_quality_index(qi, path):
    """Compute the quality index of the trial gait events detection (between 0 and 100) and produce a picture of the number surrounded by an appropriately colored circle. 
    Add quality index formula ? 

    Parameters
    ----------
        qi {int} -- quality index 
        path {str} -- file path for output fig
    """
  
    # plot qi
    max_qi=100
    xval = np.arange(0, 2*np.pi*(.05+0.90*(qi/max_qi)), 0.01)
    colormap = plt.get_cmap("RdYlGn")
    norm = mpl.colors.Normalize(0.0, 2*np.pi)
    f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4),subplot_kw=dict(projection='polar'))
    #Scatter version
    yval = np.ones_like(xval)
    ax.scatter(xval, yval, c=xval, s=300, cmap=colormap, norm=norm, linewidths=0)
  
    ax.set_axis_off()
    ax.set_ylim(0,1.5)
    if qi<10:
        ax.annotate(qi,xy=( 1.25*np.pi, .3),color=colormap(.05+0.90*(qi/max_qi)),fontsize=50)
    else :
        ax.annotate(qi,xy=( 1.18*np.pi, .5),color=colormap(.05+0.90*(qi/max_qi)),fontsize=50)
  
    plt.savefig(path, dpi=80,
                    transparent=True, bbox_inches="tight")
