import os 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap


def plot_quality_index(q, output):
    """Compute the quality index of the trial gait events detection (between 0 and 100) and produce a picture of the number surrounded by an appropriately colored circle. 

    Parameters
    ----------
        q {float} -- alternation quality index.
        output {str} -- folder path for output figure. 
    """

    # plot qi
    max_q = 100
    xval = np.arange(0, 2*np.pi*(.05+0.90*(q/max_q)), 0.01)
    norm = mpl.colors.Normalize(0.0, 2*np.pi)   
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 4), subplot_kw=dict(projection='polar'))

    # Définir les couleurs en fonction de la valeur de q
    if q < 50:
        colormap = plt.get_cmap("Reds")
    elif 50 <= q < 75:
        colormap = plt.get_cmap("Yellows")
    else:
        colormap = plt.get_cmap("Greens")

    #Scatter version    
    yval = np.ones_like(xval)    
    ax.scatter(xval, yval, c=xval, s=150, cmap=colormap, norm=norm, linewidths=1)      

    # Remplir le cercle avec un dégradé
    # theta = np.linspace(0, 2 * np.pi * (q / max_q), 100)
    # radius = np.ones_like(theta)
    # bars = ax.bar(theta, radius, width=0.05, bottom=0.0, color=cmap(theta / (2 * np.pi)), edgecolor='none')

    ax.set_axis_off()
    ax.set_ylim(0, 1.5)

    # Positionnement du texte en fonction de la valeur de q
    if q<10:
        ax.annotate(q, xy=(1.25*np.pi, .3), color=colormap(.05+0.90*(q/max_q)), fontsize=50)    
    else :
        if  q == 100:
            ax.annotate(q, xy=(1.11*np.pi, .7), color=colormap(.05+0.90*(q/max_q)), fontsize=50)
        else:
            ax.annotate(q, xy=(1.18*np.pi, .5), color=colormap(.05+0.90*(q/max_q)), fontsize=50)
    
    fig.suptitle('Quality score', fontsize = 14, fontweight='bold')
    path = os.path.join(output, "quality_index.svg")

    plt.savefig(path, dpi=300, transparent=True, bbox_inches="tight")
    

def compute_quality(steps_lim):
  """Compute the quality index of the trial gait events detection (between 0 and 100) referring to extrinsic quality: right-left step alternation

    Parameters
    ----------
        steps_lim {dataframe} -- pandas dataframe with gait events.

    Returns
    -------
        q {int} -- quality index. 
    """
  
  # estimation of stride alternation
  steps_lim_sort = steps_lim.sort_values(by = ['HS', 'TO'])
  foot_list = steps_lim_sort['Foot'].tolist()
  i = 0
  for k in range(len(foot_list)-1):
      i = i + abs(foot_list[k+1]-foot_list[k])
  q = round(100*i/(len(foot_list) -1))

  return q
