import matplotlib.pyplot as plt
import numpy as np

def boxoff():
    """
    A Matlab-like boxoff to remove top & right border of plots
    """
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

def twinboxoff():
    """
    A Matlab-like boxoff for twinx plots
    """
    plt.gca().spines['top'].set_visible(False)

def thinspines(lw=0.5):
    for axis in ['top','bottom','left','right']:
      plt.gca().spines[axis].set_linewidth(lw)
      plt.gca().tick_params(width=lw)

def find_nearest(array,value):
    """
    Find nearest value in numpy array
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx

def middles(edges):
    """
    Make middles vector from edges vector
    """
    diffs = np.median(np.diff(edges));
    edgestart = edges[0] + diffs/2.;
    edgeend = edges[-1] - diffs/2.;
    mid = np.arange(edgestart, edgeend + diffs, diffs)

    return mid

def set_fontsize(fig,fontsize):
    """
    For each text object of a figure fig, set the font size to fontsize
    http://stackoverflow.com/questions/7082597/in-matplotlib-how-do-you-change-the-fontsize-of-a-single-figure
    """
    def match(artist):
        return artist.__module__ == "matplotlib.text"

    for textobj in fig.findobj(match=match):
        textobj.set_fontsize(fontsize)

def nextcolor():
    """
    Get the next color in the default matplotlib color order
    """
    next(plt.gca()._get_lines.prop_cycler)['color']
