import matplotlib.pyplot as plt
import numpy as np

def boxoff():
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

def find_nearest(array,value):
    # http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    idx = (np.abs(array-value)).argmin()
    return idx

def middles(edges):
    # make middles vector from edges vector, for use with binit
    diffs = np.median(np.diff(edges));
    edgestart = edges[0] + diffs/2.;
    edgeend = edges[-1] - diffs/2.;
    mid = np.arange(edgestart, edgeend, diffs)
