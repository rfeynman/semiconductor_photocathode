import matplotlib as mp
from commands import getoutput as go
import matplotlib.pyplot as plt
from numpy import *
def mhist1D(x_val, scale_x = 1.0, n_bins = 10, 
            facecolor = 'red', x_label = '', y_label = '',
            plot_legend = '', grid_flag = False, backend = 'png',
            root_file_name = 'expM1D', open_flag = False) :
    fig_width_pt = 360. # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean*1.3  # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': backend, 
              'text.usetex': True,
              'axes.labelsize': 11,
              'font.size': 10,
              'font.weight' : 'bold',
              'legend.fontsize': 10,
              #'xtick.fontweight' : 'bold',
              'xtick.labelsize': 10, 
              'ytick.labelsize': 10,
              'figure.figsize': fig_size}
    x = array(x_val)*scale_x
    mp.rc('axes', linewidth=2)
    mp.rcParams.update(params)
    plt.figure(1)
    #plt.plot(x, y, markershape)
    n, bins, patches = plt.hist(x, bins = n_bins, facecolor=facecolor, alpha=0.75)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if len(plot_legend) > 0 :
        plt.legend(plot_legend)
    plt.grid(grid_flag)
    if open_flag :
      file_name = root_file_name + '.' + backend
      plt.savefig(file_name)
      plt.close(1)
      go('open ' + file_name)
    return

