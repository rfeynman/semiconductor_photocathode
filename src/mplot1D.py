import matplotlib as mp
from commands import getoutput as go
import matplotlib.pyplot as plt
from numpy import *
def mplot1D(x_val, y_val, scale_x = 1.0, scale_y = 1.0,
            markershape = 'rD-', fillstyle='full', x_label = '', y_label = '',
            xlimt = (), ylimt = (), slogx=False, slogy=False, mlinewidth = 2.8,
            plot_legend = '', plot_location = 4, plot_title = '', plot_text = '', text_x = 0.83, text_y = 0.8, 
            grid_flag = False, backend = 'pdf',
            root_file_name = 'expM1D', open_flag = False) :
    fig_width_pt = 550. # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean #*1.43  # height in inches
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
    y = array(y_val)*scale_y
    mp.rc('axes', linewidth=2)
    mp.rcParams.update(params)
    plt.figure(1)
    if (slogx and slogy) :
        plt.loglog(x, y, markershape, fillstyle=fillstyle, linewidth=mlinewidth)
    elif slogx :
        plt.semilogx(x, y, markershape, fillstyle=fillstyle, linewidth=mlinewidth)
    elif slogy :
        plt.semilogy(x, y, markershape, fillstyle=fillstyle, linewidth=mlinewidth)
    else : 
        plt.plot(x, y, markershape, fillstyle=fillstyle, linewidth=mlinewidth)
    if xlimt != () :
        plt.xlim(xlimt)
    if ylimt != () :
        plt.ylim(ylimt)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if len(plot_legend) > 0 :
        plt.legend(plot_legend, loc=plot_location)
    plt.grid(grid_flag)
    if len(plot_title) > 0 : 
        plt.title(plot_title)
    if (len(plot_text) > 0) :
        plt.text(text_x, text_y, plot_text, fontsize = 13)
    if open_flag :
      file_name = root_file_name + '.' + backend
      plt.savefig(file_name)
      plt.close(1)
      go('open ' + file_name)
    return plt

