import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.5

matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 0.25
matplotlib.rcParams['xtick.minor.size'] = 10
matplotlib.rcParams['xtick.minor.width'] = 0.25


from plotting.plotting_time_series import plot_time_series
from plotting.plotting_fft import plot_fft_all, plot_fft
from plotting.plotting_fpt import plot_FPT_all, plot_FPT
from plotting.plotting_pulses import  plot_pulses,plot_quantifiers_histograms,plot_activity, plot_consecutiveness