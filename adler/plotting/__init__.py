import matplotlib

#matplotlib.rcParams['lines.markeredgecolor'] = 'black' #matplotlib.rcParams.keys()
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.5

matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 0.25
matplotlib.rcParams['xtick.minor.size'] = 10
matplotlib.rcParams['xtick.minor.width'] = 0.25


from adler.plotting.plotting_time_series import plot_time_series, plot_time_series_square,plot_time_series_square_ou
from adler.plotting.plotting_time_series import plot_time_series_square_dist


from adler.plotting.plotting_fft import plot_fft_all, plot_fft
from adler.plotting.plotting_fpt import plot_FPT_all, plot_FPT, plot_FPT_square,plot_FPT_square_ou

from adler.plotting.plotting_pulses import  plot_pulses_square
from adler.plotting.plotting_pulses import  plot_dt_square,plot_ipi_square,plot_dm_square,plot_joint_duration_square
from adler.plotting.plotting_pulses import  plot_pulses_square_ou, plot_dt_square_ou,plot_ipi_square_ou,plot_dm_square_ou,plot_joint_duration_square_ou
#plot_pulses,plot_activity, plot_consecutiveness, plot_2d_mean_activity,plot_quantifiers_histograms

from adler.plotting.plotting_pulses import plot_2d_quantifiers,plot_2d_superposition
from adler.plotting.plotting_consecutive import plot_consecutiveness_activity,plot_consecutiveness_activity_dist


from adler.plotting.plotting_dt_alpha import plot_simulated_dt_alpha,plot_theta_alpha, plot_simulated_fixed_dt_alpha,plot_theta_fixed_alpha
from adler.plotting.plotting_simulated_dt_noise import plot_epsilon_plus,plot_epsilon_plus_in_x_minus,plot_t_plus,plot_t_plus_in_x_minus