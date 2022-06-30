from adler.plotting.dyncode_suorces import load_file,consecutive_cumulative,consecutive_non_cumulative
import numpy as np



#green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
#colors =  [green,sns.color_palette()[1],sns.color_palette()[3]]


#consecutivenes plot
def get_conc_data(dyncode_file_name):
    return load_file(dyncode_file_name)

def get_consecutive_data_dyncode(dyncode_file_name):
    ''' para el plot de consecutividad'''
    df_consecutive = get_conc_data(dyncode_file_name)[ 'an_WT_ESL']
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    norm = len(df_consecutive.index.get_level_values(0).unique())
    return np.arange(1,len(box_plot_consecutive_cumulative)+1),[i/norm for i in box_plot_consecutive_cumulative]


def get_activity_data_dyncode(dyncode_file_name):
    ''' para el plot de population activity'''
    df_consecutive = get_conc_data(dyncode_file_name)['an_WT_ESL']
    activity_experiment = df_consecutive['dt_peaks'].groupby(level='cell').sum() / df_consecutive['FRAME'].groupby(level='cell').count() *  100   
    activity_experiment_index = np.argsort(activity_experiment.values)[::-1]
    activity_experiment = [activity_experiment[j] for j in activity_experiment_index]
    
    silent_experiment = np.ones(len(activity_experiment)) * 100 - activity_experiment
    return np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),activity_experiment,silent_experiment


def get_exp_N_total_isolated_consecutive(dyncode_file_name):
    ''' para el boxplot de consecutividad'''
    df_consecutive = get_conc_data(dyncode_file_name)[ 'an_WT_ESL']
    consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
    
    isolated_N = sum(consecutive_non_cumulative_obj.is_isolated_box())
    consecutive_N = sum(consecutive_non_cumulative_obj.count_consecutive_pulses_number())
    total_N = sum(consecutive_non_cumulative_obj.get_number_of_pulses())
    return total_N,isolated_N,consecutive_N
    
def get_dyncode_pulse_rate_st(dyncode_file_name):
    '''esto es para los dos d plots'''
    df = get_conc_data(dyncode_file_name)['an_WT_ESL']
    pulse_density = []
    for cell,data in df.groupby(level="cell"):
        pulse_density.append(data.amp_peaks.count() / len(data.FRAME))
        np.median(pulse_density)
    return np.percentile(pulse_density,25), np.percentile(pulse_density,75)
