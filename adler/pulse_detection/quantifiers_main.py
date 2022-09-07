import numpy as np
import multiprocessing as mp
from math import ceil
from functools import partial 
import pandas as pd

#from adler.pulse_detection.pulse_detection_main import test_pulses_sine
from adler.data_managing_functions import check_file, save_data, download_data,points_to_time

#%% 
# =============================================================================
# =============================================================================
# =============================================================================
#               COMPUTING PULSES QUANTIFIERS MODULE
# =============================================================================
# =============================================================================
# =============================================================================


def test_pulses_quantifiers(dt,IPI,dm,joint_duration,MAX):
    
    assert len(dm) == len(IPI)
    assert len(joint_duration) == len(IPI)
    assert len(dt) == len(MAX)
    
    assert all([(i + j == k)*1 for (i,j,k) in zip(dm,joint_duration,IPI)]), 'dm:'+str(dm) + ' joint_duration: ' + str(joint_duration) + ' IPI: '+str(IPI)
    
    
    assert all([i >= 0 for i in dm]),dm
    assert all([i > 0 for i in dt])
    assert all([i > 0 for i in IPI])
    assert all([i > 0 for i in joint_duration])


def get_pulses_quantifiers(left_minima,right_minima,MAX):
    
    ''' Falta descripcion
    calcula a los quantifiers en indices'''
    IPI = []; dt = []; dm = []; joint_duration = []
    
    #pulse duration
    for left,right in zip(left_minima,right_minima):
        assert (right-left >= 0), 'left_minima: ' +str(left_minima) + ' right_minima : '+str(right_minima)
        dt.append(right - left)
    
    for M1,M2 in zip(MAX[:-1],MAX[1:]):
        # interpulse interval
        assert (M2-M1 >= 0)
        IPI.append(M2-M1)
        
        right_filter = list(filter(lambda x: (x >= M1 and x<=M2), right_minima))[0]
        left_filter = list(filter(lambda x: (x <= M2 and x>=M1), left_minima))[-1]
        
        
        #joint duration
        assert ((right_filter-M1) + (M2 -left_filter) <= M2-M1), str(M1)+ ' --- ' + str(M2) + ' --- ' + str(left_filter) + ' --- '+str(right_filter)
        assert ((right_filter-M1) + (M2 -left_filter) > 0)
        joint_duration.append((right_filter-M1) + (M2 -left_filter)) 
    
    assert all([(i<j)*1 for (i,j) in zip(left_minima,right_minima)])  
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]) 
    
    for right,left in zip(right_minima[:-1],left_minima[1:]):
        
        #silent interval
        assert (left - right >= 0),str(left) + ' ' +str(right)+'left_minima: ' +str(left_minima) + ' right_minima : '+str(right_minima)
        dm.append(left - right)
    
    test_pulses_quantifiers(dt,IPI,dm,joint_duration,MAX)
    return(dt,IPI,dm,joint_duration)
    
    
def get_pulses_quantifiers_(data_folder,save_path_name,tuple_):
    '''
    data folder: donde están los pulsos
    '''

    if 'alpha' in tuple_[1].columns:
        (i,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        delta = np.round(i/omega,4)  
        file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    elif 'alpha0' in tuple_[1].columns:
        (i,_,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        delta = np.round(i/omega,4)  
        file_name =  str(int(row.number.values[0]))+'_'+str(int(order))+'.pkl'
 
    
    if (check_file('max_'+file_name,data_folder)):
        
        print('running pulses quantifiers computation')      
        
        MAX = download_data(data_folder + 'max_'+file_name) 
        left_minima = download_data(data_folder + 'left_minima_'+ file_name) 
        right_minima = download_data(data_folder + 'right_minima_'+ file_name) 
        
        dt,IPI,dm,joint_duration = get_pulses_quantifiers(left_minima,right_minima,MAX)
        
        print('Ready! Saving files --- ')      
        save_data(dt,save_path_name+'dt_'+file_name)
        save_data(IPI,save_path_name+'IPI_'+file_name)
        save_data(dm,save_path_name+'dm_'+file_name)
        save_data(joint_duration,save_path_name+'joint_duration_'+file_name)
        print('pulses quantifiers computation finished :)')
    else:
        print(file_name,'maxima file not available')

    return(1)
   
    #%%
    
    
def compute_pulses_quantifiers(description_file,data_folder,save_path_name):
    '''
    data folder: donde están los pulsos
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    if 'alpha' in ref.keys(): tuple_ = ref.groupby(['alpha','D','order'])
    if 'alpha0' in ref.keys(): tuple_ = ref.groupby(['alpha0', 'sigma', 'tau','order'])
    get_pulses_quantifiers__ = partial(get_pulses_quantifiers_,data_folder,save_path_name)
    pool.map(get_pulses_quantifiers__,tuple_)
    pool.close()
    pool.join()
    return (2)



#%% calcula el pulse rate en pulsos sobre minutos! PARTE LA SERIE TEMPORAL ES PARA EL DOS D PLOT



def compute_pulse_rate(T,dt,d,description_file,data_folder,save_path_name):
    ''' split_flag: es para tener una serie temporal larga, y partirla en pedazos'''
    
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    if 'alpha' in ref.keys(): tuple_ = ref.groupby(['alpha','D','order'])
    if 'alpha0' in ref.keys(): tuple_ = ref.groupby(['alpha0', 'sigma', 'tau','order'])
    
    get_pulse_rate_ = partial(get_pulse_rate,T,dt,d,data_folder,save_path_name)
    pool.map(get_pulse_rate_,tuple_)
    pool.close()
    pool.join()

    return (2)

def get_pulse_rate(T,dt,d,data_folder,save_path_name,tuple_):
    '''
    data folder: donde están los pulsos
    '''
    N = 50
    if 'alpha' in tuple_[1].columns:
        (_,_,order),row = tuple_[0],tuple_[1]
        file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    elif 'alpha0' in tuple_[1].columns:
        (_,_,_,order),row = tuple_[0],tuple_[1]
        file_name =  str(int(row.number.values[0]))+'_'+str(int(order))+'.pkl'
 
    #esto es para tener una serie temporal larga, y partirla en pedazos, en principio en 50
   
    if (check_file('max_'+file_name,data_folder)):
        
        print('running pulse rate st computation')      
        pulse_rate = pulse_rate_statistics(download_data(data_folder+'max_'+file_name),np.arange(ceil(int(T/dt)/d)),int(ceil(int(T/dt)/d)/N),dt,d) 
        save_data(pulse_rate,save_path_name+'pr_'+file_name)
        print('pulse rate st computation finished :) len pulse rate',len(pulse_rate))
    else:
        pulse_rate = np.zeros(N)
        save_data(pulse_rate,save_path_name+'pr_'+file_name)
        print(file_name,'maxima file not available')

    return(1)
    
    
    
def split_len_N(ix,N):
    '''te parte ix en elementos de N elementos. El ultimo tiene falopa'''
    aux = []; i = 0
    while i < len(ix):
        aux.append(i)
        i = i + N
    return  np.split(ix,(aux[1:]))


def pulse_rate_statistics(MAX,ix,N,dt,d):
    '''ix es una lista de indices de theta, N esta en puntos
    Te da el pulse rate de cada parte de N puntos de ix (indices) de la TS'''
    pulse_rate_aux = []
        
    for ix_i in split_len_N(ix,N):    
        pulses =  list(filter(lambda x : x in MAX,ix_i))
        t_i = points_to_time(ix_i,dt,d)
        if len(ix_i) > 1: #con esto anda, pero no sería lo más preciso del mundo...
            pulse_rate_aux.append(len(pulses)/(t_i[-1]-t_i[0]))
    return pulse_rate_aux




