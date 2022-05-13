import numpy as np
import multiprocessing as mp
from math import ceil
from functools import partial 
import pandas as pd

#from adler.pulse_detection.pulse_detection_main import test_pulses_sine
from adler.data_managing_functions import check_file, save_data, download_data

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


def get_pulses_quantifiers(left_minima,right_minima,MAX):
    
    ''' Falta descripcion'''
    IPI = []; dt = []; dm = []; joint_duration = []
    
    for left,right in zip(left_minima,right_minima):
        assert (right-left >= 0), 'left_minima: ' +str(left_minima) + ' right_minima : '+str(right_minima)
        dt.append(right - left)
        
    for M1,M2 in zip(MAX[:-1],MAX[1:]):
        assert (M2-M1 >= 0)
        IPI.append(M2-M1)
        
        right_filter = list(filter(lambda x: (x >= M1 and x<=M2), right_minima))[0]
        left_filter = list(filter(lambda x: (x <= M2 and x>=M1), left_minima))[-1]
        
        ##### PROBLEMA
        #lo que está pasando es que el maximo de la derecha es tambien un minimo (el left)
        # eso lo vimos como list(filter(lambda x: (x > M1 and x < M2), left_minima))
        
        assert ((right_filter-M1) + (M2 -left_filter) <= M2-M1), str(M1)+ ' --- ' + str(M2) + ' --- ' + str(left_filter) + ' --- '+str(right_filter)
        assert ((right_filter-M1) + (M2 -left_filter) > 0)
        joint_duration.append((right_filter-M1) + (M2 -left_filter)) 
    
    assert all([(i<j)*1 for (i,j) in zip(left_minima,right_minima)])  
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]) 
    for right,left in zip(right_minima[:-1],left_minima[1:]):
        assert (left - right >= 0),str(left) + ' ' +str(right)+'left_minima: ' +str(left_minima) + ' right_minima : '+str(right_minima)

        dm.append(left - right)
    
    test_pulses_quantifiers(dt,IPI,dm,joint_duration,MAX)
    return(dt,IPI,dm,joint_duration)
    
    
def get_pulses_quantifiers_(data_folder,save_path_name,tuple_):
    '''
    data folder: donde están los pulsos
    '''

    (i,D,order),row = tuple_[0],tuple_[1]
    omega =  row.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    if (check_file('max_xf_'+file_name,data_folder)):
        
        print('running pulses quantifiers computation alpha, D : ',alpha,D)      
        
        MAX = download_data(data_folder + 'max_xf_'+file_name) 
        #MIN = download_data(data_folder + 'min_xf_'+ file_name) 
        left_minima = download_data(data_folder + 'left_minima_'+ file_name) 
        right_minima = download_data(data_folder + 'right_minima_'+ file_name) 
        test_pulses_sine(left_minima,right_minima,MAX)
        
        dt,IPI,dm,joint_duration = get_pulses_quantifiers(left_minima,right_minima,MAX)
        
        print('Ready! Saving files --- alpha, d : ',alpha,D)      
        save_data(dt,save_path_name+'dt_xf_'+file_name)
        save_data(IPI,save_path_name+'IPI_xf_'+file_name)
        save_data(dm,save_path_name+'dm_xf_'+file_name)
        save_data(joint_duration,save_path_name+'joint_duration_xf_'+file_name)
        print(alpha,D,'pulses quantifiers computation finished :)')
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
    
    tuple_ = ref.groupby(['alpha','D','order'])
    get_pulses_quantifiers__ = partial(get_pulses_quantifiers_,data_folder,save_path_name)
    pool.map(get_pulses_quantifiers__,tuple_)
    pool.close()
    pool.join()
    return (2)




