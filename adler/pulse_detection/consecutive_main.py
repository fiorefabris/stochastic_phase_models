import numpy as np
from adler.data_managing_functions import check_file, save_data, download_data
from itertools import groupby
import pandas as pd
from math import ceil
import multiprocessing as mp
from functools import partial 


def count_iterable(iterator):
    return sum(1 for i in iterator)

#%%
################################################
#### Module for computing the consecutiveness
################################################  

def get_consecutive(description_file,data_folder,save_folder):
    '''
    data folder: donde están los dt
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha','D'])
    get_consecutive_trial_ = partial(get_consecutive_trial,data_folder,save_folder)
    pool.map(get_consecutive_trial_,tuple_)
    pool.close()
    pool.join()
    return (2)

def get_consecutive_trial(data_folder,save_folder,tuple_):
    
    (i,D),row_ = tuple_[0],tuple_[1]
    
    omega =  row_.omega.unique()[0]
    delta = np.round(i/omega,4)  
    
    for (order,row) in row_.groupby(['order']):
        number      = int(row.number)
        file_name   =  str(number)+'_'+str(order)+'.pkl'
            
        if (check_file('max_'+file_name,data_folder)):       
            consecutive_trial_st_ = consecutive_trial_st(number,order,data_folder)
            isolated_pulses, consecutive_trial = consecutive_trial_st_.get_consecutive_trains_of_pulses()#
                                    
            save_data(consecutive_trial,save_folder+'c_'+file_name)
            save_data(isolated_pulses,save_folder+'i_'+file_name)
        else:
            pass
    return(0)



#%%
class consecutive_trial_st:
    '''
    esto es para los boxplots. Quiero la cantidad de pulsos consecutivos, aislados y totales de un dado traza
    '''
    def __init__(self,number,order,data_folder):
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        self.joint_duration = download_data(data_folder+'joint_duration_'+file_name)
        self.max = download_data(data_folder+'max_'+file_name)
        self.IPI = download_data(data_folder+'IPI_'+file_name)
        self.dm = download_data(data_folder+'dm_'+file_name)

        assert len(self.joint_duration) == len(self.IPI)
        assert len(self.joint_duration) == len(self.dm)
        
    def is_consecutive_trial(self):
        '''
        Output: binary list. One represents pair of consecutive pulses, and 0 pair of non-consecutive pulses
        on a certain trial (or time series)
        '''
        #NOTA: esto es para consecutive cumulative y non cumulative! 
        
        joint_duration,dm = self.joint_duration,self.dm   
        
        assert len(joint_duration) == len(dm)
        assert self.dm == [IPI - jd for IPI,jd in zip(self.IPI,self.joint_duration)]
        
        consecutive = []
        for i in range(len(joint_duration)):
            if joint_duration[i]*0.5 >= dm[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        

    def is_isolated_trial(self):
        '''computes the amount of isolated pulses
        substracts the amount of non-isolated pulses to the total pulses
        '''
        consecutive = self.is_consecutive_trial()
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(len(self.max) - count)
            

    def raise_order_consecutiveness(self,consecutive_bef):
        '''computes a higher level of consecutiveness from a binary vector of pulse trains
        (initialize with is_consecutive_trial() )
        
        Returns an empty list when there are no more sequences of pulses of a certain length
        '''
        consecutive = []
        for i,j in zip(consecutive_bef[:-1],consecutive_bef[1:]):
            consecutive.append(i*j)
        assert len(consecutive_bef) -1 == len(consecutive)
        return(consecutive)

            

    def get_consecutive_trains_of_pulses(self):
        '''
        1st output: amount of isolated (non-consecutive) pulses
        2nd output:Given one trial, returns a list with the number of pulses of length n, for the n place of the list
        '''
       
        isolated_pulses = self.is_isolated_trial() # pulsos aislados
        consecutive_bef = self.is_consecutive_trial() #pares de pulsos consecutivos
        
        box_consecutive_trial = [len(self.max),np.sum(consecutive_bef)]
        
        while box_consecutive_trial[-1] > 0:
            consecutive = self.raise_order_consecutiveness(consecutive_bef)            
            box_consecutive_trial.append(sum(consecutive))
            consecutive_bef = consecutive
        return(isolated_pulses,box_consecutive_trial[:-1])# el último elemento es el que se va a cero, por eso lo sacamos

#%%

class consecutive_trial_st_exp:
    '''
    esto es para los boxplots. Quiero la cantidad de pulsos consecutivos, aislados y totales de un dado traza
    '''
    def __init__(self,joint_duration,max_,IPI,dm):

        self.joint_duration = joint_duration
        self.max = max_
        self.IPI = IPI
        self.dm = dm

        assert len(self.joint_duration) == len(self.IPI)
        assert len(self.joint_duration) == len(self.dm)
        
    def is_consecutive_trial(self):
        '''
        Output: binary list. One represents pair of consecutive pulses, and 0 pair of non-consecutive pulses
        on a certain trial (or time series)
        '''
        #NOTA: esto es para consecutive cumulative y non cumulative! 
        
        joint_duration,dm = self.joint_duration,self.dm   
        
        assert len(joint_duration) == len(dm)
        assert self.dm == [IPI - jd for IPI,jd in zip(self.IPI,self.joint_duration)]
        
        consecutive = []
        for i in range(len(joint_duration)):
            if joint_duration[i]*0.5 >= dm[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        

    def is_isolated_trial(self):
        '''computes the amount of isolated pulses
        substracts the amount of non-isolated pulses to the total pulses
        '''
        consecutive = self.is_consecutive_trial()
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(len(self.max) - count)
            

    def raise_order_consecutiveness(self,consecutive_bef):
        '''computes a higher level of consecutiveness from a binary vector of pulse trains
        (initialize with is_consecutive_trial() )
        
        Returns an empty list when there are no more sequences of pulses of a certain length
        '''
        consecutive = []
        for i,j in zip(consecutive_bef[:-1],consecutive_bef[1:]):
            consecutive.append(i*j)
        assert len(consecutive_bef) -1 == len(consecutive)
        return(consecutive)

            

    def get_consecutive_trains_of_pulses(self):
        '''
        1st output: amount of isolated (non-consecutive) pulses
        2nd output:Given one trial, returns a list with the number of pulses of length n, for the n place of the list
        '''
       
        isolated_pulses = self.is_isolated_trial() # pulsos aislados
        consecutive_bef = self.is_consecutive_trial() #pares de pulsos consecutivos
        
        box_consecutive_trial = [len(self.max),np.sum(consecutive_bef)]
        
        while box_consecutive_trial[-1] > 0:
            consecutive = self.raise_order_consecutiveness(consecutive_bef)            
            box_consecutive_trial.append(sum(consecutive))
            consecutive_bef = consecutive
        return(isolated_pulses,box_consecutive_trial[:-1])# el último elemento es el que se va a cero, por eso lo sacamos

    

