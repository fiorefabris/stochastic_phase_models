import numpy as np
from plotting.plotting_main import download_data
from itertools import groupby

def count_iterable(iterator):
    return sum(1 for i in iterator)
#%%

class consecutive_cumulative:
    
    def __init__(self,number,order,data_folder):
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        self.joint_duration = download_data(data_folder+'joint_duration_xf_'+file_name)
        self.min = download_data(data_folder+'min_xf_'+file_name)
        self.max = download_data(data_folder+'max_xf_'+file_name)
        self.IPI = download_data(data_folder+'IPI_xf_'+file_name)
        self.dm = download_data(data_folder+'dm_xf_'+file_name)

        assert len(self.joint_duration) == len(self.IPI)
        assert len(self.joint_duration) == len(self.dm)
        assert len(self.max) == len(self.min)
        
    def is_consecutive_trial(self):
        #te devuelve par de pulsos consecutivos para cada serie temporal o trial
        y = [IPI - jd for IPI,jd in zip(self.IPI,self.joint_duration)] #silence 
        assert self.dm == y
        x = self.joint_duration
       
        consecutive = []
        for i in range(len(x)):
            if x[i]*0.5 >= y[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        
   
    def raise_order_consecutiveness(self,consecutive_bef):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad para un dado binario vector de trenes de pulsos
        consecutive = []
        for i,j in zip(consecutive_bef[:-1],consecutive_bef[1:]):
            consecutive.append(i*j)
        return(consecutive)
    
    def get_consecutive_trains_of_pulses(self):
        # devuelve una lista con la suma de la cantidad de trenes de pulsos consecutivos de longitud n dado un trial o serie temporal
        #cada lugarcito es la suma en una serie temporal o trial
        consecutive_bef = self.is_consecutive_trial() #pares de pulsos consecutivos
        box_consecutive_trial = [len(self.max),np.sum(consecutive_bef)]
        
        while box_consecutive_trial[-1] > 0:
            consecutive = self.raise_order_consecutiveness(consecutive_bef)            
            box_consecutive_trial.append(sum(consecutive))
            consecutive_bef = consecutive
        return(box_consecutive_trial[:-1])# el último elemento es el que se va a cero, por eso lo sacamos
    


#%%
class consecutive_non_cumulative:

    def __init__(self,number,order,data_folder):
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        self.joint_duration = download_data(data_folder+'joint_duration_xf_'+file_name)
        self.min = download_data(data_folder+'min_xf_'+file_name)
        self.max = download_data(data_folder+'max_xf_'+file_name)
        self.IPI = download_data(data_folder+'IPI_xf_'+file_name)
        self.dm = download_data(data_folder+'dm_xf_'+file_name)

        assert len(self.joint_duration) == len(self.IPI)
        assert len(self.joint_duration) == len(self.dm)
        assert len(self.max) == len(self.min)
        
    def is_consecutive_trial(self,individual_pulses = False):
        # te devuelve una lista con 1 representando pares de pulsos consecutivos y
        # 0 pares de pulsos no consecutivos para un dado trial o serie temporal
        joint_duration,dm = self.joint_duration,self.dm   
        assert len(joint_duration) == len(dm)
        
        consecutive = []
        for i in range(len(joint_duration)):
            if joint_duration[i]*0.5 >= dm[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)


    def is_isolated_trial(self):
        # le resta el total de pulsos consecutivos al total de pulsos
        consecutive = self.is_consecutive_trial()
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(len(self.max) - count)
            
    def raise_order_consecutiveness(self,consecutive_):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad, para dado una lista binaria de trenes de pulsos consecutivos o no cons
        new_consecutive_ = []
        for i,j in zip(consecutive_[:-1],consecutive_[1:]):
            new_consecutive_.append(i*j)
        assert len(consecutive_) -1 == len(new_consecutive_)
        return(new_consecutive_)
    
    
    def count_consecutive_pulses(self,consecutive):
        #cuenta la cantidad de unos aislados que hay en un vector 
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                if count_iterable(group) == 1:
                    count = count + 1
                else:
                    pass
        return(count) 
        
    
    def get_consecutive_trains_of_pulses(self):
        "cuenta cuantos pulsos de longitud n hay en trial o serie temporal"
        isolated_pulses = self.is_isolated_trial()
        consecutive = self.is_consecutive_trial()       
        
        consecutive_ = [isolated_pulses]

        while sum(consecutive) > 0:
            consecutive_.append(self.count_consecutive_pulses(consecutive))
            consecutive = self.raise_order_consecutiveness(consecutive)
        return(consecutive_) 
            
#%%
        
def get_consecutive_trains_of_pulses_old(self):
    #cada lugarcito es la suma en una celula
    consecutive_bef = self.is_consecutive_trial()
    box_consecutive_trial = [len(self.max),np.sum(consecutive_bef)]
    
    n = 0
    while box_consecutive_trial[-1] > 0:
        consecutive = self.raise_order_consecutiveness(consecutive_bef)
        n = n + 1 
        
        box_consecutive_trial.append(self.count_consecutive_pulses(consecutive,n))
        consecutive_bef = consecutive
    return(box_consecutive_trial[:-1])


        

def count_consecutive_pulses_old(self,consecutive,n):
   # n empieza en cero para pulsos aislados
    #La idea acá es que en cada tren de pulsos finito, hay que sumarle n = #de pulsos del tren - 1
    # para dar con el total de pulsos
    
    count = 0
    for i, group in groupby(consecutive):
        if i == 1:
            count = count + count_iterable(group) + n
    return(count)  
        
    
    
def count_consecutive_pulses(self,consecutive_):
    #cuenta la cantidad de unos aislados que hay en un vector 
    count = 0
    for j,n in enumerate(consecutive_):
        if n == 1: #si tenemos un conjunto de pulsos consecutivos
            
            # para el primer conjunto de pulsos
            #si hay mas conjuntos de pulsos y el siguiente no es consecutivo, contalo
            if j == 0: 
                if len(consecutive_) > 1: 
                    if consecutive_[j+1] == 0: 
                        count = count + 1 
            # para el primer conjunto de pulsos
            #si hay un solo conjunto de pulsos, contalo
                else:
                    count = count + 1 
            # si el conjunto de pulsos es el ultimo y el anteultimo no es consecutivo, contalo
            elif j == len(consecutive_)-1 :
                if consecutive_[j-1] == 0:
                    count = count + 1
            else:
            #si no es el primer ni el ultimo conjunto de pulsos, y los de alrededor no tienen pusos, contalo
                if (consecutive_[j-1] + consecutive_[j+1]) == 0:
                    count = count + 1
        #else: 
         #   assert n ==0
    return(count) 

        

