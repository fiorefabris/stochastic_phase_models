import numpy as np
from adler.data_managing_functions import download_data
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
#%%
    
    


#%%
# =============================================================================
# eSTE MODULO ES PARA COMPUTAR CONSECUTIVIDAD
# =============================================================================
''' module for computing consecutive trains of pulses 
    como concepto general, la box tiene en cada lugar una celula
'''


def is_isolated_cell(MAX,joint_duration,dm):
    # solo para una celula
    # calcula el numero de picos no consecutivis (isolated)
    # tiene en cuenta las celulas de un solo pulso
    # te devuelve la suma del total de picos menos los que son parte de intervalos consecutivos de pulsos

    return(np.sum(MAX) - sum(is_consecutive_cell(joint_duration,dm,True)))
    
def is_isolated_box(tuple_,data_folder):
    #Para cada conjunto de celulas, te da la estadistica de pulsos aislados
    box = []
    for row in tuple_[1]:
        file_name =  str(int(row.number))+'_'+str(int(row.order))+'.pkl'
        if (check_file('max_xf_'+file_name,data_folder)):
            box.append(is_isolated_cell(download_data(data_folder + 'max_xf_'+file_name),download_data(data_folder + 'joint_duration_xf_'+file_name),download_data(data_folder + 'dm_xf_'+file_name)))
    return(box)


def is_consecutive_cell(joint_duration,dm,number_of_pulses = False):
    # para cada celula (data), te devuelve un array con 1 representando pares de pulsos consecutivos y
    # 0 pares de pulsos no consecutivos
    #con number of pulses true te devuelve la cantidad de pulsos consecutivos . sino, te devuelve pares de pulsos
    consecutive = []
    for i in range(len(joint_duration)):
        if joint_duration[i]*0.5 >= dm[i]:
            if (number_of_pulses and (len(consecutive) ==0 or consecutive[-1] ==0)): 
                consecutive.append(1)
            consecutive.append(1)
        else:
            consecutive.append(0) 
    return(consecutive)
    

    
def is_consecutive_box(tuple_,data_folder):
    #Para cada conjunto de celulas, te da la estadistica de pulsos consecutivos en una lista
    box = []
    for row in tuple_[1]:
        file_name =  str(int(row.number))+'_'+str(int(row.order))+'.pkl'
        if (check_file('max_xf_'+file_name,data_folder)):
            box.append(is_consecutive_cell(download_data(data_folder + 'joint_duration_xf_'+file_name),download_data(data_folder + 'dm_xf_'+file_name)))
    return(box)

def raise_order_consecutiveness(box):
    #Te da lista vacia cuando no hay mas pares de pulsos. 
    # Calcula un nivel más de consecutividad
    new_box = []
    for consecutive_ in box:
        aux_consecutive_ = []
        for i,j in zip(consecutive_[:-1],consecutive_[1:]):
            aux_consecutive_.append(i*j)
        new_box.append(aux_consecutive_)
    return(new_box)

def count_consecutive_pulses(box, population = False):
    #cuenta la cantidad de unos aislados que hay en un vector 
    #si population es True , te devuelve la suma directamente sobretoda la poblacion
    count_box = []
    for consecutive_ in box:
        count = 0
        for j,n in enumerate(consecutive_):
            if n == 1:
                if j == 0:
                    if len(consecutive_) > 1:
                        if consecutive_[j+1] == 0:
                            count = count + 1
                    else:
                        count = count + 1 #este else no estoy segura! que pasa cuando hay solo un uno en el array?
                elif j == len(consecutive_)-1 :
                    if consecutive_[j-1] == 0:
                        count = count + 1
                else:
                    if (consecutive_[j-1] + consecutive_[j+1]) == 0:
                        count = count + 1
            else: 
                pass
        count_box.append(count)
    if population:
        return(np.sum(count_box)) #te tira el total sobre todo el dataset
    else:
        return(count_box) # te lo tira por celula


#%%
    
def get_consecutive_trains_of_pulses(tuple_,data_folder):
    #te da la suma sobre todo el dataset
    
    isolated = is_isolated_box(tuple_,data_folder)
    box = is_consecutive_box(tuple_,data_folder)
    
    box_plot_consecutive = [np.sum(isolated)]
    
    while np.sum([np.sum(l) for l in box]) > 0:
        box_plot_consecutive.append(count_consecutive_pulses(box,True))
        box = raise_order_consecutiveness(box)
        
    return(box_plot_consecutive)


def get_consecutive_trains_of_pulses_cells(tuple_,data_folder):
    #cada lugarcito es la suma en una celula
    isolated = is_isolated_box(tuple_,data_folder)
    box = is_consecutive_box(tuple_,data_folder)
    box_plot_consecutive = [isolated]
    
    while np.sum([np.sum(l) for l in box]) > 0:
        box_plot_consecutive.append(count_consecutive_pulses(box,False))
        box = raise_order_consecutiveness(box)
    return(box_plot_consecutive)
    
    
#%%
################################################
#### FFT FDT  Module for computing the consecutiveness
################################################  
#from adler.pulse_detection.consecutive_main import consecutive_cumulative,consecutive_non_cumulative

def get_consecutive(description_file,data_folder,save_folder):
    '''
    data folder: donde están los dt
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha','D'])
    get_consecutive_trial_D_ = partial(get_consecutive_trial,data_folder,save_folder)
    pool.map(get_consecutive_trial_D_,tuple_)
    pool.close()
    pool.join()
    return (2)

def get_consecutive_trial(data_folder,save_folder,tuple_):
    
    (i,D),row_ = tuple_[0],tuple_[1]
    omega =  row_.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    for (order,row) in row_.groupby(['order']):
        number      = int(row.number)
        file_name   =  str(number)+'_'+str(order)+'.pkl'
            
        if (check_file('min_xf_'+file_name,data_folder)):       
            consecutive_cumulative_ = consecutive_cumulative(number,order,data_folder)
            box_consecutive_cumulative_trial = consecutive_cumulative_.get_consecutive_trains_of_pulses()#
            
            consecutive_non_cumulative_ = consecutive_non_cumulative(number,order,data_folder)
            box_consecutive_non_cumulative_trial = consecutive_non_cumulative_.get_consecutive_trains_of_pulses()
            
            assert  len(box_consecutive_cumulative_trial) == len(box_consecutive_non_cumulative_trial)
            assert np.sum([(i+1)*j for i,j in enumerate(box_consecutive_non_cumulative_trial)]) == box_consecutive_cumulative_trial[0]
            
            save_data(box_consecutive_non_cumulative_trial,save_folder+'c_nc_'+file_name)
            save_data(box_consecutive_cumulative_trial,save_folder+'c_c_'+file_name)
        else:
            pass
    return(0)


