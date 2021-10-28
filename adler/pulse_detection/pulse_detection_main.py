'''
Module for pulse detection in Adler noisy time series 
'''

import numpy as np
import os
import multiprocessing as mp
from math import ceil
from functools import partial 
import pandas as pd
import time

from adler.plotting.plotting_main import check_file, save_data, download_data


#%%

def search_extremes(X,TH,W):
    '''
    search_extremes(X,TH,W)
    Search for maxima and minima on a time series X. 
    
    
    Searchs for the local maxima and minima on the time series x that are grater (maxima) than a treshold or 
    smaller (minima) than - threshold. 
    The extremes have at least W separation. 
    
    
    Parameters
    ----------
    X : list
        amplitude values of the time series
    TH : float
        threshold value
    W : minimum distance between maxima / minima
    
    
    Returns
    -------
    MAX : list 
        position of the maxima in X
    MIN : list
        position of the minima  in X
    '''
    


#NO INCLUYO los bordes PORQUE NO TENGO ESTADISTICA DE PICOS AHI
    
### Initializing variables
    MAX = [0] ; MIN = [0]
    M = 0; m = 0; 

    
    for j, x_j in zip(np.arange(1,len(X)+1),X[1:-1]):
        #el zip corta en el array mas pequeno
        #Maxima 
        if x_j > TH : 
            #si x_j supera el treshold
            for i in range(1,W):
                
                x_p = X[j+1]
                x_m = X[j-1]
            
            if (x_j >= x_p) and (x_j >= x_m ) : #Candidato a máximo
                if M == 0 or (M >=1 and ((j-MAX[M]) > W)) : 
                    #Si es el primer máximo o tiene una distancia más grande que W con
                    #el maximo anterior, anotalo
                    M = M + 1 ; MAX.append(j) 
                else: 
                    
                    if (x_j < X[MAX[M]]) : 
                        #reject new maximum
                        pass 
                    else: 
                        # overwrite the old maxima by leaving M unchanged
                        MAX[M] = j 
        #Minima
        if x_j < (-TH) :
            
            x_p = X[j+1]
            x_m = X[j-1]
        
            if (x_j <= x_p) and (x_j <= x_m) : #Candidato a mínimo
            
                #Now test for adjacent minimum, and delete the higher one.   
                if m == 0 or (m >=1 and ((j-MIN[m]) > W)):
                    
                    m = m + 1; MIN.append(j) 
    
                else: #old minimum exists
                    if x_j > X[MIN[m]] : #reject new minimum
                        pass
                    else:  # Overwrite old one by leaving m unchanged
                        MIN[m] = j
    MAX.pop(0);MIN.pop(0)
    return(MAX,MIN)
    #%%
def filter_extremes_aux(MAX,MIN,X):
    ''' filter_extremes_aux(MAX,MIN,X)
    Ensures that in between two maxima there is only one minimum, and in between 
    two minima there is a maximum. 
    
    
    This auxiliary function ensures that in between two maxima in MAX there is one 
    and only one minimum in MIN, and in between two minima in MIN there is one and 
    only one maximum in MAX.
    
    Esto está pensado para buscar los máximos del coseno. Sin embargo, 
    los usamos como referencia para encontrar los PFI y PFE de cada pulso. 
    Parameters
    -------
    X : list
        amplitude values of the time series
     MAX : list 
        position of the maxima in X
    MIN : list
        position of the minima in X
    
    Returns
    -------
    MAX : list 
        position of the maxima
    MIN : list
        position of the minima 
    
    '''
    flag_m1 = False
    # ensures that the time series quantification starts from the maxima. 
    while MAX[0] > MIN[0]:
        MIN.remove(MIN[0])
    while MAX[-1] > MIN[-1]:
        MAX.remove(MAX[-1])
 
    
    list_remove = []; flag_m2 = False; M1_ant = MAX[0]
    for M1,M2 in zip(MAX[:-1],MAX[1:]):
        
        if flag_m2: 
            #The previous M2 was the candidate to remove 
            M1 = M1_ant
            flag_m2 = False
        
        aux_list = list(filter(lambda x: (x > M1 and x < M2), MIN))
        
        while len(aux_list) > 1:
            # Removes all the minima till there is only one in aux_list
            # The minima that removes are the bigger ones. 
            aux_list_values = [X[i] for i in aux_list]
            MIN.remove(aux_list[np.argmax(aux_list_values)])
            aux_list = list(filter(lambda x: (x > M1 and x < M2), MIN))

        if len(aux_list) == 0:
            # If there is not a maxima between two minimas, remove one maximum
            # The maxima that removes are the smaller ones. 
            aux_list_values = [X[i]  for i in [M1,M2]] 
            remove_candidate = [M1,M2][np.argmin(aux_list_values)]
            list_remove.append(remove_candidate)
            
            if remove_candidate == M2:
               flag_m2 = True
               M1_ant = M1
                              
            if remove_candidate == M1:
                flag_m1 = True

        M1_ant = M1
        if flag_m1 : 
            # Tiene que volver a empezar porque sino le quedan dos minimos consecutivos
            break 
                
                
    if len(list_remove) > 0: 
        #print(list_remove)
        for NM in set(list_remove):
            MAX.remove(NM)

    return(MAX,MIN,flag_m1)
    
def filter_last_minima(MAX,MIN):
    tail = [m for m in MIN if m > MAX[-1]]
    while len(tail) > 1:
        MIN.remove(MIN[-1])
        tail = [m for m in MIN if m > MAX[-1]]
    return(MIN)

def filter_extremes(MAX,MIN,X):
    ''' Se asegura que entre dos máx haya dos mín y viceversa'''
    MAX,MIN,flag_m1 = filter_extremes_aux(MAX,MIN,X)

    while flag_m1:
        #print('flag_m1')
        MAX,MIN,flag_m1 = filter_extremes_aux(MAX,MIN,X)
    MIN = filter_last_minima(MAX,MIN)
    
    assert len(MAX) == len(MIN)
    return(MAX,MIN)
    


#%%#
def get_fixed_points(alpha):
    '''
    get_fixed_points(alpha)
    returns the fixed points for an adler dynamical system. 
    If the system is oscillatory and not excitable, it returns the point called ghost of less veoloty. 
    
    Parameters
    ----------
    alpha : float
        parameter of the adler equation. In omega units
    
    Returns
    -------
    PFE : float
        stable fixed point - angular . On the  3rd cuadrant.
    PFI : float
        unstable fixed point - angular . On the 4th cuadrant
    '''
    if alpha >= 1:
        res = np.arcsin(-1/alpha)
        PFE = -res + np.pi  #np.sin(PFE)
        PFI = (2*np.pi + res)  #np.sin(PFI)
    else: 
        PFE = -np.pi/2 + 2*np.pi
        PFI = -np.pi/2 + 2*np.pi
    return(PFE,PFI)


#%%

def pop_list(list_,remove_list_index):
    new_list = []
    for (index,value) in enumerate(list_):
        if index in remove_list_index:
            pass
        else:
            new_list.append(value)
    assert (len(list_) == len(remove_list_index) + len(new_list))
    return(new_list)
    
def get_left_minima(MAX,MIN,PFE,PFI,theta):
    '''
    get_left_minima(MAX,PFE,PFI,theta)
   from a given maxima, gets the beginning of a pulse
    
    gets the beginning of a pulse as the following: moving backwards in time from a given maxima,
    searchs from the first time point in wich the dynamical systems crosses the unstable fixed point
    

    
    Parameters
    ----------
    MAX : list 
        position of the maxima of X time series
    PFE : float
        stable fixed point - angular 
    PFI : float
        unstable fixed point - angular 
    theta : list
        amplitude values of the time series given in angles
    
    
    Returns
    -------
    left_minima : list 
        position of the starting point of each pulses 
    '''
    left_minima = []; remove_MAX = []; 
    if len(MAX) > 0:
       
        n0 = theta[MAX[0]] // (2*np.pi); #M_ant = None
        for (n,(M_id,M)) in zip(np.arange(n0,len(MAX)+n0),enumerate(MAX)):
            if M_id == 0:
                m = None 
            else:
                m = MIN[M_id-1]#este es el minimo que precede al maximo
            # On the following, ensures to use the correct quadrant
            while PFI > (theta[M] - n*2*np.pi) : # cota inferior de una función creciente de -n
                n = n - 1
            while  (theta[M] - n*2*np.pi) > 5*np.pi/2:
                n = n + 1
                
            assert (PFI < (theta[M] - n*2*np.pi) < 5*np.pi/2)                
            theta_ = theta - n*2*np.pi
            assert (len(theta_[m:M]) > 0)
           
            # searches for the left minima of the maximum M
            L = len(left_minima)
            for i, ang in enumerate(theta_[m:M][::-1]): #empieza uno después del maximo
                if ( ang - PFI <= 0 ):
                    left_minima.append(M-i); 
                    break
                else: 
                    pass
            
            if len(left_minima) - L == 1:
                pass
            else:
                remove_MAX.append(M_id)   
            #M_ant = M
            
        MAX = pop_list(MAX,remove_MAX)
        MIN = pop_list(MIN,remove_MAX)
        if (len(MAX) > 0) : assert (len(left_minima) > 0) 
    else:
        print('The trace has no maxima values')
    
    
    assert len(left_minima) == len(MAX), ' left_minima:' + str(left_minima) +' MAX:'+ str(MAX)
    
    return(left_minima, MAX, MIN)
    
def get_rigth_minima(left_minima,MAX,MIN,PFE,PFI,theta):
    '''
    get_rigth_minima(MIN,PFE,PFI,theta)
    gets the ending point of a pulse
    
    gets the ending point of a pulse as the following: moving forward in time from a given minima,
    searchs from the first time point in wich the dynamical systems crosses the stable fixed point
    
    Parameters
    ----------
    MIN : list 
        position of the minima of X time series
    PFE : float
        stable fixed point - angular 
    PFI : float
        unstable fixed point - angular 
    theta : list
        amplitude values of the time series given in angles
    
    
    Returns
    -------
    right_minima : list 
        position of the ending point of each pulses 
    '''
    right_minima = [];remove_MIN = []
    if len(MIN) > 0:
        n0 = theta[MIN[0]] // (2*np.pi)
        for (n,(m_id,m)) in zip(np.arange(n0,len(MIN)+n0),enumerate(MIN)):
            
            if m == MIN[-1]:
                M = None
            else:
                M = MAX[m_id+1] #este es el maximo que viene después del mínimo            
            while np.pi/2 > (theta[m] - n*2*np.pi) :
                n = n - 1
            while (theta[m] - n*2*np.pi) > PFE :
                n = n + 1
            
            assert (np.pi/2 < (theta[m] - n*2*np.pi) < PFE)
            theta_ = theta - n*2*np.pi
            assert (len(theta_[m:M]) > 0)
            
            L = len(right_minima)
            for i, ang in enumerate(theta_[m:M]): #empieza uno después del maximo
                #print(ang-PFE)
                if ( ang - PFE >= 0 ):
                    right_minima.append(m+i); 
                    break
                else: 
                    pass
            if len(right_minima) - L == 1:
                pass
            else: remove_MIN.append(m_id)
            
        MAX = pop_list(MAX,remove_MIN)
        MIN = pop_list(MIN,remove_MIN)
        left_minima = pop_list(left_minima,remove_MIN)
        if (len(MIN) > 0) : assert (len(right_minima) > 0)
        
    else: 
        print('The trace has no maxima values')
    assert (len(right_minima) == len(MIN)),'right_minima: '+ str(right_minima) + ' MIN:' + str(MIN)
    
    return(left_minima,right_minima,MAX,MIN )
#%%
    
    
    
#%%
def filter_maxima_sine(left_minima,right_minima,MAX):
    ''' remueve los bordes de los maximos de seno para que queden entre los minimos de pulsos, y descarta el resto'''
    
    
    remove_list_index = []
    
    for i,M in enumerate(MAX):
        if M < left_minima[0]:
            remove_list_index.append(i)
        if M > right_minima[-1]:
            remove_list_index.append(i)
    
    return(pop_list(MAX,remove_list_index))


#%%
    
def get_pulses(theta,TH,W,PFE,PFI,alpha,D):
    # Ddetecta pulsos : primero agarra los maximos y minimos, luego los filtra. Luego, calcula los comienzpy y finales de los pulsos.
    # respues testea cosas: todos los pulsos tienen unpo y solo uno proincipio, may, min y final (PONER EXACTAMENTE BIEN QUE TESTEA)
    t0 = time.time() ; print('start searching for pulses: 5 steps pending')
    
    MAX_cos,MIN_cos = search_extremes(np.cos(theta),TH,W)
    t1 = time.time() - t0
    print('step 1/5 finished: cos extremes detected',t1, 'sec')
    
    if (len(MAX_cos)>0) and (len(MIN_cos)>0): 
        MAX_cos, MIN_cos =  filter_extremes(MAX_cos,MIN_cos,np.cos(theta))
        t2 = time.time() - t1
        print('step 2/5 finished: 1 minimum between 2 cos maxima',t2, 'sec')
    
    if (len(MAX_cos)>0) and (len(MIN_cos)>0): 
        left_minima, MAX_cos ,MIN_cos = get_left_minima(MAX_cos,MIN_cos,PFE,PFI,theta)
        t3 = time.time() - t2
        print('step 3/5 finished: left pulse cos minima detected',t3, 'sec')
        
        left_minima, right_minima, MAX_cos, MIN_cos = get_rigth_minima(left_minima,MAX_cos,MIN_cos,PFE,PFI,theta)
        t4 = time.time() - t3
        print('step 4/5 finished: right pulse cos minima detected',t4, 'sec')
        
        print('testing alpha,d: ',alpha,D)
        test_pulses(left_minima,right_minima,MAX_cos,MIN_cos)

        MAX_sine , MiN_sine = search_extremes(np.sin(theta),TH,W)
        MAX_sine =  filter_maxima_sine(left_minima,right_minima,MAX_sine)
        t5 = time.time() - t4
        print('step 5/5 finished: sine extremes detected',t5, 'sec')
        test_pulses_sine(left_minima,right_minima,MAX_sine)

        
        return left_minima,right_minima,MAX_sine
    else:
        print('The trace has no maxima or minima values available')
        return [],[],[]



def main_pulse_detection(theta,alpha,omega,D,save_path_name,file_name):
    # calcula y guarda los pulsos de theta en el archivoPONER NOMBRES
    if len(theta) != 0:
        print('running pulse detection',alpha,D)      
        TH = 0.90;W = 100
        PFE , PFI = get_fixed_points(alpha)
        left_minima,right_minima,MAX = get_pulses(theta,TH,W,PFE,PFI,alpha,D)
        print('pulse detection ended')
        
        if (len(MAX)>0): #and (len(MIN)>0): 
            print('saving pulse detection results')
            save_data(MAX,save_path_name+'max_xf_'+file_name)
            save_data(left_minima,save_path_name+'left_minima_'+file_name)
            save_data(right_minima,save_path_name+'right_minima_'+file_name)
            print(file_name,'saving finished')
        else:
            print(alpha,D,'no pulses on this condition -> not saving')      
    else:
        print('ERROR: theta with length zero',omega,alpha,D)
    return(0)


    
def main_pulse_detection_(data_folder,save_path_name,tuple_,overwrite_flag = True):
   #ES UNA funcion auxiliar para paralelizar

    (i,D,order),row = tuple_[0],tuple_[1]
    omega =  row.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    if check_file(file_name,data_folder):   
        if check_file('right_minima_'+file_name,data_folder) and (not overwrite_flag):
            print('pulse detection allready done and nor overwitting it',file_name)
            pass
        else:
            print('file name:',file_name)
            theta = download_theta(file_name,data_folder)
            main_pulse_detection(theta,alpha,omega,D,save_path_name,file_name)
    else:
        print('ERROR: file not available',file_name)
    return(1)


def compute_pulse_detection(description_file,data_folder,save_path_name):
    #esta es la funcion que le calcula a cada seite temporal suspulsos
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha','D','order'])
    main_pulse_detection__ = partial(main_pulse_detection_,data_folder,save_path_name)
    pool.map(main_pulse_detection__,tuple_)
    pool.close()
    pool.join()
    return (2)


#%%
def download_theta(file_name,data_folder):
    #esta funcion es para leer la serie temporal
    list_data_folder = os.listdir(data_folder)
    if file_name in list_data_folder:
        print(file_name,' available')
        if check_file(file_name,data_folder):           
            return(download_data(data_folder + file_name))
        else:
            print(file_name,'filename empty')
            return []
    else:
        print(file_name,'filename not available')
        return []

#%%
def test_pulses(left_minima,right_minima,MAX,MIN):
    ''' testea que todos los vectores sean de la misma longitud, y que siempre venga
    primero un left minimum, después un max, despues un min, despues un rigth min'''
    
    assert len(MIN) == len(MAX)
    assert len(left_minima) == len(MAX), str(len(left_minima))+' is not equal to ' + str(len(MAX))
    assert len(right_minima) == len(MAX) , str(len(right_minima))+' is not equal to ' + str(len(MAX))
    
    #los proximos 3 iguales estan mal
    assert all([(i<=j)*1 for (i,j) in zip(left_minima,MAX)]), '----- ' + str([(j-i)*1 for (i,j) in zip(left_minima,MAX)])
    assert all([(i<=j)*1 for (i,j) in zip(MAX,MIN)])
    assert all([(i<=j)*1 for (i,j) in zip(MIN,right_minima)])
    
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]),str([(i-j) for (i,j) in zip(left_minima[1:],right_minima[:-1])]) 

def test_pulses_sine(left_minima,right_minima,MAX):
    ''' testea que todos los vectores sean de la misma longitud, y que siempre venga
    primero un left minimum, después un max, despues un min, despues un rigth min'''
    
    #assert len(MIN) == len(MAX)
    assert len(left_minima) == len(MAX), str(len(left_minima))+' (sine checking) is not equal to ' + str(len(MAX))
    assert len(right_minima) == len(MAX) , str(len(right_minima))+' (sine checking) is not equal to ' + str(len(MAX))
    
    #los proximos 3 iguales estan mal
    assert all([(i<=j)*1 for (i,j) in zip(left_minima,MAX)]), '--(sine checking)--- ' + str([(j-i)*1 for (i,j) in zip(left_minima,MAX)])
    assert all([(i<=j)*1 for (i,j) in zip(MAX,right_minima)])
    
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]),str([(i-j) for (i,j) in zip(left_minima[1:],right_minima[:-1])])

#%% 
# =============================================================================
# =============================================================================
# =============================================================================
#  ESTE MODULO ES PARA COMPUTAR LOS CUANTIFICADORES DE PULSOS
# =============================================================================
# =============================================================================
# =============================================================================


def test_pulses_quantifiers(dt,IPI,dm,joint_duration,MAX):
    
    assert len(dm) == len(IPI)
    assert len(joint_duration) == len(IPI)
    assert len(dt) == len(MAX)
    
    assert all([(i + j == k)*1 for (i,j,k) in zip(dm,joint_duration,IPI)]), 'dm:'+str(dm) + ' joint_duration: ' + str(joint_duration) + ' IPI: '+str(IPI)


def get_pulses_quantifiers(left_minima,right_minima,MAX,MIN):
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
        
        print('running pulses quantifiers computation',alpha,D)      
        
        MAX = download_data(data_folder + 'max_xf_'+file_name) 
        MIN = download_data(data_folder + 'min_xf_'+ file_name) 
        left_minima = download_data(data_folder + 'left_minima_'+ file_name) 
        right_minima = download_data(data_folder + 'right_minima_'+ file_name) 
        test_pulses(left_minima,right_minima,MAX,MIN)
        
        dt,IPI,dm,joint_duration = get_pulses_quantifiers(left_minima,right_minima,MAX,MIN)
        
        print('Ready! Saving files --- ',alpha,D)      
        save_data(dt,save_path_name+'dt_xf_'+file_name)
        save_data(IPI,save_path_name+'IPI_xf_'+file_name)
        save_data(dm,save_path_name+'dm_xf_'+file_name)
        save_data(joint_duration,save_path_name+'joint_duration_xf_'+file_name)
        print(alpha,D,'pulses quantifiers computation finished')
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



