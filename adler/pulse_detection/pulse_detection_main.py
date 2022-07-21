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

from adler.data_managing_functions import check_file, save_data, download_data, get_fixed_points
from adler.pulse_detection.quantifiers_main import get_pulses_quantifiers
from adler.pulse_detection.consecutive_main import consecutive_trial_st_exp

def pop_list(list_,remove_list_index):
    new_list = []
    for (index,value) in enumerate(list_):
        if index in remove_list_index:
            pass
        else:
            new_list.append(value)
    assert (len(list_) == len(remove_list_index) + len(new_list))
    return(new_list)

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


def test_points_max(x_j,j,X,W):
    #evalúa si x_j, contenido en x en el lugar j, 
    #es mayor o igual a los puntos dentro del entorno W
    
    test = True
    
    for i in range(1,W+1):
        #empieza en 1 y se corta en W
        x_p = X[j+i]
        x_m = X[j-i]
        if (x_j >= x_p) and (x_j >= x_m):
            pass
        else:
            test = False
            break
    return test


def test_points_min(x_j,j,X,W):
    #evalúa si x_j, contenido en x en el lugar j, 
    #es menor o igual a los puntos dentro del entorno W
    
    test = True
    
    for i in range(1,W+1):
        #empieza en 1 y se corta en W
        x_p = X[j+i]
        x_m = X[j-i]
        if (x_j <= x_p) and (x_j <= x_m):
            pass
        else:
            test = False
            break
    return test
#%%
# =============================================================================
# =============================================================================
# =============================================================================
#                   PULSE DETECTION MODULE
# =============================================================================
# =============================================================================
# =============================================================================


def search_extremes(X,TH,W):
    '''
    search_extremes(X,TH,W)
    Search for maxima and minima on a time series X. 
    
    
    Searchs for the local maxima and minima on the time series x that are grater (maxima) than a treshold or 
    smaller (minima) than these same threshold, multiplied to -1. 
    The extremes have at least W separation. 
    
    The borders are not included, because we dont have enaugh statistics there.
    #NO INCLUYO los bordes PORQUE NO TENGO ESTADISTICA DE PICOS AHI
    
    
    Parameters
    ----------
    X : list
        amplitude values of the time series
    TH : float
        threshold value
    W : float
        minimum distance between maxima / minima, given in indexes.
    
    
    Returns
    -------
    MAX : list 
        position of the maxima in X
    MIN : list
        position of the minima  in X
    '''


### Initializing variables
    MAX = [0] ; MIN = [0]
    M = 0; m = 0; 

    
    for j, x_j in zip(np.arange(W,len(X)-W),X[W:len(X)-W]):
        #quiero que empiece en W, así puedo compararlo con hasta el índice cero
        #quiero que termine en len(X)-W-1, asi puedo compararlo hasta el índice -1
        
        #Maxima 
        if x_j > TH : 
            #si x_j supera el treshold
            
            if test_points_max(x_j,j,X,W):
                #si es el maximo en todo su entorno W => candidato a maximo
                
                if M == 0 or (M >=1 and ((j-MAX[M]) > W)): 
                    #Si es el primer máximo o tiene una distancia más grande que W con
                    #el maximo anterior, anotalo
                    
                    M = M + 1 ; MAX.append(j) 
                
                else: 
                    # overwrite the old maxima by leaving M unchanged (random)
                    if  np.random.rand() > 0.5 :
                        MAX[M] = j 
                    else:
                        pass
                    
        #Minima
        elif x_j < (-TH) :
            
            if test_points_min(x_j,j,X,W):
            
                #Now test for adjacent minimum, and delete the higher one.   
                if m == 0 or (m >=1 and ((j-MIN[m]) > W)):
                
                    m = m + 1; MIN.append(j) 

                else: 
                    # overwrite the old maxima by leaving M unchanged (random)
                    if  np.random.rand() > 0.5 :
                        MIN[m] = j 
                    else:
                        pass
        else:
            #no es candidato a maximo ni minimo
            pass
        
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
    
    This is disigned to look for the may and min of cos(theta).Nevertheless, we use it as a reference 
    to find the PFE and PFI that defines a pulse. 
    
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
    if len(MAX) * len(MIN) > 0:
        while MAX[0] > MIN[0]:
            MIN.remove(MIN[0])
            if len(MAX) * len(MIN) == 0:
                MAX,MIN = [],[]
                break
    if len(MAX) * len(MIN) > 0:
        while MAX[-1] > MIN[-1]:
            MAX.remove(MAX[-1])
            if len(MAX) * len(MIN) == 0:
                MAX,MIN = [],[]
                break
    
    if len(MAX) * len(MIN) > 0:
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
    else:
        MIN,MAX = [],[]

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
    
    if len(MIN)*len(MAX) > 0:
        MIN = filter_last_minima(MAX,MIN)
    
    assert len(MAX) == len(MIN)
    assert all([M-m > 0 for M,m in zip(MIN,MAX)])
    return(MAX,MIN)
    


#%%

def read_string_max(left_minima,theta_,PFI,m,M):
    
    for i, ang in enumerate(theta_[m:M][::-1]): #empieza uno después del maximo
        test = False                        
        
        if ( ang - PFI <= 0 ):
            left_minima.append(M-i); 
            test = True
            break
        
        else: 
            pass
    return(test,left_minima)

def append_remove_MAX(test,remove_MAX,M_id):
    #test es True cuando append algo al left_minima
    if test:
        pass
    else:
        remove_MAX.append(M_id)  
    return remove_MAX


def read_string_min(right_minima,theta_,PFE,m,M):
    for i, ang in enumerate(theta_[M:m]): #empieza uno después del maximo
        test = False
        if ( ang - PFE >= 0 ):
            right_minima.append(M+i)
            test = True
            break
        else: 
            pass
    return(test,right_minima)


#%%

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
        stable fixed point - angular (4th cuadrant)
    PFI : float
        unstable fixed point - angular (3rd cuadrant)
    theta : list
        amplitude values of the time series given in angles
    
    
    Returns
    -------
    left_minima : list 
        position of the starting point of each pulses 
    '''
    #recorro del maximo al minimo anterior. Para el seno, eso es 
    # del cuadrante 1 al cuadrante -1
    
    creciente = []
    assert all([M < m for M,m in zip(MAX,MIN)])
    PFI = PFI - 2*np.pi ; assert (-np.pi< PFI < 0)
    left_minima = []; remove_MAX = []; 
    
    
    if len(MAX) <= 0:
        print('The trace has no maxima values')
    
    else:
       
        n0 = theta[MAX[0]] // (2*np.pi); 
        
        for (n,(M_id,M)) in zip(np.arange(n0,len(MAX)+n0),enumerate(MAX)):
            
            if M_id == 0:
                m = None 
            else:
                m = MIN[M_id-1] #este es el minimo que precede al maximo
                assert (M-m > 0)
            
            
            #si la fase del mínimo es más chica que la del máximo, decreciente de maximo a minimo
            if (m is None) or (theta[m] <  theta[M]):              
                # On the following, ensures to use the correct quadrant
                # el máximo tiene que estar entre el cuadrante -1 y el 3er cuadrante
                
                while PFI > (theta[M] - n*2*np.pi) : 
                    n = n - 1
                while  (theta[M] - n*2*np.pi) > np.pi:
                    n = n + 1
         
                theta_ = theta - n*2*np.pi
                assert (PFI < theta_[M] < np.pi)                
                assert (len(theta_[m:M]) > 0)
                
                
                # searches for the left minima of the maximum M
                test,left_minima = read_string_max(left_minima,theta_,PFI,m,M)
                
                if (not test) and (PFI == -np.pi/2):
                    # no llego al PFI porque el minimo está ahi cerquita. Es para el caso oscilatorio
                    if m is None:
                        remove_MAX = append_remove_MAX(test,remove_MAX,M_id)
                    else:
                        
                        if M_id <= 0: #si es el primer maximo 
                            m = None
                        else:
                            m = MAX[M_id-1] #maximo anterior
                            
                        test,left_minima = read_string_max(left_minima,theta_,PFI,m,M)
                        remove_MAX = append_remove_MAX(test,remove_MAX,M_id)
                else:
                    remove_MAX = append_remove_MAX(test,remove_MAX,M_id)
            else:
                #si la fase del mínimo es más grande que la del máximo (theta[m] >  theta[M])
                #print('creciente de maximo a minimo')
                creciente.append(MAX[M_id])
                remove_MAX.append(M_id)   


        #print(creciente/len(MAX))
        MAX = pop_list(MAX,remove_MAX)
        MIN = pop_list(MIN,remove_MAX) #remuevo el minimo anterior al maximo, total despues eso es lo que voy a usar
        if (len(MAX) > 0) : assert (len(left_minima) > 0) 
    
    
    assert len(left_minima) == len(MAX), ' left_minima:' + str(left_minima) +' MAX:'+ str(MAX)
    return(creciente,left_minima, MAX, MIN)
#%%    
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
    #recorro entre el maximo y el proximo minimo
    #para el seno, esto es del cuadrante 1 al 3.

    decreciente = []
    assert all([M < m for M,m in zip(MAX,MIN)])
    right_minima = []; remove_MAX = []; 

    if len(MAX) <= 0:
        print('The trace has no maxima values')
    else:
        n0 = theta[MAX[0]] // (2*np.pi);
    
        for (n,(M_id,M)) in zip(np.arange(n0,len(MAX)+n0),enumerate(MAX)):
            
            if M == MAX[-1]:
                m = None 
            else:
                m = MIN[M_id]#este es el minimo que sucede al maximo
                assert (m-M > 0)
            
            if (m is None) or (theta[m] >  theta[M]):
            #si la fase del mínimo es más grande que la del máximo, es creciente de máximo a mínimo

                # esto tiene que estar entre el cuadrante 1 y 3
                while 0 > (theta[M] - n*2*np.pi) :
                    n = n - 1
                while (theta[M] - n*2*np.pi) > PFE :
                    n = n + 1
                
                theta_ = theta - n*2*np.pi
                assert (0 <  theta_[M] < PFE)
                assert (len(theta_[M:m]) > 0)
            
            
                test,right_minima = read_string_min(right_minima,theta_,PFE,m,M)
                
                if (not test) and (PFE == 3/2*np.pi):
                    if m is None:
                        remove_MAX = append_remove_MAX(test,remove_MAX,M_id)
                    else:
                        
                        if M_id >= len(MAX)-1: #es el último máximo
                            m = None
                        else:
                            m = MAX[M_id+1]
                            
                        test,right_minima = read_string_min(right_minima,theta_,PFE,m,M)
                        remove_MAX = append_remove_MAX(test,remove_MAX,M_id)
                else:
                    remove_MAX = append_remove_MAX(test,remove_MAX,M_id)

            
            else:
                #si la fase del mínimo es chica grande que la del máximo 
                decreciente.append(MAX[M_id])
                remove_MAX.append(M_id)   
             
            
        MAX = pop_list(MAX,remove_MAX)
        MIN = pop_list(MIN,remove_MAX)
        left_minima = pop_list(left_minima,remove_MAX)
        if (len(MIN) > 0) : assert (len(right_minima) > 0)
        

    assert (len(right_minima) == len(MIN)),'right_minima: '+ str(right_minima) + ' MIN:' + str(MIN)   
    return(decreciente,left_minima,right_minima,MAX,MIN )


#%%
    
def get_pulses(theta,TH,W,PFE,PFI):
    '''
    get_pulses(theta,TH,W,PFE,PFI,alpha,D)
    Find the pulses of the time series theta for our definition
    
    Find the pulses of the time series theta, defined as starting in the PFI and ending in the PFI. 
    We consider that the signal is the sine of theta. This definition is for the Adler dynamical system.
    If there are no pulses, it returns empty lists.
    
    STEP 1 & 2 : detects the cosine maxima and minima, that approach the PFI and PFE. 
    STEP 3 & 4 : Computhes the beginning  and end of the pulses. 
    TEST       : verifies that all pulses have only one beggining, one maxima, one minima and one end in terms of the cosine of theta (**poner mejor que testea**). 
    STEP 5     : detects the sine maxima. 
    TEST       : TBD  
    
    Parameters
    ----------
    theta : list
        amplitude values of the time series given in angles
    TH : float
        threshold value
    W : float
        minimum distance between maxima / minima, given in indexes.
    PFE : float
        stable fixed point - angular 
    PFI : float
        unstable fixed point - angular 
    alpha : float
        parameter of the Adler deterministic equation
    D : float
        noise strengh of the dynamical system of this project

    
    Returns
    -------
    left_minima: list
        position of the starting point of each pulse 
    right_minima : list 
        position of the ending point of each pulse 
    MAX_sine : list 
        position of the maxima of a pulse of a given time series (in genral, sine of theta)

    
    '''
    

    t0 = time.time() ; #print('start searching for pulses: 5 steps pending')
    
    MAX,MIN = search_extremes(np.sin(theta),TH,W)
    t1 = time.time() - t0
    #print('finished searching local extremes',t1, 'sec')
    
    if (len(MAX)>0) and (len(MIN)>0): 
        MAX, MIN =  filter_extremes(MAX,MIN,np.sin(theta))
        t2 = time.time() - t1
       # print('finished setting 1 minimum between 2 maxima',t2, 'sec')
    
    if (len(MAX)>0) and (len(MIN)>0): 
        crec, left_minima, MAX ,MIN = get_left_minima(MAX ,MIN,PFE,PFI,theta)
        t3 = time.time() - t2
        #print('finished detecting left pulse minima ',t3, 'sec')
        
        dec, left_minima, right_minima,MAX ,MIN = get_rigth_minima(left_minima,MAX ,MIN,PFE,PFI,theta)
        t4 = time.time() - t3
        #print('finished detecting right pulse minima',t4, 'sec')
        
        #print(crec,dec)
        
        #print('testing parameters...')
        test_pulses(left_minima,right_minima,MAX)

        #MAX_sine , MIN_sine = search_extremes(np.sin(theta),TH,W)
        #MAX_sine, MIN_sine =  filter_extremes(MAX_sine,MIN_sine,np.sin(theta))        #en realidad, creo que si pondriamos otros parametros de max/min detection no haria falta hacer este paso en este lugar
        #left_minima,right_minima,MAX_sine =  filter_maxima_sine(left_minima,right_minima,MAX_sine)
        #t5 = time.time() - t4
       # print('step 5/5 finished: sine extremes detected',t5, 'sec')
        #test_pulses_sine(left_minima,right_minima,MAX_sine)
#
        
        return left_minima,right_minima,MAX
    else:
        print('The trace has no maxima or minima values available')
        return [],[],[]



def main_pulse_detection(theta,delta,omega,save_path_name,file_name):
    # calcula y guarda los pulsos de theta en el archivoPONER NOMBRES
    if len(theta) != 0:
        print('running pulse detection',delta)      
        TH = 0.90;W = 100
        PFE , PFI = get_fixed_points(delta)
        #print(theta[0:100],TH,W,PFE,PFI)
        left_minima,right_minima,MAX = get_pulses(theta,TH,W,PFE,PFI)
        print('pulse detection ended')
        
        if (len(MAX)>0): #and (len(MIN)>0): 
            #print('saving pulse detection results MAX:',MAX)
            save_data(MAX,save_path_name+'max_'+file_name)
            save_data(left_minima,save_path_name+'left_minima_'+file_name)
            save_data(right_minima,save_path_name+'right_minima_'+file_name)
            #print(file_name,'saving finished')
        else:
            print(delta,'no pulses on this condition -> not saving')      
    else:
        print('ERROR: theta with length zero',omega,delta)
    return(0)


    
def main_pulse_detection_(data_folder,save_path_name,tuple_,overwrite_flag = True):
   #ES UNA funcion auxiliar para paralelizar
   #definir delta es lo importante

    if 'alpha0' in tuple_[1].columns:
        (i,_,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        delta = np.round(i/omega,4)  
        file_name =  str(int(row.number.values[0]))+'_'+str(int(order))+'.pkl'
    
    else:
        (i,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        if 'alpha' in row.keys(): delta = np.round(i/omega,4)  
        if 'delta' in row.keys(): delta = i  
        file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    if check_file(file_name,data_folder):   
        if check_file('right_minima_'+file_name,data_folder) and (not overwrite_flag):
            print('pulse detection allready done and nor overwitting it',file_name)
            pass
        else:
            print('file name:',file_name)
            theta = download_theta(file_name,data_folder)
            main_pulse_detection(theta,delta,omega,save_path_name,file_name)
    else:
        print('ERROR: file not available',file_name)
    return(1)


def compute_pulse_detection(description_file,data_folder,save_path_name):
    #esta es la funcion que le calcula a cada seite temporal suspulsos
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    if 'alpha' in ref.keys() : tuple_ = ref.groupby(['alpha','D','order'])
    if 'delta' in ref.keys() : tuple_ = ref.groupby(['delta','D','order'])
    if 'alpha0' in ref.keys(): tuple_ = ref.groupby(['alpha0', 'sigma', 'tau','order'])
    main_pulse_detection__ = partial(main_pulse_detection_,data_folder,save_path_name)
    pool.map(main_pulse_detection__,tuple_)
    pool.close()
    pool.join()
    return (2)

#%%


def main_pulse_detection_exp(theta,delta,omega,save_path_name,file_name):
    # calcula y guarda los pulsos de theta en el archivoPONER NOMBRES
    #es para longitud de series experimentales! theta recortado
    if len(theta) != 0:
        print('running pulse detection',delta)      
        TH = 0.90;W = 100
        PFE , PFI = get_fixed_points(delta)

        left_minima,right_minima,MAX = get_pulses(theta,TH,W,PFE,PFI)
        print('pulse detection ended')
        
        if (len(MAX)>0): 
            dt,IPI,dm,joint_duration = get_pulses_quantifiers(left_minima,right_minima,MAX)
            consecutive_trial_st_ = consecutive_trial_st_exp(joint_duration,MAX,IPI,dm) 
            isolated_pulses, consecutive_trial = consecutive_trial_st_.get_consecutive_trains_of_pulses()#
            print('---------------------------------------------------',consecutive_trial,isolated_pulses)
            save_data(consecutive_trial,save_path_name+'exp_c_'+file_name)
        else:
            print(delta,'no pulses on this condition -> not saving')      
    else:
        print('ERROR: theta with length zero',omega,delta)
    return(0)

def main_pulse_detection_exp_(final_exp,data_folder,save_path_name,tuple_): 
    #ES UNA funcion auxiliar para paralelizar
    #definir delta es lo importante
    #es para longitud de series experimentales! theta recortado

    if 'alpha0' in tuple_[1].columns:
        (i,_,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        delta = np.round(i/omega,4)  
        file_name =  str(int(row.number.values[0]))+'_'+str(int(order))+'.pkl'
    
    else:
        (i,_,order),row = tuple_[0],tuple_[1]
        omega =  row.omega.unique()[0]
        if 'alpha' in row.keys(): delta = np.round(i/omega,4)  
        if 'delta' in row.keys(): delta = i  
        file_name =  str(int(row.number))+'_'+str(int(order))+'.pkl'
    
    if (check_file(file_name,data_folder) and order < len(final_exp)):   
            final_cell = final_exp[int(order)]
            theta = [i for i in download_theta(file_name,data_folder) if i < final_cell]
            main_pulse_detection_exp(theta,delta,omega,save_path_name,file_name)
    else:
        print('ERROR: file not available',file_name)
    return(1)



def compute_pulse_detection_exp(final_exp,description_file,data_folder,save_path_name):
    #esta es la funcion que le calcula a cada seite temporal suspulsos
    #es para la longitud final_exp
    #es para longitud de series experimentales! theta recortado

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    if 'alpha' in ref.keys() : tuple_ = ref.groupby(['alpha','D','order'])
    if 'delta' in ref.keys() : tuple_ = ref.groupby(['delta','D','order'])
    if 'alpha0' in ref.keys(): tuple_ = ref.groupby(['alpha0', 'sigma', 'tau','order'])
    main_pulse_detection_exp__ = partial(main_pulse_detection_exp_,final_exp,data_folder,save_path_name)
    pool.map(main_pulse_detection_exp__,tuple_)
    pool.close()
    pool.join()
    return (2)


#%%
#def test_pulses_OLD(left_minima,right_minima,MAX,MIN):
#    ''' 
#    test_pulses(left_minima,right_minima,MAX,MIN)
#    testing function for cosine outcome of the pulse detection
#    
#    Tests (1) that every vector has the same length, (2) that allways comes first the starting of the pulse,
#    then the maxima, then the minima, and then the end of the pulse.
#    
#    Parameters
#    ----------
#    left_minima: list
#        position of the starting point of each pulse for our definition
#    right_minima : list 
#        position of the ending point of each pulses for our definition
#    MAX : list 
#        position of the maximum of a pulse of a given time series (in genral, cosine of theta)
#    MIN : list 
#        position of the minimum of a pulse of a given time series (in genral, cosine of theta)   
#    
#    See also
#    ----------
#    - test_pulses_sine(left_minima,right_minima,MAX)
#
#    '''
#    
#    assert len(MIN) == len(MAX)
#    assert len(left_minima) == len(MAX), str(len(left_minima))+' is not equal to ' + str(len(MAX))
#    assert len(right_minima) == len(MAX) , str(len(right_minima))+' is not equal to ' + str(len(MAX))
#    
#    #los proximos 3 iguales estan mal, pero los ponemos porque fallan cuando hay mucho ruido
#    assert all([(i<=j)*1 for (i,j) in zip(left_minima,MAX)]), '----- ' + str([(j-i)*1 for (i,j) in zip(left_minima,MAX)])
#    assert all([(i<=j)*1 for (i,j) in zip(MAX,MIN)])
#    assert all([(i<=j)*1 for (i,j) in zip(MIN,right_minima)])
#    
#    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]),str([(i-j) for (i,j) in zip(left_minima[1:],right_minima[:-1])]) 

def test_pulses(left_minima,right_minima,MAX):
    ''' 
    test_pulses(left_minima,right_minima,MAX,MIN)
    testing function for cosine outcome of the pulse detection
    
    Tests (1) that every vector has the same length, (2) that allways comes first the starting of the pulse,
    then the maxima, and then the end of the pulse.
    
    Parameters
    ----------
    left_minima: list
        position of the starting point of each pulse for our definition
    right_minima : list 
        position of the ending point of each pulses for our definition
    MAX : list 
        position of the maximum of a pulse of a given time series (in genral, sine of theta)
 

    See also
    ----------
    - test_pulses(left_minima,right_minima,MAX,MIN)

    '''
    #assert len(MIN) == len(MAX)
    assert len(left_minima) == len(MAX), str(len(left_minima))+' (sine checking) is not equal to ' + str(len(MAX))
    assert len(right_minima) == len(MAX) , str(len(right_minima))+' (sine checking) is not equal to ' + str(len(MAX))
    
    assert all([(i<j)*1 for (i,j) in zip(left_minima,MAX)]), '--(sine checking)--- ' + str([(j-i)*1 for (i,j) in zip(left_minima,MAX)])
    assert all([(i<j)*1 for (i,j) in zip(MAX,right_minima)])
    
    #the equal is because they can share minima 
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]),str([(i-j) for (i,j) in zip(left_minima[1:],right_minima[:-1])])






