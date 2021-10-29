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
def filter_maxima_sine(left_minima,right_minima,MAX):
    ''' filter_maxima_sine(left_minima,right_minima,MAX)
    Removes the sine maxima that are not in between two minima that define a pulse, 
    and the minima that do not have a maxima in between the, 
    
    Removes the sine maxima that are not in between two minima that define a pulse.
    Remove the edges of the sine maxima so that they are between the pulse minima, and discard the rest.
    Remove the case in which two minima do not have a maxima in between, case that can happen if because of noise the system reaches the
    fixed points but the turn back. 
    
    Parameters
    ----------
    left_minima: list
        position of the starting point of each pulses 
    right_minima : list 
        position of the ending point of each pulses 
    MAX : list 
        position of the maxima of a pulse of a given time series (in genral, sine of theta)
    
    Returns
    -------
    left_minima: list
        filtered position of the starting point of each pulses 
    right_minima : list 
        filtered position of the ending point of each pulses 
    MAX : list 
        filtered position of the maxima of a pulse of a given time series (in genral, sine of theta)
        
    '''
    
    if len(left_minima) > 0 and len(right_minima) > 0 and len(MAX) > 0:
        MAX_remove_list_index = []
        
        for i,M in enumerate(MAX):
            if M < left_minima[0]:
                MAX_remove_list_index.append(i)
            if M > right_minima[-1]:
                MAX_remove_list_index.append(i)
        print('MAX_remove_list_index:', MAX_remove_list_index)
        
        #remuev el caso en que nunca tocó el maximo del seno, pero por ruido llego al del coseno
        NEW_left_minima, NEW_right_minima  = [], []
    
        for M in MAX:
            for l,r in zip(left_minima,right_minima):
                assert l < r    

                if (l < M) and (r > M):
                    #hago esto y no los elimino de la lista por si existen minimos compartidos
                    NEW_left_minima.append(l)
                    NEW_right_minima.append(r)
                    break
                else:
                    pass
        
        return(NEW_left_minima,NEW_right_minima,pop_list(MAX,MAX_remove_list_index))
    else:
        return ([],[],[])


#%%
    
def get_pulses(theta,TH,W,PFE,PFI,alpha,D):
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
        
        print('testing cosine parameters alpha, D: ',alpha,D)
        test_pulses(left_minima,right_minima,MAX_cos,MIN_cos)

        MAX_sine , MIN_sine = search_extremes(np.sin(theta),TH,W)
        MAX_sine, MIN_sine =  filter_extremes(MAX_sine,MIN_sine,np.sin(theta))        #en realidad, creo que si pondriamos otros parametros de max/min detection no haria falta hacer este paso en este lugar
        left_minima,right_minima,MAX_sine =  filter_maxima_sine(left_minima,right_minima,MAX_sine)
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
def test_pulses(left_minima,right_minima,MAX,MIN):
    ''' 
    test_pulses(left_minima,right_minima,MAX,MIN)
    testing function for cosine outcome of the pulse detection
    
    Tests (1) that every vector has the same length, (2) that allways comes first the starting of the pulse,
    then the maxima, then the minima, and then the end of the pulse.
    
    Parameters
    ----------
    left_minima: list
        position of the starting point of each pulse for our definition
    right_minima : list 
        position of the ending point of each pulses for our definition
    MAX : list 
        position of the maximum of a pulse of a given time series (in genral, cosine of theta)
    MIN : list 
        position of the minimum of a pulse of a given time series (in genral, cosine of theta)   
    
    See also
    ----------
    - test_pulses_sine(left_minima,right_minima,MAX)

    '''
    
    assert len(MIN) == len(MAX)
    assert len(left_minima) == len(MAX), str(len(left_minima))+' is not equal to ' + str(len(MAX))
    assert len(right_minima) == len(MAX) , str(len(right_minima))+' is not equal to ' + str(len(MAX))
    
    #los proximos 3 iguales estan mal, pero los ponemos porque fallan cuando hay mucho ruido
    assert all([(i<=j)*1 for (i,j) in zip(left_minima,MAX)]), '----- ' + str([(j-i)*1 for (i,j) in zip(left_minima,MAX)])
    assert all([(i<=j)*1 for (i,j) in zip(MAX,MIN)])
    assert all([(i<=j)*1 for (i,j) in zip(MIN,right_minima)])
    
    assert all([(i>=j)*1 for (i,j) in zip(left_minima[1:],right_minima[:-1])]),str([(i-j) for (i,j) in zip(left_minima[1:],right_minima[:-1])]) 

def test_pulses_sine(left_minima,right_minima,MAX):
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






