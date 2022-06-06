'''
@ Fiorella Fabris
September 22, 2020.
This are the main functions for simulating time series of Adler phase model with white noise. 

'''
import multiprocessing as mp
import numpy as np
from itertools import product
from functools import partial
import pickle
import pandas as pd
from math import ceil
import time
from adler.data_managing_functions import compute_theoretical_omega


#%%

def compute_time_series(save_path_name_description,save_path_name_data,params,dt,T,d,N):
    # =============================================================================
    # Generates the description
    # =============================================================================
    
    #save_path_name_file_name = '/home/fiore/running_25_10_2021/coding/description.xlsx'
    get_file_references(params,N,save_path_name_description)
    
    # =============================================================================
    #   Generates the trazes  
    # =============================================================================
    
    #save_path_name = '/home/fiore/running_25_10_2021/data/'
    explore_param_space_time_evolution_(params, save_path_name_data,dt,T,d,N)

#%%
# Time evolution funtion 
def time_evolution_(dt,T,d,omega,alpha,D):
    ''' Simulated data of the Adler phase model with gaussian white noise. 
    
    ------------------------------------------------------------
    INPUTS:
        - omega,alpha (real numbers): parameters of the Adler deterministic
         equation
        - D (real number): noise strengh
        - dt (positive real number, default dt=0.0001): time steps of
        the simulation
        - T (positive real number > dt, default T = 10000): Total time 
        of the simulation
        - d (integer number, default d = 5): decimation factor

    ------------------------------------------------------------
    OUTPUTS:
        - t : time
        - theta : phase simulated variables

'''

    n     = int(T/dt) # total number of steps
    
    #variables
    theta = np.zeros(ceil(n/d))
     

    #### Initial conditions ############################
    ####################################################
    np.random.seed()
    theta_past = np.random.uniform(0,2)*np.pi
    
    #### Time evolution ################################
    ####################################################

    for i in range(n-1):
        u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
        k = dt * (omega + alpha * np.sin(theta_past))
        l = np.sqrt(dt * 2 * D) * u

        theta_present = theta_past + dt/2 * (2*omega + alpha * (np.sin(theta_past) + np.sin(theta_past + l + k) )) + np.sqrt(dt * 2 * D) * u
        if i % d == 0:
            theta[i // d] = theta_past
                    
        theta_past = theta_present


    return(theta)

#%%
# ITERATORS    

def all_combinations(params):
    ''' An iterable containing all the possible 
    combinations of params values
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    return product(*params.values())

def repeat(iterable, N=1):
    ''' A generator containing N repetitions of the iterable
    '''
    for item in iterable:
        for _ in range(N):
            yield item
            


def repeat_all_combinations(params,N=1):
    ''' A generator containing N repetitions of all the possible 
    combinations of params values
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    return repeat(all_combinations(params),N)

def params_names_and_values(params):
    ''' A tuple generator: containing 1 repetitions of all the possible 
    combinations of params values (1) and the order value of each pair of parameters (0)
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    for param_count,param in enumerate(all_combinations(params)):
        yield (param_count,param)

#%% repeat and list iterators
    
def repeat_and_list(iterable, N=1):
    ''' A generator containing N repetitions of the iterable, 
    the order number of each combination parameter,
    the order parameter of the N values,
    '''
    name=0
    for item in iterable:
        name = name+1
        for order in range(N):
            yield item,name,order,

def repeat_and_list_all_combinations(params,N=1):
    ''' A generator containing N repetitions of all the possible 
    combinations of params values, the order parameter  
    and the order number of each combination parameter
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    return repeat_and_list(all_combinations(params),N)

#%%

def mp_time_evolution_(param):
    # Function for calling the main function with a tuple of parameters
    time_evolution__ = partial(time_evolution_,dt,T,d)
    return time_evolution__(*param)    


def mp_time_evolution_and_list(main_file_name,dt,T,d,param):
    # Function for calling the main function with a tuple of parameters
    time_evolution_save(*param,main_file_name,dt,T,d) 
    return(0)


def time_evolution_save(param,number,order,main_file_name,dt,T,d):
    t0= time.perf_counter()
    time_evolution__ = partial(time_evolution_,dt,T,d)
    theta = time_evolution__(*param) 
    file_name =  str(number)+'_'+str(order)+'.pkl'
    save_data(theta, main_file_name + file_name)
    t1 = time.perf_counter() - t0
    print(file_name,t1)
    return(0)


        
#%% download and save files
        
def download_data(filename):
    infile =  open(filename,'rb')
    results = pickle.load(infile)
    infile.close()
    return(results)
    
def save_data(data, filename):
    outfile= open(filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()

        
        
#%%
def explore_param_space_time_evolution_(params, main_file_name,dt,T,d,N=1, nproc = mp.cpu_count()):
    ''' Finds the numerical solution of Adler equation whit additive white noise 
     for N repetitions of each possible combinations of params values. 
    
    ------------------------------------------------------------
    INPUTS:
        - params (dictionary of lists): parameters of the equation 
        
        - N (int - default = 1): number of repetitions of each combination 
        of parameters
        - nrpoc (int) : number of parallel processes. The default value
        is the number of cpus of the local server.
        
    ------------------------------------------------------------
    OUTPUT:
        - results (list of lists): (t,theta) values for each repetition
        of each parameters combination
    '''
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    mp_time_evolution_and_list_ = partial(mp_time_evolution_and_list, main_file_name,dt,T,d)
    pool.map(mp_time_evolution_and_list_,repeat_and_list_all_combinations(params,N) ) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)

#%%
    
def get_file_references(params,N,save_path_name_file_name):
    with pd.ExcelWriter(save_path_name_file_name) as writer:
        #the references of the space parameter we are simulating
        description = pd.DataFrame(index=np.arange(len(list(all_combinations(params)))),columns=params.keys()) 
        for i,param in enumerate(all_combinations(params)):
            for c,j in enumerate(params.keys()):
                description.iloc[i][j] = param[c]  
        description.to_excel(writer, sheet_name='Plotting_parameters')


        description_file = pd.DataFrame(index=np.arange(len(list(repeat_and_list_all_combinations(params,N)))),columns=list(params.keys())+['number','order']) 
        for i,param in enumerate(repeat_and_list_all_combinations(params,N)):
            for c,j in enumerate(params.keys()):
                description_file.iloc[i][j] = param[0][c]
            description_file.iloc[i][-1] = param[2] 
            description_file.iloc[i][-2] = param[1] 
        description_file.to_excel(writer, sheet_name='File_references')
        
    writer.save()
    print('ref done!')
    return(0)

#%%
# =============================================================================
# =============================================================================
# =============================================================================
# Module for computing the time series-time invariant (non-noise hypothesis NNH)
# =============================================================================
# =============================================================================
# =============================================================================

def time_invariant_NNH_time_evolution_(dt,T,d,T_0,delta,D):
    ''' Simulated data of the Adler phase model with gaussian white noise. 
    
    ------------------------------------------------------------
    INPUTS:
        - delta (real numbers): parameters of the Adler deterministic
         equation (alpha/omega)
        - D (real number): noise strengh
        - dt (positive real number, default dt=0.0001): time steps of
        the simulation
        - T (positive real number > dt, default T = 10000): Total time 
        of the simulation
        - d (integer number, default d = 5): decimation factor
        - T_0 (integer number): characteristic time scale of pulsing
        
    ------------------------------------------------------------
    OUTPUTS:
        - t : time
        - theta : phase simulated variables

'''

    
    #variables
    n     = int(T/dt) # total number of steps
    theta = np.zeros(ceil(n/d))
    omega = 2 * np.pi/T_0 *compute_theoretical_omega(10e-5,delta)

    #### Initial conditions ############################
    ####################################################
    np.random.seed()
    theta_past = np.random.uniform(0,2)*np.pi
    
    #### Time evolution ################################
    ####################################################

    for i in range(n-1):
        u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
        k = dt * omega * (1 + delta * np.sin(theta_past))
        l = np.sqrt(dt * 2 * D) * u

        theta_present = theta_past + dt/2 * omega* (2*1 + delta * (np.sin(theta_past) + np.sin(theta_past + l + k) )) + np.sqrt(dt * 2 * D) * u
        
        if i % d == 0:
            theta[i // d] = theta_past
                    
        theta_past = theta_present


    return(theta)
#%%
def compute_time_series_NNH(save_path_name_description,save_path_name_data,params,dt,T,d,N):
    
    get_file_references(params,N,save_path_name_description)
    
    explore_param_space_time_evolution_NNH_(params, save_path_name_data,dt,T,d,N)

#%%
def explore_param_space_time_evolution_NNH_(params, main_file_name,dt,T,d,N=1, nproc = mp.cpu_count()):
    ''' 
    '''
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    mp_time_evolution_and_list_ = partial(mp_time_evolution_and_list_NNH, main_file_name,dt,T,d)
    pool.map(mp_time_evolution_and_list_,repeat_and_list_all_combinations(params,N) ) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)

#%%



def mp_time_evolution_and_list_NNH(main_file_name,dt,T,d,param):
    # Function for calling the main function with a tuple of parameters
    time_evolution_save_NNH(*param,main_file_name,dt,T,d) 
    return(0)


def time_evolution_save_NNH(param,number,order,main_file_name,dt,T,d):
    t0= time.perf_counter()
    time_evolution__ = partial(time_invariant_NNH_time_evolution_,dt,T,d)
    theta = time_evolution__(*param) 
    file_name =  str(number)+'_'+str(order)+'_NNH.pkl'
    save_data(theta, main_file_name + file_name)
    t1 = time.perf_counter() - t0
    print(file_name,t1)
    return(0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#%%
# =============================================================================
# =============================================================================
# =============================================================================
# Module for computing the time series with colored noise in alpha
# =============================================================================
# =============================================================================
# =============================================================================
    
    
# Time evolution funtion 
def time_evolution_ou(dt,T,d,omega,alpha0,sigma,tau):
    ''' Simulated data of the Adler phase model with colopr noise in alpha. 
    
    ------------------------------------------------------------
    INPUTS:

    ------------------------------------------------------------
    OUTPUTS:
        - t : time
        - theta : phase simulated variables

'''

    n     = int(T/dt) # total number of steps
    
    #variables
    theta = np.zeros(ceil(n/d))
     

    #### Initial conditions ############################
    ####################################################
    np.random.seed()
    theta_past = np.random.uniform(0,2)*np.pi
    alpha_past = alpha0
    
    #### Time evolution ################################
    ####################################################

    for i in range(n-1):
        theta_present,alpha_present = get_theta_present(theta_past,alpha_past,omega,alpha0,sigma,tau,dt)
        if i % d == 0:
            theta[i // d] = theta_past
                    
        theta_past = theta_present
        alpha_past = alpha_present


    return(theta)
    

def get_ou_process(alpha_past, alpha0,sigma,tau,dt):
    u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
    k = dt/tau * (alpha0 - alpha_past)
    l = np.sqrt(dt) * (np.sqrt(2/tau)*sigma) * u 

    alpha_present= alpha_past + dt/(2*tau) *((alpha0 - alpha_past) + (alpha0-(alpha_past + k + l))) + np.sqrt(dt) * u * (np.sqrt(2/tau)*sigma)
    return alpha_present

def get_theta_present(theta_past,alpha_past,omega,alpha0,sigma,tau,dt):
    alpha_present = get_ou_process(alpha_past, alpha0,sigma,tau,dt) 
    k = dt * omega
    l = (dt* alpha_past) * np.sin(theta_past)
    theta_present = theta_past + dt * omega +  (dt* alpha_past)/2 * (np.sin(theta_past) + np.sin(theta_past+k+l))
    return(theta_present,alpha_present)
    
#%%
def compute_time_series_ou(save_path_name_description,save_path_name_data,params,dt,T,d,N):
    # =============================================================================
    # Generates the description
    # =============================================================================
    
    #save_path_name_file_name = '/home/fiore/running_25_10_2021/coding/description.xlsx'
    get_file_references(params,N,save_path_name_description)
    
    # =============================================================================
    #   Generates the trazes  
    # =============================================================================
    
    #save_path_name = '/home/fiore/running_25_10_2021/data/'
    explore_param_space_time_evolution_ou(params, save_path_name_data,dt,T,d,N)


def explore_param_space_time_evolution_ou(params, main_file_name,dt,T,d,N=1, nproc = mp.cpu_count()):
    ''' hhh
    '''
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    mp_time_evolution_and_list_ = partial(mp_time_evolution_and_list_ou, main_file_name,dt,T,d)
    pool.map(mp_time_evolution_and_list_,repeat_and_list_all_combinations(params,N) ) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)


def mp_time_evolution_and_list_ou(main_file_name,dt,T,d,param):
    # Function for calling the main function with a tuple of parameters
    time_evolution_save_ou(*param,main_file_name,dt,T,d) 
    return(0)


def time_evolution_save_ou(param,number,order,main_file_name,dt,T,d):
    t0= time.perf_counter()
    time_evolution__ = partial(time_evolution_ou,dt,T,d)
    theta = time_evolution__(*param) 
    file_name =  str(number)+'_'+str(order)+'.pkl'
    save_data(theta, main_file_name + file_name)
    t1 = time.perf_counter() - t0
    print(file_name,t1)
    return(0)