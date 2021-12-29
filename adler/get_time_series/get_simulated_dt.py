import numpy as np
from adler.data_managing_functions import save_data, get_fixed_points,compute_theoretical_omega
from functools import partial
from itertools import product
import time
import multiprocessing as mp


#%%
# =============================================================================
# =============================================================================
# =============================================================================
#                   GET SIMULATED DT MODULE
# =============================================================================
# =============================================================================
# =============================================================================
#%%
def get_simulated_dt(dt,omega,alpha,epsilon):
    '''
    Gets the pulse duration of the Adler phase model with an epsilon on the fixews point. 
    
    ------------------------------------------------------------
    INPUTS:
        - omega,alpha :  float
                     parameters of the Adler deterministic equation
        - epsilon : float
            delta a partir del punto fijo strengh
        - dt : positive float 
            time steps of
        the simulation

    ------------------------------------------------------------
    OUTPUTS:
        - theta : list 
                a pulse in angular variables
        - DT : float
            pulse duration (minutes)
    '''
    
    #### Initial conditions ############################
    ####################################################
    
    theta_end,theta_beg = get_fixed_points(alpha/omega)
    theta_beg = theta_beg - 2* np.pi + epsilon
    if theta_beg >= theta_end:
        print('wrong epsilon : ', epsilon)
        return([],0)
    else:
        D = 0
        np.random.seed()
        theta = []
        
        #### Time evolution ################################
        ####################################################
        steps = 0
        theta_past = theta_beg
        
        
        while theta_past < (theta_end-epsilon):
            
            
            u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
            k = dt * (omega + alpha * np.sin(theta_past))
            l = np.sqrt(dt * 2 * D) * u
    
            theta_present = theta_past + dt/2 * (2*omega + alpha * (np.sin(theta_past) + np.sin(theta_past + l + k) )) + np.sqrt(dt * 2 * D) * u
    
            
            theta.append(theta_present)            
            theta_past = theta_present
            
            if theta_past <= theta_beg:
                steps = 0
                print(theta_past,theta_end)
            else:
                steps = steps + 1
        print('eps : ',epsilon,' DT :', steps*dt)
        return(theta,steps*dt)

def get_simulated_dt_alpha(main_filename,N,dt,p):
    (omega,epsilon),alphas = p
    DT = [] ; theta_example = []
    for alpha in alphas:
        DT_aux = []
        for i in range(N):
            theta,dt_aux = get_simulated_dt(dt,omega,alpha,epsilon)
            DT_aux.append(dt_aux)
        DT.append(DT_aux)
        theta_example.append(theta[::100])
    
    
    DT_filename = main_filename + 'simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(epsilon)+'.pkl'
    theta_filename = main_filename +  'theta_example_simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(epsilon)+'.pkl'

    save_data(DT, DT_filename)
    save_data(theta_example, theta_filename)
    
    return(0)
    
def get_simulated_dt_mp(main_filename,dt,N,params,nproc = mp.cpu_count()):
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    get_simulated_dt_alpha_ = partial(get_simulated_dt_alpha,main_filename,N,dt)
    
    pool.map(get_simulated_dt_alpha_,get_all_combinations_alphas(params) ) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)
    

def get_all_combinations_alphas(params):
    ''' An iterable containing all the possible 
    combinations of D and omega values, jointly with alphas list.
    Returns an tuple
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    for p in product(*(params['omega'],params['epsilon'])):
        yield p,params['alpha']





#%%
# =============================================================================
# =============================================================================
# =============================================================================
#                   GET TIME INVARIANT DT MODULE
# =============================================================================
# =============================================================================
# =============================================================================
#%%




#%%
def get_simulated_dt_fixed(dt,delta,T_0,epsilon):
    '''
    Gets the pulse duration of the Adler phase model with an epsilon on the fixews point. 
    
    ------------------------------------------------------------
    INPUTS:
        - omega,alpha :  float
                     parameters of the Adler deterministic equation
        - epsilon : float
            delta a partir del punto fijo strengh
        - dt : positive float 
            time steps of the simulation

    ------------------------------------------------------------
    OUTPUTS:
        - theta : list 
                a pulse in angular variables
        - DT : float
            pulse duration (minutes)
    '''
    
    #### Initial conditions ############################
    ####################################################
    
    theta_end,theta_beg = get_fixed_points(delta)
    theta_beg = theta_beg - 2* np.pi + epsilon
    if theta_beg >= theta_end:
        print('wrong epsilon : ', epsilon)
        return([],0)
    else:
        D = 0
        np.random.seed()
        theta = []
        
        #### Time evolution ################################
        ####################################################
        steps = 0
        theta_past = theta_beg
        
        omega = 2 * np.pi/T_0 *compute_theoretical_omega(epsilon,delta)
        
        while theta_past < (theta_end-epsilon):
            
            
            u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
            k = dt * omega*(1 + delta * np.sin(theta_past))
            l = np.sqrt(dt * 2 * D) * u
    
            theta_present = theta_past + dt/2 * omega*(2*1 + delta * (np.sin(theta_past) + np.sin(theta_past + l + k) )) + np.sqrt(dt * 2 * D) * u
    
            
            theta.append(theta_present)            
            theta_past = theta_present
            
            if theta_past <= theta_beg:
                steps = 0
                print(theta_past,theta_end)
            else:
                steps = steps + 1
        print('eps : ',epsilon,' DT :', steps*dt)
        return(theta,steps*dt)

#%%

def get_simulated_dt_fixed_alpha(main_filename,N,dt,p):
    (T_0,epsilon),deltas = p

    DT = [] ; theta_example = []
    for delta in deltas:
        DT_aux = []
        
        for i in range(N):
            theta,dt_aux = get_simulated_dt_fixed(dt,delta,T_0,epsilon)
            DT_aux.append(dt_aux)
       
        DT.append(DT_aux)
        theta_example.append(theta[::100])
    
    
    DT_filename = main_filename + 'simulated_fixed_dt_'+str(np.round(T_0,2))+'_eps_'+str(epsilon)+'.pkl'
    theta_filename = main_filename +  'theta_example_simulated_fixed_dt_'+str(np.round(T_0,2))+'_eps_'+str(epsilon)+'.pkl'

    save_data(DT, DT_filename)
    save_data(theta_example, theta_filename)
    
    return(0)
    
def get_simulated_dt_fixed_mp(main_filename,dt,N,params,nproc = mp.cpu_count()):
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    get_simulated_dt_alpha_ = partial(get_simulated_dt_fixed_alpha,main_filename,N,dt)
    
    pool.map(get_simulated_dt_alpha_,get_all_combinations_fixed_alphas(params) ) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)
    

def get_all_combinations_fixed_alphas(params):
    ''' An iterable containing all the possible 
    combinations of epsilon and T_o values, jointly with delta list.
    Returns an tuple
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    for p in product(*(params['T0'],params['epsilon'])):
        yield p,params['delta']