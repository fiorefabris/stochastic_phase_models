import numpy as np
from adler.data_managing_functions import save_data, get_fixed_points
from functools import partial
import time
import multiprocessing as mp
from itertools import product

def all_combinations(params):
    '''Is the same function as in the main module
    
    -----------------
    INPUT:
        -params: a dictionary of lists
    '''
    return product(*params.values())


def check_interval(x,PFE,PFI):
    if x <= PFI:
        return 0
    elif x>= PFE:
        return 2
    else:
        return 1
    
def get_geomspace(PFI,PFE,n):
    """
    n : cantidad de elementos que tiene la lista que te devuelve
    """
    aux = []
    for i in np.geomspace(1,PFE-PFI+1,n):
        aux.append(i-(-PFI+1))
    return (aux)
#%%
# =============================================================================
# =============================================================================
# =============================================================================
#              CONDITIONAL FIRST PASSAGE TIME MODULE
# =============================================================================
# =============================================================================
# =============================================================================
#%%

def get_epsilon_plus(init,dt,T,omega,alpha,D):
    ''' Axilliary function for computing the conditional First PAssage time probability epsilon plus, and
    measures the conditional first passage time (in steps, no time)
    
    
    The starting point is on the PFI on the -1 cuadrant, 
    and the ending point PFE is in the 3rd cuadrant.
    
    Calcula una trayectoria que empieza en init y termina en PFI (test = 0, unsuccesfull) o PFE (test = 2, succesfull)
    
    ------------------------------------------------------------
    INPUTS:
        - init (real) : initial condition. Is an angle between PFI (cuadrant -1) and PFE (cuadrant 3)
        - dt (positive real number, default dt=0.0001): time steps of
        the simulation
        - T (positive real number > dt, default T = 10000): Total time 
        of the simulation
        - omega, alpha, D: parameters of the noisy adler equation

    ------------------------------------------------------------
    OUTPUTS:
    test,theta,i+1
        - test : si el camino fue satisfactorio (2) o no satisfactorio (0)
        - theta : phase simulated variables
        - i + 1: total steps of the simulation

'''

    n     = int(T/dt) # total number of steps
    
    #variables
    theta = []
     

    #### Initial & final conditions ####################
    ####################################################
    np.random.seed()
    PFE,PFI = get_fixed_points(alpha/omega)
    theta_past = init
    assert (theta_past >= -np.pi/2) and (theta_past <= PFE)


    
    #### Time evolution ################################
    ####################################################
#
    for i in range(n-1):
#        if i % d == 0:
#            theta.append(theta_past)
    
        theta.append(theta_past)
        u = np.sqrt(-2* np.log(np.random.uniform(0,1))) * np.cos(2*np.pi*np.random.uniform(0,1))
        k = dt * (omega + alpha * np.sin(theta_past))
        l = np.sqrt(dt * 2 * D) * u

        theta_present = theta_past + dt/2 * (2*omega + alpha * (np.sin(theta_past) + np.sin(theta_past + l + k) )) + np.sqrt(dt * 2 * D) * u
        theta_past = theta_present


        test = check_interval(theta_past,PFE,PFI - 2* np.pi )
        
        if test == True:
            pass
        else:
            theta.append(theta_past)
            break

    assert (test == 0) or (test == 2)
    #plt.plot(np.sin(theta))
    return(test,theta,i+1)

#%%

def get_epsilon_plus_pop(main_filename,dt,T,p):
    """
    
    mide en steps solamente los que son success
    
    para cada condicion inicial corre get_epsilon_plus total veces
    """
    omega,alpha,D = p
   # print(omega)
    total = 500;     t0= time.perf_counter(); 
    PFE,PFI = get_fixed_points(alpha/omega)
    PFI = PFI - 2* np.pi
    #print(alpha/omega,PFI/np.pi*2,PFE/np.pi*2)
    initial_conditions = get_geomspace(PFI,PFE,100)
    
    cond_prob = []
    steps_plus = []
    
    for init in initial_conditions:        
        suc = 0
        steps_plus_aux = []
        
        for j in range(total):
                print("j : ",j, " time elapsed : ", time.perf_counter() - t0)
                test, _ , steps = get_epsilon_plus(init,dt,T,omega,alpha,D)
                if test == 2:
                    suc = suc+1
                    steps_plus_aux.append(steps)
           # print(suc/total) #numero que va entre 0 y 1
        cond_prob.append(suc)
        steps_plus.append(steps_plus_aux)

    cond_prob_filename = main_filename +  'cond_prob_omega_'+str(omega)+'_alpha_'+str(alpha/omega)+'_D_'+str(D)+'.pkl'    
    save_data((initial_conditions,cond_prob), cond_prob_filename)
    step_plus_filename = main_filename +  'step_plus_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    save_data((initial_conditions,steps_plus), step_plus_filename)

    return(0)
        
#    cond_prob_filename = main_filename +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
#    initial_conditions_filename = main_filename +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
#    save_data(cond_prob, cond_prob_filename)
#    save_data(initial_conditions, initial_conditions_filename)
#    step_plus_filename = main_filename +  'step_plus_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
#    save_data(steps_plus, step_plus_filename)


#return(initial_conditions,cond_prob)

#%%

    
def get_epsilon_plus_pop_mp(main_filename,dt,T,params,nproc = mp.cpu_count()):
    
    
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    
    get_epsilon_plus_pop_ = partial(get_epsilon_plus_pop,main_filename,dt,T)
    
    pool.map(get_epsilon_plus_pop_,all_combinations(params)) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)
    

#%%
# =============================================================================
# ESTO ES PARA CALCULAR MUCHA ESTADISTICA EN X MENOS
# =============================================================================

def get_epsilon_plus_in_x_minus_mp(main_filename,dt,T,params,nproc = mp.cpu_count()):
    
    
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    
    get_epsilon_plus_in_x_minus_ = partial(get_epsilon_plus_in_x_minus,main_filename,dt,T)
    
    pool.map(get_epsilon_plus_in_x_minus_,all_combinations(params)) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)



def get_epsilon_plus_in_x_minus(main_filename,dt,T,p):
    """
    
    la condicion inicial solamente es el PFI
    
    """
    omega,alpha,D = p
   # print(omega)
    total = 500;     t0= time.perf_counter(); 
    PFE,PFI = get_fixed_points(alpha/omega)
    PFI = PFI - 2* np.pi
    #print(alpha/omega,PFI/np.pi*2,PFE/np.pi*2)
    initial_conditions = [PFI]
    
    cond_prob = []
    steps_plus = []
    
    for init in initial_conditions:        
        suc = 0
        steps_plus_aux = []
        
        for j in range(total):
                print("j : ",j, " time elapsed : ", time.perf_counter() - t0)
                test, _ , steps = get_epsilon_plus(init,dt,T,omega,alpha,D)
                if test == 2:
                    suc = suc+1
                    steps_plus_aux.append(steps)
           # print(suc/total) #numero que va entre 0 y 1
        cond_prob.append(suc)
        steps_plus.append(steps_plus_aux)
        
    cond_prob_filename = main_filename +  'x_minus_cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    #initial_conditions_filename = main_filename +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    save_data(cond_prob, cond_prob_filename)
    #save_data(initial_conditions, initial_conditions_filename)
    step_plus_filename = main_filename +  'x_minus_step_plus_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    save_data(steps_plus, step_plus_filename)
