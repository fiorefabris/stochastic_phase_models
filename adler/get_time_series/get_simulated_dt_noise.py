import numpy as np
from adler.data_managing_functions import save_data, get_fixed_points
from functools import partial
import time
import multiprocessing as mp
from itertools import product

def all_combinations(params):
    ''' ES EL MISMO QUE EL MAIN PERO NO SE COMOIMPORTARLO Y ES TARDE! CAMBIAR Y EXPORTAR
    
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
#%%
# =============================================================================
# =============================================================================
# =============================================================================
#              CONDITIONAL FIRST PASSAGE TIME MODULE
# =============================================================================
# =============================================================================
# =============================================================================
#%%

def get_cond_prob(init,dt,T,omega,alpha,D):
    ''' Simulated data of the Adler phase model with gaussian white noise. 
    CAAAMMMBIIARR
    init entre el PFI (cuadrante -1) y PFE(cuadrante 3)
    
    ------------------------------------------------------------
    INPUTS:
        - omega,alpha (real numbers): parameters of the Adler deterministic
         equation
        - D (real number): noise strengh
        - dt (positive real number, default dt=0.0001): time steps of
        the simulation
        - T (positive real number > dt, default T = 10000): Total time 
        of the simulation

    ------------------------------------------------------------
    OUTPUTS:
        - t : time
        - theta : phase simulated variables

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
    return(test,theta)

#%%

def get_cond_prob_pop(main_filename,dt,T,p):
    omega,alpha,D = p
    print(omega)
    total = 100000
    PFE,PFI = get_fixed_points(alpha/omega)
    PFI = PFI - 2* np.pi

    initial_conditions = np.linspace(PFI,PFE,1000)
    cond_prob = []

    for init in initial_conditions:        
        suc = 0
        for j in range(total):
            test, _ = get_cond_prob(init,dt,T,omega,alpha,D)
            if test == 2:
                suc = suc+1
       # print(suc/total) #numero que va entre 0 y 1
        cond_prob.append(suc)
    
    cond_prob_filename = main_filename +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    initial_conditions_filename = main_filename +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
    save_data(cond_prob, cond_prob_filename)
    save_data(initial_conditions, initial_conditions_filename)

#return(initial_conditions,cond_prob)

#%%

    
def get_cond_prob_pop_mp(main_filename,dt,T,params,nproc = mp.cpu_count()):
    
    
    t0= time.perf_counter(); print('starting...')
    pool = mp.Pool(processes= nproc)
    
    get_cond_prob_pop_ = partial(get_cond_prob_pop,main_filename,dt,T)
    
    pool.map(get_cond_prob_pop_,all_combinations(params)) 
    pool.close() 
    pool.join()
    t1 = time.perf_counter() - t0
    print("time elapsed: ", t1)
    return(0)
    





