import pickle
import os
import numpy as np
from math import ceil


def download_data(filename):
    with open(filename,'rb') as infile:
        results = pickle.load(infile)
    return(results)

def save_data(data, filename):
    with open(filename,'wb') as outfile:
        pickle.dump(data,outfile,protocol=-1)
    

def check_file(file_name,data_folder):
    flag = False
    if os.path.isfile(data_folder + file_name):
        if os.path.getsize(data_folder + file_name) > 0:
            flag = True
    return(flag)


def time(dt,T,d):
    '''
    time(dt,T,d)
    Computes the time vector of a desired time series.
    
    Parameters
    -----------
        dt : float
            time resolution of the time series. 
        T : float
            total time lenght of the time series.
        d : integer
            the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    
    Returns
    -------
        out : numpy array
            time vector of the desired time series.
            
    Notes
    ------
    dt * d is the sampling rate.
    '''
    return(np.linspace(0,T,ceil(int(T/dt)/d)) )
    
def get_fixed_points(delta):
    '''
    get_fixed_points(delta)
    returns the fixed points for an adler dynamical system. 
    If the system is oscillatory and not excitable, it returns the point called ghost of less velocity. 
    
    Parameters
    ----------
    delta : float
        parameter of the adler equation (alpha/omega)
    
    Returns
    -------
    PFE : float
        stable fixed point - angular . On the  3rd cuadrant.
    PFI : float
        unstable fixed point - angular . On the 4th cuadrant
    '''
    if delta >= 1:
        res = np.arcsin(-1/delta)
        PFE = -res + np.pi  #np.sin(PFE)
        PFI = (2*np.pi + res)  #np.sin(PFI)
    else: 
        PFE = -np.pi/2 + 2*np.pi
        PFI = -np.pi/2 + 2*np.pi
    return(PFE,PFI)

def points_to_time(arr,dt,d):
    return np.array(arr)*dt*d #[i*dt*d for i in arr]


#%%
"""
Functions for the time-duration independent time series (no noise)

"""
def cociente(delta):
    return np.sqrt((delta + 1)/(delta - 1))

def f_(epsilon,delta):
    return (cociente(delta)-1/np.tan(1/2*(np.arccos(1/delta) + epsilon)))/ (cociente(delta)+1/np.tan(1/2*(np.arccos(1/delta) + epsilon)))
    
def compute_theoretical_omega(epsilon,delta):
    """ Esto multipicado por dos pi y dividido por T_0 te da el omega efectivo"""
    return -2 / (2*np.pi * np.sqrt(delta**2 - 1)) * np.log(f_(epsilon,delta))
    
def compute_theoretical_dt(omega,epsilon,delta):
    return -2 / (omega * np.sqrt(delta**2 - 1)) * np.log(f_(epsilon,delta))
    