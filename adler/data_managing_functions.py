import pickle
import os
import numpy as np
from math import ceil


def download_data(filename):
    infile =  open(filename,'rb')
    results = pickle.load(infile)
    infile.close()
    return(results)

def save_data(data, filename):
    outfile= open(filename,'wb')
    pickle.dump(data,outfile,protocol=-1)
    outfile.close()

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