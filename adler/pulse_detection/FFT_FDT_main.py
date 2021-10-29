'''
Data analysis quantifiers module

Author: Fiorella Fabris
Date  : 01/2020 

'''

import multiprocessing as mp
import numpy as np
import pandas as pd
from math import ceil
from functools import partial 
import os
from adler.plotting.plotting_main import check_file, save_data, download_data
from adler.pulse_detection.consecutive_main import consecutive_cumulative,consecutive_non_cumulative
#%%
################################################
#### Fast Fourier Transform computation
################################################  


def norm_fft_yf(data, Dt, max_freq = None):
    '''
    norm_fft_yf(data, Dt, max_freq = None)
    Computes a one-dimentional discrete Fourier Transform. 
    
    Computes the one-dimensional discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform (FFT)
    algorithm. Returns the FFT on a orthonomal basis so Parseval equality is satisfied. 
    
    Parameters
    -----------
    data : list of floats
        amplitude time serie.
    Dt : float
        sampling rate.
    max_freq: float, optional --- actually is not working
        maximum frequency of the FFT. 
        
    Returns
    --------
        yf : complex ndarray 
        The Fast Fourier Transfrom values, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_xf: for obtenining the frequencies involved on the FFT transform.

    '''
    
    N = data.shape[0]
    Nf = N // 2 #if max_freq is None else int(max_freq * T)
    yf =  np.fft.fft(data,norm='ortho')
    return (yf[:Nf])[1:][0::]

def norm_fft_xf(data, Dt, max_freq = None):
    '''
    norm_fft_xf(data, Dt, max_freq = None)
    Gives the frequencies involved on the one-dimentional discrete Fourier Transform. 

    
    Parameters
    -----------
    data : list of floats
        amplitude time serie.
    Dt : float
        sampling rate.
    max_freq: float, optional --- actually is not working
        maximum frequency of the FFT. 
        
    Returns
    --------
        xf : complex ndarray 
        The angular frequencies involved on the one-dimentional discrete Fourier Transfor, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_yf: for computing he one-dimentional discrete Fourier Transform. 

    '''
    N = data.shape[0]
    Nf = N // 2 #if max_freq is None else int(max_freq * T)
    xf = np.linspace(0.0, 0.5 / Dt, N // 2)
    return (xf[:Nf]*(2*np.pi))[1:][0::]


def norm_fft_statistics_aux_yf(data_folder,dt,d,ix_data):
    
    '''
    norm_fft_statistics_aux_yf(data_folder,dt,d,ix_data)
    Auxiliary function for computing the one-dimentional discrete Fourier Transform on a N repetition experiments.

    
    Parameters
    -----------
    data_folder : string
        path of where the time series are stored. 
    dt : float
        time resolution of the time series. 
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    ix_data: generator
        datasheet with the time series parameters and saving names.
                
    Returns
    --------
        out : list of complex ndarrays 
        For each time serie, a ndarray containing the Fast Fourier Transfrom values, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_statistics_aux_xf : for obtenining the frequencies involved on the FFT transform.
    norm_fft_yf : for computing he one-dimentional discrete Fourier Transform for an individual time series. 


    '''
    order = int(ix_data[1].order);number = int(ix_data[1].number)
    file_name =  str(number)+'_'+str(order)+'.pkl'
    
    if check_file(file_name,data_folder):
        return norm_fft_yf(np.cos(download_data(data_folder + file_name) ),dt * d) 
    else:
       return []
    
def norm_fft_statistics_aux_xf(data_folder,dt,d,ix_data):
    '''
    norm_fft_statistics_aux_xf(data_folder,dt,d,ix_data)
    Auxiliary function for obtaining the frequencies involved on the one-dimentional discrete Fourier Transform on a N repetition experiments.

    
    Parameters
    -----------
    data_folder : string
        path of where the time series are stored. 
    dt : float
        time resolution of the time series. 
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units)
    ix_data: generator
        datasheet with the time series parameters and saving names.
                
    Returns
    --------

        For each time serie, a ndarray containing the angular frequencies involved on the one-dimentional discrete Fourier Transform, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_statistics_aux_yf :  for computing the one-dimentional discrete Fourier Transforms. 
    norm_fft_xf: for obtenining the frequencies involved on the FFT transform.

    '''
    order = int(ix_data[1].order);number = int(ix_data[1].number)
    file_name =  str(number)+'_'+str(order)+'.pkl'
    
    if check_file(file_name,data_folder):
        return norm_fft_xf(np.cos(download_data(data_folder + file_name)),dt * d) 
    else:
        return []


def filter_nans(list_):
    ''' Filters the empty elements of a list of elements. 
    
    Parameters
    -----------
    list_: list
    list of arrays or lists
 
    Returns
    --------
    out: list
        same list as list_ without the empty elements.

 
    '''
    for arr in list_:
        aux = []
        if len(arr) == 0:
            pass
        else:
            aux.append(arr)
    return(aux)

def norm_fft_statistics(row,data_folder,dt,d):   
    '''
    norm_fft_statistics(row,data_folder,dt,d)
    Function for computing on parallel the one-dimentional discrete Fourier Transform on a N repetition experiments.
    
    Parameters
    -----------
    row: generator
    datasheet with the time series parameters and saving names.
    data_folder : string
        path of where the time series are stored. 
    dt : float
        time resolution of the time series. 
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
                
    Returns
    --------
        YF : list of floats (creo)
        The mean value of the squared absolute value of the Fourier Amplitudes. The mean is made over all the Fourier Ampolitudes of the time series.
        XF : list of floats (creo)
        The mean value of all the Fourier frequencies. The mean is made over all the Fourier frequencies related related to the time series.

 
        
    See Also
    ---------
    norm_fft_statistics_aux_xf : for obtenining the frequencies involved on the FFT transform.
    norm_fft_statistics_aux_yf :  for computing the one-dimentional discrete Fourier Transforms. 
    norm_fft_yf : for computing he one-dimentional discrete Fourier Transform for an individual time series. 
    norm_fft_xf: for obtenining the frequencies involved on the FFT transform of individual time series.
    
    Notes
    -----
    Computes the one-dimensional discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform (FFT)
    algorithm. 
    The Fourier amplitudes are on a orthonomal basis so Parseval equality is satisfied. 
    Only the positive frequencies are considered and the main frequency is not considered.
    The maximum frequency resolved is the nyquist frequency.
    The frequency resolution is 2*Dt / T, where T is the time series lenght. 
    This function works even if the simulation is not finished, or some time series files are missing or brocken. For these cases, returns empty YF and XF lists.

    '''
     
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    norm_fft_statistics_aux_yf_=partial(norm_fft_statistics_aux_yf,data_folder,dt,d)
    YF = pool.map(norm_fft_statistics_aux_yf_,row.iterrows())
    pool.close()
    pool.join()
    print('YF finished')
    
    pool2 = mp.Pool(processes= ceil(mp.cpu_count()))
    norm_fft_statistics_aux_xf_=partial(norm_fft_statistics_aux_xf,data_folder,dt,d)
    XF = pool2.map(norm_fft_statistics_aux_xf_,row.iterrows())
    pool2.close()
    pool2.join()
    print('XF finished')
    
    YF = filter_nans(YF); XF = filter_nans(XF)

    if len(YF)*len(XF)>0:
        print('len >>>')  
        return(np.mean(XF,axis=0),np.mean(np.abs(YF)**2,axis=0) )
    else:
        print('len <<<')    
        return ([],[])

    

def compute_fft_aux(save_path_name,data_folder,dt,d,i,D,row):
    ''' 
    Auxiliary function for running norm_fft_statistics and save the results
    
    '''
    
    omega =  row.omega.unique()[0]
    alpha = np.round(i/omega,4)
    print('running',alpha,D)      
    xf,yf = norm_fft_statistics(row,data_folder,dt,d)
    print('calculation ended')
    if len(xf)*len(yf) > 0:
        print('--- saving')
        save_data(xf,save_path_name+'fft_xf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
        save_data(yf,save_path_name+'fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
        print(alpha,D,'saving finished')
    else:
        print('not saving')
    print('calculation ended')
    return(0)

def compute_fft(description_file,data_folder,dt,d,save_path_name):    
    ''' 
    Function for running several experiments bla bla 
    Te devuelve numpy arrays
    '''
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);

    compute_fft_aux_= partial(compute_fft_aux,save_path_name,data_folder,dt,d)
    for (i,D), row in ref.groupby(['alpha','D']):
        compute_fft_aux_(i,D,row)
    print('the end')
    return (1)

#    '''
#    Power spectral density Fourier plotting function
#    '''

#%%

################################################
#### Module for computing the FPT
################################################  

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
    

def fpt(t,theta): 
    '''
    first_passage_time(t,theta)
    Computes the first passage time distribution for the (t,theta) time series.
    
    
    Parameters
    -----------    
        t : list
            time vector
        theta : list
            angular time series
    
    Returns
    --------
        out: list
            first passage time statistics on the time vector units
    
    See Also
    ---------
    fpt_statistics : for obtaining the first passage time distribution of a group of time series.
    
    Notes
    ------
    Here the first passage time is defined as the time needed for the angle variable to increase in one period 2pi.
            
    '''
    FPT_index = []; 
    init= theta[0]; FPT_index.append(0)
    for j,value in enumerate(theta):
        if value-init >= 2*np.pi:
            FPT_index.append(j)
            init = value
    return([t[f] - t[i] for i,f in zip(FPT_index[:-1],FPT_index[1:])]) 
    
def fpt_statistics(dt,T,d,data_folder,row):
    '''
    fpt_statistics(dt,T,d,data_folder,row)
    Computes the first passage time distribution of a group of time series.
    
    Parameters
    ----------
    dt : float
        time resolution of the time series. 
    T : float
        total time lenght of the time series.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    data_folder : string
        path of where the time series are stored.  
    row : generator
        datasheet with the time series parameters and saving names.


    Returns
    --------
        out: list
            First passage time distribution of the desired dataset.
            
    See Also
    ---------
    fpt : for obtaining the first passage time distribution of an individual time series.
    time : for obtaining the time vector of an individual time series.

    
    Notes
    ------
    Here the first passage time is defined as the time needed for the angle variable to increase in one period 2pi.
    '''

    list_data_folder = os.listdir(data_folder)
   
    ################################################
    #### FPT & checking
    ################################################  
    
    FPT = []
    for ix,data in row.iterrows():
        order = int(data.order); number = int(data.number)
        file_name =  str(number)+'_'+str(order)+'.pkl'
        if file_name in list_data_folder:
            print(file_name,' available')
            if check_file(file_name,data_folder):
               
                theta = download_data(data_folder + file_name) 
                t = time(dt,T,d)
                FPT.append(fpt(t,theta))
        else:
            print(file_name,' not available')         
    return([value for list_ in FPT for value in list_]) 
    



def compute_FPT_aux(save_path_name,data_folder,dt,d,T,tuple_):
    '''
    auxiliary function for computing on parallel the calculation of the first passage time for different experiments of N repetitions experiments each.
    
    See also
    ---------
    compute_FPT
    fpt_statistics
    
    *** FALTA COMPLETAR ESTA DESCRIPCION***
    '''
    (i,D),row = tuple_[0],tuple_[1]
    omega =  row.omega.unique()[0]
    alpha = np.round(i/omega,4)      

    FPT = fpt_statistics(dt,T,d,data_folder,row)
    save_data(FPT,save_path_name+'FPT_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
    print(alpha,D)
    return(0)

def compute_FPT(description_file,save_path_name,dt,T,d,data_folder):
    '''
    compute_FPT(description_file,save_path_name,dt,T,d,data_folder)
    Function for computing on parallel the calculation of the first passage time for different experiments of N repetitions experiments each.
    
    Parameters
    -----------
    description_file : string
        path of the datasheet with the time series parameters and saving names.
    save_path_name : string
        path in where the first passage time statistics is going to be saved.
    dt : float
        time resolution of the time series. 
    T : float
        total time lenght of the time series.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    data_folder : string
        path of where the time series are stored.  

                
    Returns
    --------
        Saves the first passage time of an experiment on a pickle on the same file.
        
    
    See also
    --------
    compute_FPT_aux
    
    *** FALTA COMPLETAR ESTA DESCRIPCION***
    '''
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    tuple_ = ref.groupby(['alpha','D'])

    compute_FPT_aux_ = partial(compute_FPT_aux,save_path_name,data_folder,dt,d,T)
    pool.map(compute_FPT_aux_,tuple_)
    pool.close()
    pool.join()
    return (1)


