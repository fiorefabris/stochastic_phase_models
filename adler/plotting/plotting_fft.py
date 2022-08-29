import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import ceil
from functools import partial 
import seaborn as sns

from adler.data_managing_functions import download_data,check_file


#%%

###############################################################################
### FFT plotting module
###############################################################################
   
###############################################################################
### Function for paralelizing the FFT plotting 
###############################################################################

def plot_fft(description_file,data_folder,dt,d,save_path_name):

    
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_fft_alpha_= partial(plot_fft_alpha,save_path_name,data_folder,dt,d)
    pool.map(plot_fft_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)
   
    
        
def moving_average(data, window_size):
    ''' Function that returns rolling averaged window smoothed data
    
    Parameters
    -----------
        data: numpy array
            one-dimentional time series that we want to smooth.
        window size: int
            the size of the rolling window.
    
    Returns
    -------
        out: numpy array
            Smoothed data.
            
    Notes
    ------
    How is it working at the beggining and end of the TS?
            
    
    '''
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(data, window, 'same')
###############################################################################
### Function for plotting the FFT time series 
###############################################################################


def plot_fft_alpha(save_path_name,data_folder,dt,d,tuple_):
    '''
    plot_fft_alpha(save_path_name,data_folder,dt,d,tuple_)
    Power spectral density Fourier plotting function.
    
    For each alpha, plots, for different values of D on a column, the mean value of the square module of the FFT.
    The x-axis are angular frequencies. 
    The red line corresponds to the smoothed data using a rolling averaged window of size 100. 
    The gray dashed line corresponds to the frequency of the deterministic dynamical system.
    
    Parameters
    ----------
    save_path_name : string
        path for the figure to be saved.
    data_folder : string
        path where the fourier coefficients and frequencies are stored. 
    dt : float
        time resolution of the time series. 
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    tuple_ : tuple generator
        the alpha coefficient, and the datasheet with the experiments parameters. 
    
    '''
    i,rows = tuple_[0],tuple_[1]
    
###############################################################################
### Parameters
###############################################################################
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    
###############################################################################
### Plotting parameters
###############################################################################    
    Cols = 1; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    ylim = [0,5e2] ; xlim = [0,4*omega]
    xticks = np.arange(0,5,1)*omega ;
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();
    text = r'$\omega = \frac{2\pi}{7 min}$' + ' ~ ' + r'$ \alpha = $'+str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    axs[0].text(0.2,1.3, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=12)

    
    for k,(D,row) in enumerate(rows.groupby(['D'])):
        ax = axs[k]; ax.grid(False);
        
        if check_file('fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl',save_path_name):
            print(alpha,np.round(D,5))
            
            xf = download_data(save_path_name+'fft_xf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
            yf = download_data(save_path_name+'fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
            
            #Plotting 
            print(yf[0])
            ax.plot(xf,yf, linewidth=0.5)
           # ax.set_yscale('log')
            #if len(xf) > 10: ax.plot(xf,moving_average(yf,100), linewidth=0.5,color ='red')
            if alpha <= 1: ax.axvline(np.sqrt(omega**2-i**2),ls = '--', color = 'gray',linewidth=0.5)
            
        else:
            print(alpha,D,'empty')


        ax.set_xlim(xlim)
        #ax.set_xscale('log')        
        #ax.set_yscale('log');
        ax.set_ylim(ylim);
        ax.set_xticks(xticks);
        ax.set_yticks(ylim)
        ax.tick_params(labelsize=10)

        
        text = 'D = ' + str(np.round(D,5))
        ax.text(1.10, 0.7, text , ha='center', va='center', transform=ax.transAxes, fontsize=12)

        
    ax.set_ylabel('power spectrum', fontsize=15);
    ax.set_xlabel('angular frequency', fontsize=15)
    ax.xaxis.set_label_coords(0.5, -0.30);
    ax.yaxis.set_label_coords(-0.10, 1)
    ax.set_xticklabels([str(int(x/omega)) + r'$\omega$' for x in xticks], fontsize=10); 
        
    #for m in range(Cols*Rows - k):
    #    fig.delaxes(axs[-m-1])

    plt.savefig(save_path_name + 'fft_alpha_'+str(alpha)+'.pdf', format='pdf')

    return(0)
#%%

#%%
from scipy.signal import find_peaks,peak_widths

#x_signal = np.linspace(0,100,200)
#signal = np.sin(x_signal) + 1
def get_quality_factor_OLD(x_signal, signal):
    #peaks, h = find_peaks(signal,height=0) #encuentra todos los picos, y el alto lo mide desde x = 0
    w0_ix = np.argmax(signal)
    
   # if len(peaks) > 0:
    #w0_ix = peaks[np.argmax(signal[peaks])] #el pico más alto
    S_w0 = signal[w0_ix] #la potencia del pico más alto
    widthx,widthy,left_ips,rigth_ips = peak_widths(signal, [w0_ix], rel_height=(1-1/np.sqrt(np.e))) # ancho (en índice), alto, donde empieza y donde termina (son puntos interpolados)
    
    
    w0 = x_signal[w0_ix] #frecuencia fundamental 
    print(w0,S_w0,widthy[0])
    beta = w0 * S_w0 / widthy[0]
    return(beta)
    #else: 
     #   return None

#%%
     

def get_quality_factor(x_signal, signal):

    w0_ix = np.argmax(signal) # el índice de la fundamental 
    S_w0 = signal[w0_ix] #la potencia del pico más alto
    
    
    rel_height=S_w0/np.sqrt(np.e)
    for ix,s in enumerate(signal[w0_ix:]):
        if s > rel_height:
            pass
        else:
            break
    
    widthx = 2*ix * np.diff(x_signal)[1] #ancho que tiene la cosa cuando la altura es rel_height
    
    
    w0 = x_signal[w0_ix] #frecuencia fundamental 
    print(w0,S_w0,widthx)
    beta = w0 * S_w0 / widthx
    return(w0, S_w0, ix+w0_ix,rel_height, beta)


def plot_fft_all(description_file,data_folder,dt,d,save_path_name):
    #para un omega, plotea los alphas uno encima del otro
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_fft_alpha_= partial(plot_fft_alpha_all,save_path_name,data_folder,dt,d)
    pool.map(plot_fft_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)

def plot_fft_alpha_all(save_path_name,data_folder,dt,d,tuple_):
    '''
    plot_fft_alpha(save_path_name,data_folder,dt,d,tuple_)
    Power spectral density Fourier plotting function.
    
    For each alpha, plots, for different values of D on a column, the mean value of the square module of the FFT.
    The x-axis are angular frequencies. 
    The red line corresponds to the smoothed data using a rolling averaged window of size 100. 
    The gray dashed line corresponds to the frequency of the deterministic dynamical system.
    
    Parameters
    ----------
    save_path_name : string
        path for the figure to be saved.
    data_folder : string
        path where the fourier coefficients and frequencies are stored. 
    dt : float
        time resolution of the time series. 
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    tuple_ : tuple generator
        the alpha coefficient, and the datasheet with the experiments parameters. 
    
    '''
    i,rows = tuple_[0],tuple_[1]
    
###############################################################################
### Parameters
###############################################################################
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    
###############################################################################
### Plotting parameters
###############################################################################    
    Cols = 1; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    colors =  sns.color_palette(sns.color_palette("viridis",Rows*1))
    colors =  colors[::1]
    ylim = [0,5e2] ; xlim = [0,4*omega]
    xticks = np.arange(0,5,1)*omega ;
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();
    text = r'$\omega = \frac{2\pi}{7 min}$' + ' ~ ' + r'$ \alpha = $'+str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    axs[0].text(0.2,1.05, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=12)

    BETA = []; D_ = []
    for k,(D,row) in enumerate(rows.groupby(['D'])):
        ax = axs[1]; ax.grid(False);
        assert check_file('fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl',data_folder)
            
        xf = download_data(data_folder+'fft_xf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
        yf = download_data(data_folder+'fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')

        w0, S_w0,wdt_ix,S_wdt, beta = get_quality_factor(xf, yf)
        if beta is not None:
            BETA.append(beta)
            D_.append(D)

            
        #Plotting 
        if D == 0: 
            ax.plot(xf,yf,linewidth=1,color =colors[k],alpha = 0.6,label = str(D))
        else: 
            ax.plot(xf,yf, linewidth=1,color =colors[k],alpha = 1,label = str(D))
            ax.plot(w0, S_w0, 'o', color = 'black',alpha = 1)
            ax.plot(xf[wdt_ix],S_wdt, 'o', color = 'black',alpha = 1)
            #if alpha <= 1: ax.axvline(np.sqrt(omega**2-i**2),ls = '--', color = 'gray',linewidth=0.5)
            


        ax.set_xlim(xlim)
        #ax.set_xscale('log')        
        #ax.set_yscale('log');
        ax.set_ylim(ylim);
    ax.set_xticks(xticks);
    ax.set_yticks(ylim)
    ax.tick_params(labelsize=10)
    axs[0].plot(D_[1:],BETA[1:],'-o') ; axs[0].set_xscale('log')

        
        #text = 'D = ' + str(np.round(D,5))
        #ax.text(1.10, 0.7, text , ha='center', va='center', transform=ax.transAxes, fontsize=12)

        
    ax.set_ylabel('power spectrum', fontsize=15);
    ax.set_xlabel('angular frequency', fontsize=15)
    ax.xaxis.set_label_coords(0.5, -0.1);
    ax.yaxis.set_label_coords(-0.10,0.5)
    ax.set_xticklabels([str(int(x/omega)) + r'$\omega$' for x in xticks], fontsize=10); 
    ax.legend(fontsize=7, ncol=1,framealpha=0, fancybox=True)

    #for m in range(Cols*Rows - k):
    #    fig.delaxes(axs[-m-1])

    plt.savefig(save_path_name + 'all_fft_alpha_'+str(alpha)+'.pdf', format='pdf')

    return(0)
    
