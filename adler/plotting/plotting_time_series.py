import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import ceil
from functools import partial 
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from adler.data_managing_functions import download_data,check_file,time
from adler.plotting.plotting_main import set_scale

###############################################################################
### Time series plotting module
###############################################################################

def plot_time_series(dt,T,d,N,delta,description_file,data_folder,save_path_name):
    '''
    plot_time_series(description_file,data_folder,save_path_name)
    Auxiliary funtion for parallelizing the time Series plotting function.
    
    Parameters
    ----------
    dt : float
        time resolution of the time series. 
    T : float
        total time lenght of the time series. Is the interval you want to plot, starting from zero.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    N : int 
        number of repetitions of the time series with same parameters.
    delta : integer
        time series plotting resolution (i.e. how often do you want to plot your time series in datapoint units)
    description_file : string
        path of the datasheet with the time series parameters and saving names.
    data_folder : string
        path of where the time series are stored.  
    save_path_name : string
        path in where the first passage time statistics is going to be saved.

    Returns
    --------
    Time series pdf files for each alpha value. 
    In each column, there is one simulation corresponding to a D value. 
    
    See also
    ---------
    plot_time_series_alpha : for plotting function for each file.

    '''
    
######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);

######## Parameters & files
    #dt = 0.0001 ;T = 10000; d=10;  delta = 1000 

######## Paralelization & running
    aux_iterator = ref.groupby(['alpha'])

    pool = mp.Pool(processes= mp.cpu_count())
    plot_time_series_alpha_ = partial(plot_time_series_alpha, data_folder,save_path_name,dt,T,d,N,delta)
    pool.map(plot_time_series_alpha_,aux_iterator )
    pool.close()
    pool.join()
    return (0)


def plot_time_series_alpha(data_folder,save_path_name,dt,T,d,N,delta,tuple_):
    '''
    plot_time_series_alpha(data_folder,save_path_name,dt,T,d,delta,tuple_)
    Auxiliary funtion for parallelizong the time Series plotting function.
    
    Parameters
    ----------
    data_folder : string
        path of where the time series are stored.  
    save_path_name : string
        path in where the first passage time statistics is going to be saved.
    dt : float
        time resolution of the time series. 
    T : float
        total time lenght of the time series. Is the interval you want to plot, starting from zero.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    N : int 
        number of repetitions of the time series with same parameters.
    delta : integer
        time series plotting resolution (i.e. how often do you want to plot your time series in datapoint units)

    tuple : iterator
        path of the datasheet with the time series parameters and saving names.


    Returns
    --------
    Time series pdf files for one alpha value. 
    In each column, there is one simulation corresponding to a D value. 
    
    See also
    ---------
    plot_time_series :  for auxiliary funtion for parallelizing the time Series plotting function.

    '''

    i,rows = tuple_[0],tuple_[1]
###############################################################################
### Parameters
###############################################################################

    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    print('alpha = ', alpha)

###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5,T+5] ; ylim = [-1.1,1.1] ;         
    Cols = 1 #TS; 
    Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    colors =  sns.color_palette(sns.color_palette("viridis",Rows*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.1)
    text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
    
    for k,(ix,row_) in  enumerate(rows.groupby(['D'])):
        
        counter = 0
        row = row_.iloc[counter]; order = int(row.order); number = int(row.number)
        file_name =  str(number)+'_'+str(order)+'.pkl'

        #Find the first file that is valid
        while (check_file(file_name,data_folder) == False) and (counter < N-1):
            counter = counter + 1
            row = row_.iloc[counter]; order = int(row.order); number = int(row.number)
            file_name =  str(number)+'_'+str(order)+'.pkl'
            
        print(counter,00)
        print(row.D,'plotting ',order,number);D = row.D
        ax = axs[k]; ax.grid(False);

        ################################################
        #### download data
        ################################################
        if check_file(file_name,data_folder):
            theta = download_data(data_folder + file_name) 
            t = time(dt,T,d)
            end = len(t)
            ax.plot(t[:end:delta],np.sin(theta)[:end:delta],linewidth=2,color=colors[k])
        ################################################
        #### Plotting
        ################################################
        text = 'D = ' + str(np.round(D,5))
        ax.text(1.05, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        if k == len(axs)-Cols:
            ax.set_ylabel(r'$\cos(\theta)$', fontsize=30);
            ax.set_xlabel('time', fontsize=30)
            ax.xaxis.set_label_coords(0.5, -0.1);
            ax.yaxis.set_label_coords(-0.05, 0.5)

        set_scale(ax, [0,T], [-1,1])
        ax.set_xticklabels([0,T])
        ax.tick_params(labelsize=20)


#    for m in range((Cols*Rows - (k+1))):
#        fig.delaxes(axs[-m-1])


    plt.savefig(save_path_name + 'time_series_alpha_'+str(alpha)+'.pdf', format='pdf')
    return(0)
    

#%%
def plot_time_series_square(dt,beg,T,d,N,Delta,description_file,data_folder,save_path_name):
    '''
    plot_time_series_alpha(data_folder,save_path_name,dt,T,d,delta,tuple_)
    Auxiliary funtion for parallelizong the time Series plotting function.
    
    Parameters
    ----------
    data_folder : string
        path of where the time series are stored.  
    save_path_name : string
        path in where the first passage time statistics is going to be saved.
    dt : float
        time resolution of the time series. 
    beg: principio de la TS. en tiempo
    T : float
        total time lenght of the time series. Is the interval you want to plot, starting from zero.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    N : int 
        number of repetitions of the time series with same parameters.
    Delta : integer
        time series plotting resolution (i.e. how often do you want to plot your time series in datapoint units)

    tuple : iterator
        path of the datasheet with the time series parameters and saving names.


    Returns
    --------
    Time series pdf files for one alpha value. 
    In each column, there is one simulation corresponding to a D value. 
    
    See also
    ---------
    plot_time_series :  for auxiliary funtion for parallelizing the time Series plotting function.

    '''

######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
    
    if 'omega' in ref_.keys() : iterator_T = ref_.groupby(['omega'])
    if 'T0' in ref_.keys() : iterator_T = ref_.groupby(['T0'])
    for T0,ref in iterator_T:
    ###############################################################################
    ### Plotting parameters
    ###############################################################################    
        xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
        Cols = len(ref.groupby(['D'])) ;
        if 'alpha' in ref.keys() : Rows = len(ref.groupby(['alpha'])) ; 
        if 'delta' in ref.keys() : Rows = len(ref.groupby(['delta'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
    ###############################################################################
    ### Figure
    ###############################################################################    
    
        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
        #text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' 
        #axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
        
        for col,(D,col_) in  enumerate(ref.groupby(['D'])):
            if 'alpha' in ref.keys() : iterator = col_.groupby(['alpha'])
            if 'delta' in ref.keys() : iterator = col_.groupby(['delta'])
            for row, (alpha,row_)  in  enumerate(iterator):
                if 'alpha' in ref.keys() : delta= np.round(alpha/col_.omega.unique()[0],4)  
                if 'delta' in ref.keys() : delta = alpha
                order = int(row_.order.unique()[0]); number = int(row_.number.unique()[0])
                file_name =  str(number)+'_'+str(order)+'.pkl'
                ax = axs[row,col]; ax.grid(False);
                
                ################################################
                #### download data
                ################################################
                if check_file(file_name,data_folder):            
                    
                    theta = download_data(data_folder + file_name) 
                    t = time(dt,T+beg,d)
                    end = len(t)
                    beg_ = int(beg/(dt*d))
                    ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=colors[col])
                
                ###############################################
                #### Plotting
                ################################################
                if row == 0:
                    text = 'D = ' + str(np.round(D,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'delta = ' + str(delta)
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
    
                ax.set_ylim(ylim);
                ax.set_xlim(xlim)
                
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                    ax.set_xlabel('time', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.1);
                    ax.yaxis.set_label_coords(-0.05, 0.5)
                
                set_scale(ax,[beg,T], [0,2])
                ax.set_xticklabels([beg,T])
                ax.set_yticklabels([0,2])
                ax.tick_params(labelsize=20)
    
    
    #    for m in range((Cols*Rows - (k+1))):
    #        fig.delaxes(axs[-m-1])
    
    
        plt.savefig(save_path_name + 'time_series_square_'+str(T0)+'.pdf', format='pdf')
    return(0)

#%%
def plot_time_series_square_ou(dt,beg,T,d,N,Delta,description_file,data_folder,save_path_name):
    '''
    plot_time_series_alpha(data_folder,save_path_name,dt,T,d,delta,tuple_)
    Auxiliary funtion for parallelizong the time Series plotting function.
    
    Parameters
    ----------
    data_folder : string
        path of where the time series are stored.  
    save_path_name : string
        path in where the first passage time statistics is going to be saved.
    dt : float
        time resolution of the time series. 
    beg: principio de la TS. en tiempo
    T : float
        total time lenght of the time series. Is the interval you want to plot, starting from zero.
    d : integer
        the decimation factor of the experiment (i.e. how often do you save your time series in dt units).
    N : int 
        number of repetitions of the time series with same parameters.
    Delta : integer
        time series plotting resolution (i.e. how often do you want to plot your time series in datapoint units)

    tuple : iterator
        path of the datasheet with the time series parameters and saving names.


    Returns
    --------
    Time series pdf files for one alpha value. 
    In each column, there is one simulation corresponding to a D value. 
    
    See also
    ---------
    plot_time_series :  for auxiliary funtion for parallelizing the time Series plotting function.

    '''

######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
    omega = ref_.omega.unique()
    ###############################################################################
    ### Plotting parameters
    ###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
    ###############################################################################
    ### Figure
    ###############################################################################    
    
        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):
                delta= np.round(alpha/col_.omega.unique()[0],4)  
                order = int(row_.order.values[0]); number = int(row_.number.values[0])
                file_name =  str(number)+'_'+str(order)+'.pkl'
                ax = axs[row,col]; ax.grid(False);
                
                ################################################
                #### download data
                ################################################
                if check_file(file_name,data_folder):            
                    
                    theta = download_data(data_folder + file_name)
                    delta_arr =  download_data(data_folder + 'alpha_' + file_name) / omega
                    t = time(dt,T+beg,d)
                   
                    end = len(t) #end is ix
                    beg_ = int(beg/(dt*d)) # beg_is ix
                   
                    ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=colors[col])
                    divider = make_axes_locatable(ax)
                    ax_alpha = divider.append_axes("bottom", size="100%", pad=0.1, sharex=ax)
                    ax_alpha.plot(t[beg_:end:Delta],delta_arr[beg_:end:Delta],linewidth=2,color='gray')
                    ax_alpha.plot(t[beg_:end:Delta],np.ones(len(t[beg_:end:Delta])),':',color='gray')
                ###############################################
                #### Plotting
                ################################################
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
    
                ax.set_ylim(ylim);
                ax.set_xlim(xlim)
                ax_alpha.set_xlim(xlim)
                ax_alpha.set_ylim([0.5,1.5])
                
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                    ax.yaxis.set_label_coords(-0.05, 0.5)
                    ax_alpha.set_xlabel('time', fontsize=30)
                    ax_alpha.xaxis.set_label_coords(0.5, -0.1);
                   
                
                set_scale(ax,[beg,T], [0,2])
                ax.set_xticklabels([None,None])
                ax.set_yticklabels([0,2])
                ax.tick_params(labelsize=20)
                set_scale(ax_alpha,[beg,T], [-1,2])
                ax_alpha.set_xticklabels([beg,T])
                ax_alpha.set_yticklabels([-1,2])
                ax_alpha.tick_params(labelsize=20)

                print(alpha)
    
    
        
        plt.savefig(save_path_name + 'time_series_square_ou_'+str(delta)+'.pdf', format='pdf')
    return(0)

#%%
def plot_time_series_square_dist(dt,beg,T,d,N,Delta,description_file,data_folder,save_path_name):
    '''
    
    '''

######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    plot_n = len(ref.groupby(['D']))
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.05,2.05] ;         
    colors =  sns.color_palette(sns.color_palette("viridis",plot_n*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    for color_ix,(D,ref_) in enumerate(ref.groupby(['D'])):

        fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
        axs = axs.ravel(); 
        
        for ix, (alpha,df_)  in  enumerate(ref_.groupby(['alpha'])):
            if ix < len(axs):
                ax = axs[ix]
                delta= np.round(alpha/ref_.omega.unique()[0],4)  
                print(delta,D)
            
                order = int(df_.order); number = int(df_.number)
                file_name =  str(number)+'_'+str(order)+'.pkl'
                
            
                ################################################
                #### download data
                ################################################
                if check_file(file_name,data_folder):            
                    
                    theta = download_data(data_folder + file_name) 
                    t = time(dt,T,d)
                    end = len(t)
                    beg_ = int(beg/(dt*d))
                    assert len(t) == len(theta), (len(t),len(theta))
                    print(len(t[beg_:end:Delta]),len(1+np.sin(theta)[beg_:end:Delta]))
                    ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=colors[color_ix])
                    
                ax.set_ylim(ylim);
                ax.set_xlim(xlim)
                
                ###############################################
                #### Plotting
                ################################################
                text = str(np.round(delta,5))
                ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                
                if (ix == len(axs)-6):
                    ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                    ax.set_xlabel('time (min)', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.1);
                    ax.yaxis.set_label_coords(-0.05, 0.5)
                
                set_scale(ax,[beg,T], [0,2])
                ax.set_xticklabels([beg,T])
                ax.set_yticklabels([0,2])
                ax.tick_params(labelsize=20)
    

    
        plt.savefig(save_path_name + 'time_series_dist_D'+ str(D)+'.pdf', format='pdf')
    return(0)

#%%