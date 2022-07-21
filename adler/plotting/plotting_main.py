'''
Auxiliary plotting functions module

Author: Fiorella Fabris
Date  : 01/2020 

'''

import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from adler.data_managing_functions import download_data,check_file,time,points_to_time
from math import ceil
from functools import partial
import multiprocessing as mp


def silent_ax(ax):
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

def set_scale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([]))

    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.tick_params(axis="y",direction="out", width=1,labelsize=8)
    ax.tick_params(axis="x",direction="out", width=1,labelsize=8)


#%%for plotting quantifiers histograms and 2d plots


def download_quantifiers(row_,data_folder,T,dt,d):
    ''' para pulse rate parto la serie temporal en 200 veces'''
    DT = []; IPI = []; joint_duration = []; dm = []; pulse_rate = []
    for (order,row) in row_.groupby(['order']):
    #para cada ensayo con todos los mismos parámetros
        number = int(row.number.values[0])
        file_name =  str(number)+'_'+str(order)+'.pkl'

        if (check_file('dt_'+file_name,data_folder)):        
        
            DT = DT + download_data(data_folder+'dt_'+file_name)
            IPI = IPI + download_data(data_folder+'IPI_'+file_name)
            dm = dm + download_data(data_folder+'dm_'+file_name)
            joint_duration = joint_duration + download_data(data_folder+'joint_duration_'+file_name)
            
            pulse_rate = pulse_rate + [len(download_data(data_folder+'max_'+file_name)) / T]
            #pulse_rate + download_data(data_folder+'pr_ns_'+file_name)
        else:
            pass
    return(points_to_time(DT,dt,d),points_to_time(IPI,dt,d),points_to_time(joint_duration,dt,d),points_to_time(dm,dt,d),pulse_rate)



#%%


#for plotting 2d plots
def create_df(ref,data_folder,dt,T,d):
    ''' creates dataframes where 2d alpha are indexes and D are columns'''
    omega = ref.omega.unique()[0]
    Cols = ref.D.unique()
    Rows = ref.alpha.unique()
    mean_dt_matrix,mean_ipi_matrix,mean_fpt_matrix,mean_pulses_matrix = [],[],[],[]
    
    for alpha,row_ in ref.groupby(['alpha']):
        aux_dt,aux_ipi,aux_fpt,aux_pulses = [],[],[],[]    
        
        for D,col_ in row_.groupby(['D']):
            DT,IPI,_,_ ,pulse_rate = download_quantifiers(col_,data_folder,T,dt,d)
           
            aux_dt.append(np.median(DT));aux_ipi.append(np.median(IPI))
            aux_pulses.append(np.median(pulse_rate))
            
             
            fpt_file_name = 'FPT_'+str(omega)+'_'+str(np.round(alpha/omega,4) )+'_'+str(D)+'.pkl'
            if check_file(fpt_file_name,data_folder) :
                FPT = download_data(data_folder+fpt_file_name)
            else:
                FPT = []
            aux_fpt.append(np.median(FPT))
            
            #activity,_,_ = load_activity(col_,data_folder,dt,T,d); 
           
            #assert ((np.mean(activity) <= 100) or (len(activity) == 0)),(alpha,D,activity)
           # if len(activity) > 0: aux_act.append(np.mean(activity))
           # else:aux_act.append(0)
        mean_dt_matrix.append(aux_dt);mean_ipi_matrix.append(aux_ipi),mean_fpt_matrix.append(aux_fpt),mean_pulses_matrix.append(aux_pulses)
    return (pd.DataFrame(mean_dt_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_ipi_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_fpt_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_pulses_matrix,columns = Cols,index = Rows))


def load_activity(row_,data_folder,dt,T,d):
    '''
    load_activity(row_,dt,T,d). Calcula la activity para cada par D, alpha.
    
    row_: DataFrame
        variable que viene del excel. Es el grupo de filas del excel de descripción,
        con un sólo D y alpha, pero con todos los trials. 
    '''    

    activity = []; n_cell = 0    
    for (order,row) in row_.groupby(['order']):
        
        number      = int(row.number)
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        
        if (check_file('dt_'+file_name,data_folder)):        
            
            duration_cell   = download_data(data_folder+'dt_'+file_name)     
            n_cell = n_cell + 1
            activity = activity + [sum(duration_cell) / len(time(dt,T,d)) *  100]
            
    activity = np.sort(activity)[::-1] #orden descendente
    silent = np.ones(len(activity)) * 100 - activity

    return activity,silent,n_cell

def load_activity_dist(row_,data_folder,dt,T,d):
    '''
    load_activity(row_,dt,T,d). Calcula la activity para cada D y alpha una dist.
    
    row_: DataFrame
        variable que viene del excel. Es el grupo de filas del excel de descripción,
        con un sólo D, pero con todos los trials (number y order diferentes). 
    '''    

    activity = []; n_cell = 0    
    for (number,order),row in row_.groupby(['number','order']):
        
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        
        if (check_file('dt_'+file_name,data_folder)):        
            
            duration_cell   = download_data(data_folder+'dt_'+file_name)     
            n_cell = n_cell + 1
            activity = activity + [sum(duration_cell) / len(time(dt,T,d)) *  100]
        elif  (check_file(file_name,data_folder)):
            activity = activity + [0]
            n_cell = n_cell + 1 
        else:
            print("no TS for file_name")
            
    activity = np.sort(activity)[::-1] #orden descendente
    silent = np.ones(len(activity)) * 100 - activity

    return activity,silent,n_cell
    
    
#%%
def mask_arr(beg,end,arr):
    return (np.array([i-beg for i in arr if (i  < end and i >= beg)]))
#%%
def tune_plot(ax,x_label,y_label,xlim,xscale,ylim,yscale,axlabel_fs = 10,ticklabel_fs = 8):
    ''' Esto es una funcion auxiliar para plotear. La idea es escribir una vez
    el tamano de las labels, las escalas y esas cosas'''
    
    #ax.set_xlabel(x_label, fontsize=axlabel_fs);
    #ax.set_ylabel(y_label, fontsize=axlabel_fs);
    ax.set_xlim(xlim)
    ax.set_xticklabels([np.round(item/xscale,3) for item in ax.get_xticks()],fontsize=ticklabel_fs)
    ax.set_ylim(ylim)
    ax.set_yticklabels([np.round(item*yscale,3) for item in ax.get_yticks()],fontsize=ticklabel_fs)

def compute_st_values(ax,samples,bins,scale,fs = 6):
    ''' Esto esuna función auxiliar paraplotear.
    computa la moda de samples usando bins, y los cuartiles de samples usando samples. 
    samples es el dataset que queres estudiar, y bins es el histograma que estas graficando.
    scale es por loque tengo que dividir samples para llegar a minutos. fs es el fontsize del texto. '''           
    
    if len(samples) <= 0:
        pass
    else:
        mode = (bins[1][np.argmax(bins[0])] + bins[1][np.argmax(bins[0])+1])/2 ; 
        mode = 'mode: '+str(np.round(mode/scale,2))+' min \n'
        ax.text(1, 0.75, mode, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 
        ax.text(1, 0.9,r'$Q$: '+str(np.round(np.quantile(samples,0.25)/scale,2))+' ; '+str(np.round(np.quantile(samples,0.5)/scale,2))+' ; '
                +str(np.round(np.quantile(samples,0.75)/scale,2)) , ha='right', va='center', transform=ax.transAxes, fontsize=fs)
        ax.text(1, 0.7, 'total data: ' + str(len(samples)), ha='right', va='center', transform=ax.transAxes, fontsize=fs) 
    
    #chequeamos que no haya numeros raros
        if True:
            neg = sum([1 for i in samples if i < bins[1][0]]) 
            pos =  sum([1 for i in samples if i > bins[1][-1]])
            text_borders = '< lower bound :'+str(neg) + '\n > upper bound :'+str(pos)
            print(text_borders)
            #ax.text(1, 0.75, text_borders, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 