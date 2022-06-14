'''
Auxiliary plotting functions module

Author: Fiorella Fabris
Date  : 01/2020 

'''

import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from adler.data_managing_functions import download_data,check_file,time



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
def points_to_time(arr,dt,d):
    return np.array(arr)*dt*d #[i*dt*d for i in arr]

def download_quantifiers(row_,data_folder,dt,d):
    DT = []; IPI = []; joint_duration = []; dm = []
    for (order,row) in row_.groupby(['order']):
    #para cada ensayo
        number = int(row.number.values[0])
        file_name =  str(number)+'_'+str(order)+'.pkl'

        if (check_file('dt_'+file_name,data_folder)):        
        
            DT = DT + download_data(data_folder+'dt_'+file_name)
            IPI = IPI + download_data(data_folder+'IPI_'+file_name)
            dm = dm + download_data(data_folder+'dm_'+file_name)
            joint_duration = joint_duration + download_data(data_folder+'joint_duration_'+file_name)
        else:
            pass
    return(points_to_time(DT,dt,d),points_to_time(IPI,dt,d),points_to_time(joint_duration,dt,d),points_to_time(dm,dt,d))
    
        #for plotting 2d plots
def create_df(ref,data_folder,dt,T,d):
    ''' creates dataframes where 2d alpha are indexes and D are columns'''
    omega = ref.omega.unique()[0]
    Cols = ref.D.unique()
    Rows = ref.alpha.unique()
    mean_dt_matrix,mean_ipi_matrix,mean_fpt_matrix,mean_act_matrix = [],[],[],[]
    
    for alpha,row_ in ref.groupby(['alpha']):
        aux_dt,aux_ipi,aux_fpt,aux_act = [],[],[],[]    
        
        for D,col_ in row_.groupby(['D']):
            DT,IPI,_,_ = download_quantifiers(col_,data_folder,dt,d)
            aux_dt.append(np.mean(DT));aux_ipi.append(np.mean(IPI))
            
            fpt_file_name = 'FPT_'+str(omega)+'_'+str(np.round(alpha/omega,4) )+'_'+str(D)+'.pkl'
            if check_file(fpt_file_name,data_folder) :
                FPT = download_data(data_folder+fpt_file_name)
            else:
                FPT = []
            aux_fpt.append(np.mean(FPT))
            
            activity,_,_ = load_activity(col_,data_folder,dt,T,d); 
            assert ((np.mean(activity) <= 100) or (len(activity) == 0)),(alpha,D,activity)
            if len(activity) > 0: aux_act.append(np.mean(activity))
            else:aux_act.append(0)
        mean_dt_matrix.append(aux_dt);mean_ipi_matrix.append(aux_ipi),mean_fpt_matrix.append(aux_fpt),mean_act_matrix.append(aux_act)
    return (pd.DataFrame(mean_dt_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_ipi_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_fpt_matrix,columns = Cols,index = Rows),pd.DataFrame(mean_act_matrix,columns = Cols,index = Rows))


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
    
#%%
def mask_arr(beg,end,arr):
    return (np.array([i-beg for i in arr if (i  < end and i >= beg)]))
