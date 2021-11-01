import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.plotting.plotting_main import download_data,check_file

#%%
def download_dt_alpha(row,data_folder):
    '''
    te arma una lista con todos los dt, en donde cada lugar de la lista (que es una lista) es un alpha distinto
    si un dt no esta en la carpeta, devuelve un cero
    
    Esta preparado para solo tomar un numero de orden, es decir una serie temporal por parametros.
    '''
    dt = []; 
    for (number,row_) in row.groupby(['number']):
        order = row_.order.values[0]
        file_name =  str(number)+'_'+str(order)+'.pkl'
        if (check_file('dt_xf_'+file_name,data_folder)):        
            dt.append(download_data(data_folder+'dt_xf_'+file_name))
        else:
            dt.append([0])        
    return(dt)



#%%
#tiene ax[1]
# el numero de orden order donde levantas el archivo es cero. 

def plot_dt_alpha(description_file,data_folder,save_path_name):
    
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);

    colors =  sns.color_palette(sns.color_palette("viridis",len(ref.groupby(['D']))))

###############################################################################
### Figure
###############################################################################  
    
    fig, axs = plt.subplots(2, 1, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();        
    

    for k,(D,row) in enumerate(ref.groupby(['D'])):
        
###############################################################################
### Parameters and data -- plotting
###############################################################################

        omega =  row.omega.unique()[0]
        alphas = row.alpha.values/omega
        
        dt =download_dt_alpha(row,data_folder)
        
        scale = 1000 #es por lo que tengo que dividir para llegar a minutos
        mean_dt = [np.mean(i)/scale for i in dt]
        std_dt  = [np.std(i)/scale for i in dt]
        
        axs[0].plot(alphas,mean_dt,'-o',linewidth = 1,color=colors[k],label = D)
        axs[0].fill_between(alphas,[i-j for i,j in zip(mean_dt,std_dt)],[i+j for i,j in zip(mean_dt,std_dt)],linewidth = 0, color =colors[k],alpha = 0.2)
        
        axs[1].plot(alphas,mean_dt,'-o',linewidth = 1, color=colors[k],label = D)
        axs[1].fill_between(alphas,[i-j for i,j in zip(mean_dt,std_dt)],[i+j for i,j in zip(mean_dt,std_dt)],linewidth = 0,color =colors[k],alpha = 0.2)


    axs[1].set_xscale('log',basex=2); axs[0].set_xscale('linear')
    axs[1].set_ylabel('mean duration (min)', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    #axs[0].set_xlim([0,100]);
    axs[1].set_ylim([-1,100]);
    
    axs[1].xaxis.set_label_coords(0.5, -0.1);
    axs[1].yaxis.set_label_coords(-0.1, 0.5)
    
    axs[0].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    axs[1].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    
    xticks = [np.round(i,2) for i in alphas]
    axs[1].set_xticks(xticks); axs[1].set_xticklabels(xticks)
    axs[1].tick_params(labelsize=10)

    plt.savefig(save_path_name + 'dt_alpha.pdf', format='pdf')
    return(0)

