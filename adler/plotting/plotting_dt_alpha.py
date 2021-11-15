import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.data_managing_functions import download_data,check_file,get_fixed_points
from adler.get_time_series.get_simulated_dt import get_all_combinations_alphas
from adler.plotting.plotting_main import set_scale,silent_ax

''' Este modulo es para el primer approach que tuvimos para simular dt'''

#el resultado teorico empueza en alphas[1:]
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

def cociente(delta):
    return np.sqrt((delta - 1)/(delta + 1))

def f_(epsilon,delta):
    return epsilon * cociente(delta) / (2 - epsilon * cociente(delta))
    
def compute_theoretical_dt(omega,epsilon,delta):
    return -2 / (omega * np.sqrt(delta**2 - 1)) * np.log(f_(epsilon,delta))


#%%
#compute_theoretical_dt(omega,1e-100,1.0005) #tiende a cero eventualemnte
#
#EPS   =np.linspace(2,5,5) #PASA ALGO RARO EN 5
#omega = 2* np.pi/7
#DEL = [1.001,1.005,1.01,1.03,1.05,1.1,1.3,1.5]
#
#for delta in DEL:
#    aux = [compute_theoretical_dt(omega,eps,delta) for eps in EPS]
#    #aux = [np.log(f_(eps,delta)) for eps in EPS]
#    plt.plot(EPS,aux,'-o')# el eje y es en minutos creo
#    plt.ylim([0,100])
#
#for epsilon in EPS:
#    aux = [compute_theoretical_dt(omega,epsilon,delta) for delta in DEL]
#    plt.plot(DEL,aux,'-o') 
#    plt.ylim([0,100])
#%%
#tiene ax[1]
# el numero de orden order donde levantas el archivo es cero. 

''' Este modulo es para el segundo approach que tuvimos para simular dt'''


def plot_dt_alpha(description_file,data_folder,save_path_name):
    
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);

    colors =  sns.color_palette(sns.color_palette("viridis",len(ref.groupby(['D']))))

###############################################################################
### Figure
###############################################################################  
   
    plt.clf()
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
        #axs[0].fill_between(alphas,[i-j for i,j in zip(mean_dt,std_dt)],[i+j for i,j in zip(mean_dt,std_dt)],linewidth = 0, color =colors[k],alpha = 0.2)
        
        
        axs[1].plot(alphas,mean_dt,'-o',linewidth = 1, color=colors[k],label = D)
        #axs[1].fill_between(alphas,[i-j for i,j in zip(mean_dt,std_dt)],[i+j for i,j in zip(mean_dt,std_dt)],linewidth = 0,color =colors[k],alpha = 0.2)

    EPS = np.linspace(1,10,5)
    colors_eps =  sns.color_palette(sns.color_palette("Greys",len(EPS)))

    for m,eps in enumerate(EPS):
        deltas = np.linspace(1.0005005,1.5,1000) #alphas[1:] 
        aux = [compute_theoretical_dt(omega,eps,delta) for delta in deltas]
        axs[0].plot(deltas,aux,'-',color = colors_eps[m],alpha = 0.8,label = eps)


    axs[1].set_xscale('linear'); axs[0].set_xscale('linear') #axs[1].set_xscale('log',basex=2)
    axs[1].set_ylabel('mean duration (min)', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    axs[0].set_ylim([-1,100]);axs[0].set_xlim([1,1.1]); 
    axs[1].set_ylim([-1,100]);axs[1].set_xlim([1,1.1]); 
    
    axs[1].xaxis.set_label_coords(0.5, -0.1);
    axs[1].yaxis.set_label_coords(-0.1, 0.5)
    
    axs[0].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    axs[1].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    
    #xticks = [np.round(i,2) for i in alphas]
   # axs[1].set_xticks(xticks); axs[1].set_xticklabels(xticks)
    axs[1].tick_params(labelsize=10)

    plt.savefig(save_path_name + 'dt_alpha.pdf', format='pdf')
    return(0)

#%%
# =============================================================================
# =============================================================================
# =============================================================================
#                   PLOTTING SIMULATED DT MODULE
# =============================================================================
# =============================================================================
# =============================================================================


###############################################################################
### Figure data
###############################################################################  
    
    
def plot_simulated_dt_alpha(data_folder,save_path_name,params):
    
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['epsilon'])))
    k = 0 ; 
       
    fig, axs = plt.subplots(2, 1, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel(); 
  

    for (omega,epsilon),alphas in get_all_combinations_alphas(params):
        alphas = [a/omega for a in alphas]


        dt_filename = data_folder + 'simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(np.round(epsilon,2))+'.pkl'
        if check_file('simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(np.round(epsilon,2))+'.pkl',data_folder):
            dt =  download_data(dt_filename)
            theta = data_folder + 'theta_example_simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(np.round(epsilon,2))+'.pkl'
            theta = download_data(theta)
        else:
            dt = [[] for i in alphas]

        scale = 1 #1000 #es por lo que tengo que dividir para llegar a minutos
        mean_dt = [np.mean(i)/scale for i in dt]
        print(mean_dt)
        axs[0].plot(alphas,mean_dt,'-o',linewidth = 1,color=colors[k],label = epsilon)
        axs[1].plot(alphas,mean_dt,'-o',linewidth = 1, color=colors[k],label = epsilon)
        
        k = k+1
        
#    EPS = np.linspace(1,10,5)
#    colors_eps =  sns.color_palette(sns.color_palette("Greys",len(EPS)))

#    for m,eps in enumerate(EPS):
#        deltas = np.linspace(1.0005005,1.5,1000) #alphas[1:] 
#        aux = [compute_theoretical_dt(omega,eps,delta) for delta in deltas]
#        axs[0].plot(deltas,aux,'-',color = colors_eps[m],alpha = 0.8,label = eps)


    axs[1].set_xscale('linear'); axs[0].set_xscale('linear') #axs[1].set_xscale('log',basex=2)
    axs[1].set_ylabel('mean duration (min)', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    #axs[0].set_ylim([-1,100]);axs[0].set_xlim([1,1.1]); 
    #axs[1].set_ylim([-1,100]);axs[1].set_xlim([1,1.1]); 
    
    axs[1].xaxis.set_label_coords(0.5, -0.1);
    axs[1].yaxis.set_label_coords(-0.1, 0.5)
    
    axs[0].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    axs[1].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    
    #xticks = [np.round(i,2) for i in alphas]
   # axs[1].set_xticks(xticks); axs[1].set_xticklabels(xticks)
    axs[1].tick_params(labelsize=10)

    plt.savefig(save_path_name + 'simulated_dt_alpha.pdf', format='pdf')
    return(0)
    
    
#%%
    
###############################################################################
### theta plotting
###############################################################################  
    
def plot_theta_alpha(data_folder,save_path_name,params):

    
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['epsilon'])))
    
    EPS = len(params['epsilon']); ALP = len(params['alpha'])
    fig, axs = plt.subplots(EPS,ALP, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    #axs = axs.ravel(); 
  
    #para cada fila
    for k,((omega,epsilon),alphas) in enumerate(get_all_combinations_alphas(params)):
        alphas = [a/omega for a in alphas]
        theta_filename = 'theta_example_simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(np.round(epsilon,2))+'.pkl'
        
        if check_file(theta_filename ,data_folder):
            theta = download_data(data_folder + theta_filename)
        else:
            theta = [[] for i in alphas]
        
        for l,th in enumerate(theta):
            ax = axs[k,l]
            theta_end,theta_beg = get_fixed_points(alphas[l]/omega)

            ax.plot(np.cos(th),np.sin(th),linewidth=2,color=colors[k])
            ax.plot(np.cos(theta_beg),np.sin(theta_beg),'o',color='red')
            ax.plot(np.cos(theta_end),np.sin(theta_end),'o',color='red')
            ax.set_ylim([-1.05,1.05]); ax.set_xlim([-1.05,1.05])
              
            if k == EPS-1 and l == ALP-1 :
                ax.set_ylabel(r'$\sin(\theta)$', fontsize=10);
                ax.set_xlabel(r'$\cos(\theta)$', fontsize=10)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
                set_scale(ax,[-1,1],[-1,1])
                ax.set_xticklabels([-1,1]); ax.set_yticklabels([-1,1]); ax.tick_params(labelsize=10) 

            elif k == 0:
                text = r'$\alpha = $' + str(alphas[l]) 
                ax.text(0.7, 0.8, text , ha='center', va='center', transform=ax.transAxes, fontsize=10)
                silent_ax(ax)


            else:
                silent_ax(ax)

            if l == 0 :
                text = 'eps = ' + str(epsilon)
                ax.text(0.7, 1.1, text , ha='center', va='center', transform=ax.transAxes, fontsize=10)

        
    plt.savefig(save_path_name + 'time_series_.pdf', format='pdf')
    return(0)
