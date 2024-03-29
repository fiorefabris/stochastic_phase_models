import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.data_managing_functions import download_data,check_file,get_fixed_points,compute_theoretical_dt
from adler.get_time_series.get_simulated_dt import get_all_combinations_alphas,get_all_combinations_fixed_alphas
from adler.plotting.plotting_main import set_scale,silent_ax


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
    colors_eps =  sns.color_palette(sns.color_palette("Greys",len(params['epsilon'])))
    k = 0 ; 
       
    fig, axs = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel(); 
  

    for (omega,epsilon),alphas in get_all_combinations_alphas(params):
        alphas = [a/omega for a in alphas]


        dt_filename = data_folder + 'simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(epsilon)+'.pkl'
        if check_file('simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(epsilon)+'.pkl',data_folder):
            dt =  download_data(dt_filename)
        else:
            dt = [[] for i in alphas]
            print("no file")
            
        scale = 1 #1000 #es por lo que tengo que dividir para llegar a minutos
        mean_dt = [np.mean(i)/scale for i in dt]
        print(epsilon,mean_dt)

        axs[0].plot(alphas,mean_dt,'-o',linewidth = 1,color=colors[k],label = epsilon)
        axs[1].plot(alphas,mean_dt,'-o',linewidth = 1, color=colors[k],label = epsilon)

        synth_alphas = np.linspace(1e-10,3,10000) #alphas[1:] 
        aux = [compute_theoretical_dt(omega,epsilon,delta) for delta in synth_alphas]
        axs[0].plot(synth_alphas,aux,'-',color = colors_eps[k],alpha = 0.8,label = epsilon)
        k = k+1


    axs[1].set_xscale('linear'); axs[0].set_xscale('linear') #axs[1].set_xscale('log',basex=2)
    axs[1].set_ylabel('mean duration (min)', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    axs[0].set_ylim([-1,100]);axs[0].set_xlim([0.97,2]); 
    axs[1].set_ylim([-1,100]);axs[1].set_xlim([0.97,2]); 
    
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
        theta_filename = 'theta_example_simulated_dt_'+str(np.round(omega,2))+'_eps_'+str(epsilon)+'.pkl'
        
        print(check_file(theta_filename ,data_folder))
        print(theta_filename)
        if check_file(theta_filename ,data_folder):
            theta = download_data(data_folder + theta_filename)
            
        else:
            theta = [[] for i in alphas]
            print(theta)
        for l,th in enumerate(theta):
            ax = axs[k,l]
            theta_end,theta_beg = get_fixed_points(alphas[l])

            ax.plot(np.cos(th),np.sin(th),linewidth=1,color=colors[k])
            ax.plot(np.cos(theta_beg),np.sin(theta_beg),'o',color='red')
            ax.plot(np.cos(theta_end),np.sin(theta_end),'o',color='red')
            ax.set_ylim([-1.05,1.05]); ax.set_xlim([-1.05,1.05])
              
            if k == EPS-1 and l == 0 :
                ax.set_ylabel(r'$\sin(\theta)$', fontsize=10);
                ax.set_xlabel(r'$\cos(\theta)$', fontsize=10)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
                set_scale(ax,[-1,1],[-1,1])
                ax.set_xticklabels([-1,1]); ax.set_yticklabels([-1,1]); ax.tick_params(labelsize=10) 

            elif k == 0:
                text = r'$\alpha = $' + str(alphas[l]) 
                ax.text(0.7, 0.8, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)
                silent_ax(ax)


            else:
                silent_ax(ax)

            if l == 0 :
                text = 'delta tita = ' + str(epsilon)
                ax.text(0.7, 1.1, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)

        
    plt.savefig(save_path_name + 'time_series_.pdf', format='pdf')
    return(0)




#%%
# =============================================================================
# =============================================================================
# =============================================================================
#                   PLOTTING SIMULATED FIXED  DT MODULE
# =============================================================================
# =============================================================================
# =============================================================================


###############################################################################
### Figure data
###############################################################################  
    
    
def plot_simulated_fixed_dt_alpha(data_folder,save_path_name,params):
    
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['epsilon'])))
    colors_eps =  sns.color_palette(sns.color_palette("Greys",len(params['epsilon'])))
    k = 0 ; 
       
    fig, axs = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel(); 
  

    for (T_0,epsilon),deltas in get_all_combinations_fixed_alphas(params):


        dt_filename = data_folder + 'simulated_fixed_dt_'+str(np.round(T_0,2))+'_eps_'+str(epsilon)+'.pkl'
        if check_file('simulated_fixed_dt_'+str(np.round(T_0,2))+'_eps_'+str(epsilon)+'.pkl',data_folder):
            dt =  download_data(dt_filename)
        else:
            dt = [[] for i in deltas]
            print("no file")
            
        scale = 1 #1000 #es por lo que tengo que dividir para llegar a minutos
        mean_dt = [np.mean(i)/scale for i in dt]

        axs[0].plot(deltas,mean_dt,'-o',linewidth = 1,color=colors[k],alpha = 0.6,label = epsilon)
        axs[1].plot(deltas,mean_dt,'-o',linewidth = 1, color=colors[k],alpha = 0.6,label = epsilon)

        #synth_alphas = np.linspace(1e-10,3,10000) #alphas[1:] 
        #aux = [compute_theoretical_dt(omega,epsilon,delta) for delta in synth_alphas]
        #axs[0].plot(synth_alphas,aux,'-',color = colors_eps[k],alpha = 0.8,label = epsilon)
        k = k+1


    axs[1].set_xscale('linear'); axs[0].set_xscale('linear') #axs[1].set_xscale('log',basex=2)
    axs[1].set_ylabel('mean duration (min)', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    axs[0].set_ylim([-1,10]);axs[0].set_xlim([0.97,2]); 
    axs[1].set_ylim([-1,10]);axs[1].set_xlim([0.97,2]); 
    
    axs[1].xaxis.set_label_coords(0.5, -0.1);
    axs[1].yaxis.set_label_coords(-0.1, 0.5)
    
    axs[0].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    axs[1].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    
    #xticks = [np.round(i,2) for i in alphas]
   # axs[1].set_xticks(xticks); axs[1].set_xticklabels(xticks)
    axs[1].tick_params(labelsize=10)

    plt.savefig(save_path_name + 'simulated_fixed_dt_alpha.pdf', format='pdf')
    return(0)
    
    
#%%
    
###############################################################################
### theta plotting
###############################################################################  
    
def plot_theta_fixed_alpha(data_folder,save_path_name,params):

    
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['epsilon'])))
    
    EPS = len(params['epsilon']); ALP = len(params['delta'])
    fig, axs = plt.subplots(EPS,ALP, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    #axs = axs.ravel(); 
  
    #para cada fila
    for k,((T_0,epsilon),deltas) in enumerate(get_all_combinations_fixed_alphas(params)):
        theta_filename = 'theta_example_simulated_fixed_dt_'+str(np.round(T_0,2))+'_eps_'+str(epsilon)+'.pkl'
        
        print(check_file(theta_filename ,data_folder))
        print(theta_filename)
        if check_file(theta_filename ,data_folder):
            theta = download_data(data_folder + theta_filename)
            
        else:
            theta = [[] for i in deltas]
            print(theta)

        for l,th in enumerate(theta):
            ax = axs[k,l]
            theta_end,theta_beg = get_fixed_points(deltas[l])
            
            ax.plot(np.cos(th),np.sin(th),linewidth=1,color=colors[k])
            ax.plot(np.cos(theta_beg),np.sin(theta_beg),'o',color='red',markersize=1.5)
            ax.plot(np.cos(theta_end),np.sin(theta_end),'o',color='red',markersize=1.5)
            ax.set_ylim([-1.05,1.05]); ax.set_xlim([-1.05,1.05])
              
            if k == EPS-1 and l == 0 :
                ax.set_ylabel(r'$\sin(\theta)$', fontsize=10);
                ax.set_xlabel(r'$\cos(\theta)$', fontsize=10)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
                set_scale(ax,[-1,1],[-1,1])
                ax.set_xticklabels([-1,1]); ax.set_yticklabels([-1,1]); ax.tick_params(labelsize=10) 

            elif k == 0:
                text = r'$\alpha = $' + str(deltas[l]) 
                ax.text(0.7, 0.8, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)
                silent_ax(ax)


            else:
                silent_ax(ax)

            if l == 0 :
                text = 'delta tita = ' + str(epsilon)
                ax.text(0.7, 1.1, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)

        
    plt.savefig(save_path_name + 'time_series_fixed.pdf', format='pdf')
    return(0)