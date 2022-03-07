#import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.data_managing_functions import download_data, check_file
from adler.get_time_series.main import all_combinations
from adler.plotting.plotting_main import set_scale,silent_ax
import seaborn as sns
#%%
    
###############################################################################
###  plotting
###############################################################################  
    
def plot_theta_alpha(data_folder,save_path_name,params):

    
    #colors =  sns.color_palette(sns.color_palette("viridis",len(params['alpha'])))
    
    D_ = len(params['D']); ALP = len(params['alpha'])
    
    fig, axs = plt.subplots(ALP,D_, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel(); 
  
    #para cada fila
    for k,(omega,alpha,D) in enumerate(all_combinations(params)):
        
        ax = axs[k]
        print(alpha/omega,omega,D)
        cond_prob_filename = data_folder +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        initial_conditions_filename = data_folder +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        
        if check_file(cond_prob_filename,""):
            cond_prob = download_data(cond_prob_filename)
            initial_conditions = download_data(initial_conditions_filename)
    
            ax.plot(initial_conditions,cond_prob,linewidth=1) #,color=colors[k]
            ax.plot(initial_conditions,cond_prob,'o', markersize = 2)
            ax.set_ylim([-1,110]); ax.set_xlim([-np.pi/2,3/2*np.pi])
              
        if k == (D_*ALP - D_ ):
            ax.set_ylabel("ocurrences", fontsize=10);
            ax.set_xlabel("initial conditions", fontsize=10)
            ax.xaxis.set_label_coords(0.5, -0.5);
            ax.yaxis.set_label_coords(-0.4, 0.5)
            #set_scale(ax,[-np.pi/2,3/2*np.pi],[-0.5,110])
            ax.set_xticks([-np.pi/2,3/2*np.pi]);ax.set_yticks([-1,110])
            ax.set_yticklabels([-1,110]); ax.set_xticklabels([r'$-\frac{\pi}{2}$',r'$\frac{3 \pi}{2}$']); ax.tick_params(labelsize=10) 
        else:
            silent_ax(ax)
            
        if k < D_:
            text = 'D = ' + str(D)
            ax.text(0.7, 1.1, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)


        if k % D_ == 0 :
            text = r'$\alpha = $' + str(np.round(alpha/omega,4)) 
            ax.text(0.7, 0.8, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)


        
    plt.savefig(save_path_name + 'plotting.pdf', format='pdf')
    return(0)

#%%
    
###############################################################################
###  plotting
###############################################################################  
    
def plot_first_IV(data_folder,save_path_name,params):

    
    

    
    fig, axs = plt.subplots(2,2, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    ax = axs[0,0]; 
  
    omega = params['omega'][0]
    
########################## D plot #############################################  
   
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['alpha'])))

    for k,alpha in enumerate(params['alpha']):
        
        result = []; D_aux = []
       
        for D in params['D']:                    
            print(alpha/omega,omega,D)
            cond_prob_filename = data_folder +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        #initial_conditions_filename = data_folder +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        
            if check_file(cond_prob_filename,""):
                cond_prob = download_data(cond_prob_filename)
                #initial_conditions = download_data(initial_conditions_filename)
                result.append(cond_prob[0])
                D_aux.append(D)
                        
        #ax.plot(initial_conditions,cond_prob,linewidth=1) #,color=colors[k]
        ax.plot(D_aux,result,'o', markersize = 4,color = colors[k],label = str(alpha/omega))
            
        #ax.set_ylim([-1,110]); ax.set_xlim([-np.pi/2,3/2*np.pi])
              
        ax.set_ylabel("ocurrences", fontsize=10);
        ax.set_xlabel("noise D", fontsize=10)
        ax.xaxis.set_label_coords(0.5, -0.1);
        ax.yaxis.set_label_coords(-0.15, 0.5)
        ax.legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
        #set_scale(ax,[-np.pi/2,3/2*np.pi],[-0.5,110])
        #ax.set_xticks([-np.pi/2,3/2*np.pi]);ax.set_yticks([-1,110])
        #ax.set_yticklabels([-1,110]); ax.set_xticklabels([r'$-\frac{\pi}{2}$',r'$\frac{3 \pi}{2}$']); ax.tick_params(labelsize=10) 
        
    
    
######################## alpha plot ##########################################  
    ax = axs[1,0]; 
    colors =  sns.color_palette(sns.color_palette("viridis",len(params['D'])))

    for k,D in enumerate(params['D']):
        
        result = []; ALP_aux = []
       
        for alpha in params['alpha']:                    
            print(alpha/omega,omega,D)
            cond_prob_filename = data_folder +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        #initial_conditions_filename = data_folder +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        
            if check_file(cond_prob_filename,""):
                cond_prob = download_data(cond_prob_filename)
                #initial_conditions = download_data(initial_conditions_filename)
                result.append(cond_prob[0])
                ALP_aux.append(D)
                        
        #ax.plot(initial_conditions,cond_prob,linewidth=1) #,color=colors[k]
        ax.plot(ALP_aux,result,'o', markersize = 4,color = colors[k],label = str(D))
            
        #ax.set_ylim([-1,110]); ax.set_xlim([-np.pi/2,3/2*np.pi])
              
        ax.set_ylabel("ocurrences", fontsize=10);
        ax.set_xlabel("alpha", fontsize=10)
        ax.xaxis.set_label_coords(0.5, -0.1);
        ax.yaxis.set_label_coords(-0.15, 0.5)
        ax.legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)

        

        
    plt.savefig(save_path_name + 'plotting_FIC.pdf', format='pdf')
    return(0)

