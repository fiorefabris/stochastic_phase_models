'''
Auxiliary plotting functions module

Author: Fiorella Fabris
Date  : 01/2020 

'''

import pickle
import matplotlib.ticker as ticker
import os

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

    
    

