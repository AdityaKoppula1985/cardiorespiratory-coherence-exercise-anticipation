# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:46:09 2024

@author: g.tec
"""

import neuron
import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import FormatStrFormatter


from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from tkinter import Toplevel, Button, Tk, Menu  
from matplotlib import style
import tkinter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy 
import biosppy
import seaborn as sns

import mne
import time
import os
import random
import adi
def stdz(x):
    return (x-np.mean(x))/np.std(x)
def rstdz(y,x):
    return (y*np.std(x)/np.std(y))+np.mean(x)
def rmsd(x,y):
    out= np.sqrt(np.mean((x-y)**2))
    return out




#%% simulation data
conds=['rav','ravs','rrav','rrav2','rravs','rravs2','rrv','rrv2','rrvs','rrvs2','rs','sinusoidal']
auto=['AT','LT','BT','HT','SM','SM_BT','SM_LT','SM_HT']

simcohavg_df=pd.DataFrame()
combs=[]
for i in np.arange(len(conds)):
    condsx=conds[i]
    for j in np.arange(len(auto)):
        autox=auto[j]
        fileid='cohmat_'+condsx+'_'+autox
        path='E:/crcsims/'+condsx+'/'
        dat=np.loadtxt(path+fileid+'.txt')
        freq=dat[:,0]
        cohavg=np.mean(dat[:,1:],axis=1)
        simcohavg_df['freq']=freq
        simcohavg_df[fileid]=cohavg
simcohavg_df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/crcsims/simcohavg.csv',index=None)
        
#%%    
simcohavg_df=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/crcsims/simcohavg.csv')   
cohavg_df=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/crcsims/exptcohavg.csv')

err_df=pd.DataFrame(columns=cohavg_df.columns[1:],index=simcohavg_df.columns[1:])
for i in np.arange(cohavg_df.shape[1]-1):
    sig1=np.array(cohavg_df.iloc[:,i+1])
    for j in np.arange(simcohavg_df.shape[1]-1):
        sig2=np.array(simcohavg_df.iloc[:,j+1])
        err_df.iloc[j,i]=rmsd(sig1,sig2)

        
        
base_bestfit=err_df['Base'].sort_values()[:3]
print('base_bestfit')
print(base_bestfit) 
    
h35l_bestfit=err_df['H35L'].sort_values()[:3]
print('h35l_bestfit')
print(h35l_bestfit)

h35r_bestfit=err_df['H35R'].sort_values()[:3]
print('h35r_bestfit')
print(h35r_bestfit)

h50l_bestfit=err_df['H50L'].sort_values()[:3]
print('h50l_bestfit')
print(h50l_bestfit)

h50r_bestfit=err_df['H50R'].sort_values()[:3]
print('h50r_bestfit')
print(h50r_bestfit)

be50l_bestfit=err_df['BE50L'].sort_values()[:3]
print('BE50L_bestfit')
print(be50l_bestfit)


be50r_bestfit=err_df['BE50R'].sort_values()[:3]
print('BE50R_bestfit')
print(be50r_bestfit)

be70l_bestfit=err_df['BE70L'].sort_values()[:3]
print('BE70L_bestfit')
print(be70l_bestfit)


be70r_bestfit=err_df['BE70R'].sort_values()[:3]
print('BE70R_bestfit')
print(be70r_bestfit)

be50bl_bestfit=err_df['BE50BL'].sort_values()[:3]
print('BE50BL_bestfit')
print(be50bl_bestfit)


be70bl_bestfit=err_df['BE70BL'].sort_values()[:3]
print('BE70BL_bestfit')
print(be70bl_bestfit)
