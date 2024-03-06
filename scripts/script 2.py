# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:41:03 2024

@author: g.tec
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 09:11:46 2023

@author: g.tec
"""

import neuron
import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import FormatStrFormatter

import tkinter
from tkinter import *
#from matplotlib.backends.backend_tkagg import (
#    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
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
# import spectrum
# import pyhrv
# import pyhrv.tools as tools
# from pyhrv.hrv import hrv
# import pyhrv.frequency_domain as fd
# import pyhrv.nonlinear as nl
# import sys
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

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
from scipy import signal

import pandas as pd
import numpy as np


from neuron import h
from neuron import h,rxd
from neuron.units import s,ms,mV


    
    
    
    
    

def myzeropad(arr,length):
    addzeros=length-len(arr)
    if length%2==0:
        return np.pad(arr,int(addzeros/2),mode='constant')
    else:
        return np.pad(arr,([int(addzeros//2),int((addzeros//2)+1)]),mode='constant')
   
def uwhis2(dat,lq,uq,cut):
    iqr=np.quantile(dat,uq)-np.quantile(dat,lq)
    out1=np.quantile(dat,uq)+(cut*iqr)
    return out1

def logy(x,lq,uq,cut):
    out1=uwhis2(x,lq,uq,cut)/(1+np.exp(-1*(x-np.median(x))))
    delta=np.mean(x)-np.mean(out1)
    return out1+delta
irr=['P3','P4','P5','P6','P7','P8','P12','P14','P15','P18','P20','P21','P22','P23','P24','P25']
no2=['P10','P11','P13','P16','P17','P28','P29','P30']
no=['P2','P9','P19','P26','P27']



def logity(x):
    var1=x/max(np.abs(x))
    var2=var1+np.abs(min(var1))
    var3=var2/max(var2)
    return scipy.special.logit(var3)


os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
h.nrn_load_dll('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2/nrnmech.dll')
#%%


import time

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

start = time.time()
 
os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
#h.nrn_load_dll('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2/nrnmech.dll')



dur=300000
soma=h.Section(name='soma')
soma.L=100
soma.diam=20
soma.cm=2
Vm=[-80]
ti=[0]
soma.insert('cat1g')
soma.insert('IKr')
soma.insert('IKs')
soma.insert('IK1')
soma.insert('ICaL')
soma.insert('INaf')
soma.insert('htc')
soma.insert('pas')

for seg  in soma:
    seg.pas.g=1e-6
    seg.pas.e=-80
    seg.htc.ghbar=1e-5
    seg.htc.eh=-10
    seg.IKr.Tauact=1.2 
    seg.IKr.gKr=0.015
    seg.cat1g.gbar=0.094

iclamp=h.IClamp(soma(0.5))
iclamp.delay=20000
iclamp.dur=dur
iclamp.amp=0
t1=np.arange(0,dur,0.025)
t2=h.Vector(t1)


exp=0.3

l=[0.2]*60

d1=np.sort(l)
indx=np.arange(0,len(d1)).astype(int)
#d2=np.array(random.sample(list(d1),len(d1)))
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K[:4]=[3,6,4,10]
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  

y1=sig#[1:len(t1)+1]

y1_1=y1*10

y2=h.Vector(-1.5+((y1_1)*1e-2))
#y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
#y2=h.Vector(y1*1e-2)
y2.play(iclamp._ref_amp,t2,1)

I_inj=h.Vector().record(iclamp._ref_i)

v1=h.Vector().record(soma(0.5)._ref_v)
#v2=h.Vector().record(soma(0.5)._ref_i_membrane)
tx1=h.Vector().record(h._ref_t)


h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time1=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time1[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time1)

hr=60/(1e-3*peakT_inter)
# hrmin=min(hr[800000:(290*40000)])
# hrmax=max(hr[800000:(290*40000)])
# DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
#t=np.arange(0,duration,1/fs)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
freq=coh[0]    
    
    


end = time.time()
 
# print the difference between start
# and end time in milli. secs
print("The time of execution of above program is :",
      (end-start)/60, "mins")

np.savetxt('E:/crcsims/coh.txt',coh)



now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#%% sinusoidal

###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+((y1_1)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_BT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+((y1_1)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_AT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+((y1_1)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_LT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+((y1_1)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_HT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(y1_1+y1_2)*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(y1_1+y1_2)*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_LT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(y1_1+y1_2)*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_BT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(y1_1+y1_2)*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_HT_rawsin.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
#%% sinusoidal logarithmic resp 1

###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
   
np.savetxt('E:/crcsims/comat_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    #freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=60,loc=0.15,scale=0.25)
    l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig#[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
#%% Resp freq variability resp 

###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=61,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrv_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=61,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=61,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=80,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(80)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=80,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(60)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=80,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(80)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=80,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(80)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=80,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(80)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
#%% Resp amplitude variability 
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rav_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rav_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rav_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rav_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rav_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rav_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rav_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rav_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#%% sighs
#%% respiratory rate varibility + resp amplitude variabillity
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrav_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrav_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrav_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrav_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#%% resp rate variability and sighs
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#%% resp amplitude variability + sighs
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_ravs_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_ravs_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_ravs_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_ravs_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_ravs_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_ravs_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_ravs_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    #l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_ravs_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
#%% resp rate, amplitude variability + sighs

###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.3
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#%% resp rate/freq variability 2
##basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrv2_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv2_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(1000)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv2_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/comat_rrv2_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv2_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv2_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv2_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=1*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrv2_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
#%% respiratory rate varibility + resp amplitude variabillity 2
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrav2_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrav2_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrav2_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrav2_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav2_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav2_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav2_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    #K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrav2_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)



#%%resp rate variability and sighs 2
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs2_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs2_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs2_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rrvs2_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs2_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs2_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs2_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    #K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rrvs2_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)





#%%resp rate, amplitude variability + sighs 2
###basal tone
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    freq=coh[0]
    cohmat[:,j+1]=coh[1]
    if j==0:
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs2_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
###ABSENT TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs2_AT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


###Low TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs2_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###High TONE
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
#fig,ax=plt.subplots()
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*100
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]
    
    #ax.plot(y1)

    y1_1=y1*10
    #y1_1=logity(y1)

    #y1_1=-15*np.sin(2*np.pi*0.015*t1*1e-3)
    #y1_2=30*np.sin(2*np.pi*0.08*t1*1e-3)

    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/comat_rravs2_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=0
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs2_SM.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+LT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs2_SM_LT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+BT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-1.5
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs2_SM_BT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

###SM+HT
from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
start = time.time()
cohmat=np.zeros((129,21))
at=-2
for j in np.arange(20):
    print('*******************'+'start_loop: '+str(j))
    
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
    Vm=[-80]
    ti=[0]
    soma.insert('cat1g')
    soma.insert('IKr')
    soma.insert('IKs')
    soma.insert('IK1')
    soma.insert('ICaL')
    soma.insert('INaf')
    soma.insert('htc')
    soma.insert('pas')

    for seg  in soma:
        seg.pas.g=1e-6
        seg.pas.e=-80
        seg.htc.ghbar=1e-5
        seg.htc.eh=-10
        seg.IKr.Tauact=1.2 
        seg.IKr.gKr=0.015
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)


    # fc = 0.2  
    # fm = 0.016
    # beta =0.2/200
    # modulation = np.sin(2 * np.pi * fm * t1*1e-3)    
    # fm_signal = (np.cos(2 * np.pi * (fc + beta * modulation) * t1*1e-3)+1e-2)
    # #y1=-amp[j]*fm_signal
    # #y1=-amp[j]*(4+np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal
    # y1=-(5+0.5*np.sin(2*np.pi*0.01*t1*1e-3))*fm_signal*4
    exp=0.8
    #random.seed(80)
    freqs=np.arange(0.15,0.25,0.005)
    l=scipy.stats.powerlaw.rvs(exp,size=100,loc=0.15,scale=0.25)
    #l=[0.2]*60
    #fac=(max(freqs)-min(freqs))/(max(l)-min(l))
    #d1=np.sort((l*fac) + min(freqs))
    d1=np.sort(l)
    indx=np.arange(0,len(d1)).astype(int)
    #d2=np.array(random.sample(list(d1),len(d1)))
    indx2=np.array(random.sample(list(indx),len(indx)))
    K=np.ones(100)
    K[:4]=[3,6,4,10]
    K=K+np.random.normal(0,0.3,len(K))
    sig=[]
    for i in np.arange(len(indx2)):
        ind=indx2[i]
        t=np.arange(0,1/d1[ind],1/40000)
        x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
        sig=np.append(sig,x)  

    y1=sig[1:len(t1)+1]

    y1_1=y1*10
    y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)
    y2=h.Vector(at+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
    #y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
    #y2=h.Vector(y1*1e-2)
    y2.play(iclamp._ref_amp,t2,1)

    I_inj=h.Vector().record(iclamp._ref_i)

    v1=h.Vector().record(soma(0.5)._ref_v)
    #v2=h.Vector().record(soma(0.5)._ref_i_membrane)
    tx1=h.Vector().record(h._ref_t)


    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time1=np.array(tx1)

    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time1[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time1)

    hr=60/(1e-3*peakT_inter)
    # hrmin=min(hr[800000:(290*40000)])
    # hrmax=max(hr[800000:(290*40000)])
    # DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    #t=np.arange(0,duration,1/fs)
    tr=np.arange(0,300,1/4)
    cs1=scipy.interpolate.CubicSpline(time1*1e-3, hr[:len(time1)])
    x1=cs1(tr)
    cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time1)])
    resp=cs2(tr)
    coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
    
    cohmat[:,j+1]=coh[1]
    if j==0:
        freq=coh[0]
        cohmat[:,0]=freq
    
np.savetxt('E:/crcsims/cohmat_rravs2_SM_HT.txt',cohmat)
end = time.time()
print("The time of execution of above program is :",
      (end-start)/60, "mins")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)




#%%

