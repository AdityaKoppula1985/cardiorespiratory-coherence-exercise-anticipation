# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 12:06:11 2024

@author: g.tec
"""


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib import style
import numpy as np
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
def rmsd(x,y):
    out= np.sqrt(np.mean((x-y)**2))
    return out
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
from scipy import signal
import neuron
from neuron import h
from neuron import h,rxd
from neuron.units import s,ms,mV
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
#%% simulation pipeline for breathing variability
tim1=np.arange(0,300,0.025*1e-3)

#resp 1: sinusoidal breathing
y1=np.sin(2*np.pi*0.2*tim1)
x1=y1

exp=0.8

#resp 2: respiratory rate variability
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(exp,size=70,loc=0.15,scale=0.25)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(70)
sig=[]
freqs1=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)    
    freqs1=np.append(freqs1,d1[ind])
tim=np.arange(0,len(sig)/40000,1/40000)
x2=sig    
tim2=tim

#resp 3: respiratory rate and depth variability
random.seed(70)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(exp,size=70)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(70)
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)     
tim=np.arange(0,len(sig)/40000,1/40000)    
tim3=tim
x3=sig

#resp 4: respiratory rate & depth variability and sighs
random.seed(85)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(exp,size=70)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(70)
K[:6]=[3,6,4,10,7,8]
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)     
tim=np.arange(0,len(sig)/40000,1/40000)    
tim4=tim
x4=sig


fig,ax=plt.subplots(2,2,sharex=True)
ax[0,0].plot(tim1,x1,linewidth=2)
ax[0,0].grid()
x[0,0].set_title('(a)',fontsize=12,fontweight='bold',loc='left')
ax[0,0].set_ylim(-2,2)
ax[1,0].plot(tim2,x2,linewidth=2)
ax[1,0].grid()
ax[1,0].set_title('(b)',fontsize=12,fontweight='bold',loc='left')
ax[1,0].set_xlabel('Time [sec]',fontsize=12,fontweight='bold')
ax[1,0].set_ylim(-2,2)
ax[0,1].plot(tim3,x3,linewidth=2)
ax[0,1].grid()
ax[0,1].set_title('(c)',fontsize=12,fontweight='bold',loc='left')
ax[0,1].set_ylim(-2,2)
ax[1,1].plot(tim4,x4,linewidth=2)
ax[1,1].grid()
ax[1,1].set_title('(d)',fontsize=12,fontweight='bold',loc='left')
ax[1,1].set_xlabel('Time [sec]',fontsize=12,fontweight='bold')
plt.xlim(0,300)
plt.suptitle('Simulation of Respiration',fontsize=18,fontweight='bold')

#%% pacemaker cell spontaneous output

dur=30000
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
    seg.IKr.gKr=0.01
    seg.cat1g.gbar=0.094

iclamp=h.IClamp(soma(0.5))
t1=np.arange(0,dur,0.025)
t2=h.Vector(t1)


v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)


h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)


fig,ax=plt.subplots(2,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[1].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[1].set_xlabel('Time [s]',fontsize=12,fontweight='bold')
plt.xlim(0,30)
ax[1].grid()
ax[0].set_title('(a)',fontsize=12,fontweight='bold',loc='left')
ax[1].set_title('(b)',fontsize=12,fontweight='bold',loc='left')
plt.suptitle('Output from the pacemaker model',fontsize=18,fontweight='bold')

#%% simulation of vagal tone
icl=[-1,-1.8,-2.5]
vt=['low vagal tone','medium vagal tone','high vagal tone']
c=['blue','green','red']

fig,ax=plt.subplots(2,1,sharex=True)
for i in np.arange(len(icl)):
  
    print('inj current:'+str(icl[i]))
    os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
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
        seg.IKr.gKr=0.01
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=icl[i]
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)

    I_inj=h.Vector().record(iclamp._ref_i)
    tx1=h.Vector().record(h._ref_t)
    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)

    y=np.array(v1)
    time=np.array(tx1)
    
  
    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time)
    ax[0].plot(time*1e-3,np.array(I_inj),color=c[i],label=vt[i],linewidth=3)
    ax[1].plot(time*1e-3,60/(1e-3*peakT_inter),color=c[i],label=vt[i],linewidth=3)
ax[0].legend()
ax[1].legend()      
ax[0].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[0].grid()
ax[0].set_ylim(-3,1)
ax[1].grid()
ax[1].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[1].set_xlabel('Time[s]',fontsize=12,fontweight='bold')    
ax[0].set_title('(a)',fontsize=12,fontweight='bold',loc='left')
ax[1].set_title('(b)',fontsize=12,fontweight='bold',loc='left')
plt.xlim(0,300)
plt.suptitle('Effect of Vagal tone strength',fontsize=18,fontweight='bold')
 
#%% simulation of vagal tone & respiratory vagal modulation
icl=[-0.3,-1,-2]

vt=['low vagal tone','medium vagal tone','high vagal tone']
c=['blue','green','red']
DF=pd.DataFrame(columns=['cond','hrmin','hrmax'])
fig,ax=plt.subplots(2,1,sharex=True)
for i in np.arange(len(icl)):
    
    print('inj current:'+str(icl[i]))
    os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
    dur=300000
    soma=h.Section(name='soma')
    soma.L=100
    soma.diam=20
    soma.cm=2
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
        seg.IKr.gKr=0.01
        seg.cat1g.gbar=0.094

    iclamp=h.IClamp(soma(0.5))
    iclamp.delay=20000
    iclamp.dur=dur
    iclamp.amp=0
    t1=np.arange(0,dur,0.025)
    t2=h.Vector(t1)
    y1=30*np.sin(2*np.pi*0.2*t1*1e-3)
    y2=h.Vector(icl[i]+(y1*1e-2))
    y2.play(iclamp._ref_amp,t2,1)
    I_inj=h.Vector().record(iclamp._ref_i)
    v1=h.Vector().record(soma(0.5)._ref_v)
    tx1=h.Vector().record(h._ref_t)
    h.load_file('stdrun.hoc')
    h.finitialize(-80*mV)
    h.continuerun(dur)
    y=np.array(v1)
    time=np.array(tx1)
      
    peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
    tpeaks=time[peaks1[0]]
    peakT=np.diff(tpeaks)
    cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
    peakT_inter=cs1(time)
    
    hr=60/(1e-3*peakT_inter)
    hrmin=min(hr[800000:(290*40000)])
    hrmax=max(hr[800000:(290*40000)])
    DF=DF.append({'cond':vt[i],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
    ax[0].plot(time*1e-3,np.array(I_inj),color=c[i],label=vt[i],linewidth=3)
    ax[1].plot(time*1e-3,60/(1e-3*peakT_inter),color=c[i],label=vt[i],linewidth=3)
ax[0].legend(framealpha=0.8)
ax[1].legend(framealpha=0.8)      
ax[0].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[0].grid()
ax[0].set_ylim(-3,1)
ax[1].grid()
ax[1].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[1].set_xlabel('Time[s]',fontsize=12,fontweight='bold')    
ax[0].set_title('(a)',fontsize=12,fontweight='bold',loc='left')
ax[1].set_title('(b)',fontsize=12,fontweight='bold',loc='left')
plt.xlim(0,300)
plt.suptitle('Effect of Vagal tone strength/basal heart rate on RSA',fontsize=18,fontweight='bold')
#%% simulation of vagal tone, respiratory-vagal modulation & depth of breathing

icl=[-0.6,-2]
amp=[30,60]

color=['blue','green','brown','red']

k=0
DF=pd.DataFrame(columns=['vt','vm','hrmin','hrmax'])
fig,ax=plt.subplots(2,1,sharex=True)
for i in np.arange(len(icl)):
    print('inj current:'+str(icl[i]))
    for j in np.arange(len(amp)):
        print('mod amp:'+str(amp[j]))
    
        os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
        
        dur=300000
        soma=h.Section(name='soma')
        soma.L=100
        soma.diam=20
        soma.cm=2
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
            seg.IKr.gKr=0.01
            seg.cat1g.gbar=0.094
    
        iclamp=h.IClamp(soma(0.5))
        iclamp.delay=20000
        iclamp.dur=dur
        iclamp.amp=0
        t1=np.arange(0,dur,0.025)
        t2=h.Vector(t1)
        y1=amp[j]*np.sin(2*np.pi*0.2*t1*1e-3)    
        y2=h.Vector(icl[i]+(y1*1e-2)) 
        y2.play(iclamp._ref_amp,t2,1)   
        I_inj=h.Vector().record(iclamp._ref_i)
        v1=h.Vector().record(soma(0.5)._ref_v)
        tx1=h.Vector().record(h._ref_t)    
        h.load_file('stdrun.hoc')
        h.finitialize(-80*mV)
        h.continuerun(dur)   
        y=np.array(v1)
        time=np.array(tx1)
        peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
        tpeaks=time[peaks1[0]]
        peakT=np.diff(tpeaks)
        cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
        peakT_inter=cs1(time)
        
        hr=60/(1e-3*peakT_inter)
        hrmin=min(hr[800000:(290*40000)])
        hrmax=max(hr[800000:(290*40000)])
        DF=DF.append({'vt':icl[i],'vm':amp[j],'hrmin':hrmin,'hrmax':hrmax},ignore_index=True)
        ax[0].plot(time*1e-3,np.array(I_inj),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2),color=color[k],linewidth=2)
        ax[1].plot(time*1e-3,60/(1e-3*peakT_inter),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2),color=color[k],linewidth=2)
        k=k+1
ax[0].legend(framealpha=0.8)
ax[1].legend(framealpha=0.8)      
ax[0].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[0].grid()
ax[0].set_ylim(-3,1)
ax[1].grid()
ax[1].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[1].set_xlabel('Time[s]',fontsize=12,fontweight='bold')    
plt.xlim(0,300)
plt.suptitle(' Effect of interaction between Vagal tone and vagal modulation on RSA',fontsize=18,fontweight='bold')
#%% Coherence with various patterns of patterns with/without logistic transformation

#%% sinusoidal breathing

###without logistic transformation
icl=[-1,-2]
amp=[10,20]
color=['blue','green','brown','red']

k=0
coharray=np.zeros((129,5))
fig,ax=plt.subplots(1,1)
for i in np.arange(len(icl)):
    print('inj current:'+str(icl[i]))
    for j in np.arange(len(amp)):
        print('mod amp:'+str(amp[j]))
        os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
        #h.nrn_load_dll('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2/nrnmech.dll')
        fig2,ax2=plt.subplots(3,1,sharex=True)
        dur=300000
        soma=h.Section(name='soma')
        soma.L=100
        soma.diam=20
        soma.cm=2
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
            seg.IKr.gKr=0.01
            seg.cat1g.gbar=0.094
    
        iclamp=h.IClamp(soma(0.5))
        iclamp.delay=20000
        iclamp.dur=dur
        iclamp.amp=0
        t1=np.arange(0,dur,0.025)
        t2=h.Vector(t1)
        y1=amp[j]*np.sin(2*np.pi*0.2*t1*1e-3)
        y2=h.Vector(icl[i]+(y1*1e-2))
        
        # y2=h.Vector(icl[i]+(logy(y1,0.25,0.75,5)*1e-2))
      
        y2.play(iclamp._ref_amp,t2,1)
    
        I_inj=h.Vector().record(iclamp._ref_i)
        v1=h.Vector().record(soma(0.5)._ref_v)
        tx1=h.Vector().record(h._ref_t) 
        h.load_file('stdrun.hoc')
        h.finitialize(-80*mV)
        h.continuerun(dur)    
        y=np.array(v1)
        time=np.array(tx1)
        peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
        tpeaks=time[peaks1[0]]
        peakT=np.diff(tpeaks)
        cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
        peakT_inter=cs1(time)
        
        hr=60/(1e-3*peakT_inter)

        tr=np.arange(0,300,1/4)
        cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
        x1=cs1(tr)
        cs2=scipy.interpolate.CubicSpline(time*1e-3, np.array(I_inj)[:len(time)])
        resp=cs2(tr)
        
        coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
        freq=coh[0]
        coharray[:,0]=freq
        coharray[:,k+1]=coh[1]
        ax.plot(coh[0],coh[1],label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2),color=color[k])
        ax2[0].plot(time*1e-3,y,label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
       # ax2[0].grid()
        ax2[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
        ax2[1].plot(time*1e-3,np.array(I_inj),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
        ax2[2].plot(time*1e-3,60/(1e-3*peakT_inter),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
        ax2[0].grid()
        ax2[1].grid()
        ax2[2].grid()
        ax2[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')  
        ax2[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
        ax2[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
        ax2[0].set_title('(A)',loc='left',fontsize=12,fontweight='bold')
        ax2[1].set_title('(B)',loc='left',fontsize=12,fontweight='bold')
        ax2[2].set_title('(C)',loc='left',fontsize=12,fontweight='bold')
        plt.suptitle('Sinusoidal respiratory-vagal input without logistic transformation'+' ('+'vt:'+str(icl[i])+' nA, '+'vm:'+str(amp[j]*1e-2)+' nA)',fontsize=16,fontweight='bold')
        k=k+1

ax.legend(framealpha=0.8)
ax.grid()
ax.set_xlabel('Frequency [Hz]',fontsize=12,fontweight='bold') 
ax.set_ylabel('Coherence',fontsize=12,fontweight='bold')
ax.set_ylim(0,1.1)
ax.set_xlim(0,1)
 
np.savetxt('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/coharray_sin_nolog.txt',coharray)


### with logistic transformation

icl=[-1,-2]
amp=[10,20]
color=['blue','green','brown','red']

k=0
coharray=np.zeros((129,5))
fig,ax=plt.subplots(1,1)
for i in np.arange(len(icl)):
    print('inj current:'+str(icl[i]))
    for j in np.arange(len(amp)):
        print('mod amp:'+str(amp[j]))
        os.chdir('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2')
        #h.nrn_load_dll('D:/iithLIBRARY/IIT/CRCmodeling/pacemkr/pace2/nrnmech.dll')
        fig2,ax2=plt.subplots(3,1,sharex=True)
        dur=300000
        soma=h.Section(name='soma')
        soma.L=100
        soma.diam=20
        soma.cm=2
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
            seg.IKr.gKr=0.01
            seg.cat1g.gbar=0.094
    
        iclamp=h.IClamp(soma(0.5))
        iclamp.delay=20000
        iclamp.dur=dur
        iclamp.amp=0
        t1=np.arange(0,dur,0.025)
        t2=h.Vector(t1)
        y1=amp[j]*np.sin(2*np.pi*0.2*t1*1e-3)
        #y2=h.Vector(icl[i]+(y1*1e-2))
        
        y2=h.Vector(icl[i]+(logy(y1,0.25,0.75,5)*1e-2))
      
        y2.play(iclamp._ref_amp,t2,1)
    
        I_inj=h.Vector().record(iclamp._ref_i)
        v1=h.Vector().record(soma(0.5)._ref_v)
        tx1=h.Vector().record(h._ref_t) 
        h.load_file('stdrun.hoc')
        h.finitialize(-80*mV)
        h.continuerun(dur)    
        y=np.array(v1)
        time=np.array(tx1)
        peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
        tpeaks=time[peaks1[0]]
        peakT=np.diff(tpeaks)
        cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
        peakT_inter=cs1(time)
        
        hr=60/(1e-3*peakT_inter)

        tr=np.arange(0,300,1/4)
        cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
        x1=cs1(tr)
        cs2=scipy.interpolate.CubicSpline(time*1e-3, np.array(I_inj)[:len(time)])
        resp=cs2(tr)
        
        coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)
        freq=coh[0]
        coharray[:,0]=freq
        coharray[:,k+1]=coh[1]
        ax.plot(coh[0],coh[1],label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2),color=color[k])
        ax2[0].plot(time*1e-3,y,label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
       # ax2[0].grid()
        ax2[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
        ax2[1].plot(time*1e-3,np.array(I_inj),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
        ax2[2].plot(time*1e-3,60/(1e-3*peakT_inter),label='vt:'+str(icl[i])+', '+'vm:'+str(amp[j]*1e-2))
        ax2[0].grid()
        ax2[1].grid()
        ax2[2].grid()
        ax2[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')  
        ax2[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
        ax2[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
        ax2[0].set_title('(A)',loc='left',fontsize=12,fontweight='bold')
        ax2[1].set_title('(B)',loc='left',fontsize=12,fontweight='bold')
        ax2[2].set_title('(C)',loc='left',fontsize=12,fontweight='bold')
        plt.suptitle('Sinusoidal respiratory-vagal input with logistic transformation'+' ('+'vt:'+str(icl[i])+' nA, '+'vm:'+str(amp[j]*1e-2)+' nA)',fontsize=16,fontweight='bold')
        k=k+1

ax.legend(framealpha=0.8)
ax.grid()
ax.set_xlabel('Frequency [Hz]',fontsize=12,fontweight='bold') 
ax.set_ylabel('Coherence',fontsize=12,fontweight='bold')
ax.set_ylim(0,1.1)
ax.set_xlim(0,1)
 
np.savetxt('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/coharray_sin_log.txt',coharray)
#%% Respiratory rate variability
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


random.seed(80)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
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
y2=h.Vector(-1.5+((y1_1)*1e-2))
# y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)
y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
ax[0].set_title('(A)',loc='left',fontsize=12,fontweight='bold')
ax[1].set_title('(B)',loc='left',fontsize=12,fontweight='bold')
ax[2].set_title('(C)',loc='left',fontsize=12,fontweight='bold')
#plt.suptitle('Effect of respiratory rate variability (without logistic transformation)',fontsize=16,fontweight='bold')
plt.suptitle('Effect of respiratory rate variability (without logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRV_coh_nolog.csv')


# with logistic transfromation
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


random.seed(80)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
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
# y2=h.Vector(-1.5+((y1_1)*1e-2))
y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)
y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1220],resp[80:1220],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
ax[0].set_title('(A)',loc='left',fontsize=12,fontweight='bold')
ax[1].set_title('(B)',loc='left',fontsize=12,fontweight='bold')
ax[2].set_title('(C)',loc='left',fontsize=12,fontweight='bold')
plt.suptitle('Effect of respiratory rate variability (with logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRV_coh_log.csv')
#%% Respiratory rate and depth variability

# without logistic transformation

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
random.seed(80)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  

y1=sig[1:len(t1)+1]

y1_1=y1*10
y2=h.Vector(-1.5+((y1_1)*1e-2))
# y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)

tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & amplitude variability (without logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAV_coh_nolog.csv')

# With logistic transformation
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
random.seed(80)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  

y1=sig[1:len(t1)+1]

y1_1=y1*10
# y2=h.Vector(-1.5+((y1_1)*1e-2))
y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)*1e-2))
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)

tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & amplitude variability (with logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAV_coh_log.csv')
#%% Respiratory rate. depth variability and sighs

# without logistic transformation
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

random.seed(70)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))

indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K[:4]=[3,6,4,10]
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  
   
y1=sig[:len(t1)]
y1_1=y1*10
y2=h.Vector(-1.5+((y1_1)*1e-2))
# y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5))*1e-2)
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('SAN rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & amplitude variability with sighs (without logistic transformation) ',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAVS_coh_nolog.csv')

# with logistic transformation

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

random.seed(70)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))

indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K[:4]=[3,6,4,10]
K=K+np.random.normal(0,0.3,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  
   
y1=sig[:len(t1)]
y1_1=y1*10
# y2=h.Vector(-1.5+((y1_1)*1e-2))
y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5))*1e-2)
y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)
h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('SAN rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & amplitude variability with sighs (with logistic transformation) ',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAVS_coh_log.csv')
#%% Respiratory rate. depth variability and sighs with sympathetic modulation at low/high vago-sympathetic tone

# low vagal tone
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

random.seed(70)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
#d2=np.array(random.sample(list(d1),len(d1)))
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K[:4]=[3,6,4,10]
random.seed(85)
K=K+np.random.normal(0,0.4,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  
#tim=np.arange(0,len(sig)/40000,1/40000)   
y1=sig[:len(t1)]

y1_1=y1*10
y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)

y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
# y2=h.Vector(-2+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)

y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)

h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & depth variability, sighs and sympathetic modulation (with logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAVSsymp_coh_log1.csv')

#high vagal tone
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

random.seed(70)
freqs=np.arange(0.15,0.25,0.005)
l=scipy.stats.powerlaw.rvs(0.6,size=60)
fac=(max(freqs)-min(freqs))/(max(l)-min(l))
d1=np.sort((l*fac) + min(freqs))
indx=np.arange(0,len(d1)).astype(int)
#d2=np.array(random.sample(list(d1),len(d1)))
indx2=np.array(random.sample(list(indx),len(indx)))
K=np.ones(60)
K[:4]=[3,6,4,10]
random.seed(85)
K=K+np.random.normal(0,0.4,len(K))
sig=[]
for i in np.arange(len(indx2)):
    ind=indx2[i]
    t=np.arange(0,1/d1[ind],1/40000)
    x=K[ind]*np.sin(2*np.pi*d1[ind]*t+(3*np.pi/2))
    sig=np.append(sig,x)  
#tim=np.arange(0,len(sig)/40000,1/40000)   
y1=sig[:len(t1)]

y1_1=y1*10
y1_2=15*np.sin(2*np.pi*0.08*t1*1e-3)

# y2=h.Vector(-1.5+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)
y2=h.Vector(-2+(logy(y1_1,0.25,0.75,5)+logy(y1_2,0.25,0.75,5))*1e-2)

y2.play(iclamp._ref_amp,t2,1)
I_inj=h.Vector().record(iclamp._ref_i)
v1=h.Vector().record(soma(0.5)._ref_v)
tx1=h.Vector().record(h._ref_t)

h.load_file('stdrun.hoc')
h.finitialize(-80*mV)
h.continuerun(dur)

y=np.array(v1)
time=np.array(tx1)

peaks1=scipy.signal.find_peaks(y,height=0,distance=12000)
tpeaks=time[peaks1[0]]
peakT=np.diff(tpeaks)
cs1=scipy.interpolate.CubicSpline(tpeaks[:len(tpeaks)-1],peakT)
peakT_inter=cs1(time)

hr=60/(1e-3*peakT_inter)
tr=np.arange(0,300,1/4)
cs1=scipy.interpolate.CubicSpline(time*1e-3, hr[:len(time)])
x1=cs1(tr)
cs2=scipy.interpolate.CubicSpline(t1*1e-3, y1_1[:len(time)])
resp=cs2(tr)
coh=scipy.signal.coherence(x1[80:1160],resp[80:1160],fs=4,nperseg=256)

df=pd.DataFrame()
df['freq']=coh[0]
df['coh']=coh[1]

fig,ax=plt.subplots(3,1,sharex=True)
ax[0].plot(time*1e-3,y)
ax[0].grid()
ax[0].set_ylabel('Membrane potential [mV]',fontsize=12,fontweight='bold')
ax[1].plot(time*1e-3,np.array(I_inj),color='r',label='Total Inj current')
ax[1].plot(t1*1e-3,(y1_1*1e-2),color='b',label='Resp')
ax[1].set_ylabel('Injected current (nA)',fontsize=12,fontweight='bold')
ax[1].grid()
ax[1].legend(framealpha=0.7)

ax[2].plot(time*1e-3,60/(1e-3*peakT_inter),color='red')
ax[2].set_ylabel('Pacemaker rate (spikes/min)',fontsize=12,fontweight='bold')
ax[2].set_xlabel('Time[s]',fontsize=12,fontweight='bold')
ax[1].set_ylim(-2.5,1.5)
ax[2].grid()
plt.xlim(0,300)
plt.suptitle('Effect of respiratory rate & depth variability, sighs and sympathetic modulation (with logistic transformation)',fontsize=16,fontweight='bold')
df.to_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/RRAVSsymp_coh_log2.csv')
#%% combined plot of coherence for all respiratroy variability patterns: sanity check
fig,ax=plt.subplots(2,2,sharex=True,sharey=True)

sin_nolog=np.loadtxt('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/coharray_sin_nolog.txt')

ax[0,0].plot(sin_nolog[:,0],sin_nolog[:,1],label='vt: -1 nA, vm: 0.1 nA',linewidth=2)
ax[0,0].plot(sin_nolog[:,0],sin_nolog[:,2],label='vt: -1 nA, vm: 0.2 nA',linewidth=2)
ax[0,0].plot(sin_nolog[:,0],sin_nolog[:,3],label='vt: -2 nA, vm: 0.1 nA',linewidth=2)
ax[0,0].plot(sin_nolog[:,0],sin_nolog[:,4],label='vt: -2 nA, vm: 0.2 nA',linewidth=2)
ax[0,0].legend(framealpha=0.7)
ax[0,0].grid()
ax[0,0].set_title('(a)',fontsize=12,fontweight='bold',loc='left')
ax[0,0].set_xlim(0,1)
ax[0,0].set_ylim(0,1.2)

sin_log=np.loadtxt('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/coharray_sin_log.txt')

ax[0,1].plot(sin_log[:,0],sin_log[:,1],label='vt: -1 nA, vm: 0.1 nA',linewidth=2)
ax[0,1].plot(sin_log[:,0],sin_log[:,2],label='vt: -1 nA, vm: 0.2 nA',linewidth=2)
ax[0,1].plot(sin_log[:,0],sin_log[:,3],label='vt: -2 nA, vm: 0.1 nA',linewidth=2)
ax[0,1].plot(sin_log[:,0],sin_log[:,4],label='vt: -2 nA, vm: 0.2 nA',linewidth=2)
ax[0,1].legend(framealpha=0.7)
ax[0,1].grid()
ax[0,1].set_title('(b)',fontsize=12,fontweight='bold',loc='left')

rrv_nolog=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rrv_coh_nolog.csv').iloc[:,[1,2]]
rrav_nolog=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rrav_coh_nolog.csv').iloc[:,[1,2]]
rravs_nolog=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravs_coh_nolog.csv').iloc[:,[1,2]]
rravssymp_nolog=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravssymp_coh_nolog.csv').iloc[:,[1,2]]
rravssymp_nolog2=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravssymp_coh_nolog2.csv').iloc[:,[1,2]]

ax[1,0].plot(rrv_nolog.iloc[:,0],rrv_nolog.iloc[:,1],label='RRV',linewidth=2)
ax[1,0].plot(rrav_nolog.iloc[:,0],rrav_nolog.iloc[:,1],label='RRAV',linewidth=2)
ax[1,0].plot(rravs_nolog.iloc[:,0],rravs_nolog.iloc[:,1],label='RRAVS',linewidth=2)
ax[1,0].plot(rravssymp_nolog.iloc[:,0],rravssymp_nolog.iloc[:,1],label='RRAVS Symp 1',linewidth=2)
ax[1,0].plot(rravssymp_nolog2.iloc[:,0],rravssymp_nolog2.iloc[:,1],label='RRAVS Symp 2',linewidth=2)
ax[1,0].legend(framealpha=0.7)
ax[1,0].grid()
ax[1,0].set_title('(c)',fontsize=12,fontweight='bold',loc='left')

rrv_log=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rrv_coh_log.csv').iloc[:,[1,2]]
rrav_log=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rrav_coh_log.csv').iloc[:,[1,2]]
rravs_log=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravs_coh_log.csv').iloc[:,[1,2]]
rravssymp_log=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravssymp_coh_log.csv').iloc[:,[1,2]]
rravssymp_log2=pd.read_csv('D:/iithLIBRARY/IIT/CRCmodeling/figs/thesis plots/chap 4/new data plots 2/rravssymp_coh_log2.csv').iloc[:,[1,2]]

ax[1,1].plot(rrv_log.iloc[:,0],rrv_log.iloc[:,1],label='RRV',linewidth=2)
ax[1,1].plot(rrav_log.iloc[:,0],rrav_log.iloc[:,1],label='RRAV',linewidth=2)
ax[1,1].plot(rravs_log.iloc[:,0],rravs_log.iloc[:,1],label='RRAVS',linewidth=2)
ax[1,1].plot(rravssymp_log.iloc[:,0],rravssymp_log.iloc[:,1],label='RRAVS Symp 1',linewidth=2)
ax[1,1].plot(rravssymp_log2.iloc[:,0],rravssymp_log2.iloc[:,1],label='RRAVS Symp 2',linewidth=2)
ax[1,1].legend(framealpha=0.7)
ax[1,1].grid()
ax[1,1].set_title('(d)',fontsize=12,fontweight='bold',loc='left')

plt.suptitle ('Factors affecting Model "cardiorespiratory" Coherence ',fontsize=18,fontweight='bold')

ax[1,0].set_xlabel('Frequency [Hz]',fontsize=12,fontweight='bold') 
ax[1,1].set_xlabel('Frequency [Hz]',fontsize=12,fontweight='bold') 
ax[0,0].set_ylabel('Coherence',fontsize=12,fontweight='bold') 
ax[1,0].set_ylabel('Coherence',fontsize=12,fontweight='bold')

ax[1,0].set_xticklabels(labels=ax[1,0].get_xticklabels(),fontweight='bold')
ax[1,1].set_xticklabels(labels=ax[1,1].get_xticklabels(),fontweight='bold')
ax[0,0].set_yticklabels(labels=ax[0,0].get_yticklabels(),fontweight='bold')
ax[1,0].set_yticklabels(labels=ax[1,0].get_yticklabels(),fontweight='bold')
#%%