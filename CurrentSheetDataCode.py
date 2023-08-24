# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:06:57 2023

@author: jackm
"""

#Imports/Download Data

import cdflib
import xarray
import pandas as pd
import pyspedas
import pytplot
import matplotlib.pyplot as plt
import os
import numpy as np

#date = '2023-03-16'
#t1 = '18:00:00'
#t2 = '24:00:00'
date = '2021-04-29'
t1 = '3:00:00'
t2 = '6:00:00'
tstring1 = date + '/' + t1
tstring2 = date + '/' + t2
tstring = date +','+t1.replace(':',';')+'-'+t2.replace(':',';')

pyspedas.psp.fields(trange=[tstring1,tstring2],time_clip=True)
tinit, Binit = pytplot.get_data('psp_fld_l2_mag_RTN').times, pytplot.get_data('psp_fld_l2_mag_RTN').y

#Getting PVI

os.mkdir(r'C:/Users/jackm/Desktop/SummerInternship/CurrentSheetPlots/{}'.format(tstring))


n1 = []

PVI_total = []

steps = [1,4,16,64,256,1028,1028*4]
dt = [3,3,3,3,6,10,20]
findex = []
tcal = pyspedas.time_double(date+ ' 0:00:00')
ttemp = tinit-tcal
tstep = np.min([ttemp[i+1]-ttemp[i] for i in range(len(ttemp)-1)])
t = np.arange(ttemp[0],ttemp[-1],tstep)
for a in range(7):
    step = steps[a]
    dB = []
    #tao = (tstep)*step
    #print(tao)
    
    Btemp = np.array(Binit).transpose()

    B = [np.interp(t,ttemp,Btemp[0]),np.interp(t,ttemp,Btemp[1]),np.interp(t,ttemp,Btemp[2])]


    period = int(3600/tstep)*2
    for i in range((len(B[0])-step)):
        diff = []
        for j in range(3):
            diff.append((B[j][i+step]-B[j][i]))
        dB.append(diff)
        
    dB = np.array(dB).transpose()   
    dB2 = np.array([[b**2 for b in bb] for bb in dB])
    bean = True
    start = 0
    sigma = []
    while bean:
        if start+period>len(dB2[0]):
            period = len(dB2[0])-start   
            bean = False
        if start == len(dB2[0]):
            bean = False
            break
        '''
        s = []
        for i in [0,1,2]:
            sum = 0
            av = np.average(dB[i][start:start+period])
            for j in range(start,start+period):
                sum = sum + (dB[i][j]-av)**2
            s.append(sum/len(dB[i]))
            '''
        s = [np.average(dB2[0][start:start+period]),np.average(dB2[1][start:start+period]),np.average(dB2[2][start:start+period])]
        sigma.append(s)
        start += period 
    dB_norm = []    
    period = int(3600/tstep)*2
    for i in range(len(dB2[0])):
        si = sigma[int(i/period)]
        dB_norm.append([dB2[0][i]/si[0],dB2[1][i]/si[1],dB2[2][i]/si[2]])
        
    PVI = np.array([np.sqrt(np.sum(b)) for b in dB_norm])
    PVI_total.append(PVI)
    
    interval = int(dt[a]/(tstep))
    bean = False
    count = 0
    i = 0
    floor = 5
    sub = []
    n = []


    #Collecting all High PVI Events
    while i < len(PVI):
        if PVI[i]>floor:      
            start = i
            while PVI[i+1]>floor:
                i+=1
                if i+1>=len(PVI):
                    break
            diff = i-start
            if diff<2*interval:
                extend = int(interval-diff/2)
            else:
                extend = 0
            sub.append((np.max([start-extend,0]),np.min([i+extend,len(t)-1])))
            n.append(step)
    
        i+=1
                
    #combining close events
    tolerance = -dt[a]/2
    tolerance = 0
    index = []
    count = -1
    start = True
    end = False
    for b in range(len(sub)):
        if start:
            add = 0
            index.append([])
            count += 1
            index[count].append(sub[b][0])
            start = False
            
        if b == len(sub)-1:
            end == True
        else:
            t0 = np.average(sub[b])
            t1 = np.average(sub[b+1])
            if t[sub[b+1][0]]-t[sub[b][1]] < (tolerance-add):
                end = False
                add+=dt[a]/30
            else:
                end = True
        if end:
            index[count].append(sub[b][1])
            start = True
            #removing redundancies 
            bean = True
            
            for f in findex:
                
                if index[count][0]>f[0]-step and index[count][1]<f[1]+step:
                    bean = False
                    break
                if np.abs(np.average([t[f[0]],t[f[1]]])-np.average([t[index[count][0]],t[index[count][1]]]))<dt[a]:
                    bean = False
                    break
                if np.abs(index[count][0]-f[0])<dt[a] and np.abs(index[count][1]-f[1])<dt[a]:
                    bean = False
                    break
            
            if bean:
                findex.append(index[count])
                n1.append(n[b])

findex_t = [t[findex[i]] for i in range(len(findex))]
sub_t = [(t[sub[i][0]],t[sub[i][1]]) for i in range(len(sub))]


bp  = []    
tp = []
    
#Performing MVA, Saving Plots 

for f in range(len(findex)):
    bounds = findex[f]
    tt = t[bounds[0]:bounds[1]]
    b = np.array(B).transpose()[bounds[0]:bounds[1]]
    M = []
    for i in [0,1,2]:
        M.append([])
        for j in [0,1,2]:
            sum = 0
            for z in range(len(b)):
                sum += b[z][i]*b[z][j]
            m = sum/len(b) - np.average(b[:,i])*np.average(b[:,j])
            M[i].append(m)

    w,v = np.linalg.eig(M)    
    lamda = 0    
    max = np.max(w)
    min = np.min(w)   
    for i in range(len(w)):
        if np.abs(w[i]-max) < 0.000001:
            l = v[:,i]
        elif np.abs(w[i]-min) < 0.000001:
            n = v[:,i]
        else:
            m = v[:,i]
    
    bl = []
    bn = []
    bm = []
    bean = False
    for i in range(len(b)):
        bl.append(b[i][0]*l[0]+b[i][1]*l[1]+b[i][2]*l[2])
        #bn.append(b[i][0]*n[0]+b[i][1]*n[1]+b[i][2]*n[2])
        #bm.append(b[i][0]*m[0]+b[i][1]*m[1]+b[i][2]*m[2])
        
    plt.plot(tt,bl)
    plt.title('Count:{} , Tao = {}'.format(f,n1[f]))
    plt.ylabel('Magnetic Field [nT]')
    plt.xlabel('Time [s]')
    plt.savefig(r'C:/Users/jackm/Desktop/SummerInternship/CurrentSheetPlots/{}/[{:.3f}-{:.3f}].png'.format(tstring,t[bounds[0]],t[np.min([len(t)-1,bounds[1]])]))
    plt.close()
    #bp.append([bl,bn,bm])
    #tp.append(tt)
    
print(len(findex))
#print(len(sub))
#%%
#PLOTTING ENTIRE B MEASURMENTS
tt = (t[::step][:-1] - pyspedas.time_double('2021-04-29 0:00:00'))
plt.plot(tt,B[:,0][:-1])
plt.show()
plt.plot(tt,PVI)
#%%
tlow = 0
thigh = 20
index = [i for i in range(len(PVI)) if t[i]-tcal>tlow and t[i]-tcal<thigh]
plt.plot(t[index]-tcal,PVI[index])

count = 0
i = 0
while t[i]-tcal < 20:
    if PVI[i]>5:
        count+=1
        print(i)
        
    i+=1
print(count)
