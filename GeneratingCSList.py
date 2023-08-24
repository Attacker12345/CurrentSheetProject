# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:56:31 2023

@author: jackm
"""
import cdflib
import xarray
import pandas as pd
import pyspedas
import pytplot
import matplotlib.pyplot as plt
import os
import numpy as np


f = open('/Users/kimberleymawhinney/Desktop/Internship/currentsheetdata.csv','r')
f.readline()
s = 'a;lskdfj'
dates = []
interv = []
tbounds_init = []
while len(s)>1:
    s = f.readline()
    if s == '':
        break
    i0 = s.find(',')
    i = s.find(',',i0+1)
    i1 = s.find(',',i+1)
    i2 = s.find(',',i1+1)
    date = s[1:i0]
    t1 = float(s[i+1:i1])
    t2 = float(s[i1+1:i2])
    if date == '2021-08-11':
        break
    else:
        dates.append(date)
        tbounds_init.append([t1,t2])
        interv.append(s[1:i-1])
    
f.close()

datetimes = ['2023-03-16,18;00;00-24;00;00','2023-03-17,0;00;00-2;00;00',
             '2021-04-29,8;00;00-10;00;00','2021-08-10,14;00;00-18;00;00',
             '2021-08-10,18;00;00-19;00;00','2021-08-11,6;00;00-12;00;00',
             '2021-04-29,3;00;00-6;00;00']

for d in datetimes:
#for d in ['2023-03-17,0;00;00-2;00;00']:
    date = d[:10]
    
        
    #READING CURRENT SHEET BOUNDS
    
    direct = sorted(os.listdir('/Users/kimberleymawhinney/Desktop/Internship/CurrentSheetPlots/'+d+' Final'))
    
    for st in direct:
        if st[0] == 'i':
            under = st.find('_')
            dash = st.find('-')
            t1 = float(st[under+1:dash])
            t2 = float(st[dash+1:-4])
            tbounds_init.append([t1,t2])
            dates.append(date)
            interv.append(d)

           
datetimes = ['2021-04-28,10;00;00-12;00;00','2021-04-28,12;00;00-14;00;00',
            '2021-11-22,2;45;00-6;00;00','2021-11-22,6;00;00-10;30;00',
            '2021-11-21,21;30;00-24;00;00','2021-11-22,0;00;00-0;45;00','2022-09-06,18;00;00-24;00;00'
          ,'2022-09-07,0;00;00-6;00;00','2022-09-07,6;00;00-12;00;00','2022-09-06,12;00;00-18;00;00'
          ,'2022-09-07,12;00;00-18;00;00','2021-08-10,18;00;00-19;00;00',
          '2023-03-16,18;00;00-24;00;00','2023-03-17,0;00;00-2;00;00','2021-04-29,8;00;00-10;00;00',
          '2021-08-10,14;00;00-18;00;00','2021-08-11,6;00;00-12;00;00','2021-04-29,3;00;00-6;00;00']
        

transitionbounds = []
tbounds = np.array([[]])
tbounds_init = np.array(tbounds_init)
interv = np.array(interv)
dates = []
for i in range(len(datetimes)):
#for i in range(14):
    
    date = datetimes[i][:10]

    q = np.where(interv == datetimes[i])
    #print(q)
    #trash,tbounds_temp = zip(*sorted(zip(np.transpose(tbounds_init[q])[0],tbounds_init[q])))
    q1 = np.argsort(np.transpose(tbounds_init[q])[0])
    tbounds_temp = tbounds_init[q][q1]

    
    transitionbounds_init = []
    if i < 11:
        d = datetimes[i].replace(';', '_')
        direct = os.listdir('/Users/kimberleymawhinney/Desktop/Internship/CurrentSheetPlots [ARCHIVE]/'+d+' Transition Bounds')
    else:
        d = datetimes[i]
        direct = os.listdir('/Users/kimberleymawhinney/Desktop/Internship/CurrentSheetPlots/'+d+' Transition Bounds')
    for st in direct:
        if not st == '.DS_Store':
            under = st.find('_')
            dash = st.find('-')
            t1 = float(st[under+1:dash])
            t2 = float(st[dash+1:-4])
            transitionbounds_init.append([t1,t2])
 
    transitionbounds_init = np.array(transitionbounds_init)
    q2 = np.argsort(np.transpose(transitionbounds_init)[0])
    if i == 0:
        tbounds = np.array(tbounds_temp)
        transitionbounds = np.array(transitionbounds_init[q2])
    else:
        tbounds = np.concatenate((tbounds,np.array(tbounds_temp)))
        transitionbounds = np.concatenate((transitionbounds,np.array(transitionbounds_init[q2])))
        
    for j in range(len(tbounds_temp)):
        dates.append(date)
        
    
    for i in range(len(tbounds)):
        if transitionbounds[i][0]<tbounds[i][0] or transitionbounds[i][1]>tbounds[i][1]:
            test = range(i-1,i+1)
            for j in test:
                if transitionbounds[j][0]>tbounds[i][0] or transitionbounds[j][1]<tbounds[i][1]:
                    transitionbounds[i],transitionbounds[j] =  transitionbounds[j],transitionbounds[i]
                    break
                
qf = np.argsort(dates)
tbounds = tbounds[qf]
transitionbounds = transitionbounds[qf]
dates = np.array(dates)[qf]

import csv

with open('/Users/kimberleymawhinney/Desktop/Internship/currentsheets.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerow(['Date','tL (bounds)','tR (bounds)','tL (transition)','tR (transition)'])
    for i in range(len(tbounds)):
        writer.writerow([dates[i],tbounds[i][0],tbounds[i][1],transitionbounds[i][0],transitionbounds[i][1]])
        #writer.writerow([dates[i],tbounds[i][0],tbounds[i][1]])
  
#%%
for i in range(len(tbounds)):
    if transitionbounds[i][0]<tbounds[i][0]:
        print(dates[i],tbounds[i],transitionbounds[i])
        
    if transitionbounds[i][1]>tbounds[i][1]:
        print(dates[i-1],tbounds[i-1],transitionbounds[i-1])
        print(dates[i],tbounds[i],transitionbounds[i])
        print(dates[i+1],tbounds[i+1],transitionbounds[i+1])
        


