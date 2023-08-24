#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:20:17 2023

@author: kimberleymawhinney
"""

import cdflib
import xarray
import pandas as pd
import pyspedas
import pytplot
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.optimize as fitter
import csv


d = 'Date'
datelist = []
tlist = []
intervbounds = []
with open('/Users/kimberleymawhinney/Desktop/Internship/currentsheets.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if not row[0] == 'Date':
            if d == row[0]:
                tlist.append(float(row[1]))
                tlist.append(float(row[2]))
            else:
                datelist.append(row[0])
                if not d == 'Date':
                    intervbounds.append([np.min(tlist),np.max(tlist)])
                    tlist = []
                d = row[0]
    intervbounds.append([np.min(tlist),np.max(tlist)])
    
    
#reading electron temp data [JASPER]

f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/coretpar.txt', 'r') 
f.readline()
s = 'a;lskdfj'
lines1 = []
while len(s)>1:
    s = f.readline()
    if s == '':
        break
    i = s.find('/')
    i1 = s.find('         ')
    lines1.append([s[1:i],s[i+1:i1],s[i1+9:-1]])

f.close()
        
f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/coretperp.txt', 'r') 
f.readline()
s = 'a;lskdfj'
lines2 = []
while len(s)>1:
    s = f.readline()
    if s == '':
        break
    i = s.find('/')
    i1 = s.find('         ')
    lines2.append([s[1:i],s[i+1:i1],s[i1+9:-1]])
            
f.close()

lines1 = np.array(lines1)
lines2 = np.array(lines2)

q = np.where(lines1[:,2] == '         -NaN')

lines = np.delete(np.array([[lines1[i][0],lines1[i][1],lines1[i][2],lines2[i][2]] for i in range(len(lines1))]),q,axis = 0)



def mag(v):
    return(np.sqrt(np.sum([a**2 for a in v])))

Brms_list = []
B0_list = []
betainterval = []
ebetainterval = []
alpha_list = []
tao_list = []
meanthickness_list = []
medianthickness_list = []
lowthickness_list = []

betalist = []
ebetalist = []
betainclude = []
for u in range(len(datelist)):
    Brms_list.append([])
    B0_list.append([])
    betainterval.append([])
    ebetainterval.append([])
    date = datelist[u]
    t1 = '0:00:00'
    t2 = '24:00:00'
    tstring1 = date + '/' + t1
    tstring2 = date + '/' + t2
    if date[:4] == '2023':
        year = 2023
        month = 3
        if date[-2:] == '17':
            day = 17
            st = '1700'
        else:
            day = 16
            st = '1618'
        
        epoch = cdflib.cdfepoch.compute_tt2000([year,month,day,0,0,0])

        a = cdflib.cdf_to_xarray(r'/Users/kimberleymawhinney/Desktop/Internship/CDF/Fields/psp_fld_l2_mag_RTN_202303'+st+'_v02.cdf')
        t = (a['epoch_mag_RTN'].to_numpy() - epoch)*10**-9
        B = a['psp_fld_l2_mag_RTN'].to_numpy()
    
    else:
        tcal = pytplot.time_double(date+' 0:00:00')
        pyspedas.psp.fields(trange=[tstring1,tstring2],time_clip=True)
        t, B = pytplot.get_data('psp_fld_l2_mag_RTN').times-tcal, pytplot.get_data('psp_fld_l2_mag_RTN').y   
    
    year = int(date[:4])
    month = int(date[5:7])
    day = int(date[8:])
    epoch = cdflib.cdfepoch.compute_tt2000([year,month,day,0,0,0])
    dstring = date.replace('-','')
    a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Sweap/psp_swp_spi_sf00_L3_mom_'+dstring+'_v04')
    ttemper = (a['Epoch'].to_numpy() - epoch)*10**-9
    #Obtaining Temp Data for specific day
    temp = (a['TEMP']).to_numpy()
    
    nanfinder = np.isnan(temp)
    q = np.where(nanfinder == False)
    q1 = np.where(temp[q]>0)
    temp = temp[q][q1]
    ttemper = ttemper[q][q1]

    #Obtaining Electron Density Data, Filtering Out Invalid Values
    if year == 2021:
        a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Fields/psp_fld_l3_rfs_lfr_qtn_20210109_20211201_v01')
        tdens2 = (a['Epoch'].to_numpy() - epoch)*10**-9
        dens2 =  a['N_elec'].to_numpy()
        q1 = np.where(tdens2>0) 
        q2 = np.where(tdens2<3600*24)
        q = np.intersect1d(q1, q2)
        tdens1 = tdens2[q]
        dens1 = dens2[q]
        p = np.where(dens1>0)
        dens = dens1[p]
        tdens = tdens1[p]
    else:          
        tdens = []
        dens = []
        if year == 2022 and month>=8 and month<=9:
            f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/qtn_2022-08-30_2022-09-14.txt', 'r')
            f.readline()
            
            s = 'a;lskdfj'
            while len(s)>1:
                s = f.readline()
                if s == '':
                    break
                if int(s[5:7]) == month and int(s[8:10]) == day:
                    i = s.find(',')
                    i1 = s.find(',',i+1)
                    i2 = s.find(',',i1+1)
                    timestr = s[11:i]
                    if i2>i1+1:
                        tdens.append(int(timestr[0:2])*3600+int(timestr[3:5])*60+float(timestr[6:]))
                        dens.append(float(s[i1+1:i2]))     
            f.close()
            
        if year ==2023: 
            if day==16:
                f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/Jack_E15_20230316.txt', 'r')
            else:
                f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/Jack_E15_20230317.txt', 'r')
            f.readline()
            s = 'a;lskdfj'
            while len(s)>1:
                s = f.readline()
                if s == '':
                    break
                i = s.find('\t')
                if i == -1:
                    i = s.find('	')
                if not s[i+1:-1] == 'NaN':
                    tdens.append(float(s[:i])*3600)
                    dens.append(float(s[i+1:-1]))
            f.close()
    
    dens = np.array(dens)
    
    if year == 2021 and month ==4:
        a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Electron/psp_fld_l3_sqtn_rfs_V3V4_20210419_20210505_v01_LESIA_E08.cdf')
        et_init = (a['Epoch'].to_numpy() - epoch)*10**-9
        eTemp_init=  a['electron_core_temperature'].to_numpy()
        q1 = np.where(et_init>0) 
        q2 = np.where(et_init<3600*24)
        q = np.intersect1d(q1, q2)
        eTemp_temp = eTemp_init[q]
        p = np.where(eTemp_temp>0)
        eTemp1 = eTemp_temp[p]
        et1 = et_init[q][p]
    elif year == 2023:
        f = open('/Users/kimberleymawhinney/Desktop/Internship/TXT/LFRV1-V2Ne_Tc_Thk_mimo_2023-03-07_2023-03-24.txt','r')
        f.readline()
        eTemp1 = []
        et1 = []
        s = 'a;lskdfj'
        while len(s)>1:
            s = f.readline()
            if s == '':
                break
            if int(s[5:7]) == month and int(s[8:10]) == day:
                i = s.find(',')
                i1 = s.find(',',i+1)
                i2 = s.find(',',i1+1)
                i3 = s.find(',',i2+1)
                i4 = s.find(',',i3+1)
                i5 = s.find(',',i4+1)
                timestr = s[11:i]
                if i5>i4+1:
                    et1.append(int(timestr[0:2])*3600+int(timestr[3:5])*60+float(timestr[6:]))
                    eTemp1.append(float(s[i4+1:i5]))     
        eTemp1 = np.array(eTemp1)
        et1 = np.array(et1)
        f.close()
    else:
        a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Electron/psp_fld_l3_sqtn_rfs_V1V2_'+dstring+'_v2.0.cdf')
        et_init = (a['Epoch'].to_numpy() - epoch)*10**-9
        eTemp_init=  a['electron_core_temperature'].to_numpy()
        p = np.where(eTemp_init>0)
        eTemp1 = eTemp_init[p]
        et1 = et_init[p]

    #GETTING JASPER DATA    
    q = np.where(lines[:,0] == date)
    et2 = []
    eTemp2 = []
    for l in lines[q]:
        timestr = l[1]
        et2.append(int(timestr[0:2])*3600+int(timestr[3:5])*60+float(timestr[6:]))
        eTemp2.append((float(l[2])+float(l[3]))/2)

    eTemp1 = np.array(eTemp1)
    eTemp2 = np.array(eTemp2)
    #CALCULATING BRMS FOR 10 Minute Intervals
    
    marker = intervbounds[u][0]
    step = 600
    while(True):
        if marker+step>intervbounds[u][1]:
            break
        else:
            include = True
            q1 = np.where(t>marker)
            q2 = np.where(t<marker+step)
            q = np.intersect1d(q1, q2)
            ttemp = t[q]
            Btemp = np.array(B[q])
            B0 = np.transpose([np.average(b) for b in np.transpose(Btemp)])
            dB = np.array([Btemp[i]-B0 for i in range(len(Btemp))])
            Brms = np.sqrt(np.average([np.sum(dB[i]**2) for i in range(len(dB))])) 
            marker = marker + step
            
            q1 = np.where(ttemper>marker)
            q2 = np.where(ttemper<marker+step)
            q = np.intersect1d(q1, q2) 
            if len(q)<1:
                include = False
            else:
                T = np.average(temp[q])
            
            q1 = np.where(tdens>marker)
            q2 = np.where(tdens<marker+step)
            q = np.intersect1d(q1, q2)
            if len(q)<1:
                include = False
            else:
                D = np.average(dens[q])
                
            q1 = np.where(et2>marker)
            q2 = np.where(et2<marker+step)
            q = np.intersect1d(q1, q2)
            if len(q)<1:
                q1 = np.where(et1>marker)
                q2 = np.where(et1<marker+step)
                q = np.intersect1d(q1, q2)
                if len(q)<1:
                    include = False
                else:
                    eT = np.average(eTemp1[q])
            else:
                eT = np.average(eTemp2[q])
            if include:
                beta = 0.48*D*T/(mag(B0)**2)
                ebeta = 0.48*D*eT/(mag(B0)**2)
                Brms_list[u].append(Brms)
                B0_list[u].append(np.sqrt(np.sum(B0**2)))
                betainterval[u].append(beta)
                ebetainterval[u].append(ebeta)
            
            
    marker = intervbounds[u][0]
    step = 3600
    while(True):
        if marker+step>intervbounds[u][1]:
            break
        else:
            transitionbounds1 = []
            dates = []
            with open('/Users/kimberleymawhinney/Desktop/Internship/currentsheets.csv', 'r') as file:
                reader = csv.reader(file)
                for row in reader:
                    if not row[0] == 'Date':
                        dates.append(row[0])
                        transitionbounds1.append([float(row[3]),float(row[4])])
            dates = np.array(dates)
            q = np.where(dates == date)
            transitionbounds1 = np.array(transitionbounds1)[q]
            
            transitionbounds = []
            for i in range(len(transitionbounds1)):
                if transitionbounds1[i][0]>marker and transitionbounds1[i][1]<marker+step:
                    transitionbounds.append(transitionbounds1[i])
            if len(transitionbounds)>10:
                transitionbounds = np.array(transitionbounds)
                thickness = [transition[1]-transition[0] for transition in transitionbounds]
                meanthickness_list.append(np.average(thickness)/2)
                medianthickness_list.append(np.median(thickness)/2)
                lowthickness_list.append(np.percentile(thickness,10)/2)
                        
                #FINDING TAO-ALPHA Relationship throughout the interval
                
    
                tdiff = [t[i+1]-t[i] for i in range(len(t)-1)]
                tstep = np.average(tdiff)
                tt = np.arange(marker,marker+step,tstep)
                BB = np.transpose([np.interp(tt, t, b) for b in np.transpose(B)])
                
                tao10 = np.linspace(np.log10(tstep),np.log10(30),120)
                taos = 10**tao10
                slist1 = [int(tao/tstep) for tao in taos]
                slist1[0] = 1
                slist = []
                stor = 0
                for i in range(len(slist1)):
                    if not slist1[i] == stor:
                        slist.append(slist1[i])
                        stor = slist1[i]
                        
                tao = []
                alpha = []
                #for s in range(1,upperlim,int(upperlim/100)):
                for s in slist:
                    tao.append(s*tstep)
                    alpha.append(np.average([np.arccos(np.dot(BB[i],BB[i+s])/(mag(BB[i])*mag(BB[i+s]))) for i in range(0,len(BB)-s,50)]))
                alpha_list.append(alpha)
                tao_list.append(tao)
                
                q1 = np.where(ttemper>marker)
                q2 = np.where(ttemper<marker+step)
                q = np.intersect1d(q1, q2) 
                if len(q)<1:
                    include = False
                else:   
                    T = np.average(temp[q])
                
                q1 = np.where(tdens>marker)
                q2 = np.where(tdens<marker+step)
                q = np.intersect1d(q1, q2)
                if len(q)<1:
                    include = False
                else:
                    D = np.average(dens[q])
                    
                q1 = np.where(et2>marker)
                q2 = np.where(et2<marker+step)
                q = np.intersect1d(q1, q2)
                if len(q)<1:
                    q1 = np.where(et1>marker)
                    q2 = np.where(et1<marker+step)
                    q = np.intersect1d(q1, q2)
                    if len(q)<1:
                        include = False
                    else:
                        eT = np.average(eTemp1[q])
                else:
                    eT = np.average(eTemp2[q])
                    
                if include:
                    betalist.append(0.48*D*T/(mag(B0)**2))
                    ebetalist.append(0.48*D*eT/(mag(B0)**2))
                else:
                    betalist.append(0)
                    ebetalist.append(0)
                    
                betainclude.append(include)
                
            marker = marker+step
            
    
    print(date)
    

#%%
criticalpoint = []

#for u in range(len(alpha_list)):
for u in [39]:
    plt.plot(tao_list[u],alpha_list[u])
    taolog = np.log10(tao_list[u])
    alphalog = np.log10(alpha_list[u])
    direv = np.array([(alphalog[i+1]-alphalog[i])/(taolog[i+1]-taolog[i]) for i in range(len(alphalog)-1)])
    # direv_diff = [(direv[i+5]-direv[i-5])/(taolog[i+5]-taolog[i-5]) for i in range(5,len(direv)-5)]
    # direv_diff2 = [(direv_diff[i+7]-direv_diff[i-7])/(taolog[i+7]-taolog[i-7]) for i in range(7,len(direv_diff)-7)]
    # index = np.argmax(direv_diff)
    # print(datelist[u],tao_list[u][index])
    plt.title('Alpha vs Tao')
    plt.xlabel('Tao')
    plt.ylabel('Alpha')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
    # #plt.plot(taolog,alphalog)
    # #plt.show()
    
    # q1 = np.where(taolog>0)
    # q2 = np.where(taolog<-1.5)
    
    # def model(x,a,b):
    #     return a*x+b
    
    # par0  = [1,0]
    # par, cov = fitter.curve_fit(model, taolog[q1], alphalog[q1], par0)

    # x = np.linspace(np.min(taolog),np.max(taolog),100)
    # y1 = model(x,par[0],par[1])

    # plt.plot(x,y1,color  = 'orange')
    
    
    # par0  = [1,0]
    # par, cov = fitter.curve_fit(model, taolog[q2], alphalog[q2], par0)

    # y2 = model(x,par[0],par[1])

    # plt.plot(x,y2,color  = 'orange')
    
    # plt.show()
    
    q1 = np.where(taolog[:-1]>0.3)
    q2 = np.where(taolog[:-1]<-1.5)
    
    alow = np.average(direv[q1])
    ahigh = np.max(direv[q2])
    
    diff = ahigh-alow
    
    q = np.where(direv-alow<diff/2)
    qq = np.where(taolog>-1.8)
    qincep = np.where(q[0]>np.min(qq[0]))
    index = q[0][qincep][0]
    
    criticalpoint.append(tao_list[u][index])
    

    plt.plot(taolog[:-1],direv)
    plt.title('First Derivative')
    plt.xlabel('Log(Tao)')
    plt.ylabel('Log(Alpha)')
    
    yy = np.linspace(np.min(direv),np.max(direv),100)
    xx = np.ones(100)*taolog[index]
    
    plt.plot(xx,yy)
    
    plt.plot(np.linspace(np.min(taolog),np.max(taolog),100),np.ones(100)*alow)
    plt.plot(np.linspace(np.min(taolog),np.max(taolog),100),np.ones(100)*ahigh)
    
    plt.show()
    
    # plt.plot(taolog[5:-6],direv_diff)
    # plt.title('Second Derivative')
    # plt.xlabel('Log(Tao)')
    # plt.ylabel('Log(Alpha)')
    # plt.show()
    
    # plt.plot(taolog[12:-13],direv_diff2)
    # plt.title('Third Derivative')
    # plt.xlabel('Log(Tao)')
    # plt.ylabel('Log(Alpha)')
    # plt.show()

lowthickness_list = np.array(lowthickness_list)
meanthickness_list = np.array(meanthickness_list)
medianthickness_list = np.array(medianthickness_list)
criticalpoint = np.array(criticalpoint)

    
#%%
plt.scatter(meanthickness_list,criticalpoint)
plt.title('Critical Point vs Average Current Sheet Temporal Thickness')
plt.xlabel('CS Thickness')
plt.ylabel('Critical Tao Value')

def model(x,a,b):
    return a*x**b

par0  = [1,1]
par, cov = fitter.curve_fit(model, meanthickness_list, criticalpoint, par0)

x = np.linspace(np.min(meanthickness_list),np.max(meanthickness_list),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

print(par)
plt.show()

plt.scatter(medianthickness_list,criticalpoint)
plt.title('Critical Point vs Median Current Sheet Temporal Thickness')
plt.xlabel('CS Thickness')
plt.ylabel('Critical Tao Value')


par0  = [1,1]
par, cov = fitter.curve_fit(model, medianthickness_list, criticalpoint, par0)

x = np.linspace(np.min(medianthickness_list),np.max(medianthickness_list),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

print(par)
plt.show()

ind = np.argmax(lowthickness_list)
r = np.array(range(len(lowthickness_list)))
q = np.where(np.logical_not(r == ind))


plt.scatter(lowthickness_list[q],criticalpoint[q])
plt.title('Critical Point vs 10% Thinnest Current Sheet Temporal Thickness')
plt.xlabel('CS Thickness')
plt.ylabel('Critical Tao Value')


par0  = [1,1]
par, cov = fitter.curve_fit(model, lowthickness_list[q], criticalpoint[q], par0)

x = np.linspace(np.min(lowthickness_list[q]),np.max(lowthickness_list[q]),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

print(par)
plt.show()
#%%
betainclude = np.array(betainclude)
q = np.where(betainclude)
betalist = np.array(betalist)
ebetalist = np.array(ebetalist)
betaslist = np.array([betalist,ebetalist])
betalabels = ['Proton Beta','Electron Beta']

for u in [0,1]:

    plt.scatter(betaslist[u][q],meanthickness_list[q],label = 'mean cs thickness')
    plt.scatter(betaslist[u][q],criticalpoint[q],color = 'green', label = 'criitcal point')
    
    plt.title('Tao Values vs. '+betalabels[u])
    plt.xlabel('Beta')
    plt.ylabel('Tao')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.show()


#%%
Brms_array = np.concatenate(Brms_list)
B0_array = np.concatenate(B0_list) 
beta_array = np.concatenate(betainterval)
ebeta_array = np.concatenate(ebetainterval)
Brms_norm = [Brms_array[i]/B0_array[i] for i in range(len(Brms_array))]
plt.scatter(beta_array,Brms_norm)
plt.title('Normalized Brms vs Proton Beta')
plt.xlabel('Beta')
plt.ylabel('Brms/B0')
plt.xscale('log')
plt.yscale('log')
        
def model(x,a,b):
    return a*x**b

par0  = [1,1]
par, cov = fitter.curve_fit(model, beta_array, Brms_norm, par0)

x = np.linspace(np.min(beta_array),np.max(beta_array),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

print(par)

plt.show()

plt.scatter(ebeta_array,Brms_norm)
plt.title('Normalized Brms vs Electron Beta')
plt.xlabel('Beta')
plt.ylabel('Brms/B0')
plt.xscale('log')
plt.yscale('log')
        
def model(x,a,b):
    return a*x**b

par0  = [1,1]
par, cov = fitter.curve_fit(model, ebeta_array, Brms_norm, par0)

x = np.linspace(np.min(ebeta_array),np.max(ebeta_array),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

print(par)

plt.show()
    
    
