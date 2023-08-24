# -*- coding: utf-8 -*-
"""
Created on Sun May 28 20:00:19 2023

@author: jackm
"""

#IMPORTS
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
import math


shear_angle_array = []
J0_array = []
Ja_array = []
Jth_array = []
inertial_length_array = []
lamda_array = []
lamda2_array = []
nlangle = []
Bxav_array = []
deltaBx_array = []
deltaB_array = []
Bav_array = []
deltaBmax_array = []
Jpeak_array = []
beta_array = []
ebeta_array = []
vn_array = []
v_array = []
dist_array = []
eTemp_array1 = []
eTemp_array2 = []
eTemp_array = []
eTemp_valid = []
date_array = []
temp_array = []
exempt = 0
exempt_array = []
Bzdiff_array = []
Bzav_array = []
dens_array = []
vd_array = []
cs_array = []

tbounds = []
transitionbounds = []
dates = []
with open('/Users/kimberleymawhinney/Desktop/Internship/currentsheets.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if not row[0] == 'Date':
            dates.append(row[0])
            tbounds.append([float(row[1]),float(row[2])])
            transitionbounds.append([float(row[3]),float(row[4])])
            
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

plot = False
csnum = len(dates)
date = ''
for index in range(csnum):
#for index in qnique[0]:
    if not date == dates[index]:
        date = dates[index]
        year = int(date[:4])
        month = int(date[5:7])
        day = int(date[8:])
        print(year,month,day)
        
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
            #Obtaining Ion Velocity Data
            a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Fields/psp_fld_l2_mag_RTN_202303'+st+'_v02.cdf')
            t = (a['epoch_mag_RTN'].to_numpy() - epoch)*10**-9
            B = a['psp_fld_l2_mag_RTN'].to_numpy()
        else:
            tcal = pytplot.time_double(date+' 0:00:00')
        
            pyspedas.psp.fields(trange=[tstring1,tstring2],time_clip=True)
            t, B = pytplot.get_data('psp_fld_l2_mag_RTN').times-tcal, pytplot.get_data('psp_fld_l2_mag_RTN').y
        

        epoch = cdflib.cdfepoch.compute_tt2000([year,month,day,0,0,0])
        #Obtaining Ion Velocity Data
        dstring = date.replace('-','')
        a = cdflib.cdf_to_xarray('/Users/kimberleymawhinney/Desktop/Internship/CDF/Sweap/psp_swp_spi_sf00_L3_mom_'+dstring+'_v04')
        tt = (a['Epoch'].to_numpy() - epoch)*10**-9
        
        #Obtaining Temp Data for specific day
        temp = (a['TEMP']).to_numpy()
        nanfinder = np.isnan(temp)
        q = np.where(nanfinder == False)
        q1 = np.where(temp[q]>0)
        temp = temp[q][q1]
        
        tt = tt[q][q1]
        
        vi_sun = a['VEL_RTN_SUN'].to_numpy()
        vsc_sun = a['SC_VEL_RTN_SUN'].to_numpy()
        vi = np.array(vi_sun - vsc_sun)[q][q1]
        
        
        
        dist  = a['SUN_DIST'].to_numpy()[q][q1]
        
        
        
        
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
                
        
                
        #GETTING ELECTRON TEMP 
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
  
#FINDING DATA FOR t,B and Ion velocity in Timeframe of Current Sheet
    halfbounds = (transitionbounds[index][1]-transitionbounds[index][0])/2
    if month == 9:     
        if transitionbounds[index][0]-tbounds[index][0]<halfbounds:
            tbounds[index][0] = transitionbounds[index][0]-halfbounds
        if tbounds[index][1] - transitionbounds[index][1]<halfbounds:
            tbounds[index][1] = transitionbounds[index][1]+halfbounds

    tcs = []
    Bcs = []
    for i in range(len(t)):
        if t[i]>tbounds[index][0] and t[i]<tbounds[index][1]:
            tcs.append(t[i])
            Bcs.append(B[i])
    
            
    midd = np.average(tcs)
    ttdiff = np.abs([ttt - midd for ttt in tt])
    windex = np.argmin(ttdiff)
    vics = vi[windex]
    distcs = dist[windex]/695700
    
    
    T = temp[windex]
        
    etdiff1 = np.abs([ee-midd for ee in et1])
    windex1 = np.argmin(etdiff1)
    eT1 = eTemp1[windex1]
    
        
    etdiff2 = np.abs([ee-midd for ee in et2])
    windex2 = np.argmin(etdiff2)
    eT2 = eTemp2[windex2]
    
    
    if etdiff1[windex1]<600:
        eT = eT1
        
        eTemp_v = True
    elif etdiff2[windex2]<600:
        eT = eT2
        eTemp_v = True
    elif etdiff1[windex1]<etdiff2[windex2]:
        eT = eT1
        eTemp_v = False
    else:
        eT = eT2
        eTemp_v = False
        

    
    #Defining Bl and Br on the border of the current sheet    
    Bcs = np.array(Bcs)
    tcs = np.array(tcs)
    q1 = np.where(tcs<transitionbounds[index][0])
    q2 = np.where(tcs>transitionbounds[index][1])
    Bl = [np.average(Bcs[q1,i]) for i in [0,1,2]]
    Br = [np.average(Bcs[q2,i]) for i in [0,1,2]]
    
    def mag(v):
        return(np.sqrt(np.sum([a**2 for a in v])))
    
    #Determining Normal Vector,Shear Angle, Ion Speed Across Normal
    u = np.cross(Bl,Br)
    n = u/mag(u)
    shear_angle = np.arccos(np.dot(Bl,Br)/(mag(Bl)*mag(Br)))*180/np.pi
    
    vn = vics[0]*n[0]+vics[1]*n[1]+vics[2]*n[2]
    

    #MVA to find l vector
    M = []
    for i in [0,1,2]:
        M.append([])
        for j in [0,1,2]:
            sum = 0
            for z in range(len(Bcs)):
                sum += Bcs[z][i]*Bcs[z][j]
            m = sum/len(Bcs) - np.average(Bcs[:,i])*np.average(Bcs[:,j])
            M[i].append(m)
    
    w,v = np.linalg.eig(M)     
    
    max = 0  
    for i in range(len(w)):
        if w[i]>max:
            dim = i
            max = w[i]
    l = v[:,dim]

    
    #Determing Orthonormal Vector Basi
    x1 = l - n*np.dot(l,n)
    x = x1/mag(x1)
    if np.dot(Bcs[0],x)>0:
        x = -x
    z = n
    y = np.cross(z,x)
    
    nl = np.arccos(np.dot(n,l))*180/np.pi
    
    Bx = [np.dot(b,x) for b in Bcs]
    By = [np.dot(b,y) for b in Bcs]
    Bz = [np.dot(b,z) for b in Bcs]
    
    #Bzdiff = np.max(Bz)-np.min(Bz)
    Bzdiff = np.max(np.abs(Bz))
    Bzav = np.average(np.abs(Bz))
   
    
    Bxl = np.dot(Bl,x)       
    Bxr = np.dot(Br,x)
    Bxav = (Bxl+Bxr)/2
    deltaBx = np.abs(Bxr-Bxl)
    
    assym = Bxav/deltaBx
    
    deltaB = np.abs(mag(Bl)-mag(Br))
    
    Bmag = [mag(b) for b in Bcs]
    deltaBmax = np.max(Bmag)-np.min(Bmag)
    
    Bav = (mag(Bl)+mag(Br))/2
    
    
    Bav2 = np.average(Bmag)
    
    #Finding Central Region
    fac = 0.25
    c = int(len(tcs)/2)
    count = 0
    while(True):
        if np.abs(Bx[c]-Bxav)<fac/2*deltaBx:
            break
        else:
            count+=1
            c = c + count*(-1)**count
            if c==0 or c==len(tcs):
                print('Error finding current sheet center for cs {}'.format(index))
    center = []
    cindex = []
    e = c
    while np.abs(Bx[e-1]-Bxav)<fac*deltaBx:
        e += -1
    center.append(tcs[e])
    cindex.append(e)
    e = c
    while np.abs(Bx[e+1]-Bxav)<fac*deltaBx:
        e +=1
    cindex.append(e)
    center.append(tcs[e])
    
    
    
    #Computing Current Density
    dBx = []
    dBy = []
    dBx.append((Bx[1] - Bx[0])/(tcs[1] - tcs[0]))
    dBy.append((By[1] - By[0])/(tcs[1] - tcs[0]))
    
    for i in range(1,len(Bx)-1):
        dBx.append((Bx[i+1] - Bx[i-1])/(tcs[i+1]-tcs[i-1]))
        dBy.append((By[i+1] - By[i-1])/(tcs[i+1]-tcs[i-1]))
    
    dBx.append((Bx[-1] - Bx[-2])/(tcs[-1] - tcs[-2]))
    dBy.append((By[-1] - By[-2])/(tcs[-1] - tcs[-2]))
    
    Jx = np.array(dBy)*(10**4/(4*np.pi)*1/vn)
    Jy = np.array(dBx)*(-10**4/(4*np.pi)*1/vn)
    J0 = np.abs(np.average(Jy[cindex[0]:cindex[1]]))
    
    if vn<0:
        Jpeak = np.max(Jy[cindex[0]:cindex[1]])
        if not Jpeak>0:
            Jpeak = -np.min(Jy[cindex[0]:cindex[1]])
    if vn>0:
        Jpeak = -np.min(Jy[cindex[0]:cindex[1]])
        if not Jpeak>0:
                Jpeak = np.max(Jy[cindex[0]:cindex[1]])
    
    #CURRENT SHEET THICKNESS
    lamda = 10**4/(4*np.pi)*deltaBx/(2*J0)
    
    lamda2 = np.abs((transitionbounds[index][1]-transitionbounds[index][0])*vn/2)
    
    
    
    #Spacial Profile
    min = deltaB+Bav
    for i in range(len(tcs)):
        if np.abs(Bx[i]-Bxav)<min:
            min = np.abs(Bx[i]-Bxav)
            cent = i
    
    z1cs = tcs*-vn
    zcs = z1cs - z1cs[cent]
    
    #Ion Inertail Length, Alfven current density/velocity, and Beta
    
    #Finding closest data point of density to data
    midd = np.average(tcs)
    tdensdiff = np.abs([ttt - midd for ttt in tdens])
    windex = np.argmin(tdensdiff)
    Ncs = dens[windex]
    
    vd = np.array(J0)/(0.16*Ncs)
    cs = 434*np.sqrt(eT/1000)/np.sqrt(2)
    #Found closest temp measurent with ion velocity above
    '''
    tempcs = []
    min = 
    for i in range(len(tt)):
        if tt[i]>tbounds[index][0] and tt[i]<tbounds[index][1]:
            tempcs.append(temp[i])
    T = np.average(tempcs)
    '''
    beta = 0.48*Ncs*T/(Bav**2)
    
    
    ebeta = 0.48*Ncs*eT/(Bav**2)
    
        
    
    wpi = 2*np.pi*9000/np.sqrt(1836)*np.sqrt(Ncs)
    inertial_length = 3*10**5/wpi
    
    
    Va = 434*Bav/(20*np.sqrt(Ncs))
    Ja = 0.16*np.sqrt(Ncs)*434*Bav/20
    
    
    Jth = Ja*np.sqrt(beta)
    
    #print(index, lamda/inertial_length)
    
    
    '''
    if lamda/inertial_length<0.15:
        print(index, ':', lamda/inertial_length, ',', lamda , ',', inertial_length, 'Dens =',Ncs)
        plot = True
    else:
        plot = False
    '''
   # print(beta,np.array(tbounds[index])/3600,Ncs,T)

    append = True
    
    if lamda/inertial_length >500:
        append =False
    if np.abs(assym)>0.5:
        append = False
    if np.abs(nl-90)>10:
        append = False
    if beta<0.0005:
        append = False
    if Bzdiff>deltaBx/2:
        append=False
    for i in range(index+1,np.min([index+5,csnum])):
        if np.abs(tbounds[i][0]-tbounds[index][0])<halfbounds/8 and np.abs(tbounds[i][1]-tbounds[index][1])<halfbounds/8:
            append = False
    #plot = not append
    if append:
        eTemp_array2.append(eT2)
        eTemp_array1.append(eT1)  
        dist_array.append(distcs)
        eTemp_valid.append(eTemp_v)
        eTemp_array.append(eT)
        deltaBx_array.append(deltaBx)
        Bxav_array.append(Bxav)
        nlangle.append(nl)
        v_array.append(mag(vics))
        vn_array.append(vn)
        shear_angle_array.append(shear_angle)
        Bav_array.append(Bav)
        deltaBmax_array.append(deltaBmax)
        deltaB_array.append(deltaB)
        lamda2_array.append(lamda2)
        lamda_array.append(lamda)
        Jpeak_array.append(Jpeak)
        J0_array.append(J0)
        beta_array.append(beta)
        ebeta_array.append(ebeta)
        inertial_length_array.append(inertial_length)
        Ja_array.append(Ja)
        Jth_array.append(Jth)
        date_array.append(date)
        temp_array.append(T)
        exempt_array.append(exempt)
        Bzdiff_array.append(Bzdiff/deltaBx)
        Bzav_array.append(Bzav/deltaBx)
        cs_array.append(cs)
        vd_array.append(vd)
        dens_array.append(Ncs)
    else:
        exempt+=1
    if plot:
    
        #PLOTTING CS with central Region Markers
        #plt.plot(tcs,Bx)
        #plt.show()

        
        # plt.plot(np.ones(100)*center[0],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.plot(np.ones(100)*center[1],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxav)
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxl)
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxr)
        # plt.title('Bx with Marked Central Region and Boundries, CS{}'.format(index+1))
        # plt.ylabel('B [nT]')
        # plt.xlabel('Time [s]')
        # plt.show()
        
        # plt.plot(np.ones(100)*tcs[cindex[0]],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.plot(np.ones(100)*tcs[cindex[1]],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxav)
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxl)
        # plt.plot(np.linspace(tcs[0],tcs[-1],100),np.ones(100)*Bxr)
        # plt.ylabel('B [nT]')
        # plt.xlabel('t [s]')
        # plt.title(' Bx Group{} , CS{}'.format(date,index))
        
        # plt.plot(np.ones(100)*tcs[q1[0][-1]],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.plot(np.ones(100)*tcs[q2[0][0]],np.linspace(-deltaBx/2+Bxav,deltaBx/2+Bxav,100))
        # plt.show()
        # #print(tbounds[index])
        
        # plt.plot(zcs,Bz)
        # plt.ylabel('B [nT]')
        # plt.xlabel('z [km]')
        # plt.title(' Bz Group{} , CS{}'.format(date,index))
        # plt.show()
        
        # #Plotting Current Density
        # plt.plot(tcs,Jx)
        # plt.title('Jx')
        # plt.xlabel('Time')
        # plt.show()
        
        # plt.plot(zcs,Jy)
        # plt.plot(np.ones(100)*zcs[cindex[0]],np.linspace(np.min(Jy),np.max(Jy),100))
        # plt.plot(np.ones(100)*zcs[cindex[1]],np.linspace(np.min(Jy),np.max(Jy),100))
        # plt.title('Jy Group{} , CS{}'.format(date,index))
        # plt.xlabel('z [km]')
        # plt.show()
        
    
        #Plotting Bfield with xyz coordinates, spacial profile
        plt.plot(zcs,Bx, label = 'Bx')
        #plt.plot(zcs,By, label = 'By')
        plt.plot(zcs,Bz, label = 'Bz')
        plt.title('B field in xyz coordinates, spacial profile')
        plt.ylabel('B [nT]')
        plt.xlabel('z [km]')
        plt.legend()
        plt.show()
        
csnum = csnum-exempt
betas_array = np.array([beta_array,ebeta_array])
beta_labels = ['Proton Beta', 'Electron Beta']
#%%
#MAKING PLOTS

#plot_title = ': 09/07/2022 0:00-6:00'
plot_title = ': ALL INTERVALS'

csnum = len(lamda_array)
norm_thickness = [lamda_array[i]/inertial_length_array[i] for i in range(csnum)]
plt.scatter(norm_thickness,np.array(shear_angle_array),s=10)
#plt.xlim(0,5)
plt.xscale("log")
plt.yscale("log")
plt.title('Shear Angle vs. Normalized Thickness'+plot_title)
plt.xlabel('Lamda/Lamda p ')
plt.ylabel('delta theta')
plt.show()

norm_current = [J0_array[i]/Ja_array[i] for i in range(csnum)]
plt.scatter(norm_thickness,norm_current,s=10)
plt.title('Normalized Current vs. Normalized Thickness'+plot_title)
plt.xlabel('Lamda/Lamda p ')
plt.ylabel('J0/Ja')
plt.yscale('log')
plt.xscale('log')
plt.show()

plt.scatter(inertial_length_array,lamda_array,s=3)
plt.yscale('log')
plt.xscale('log')
plt.title('Thickness vs Ion Inertial Length'+plot_title)
plt.xlabel('Lamda p')
plt.ylabel('Lamda')
plt.show()


gyroradius = [inertial_length_array[i]*np.sqrt(beta_array[i]) for i in range(csnum)]
plt.scatter(gyroradius,lamda_array,s=3)
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.1)
plt.ylim(1,100)
plt.title('Thickness vs Ion Gyroradius'+plot_title)
plt.xlabel('Roh p')
plt.ylabel('Lamda')
plt.show()


plt.hist(norm_thickness,bins=5000)
plt.xlim(0,5)
plt.title('Histogram of Normalized Thickness (Equation)'+plot_title)
plt.xlabel('Lamda/Lamda p')
plt.show()

norm_thickness_gyro = [lamda_array[i]/gyroradius[i] for i in range(csnum)]
plt.hist(norm_thickness_gyro,bins=3000)
plt.title('Histogram of Normalized Thickness with gyroradius'+plot_title)
plt.xlabel('Lamda/wp')
plt.xlim(0,20)
plt.show()

norm_thickness2 = [lamda2_array[i]/inertial_length_array[i] for i in range(csnum)]
plt.hist(norm_thickness2,bins=100)
plt.xlim(0,5)
plt.title('Histogram of Normalized Thickness (Selected Bounds)'+plot_title)
plt.xlabel('Lamda/Lamda p')
plt.show()

plt.hist(nlangle,bins=20)
plt.title('Histogram of Angle between n,l vectors'+plot_title)
plt.xlabel('Alpha_nl')
plt.show()

assymx_array = np.sort([Bxav_array[i]/deltaBx_array[i] for i in range(csnum)])
plt.subplot(2,1,1)
plt.title('Histogram/CDF of Current Sheet Assymetry, x axis'+plot_title)
plt.hist(assymx_array,bins = 20)
x = np.linspace(assymx_array[0],assymx_array[-1],100)
y = [np.where(assymx_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlabel('<Bx>/deltaBx')
plt.show()

assym_array = np.sort([deltaB_array[i]/Bav_array[i] for i in range(csnum)])
plt.subplot(2,1,1)
plt.title('Histogram of Current Sheet Assymetry, magnitude'+plot_title)
plt.hist(assym_array,bins = 200)
x = np.linspace(assym_array[1],assym_array[-1],100)
y = [np.where(assym_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlabel('deltaB/<B>')
plt.show()

assymax_array = np.sort([deltaBmax_array[i]/Bav_array[i] for i in range(csnum)])
plt.subplot(2,1,1)
plt.title('Histogram of Current Sheet Assymetry, max diff'+plot_title)
plt.hist(assymax_array,bins = 200)
x = np.linspace(assymax_array[0],assymax_array[-1],100)
y = [np.where(assymax_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlabel('deltaBmax/<B>')
plt.show()

Jpeak1_array = np.sort(Jpeak_array)
plt.subplot(2,1,1)
plt.title('Historgram/CDF of Jpeak'+plot_title)
plt.xlim(0,12000)
plt.hist(Jpeak_array,bins = 800)
x = np.linspace(Jpeak1_array[0],Jpeak1_array[-1],1000)
y = [np.where(Jpeak1_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlim(0,12000)
plt.xlabel('Jpeak')
plt.show()


plt.subplot(2,1,1)
J01_array = np.sort(J0_array)
plt.title('Historgram/CDF of J0'+plot_title)
plt.hist(J0_array,bins =600)
plt.xlim(0,10000)
x = np.linspace(J01_array[0],J01_array[-1],1000)
y = [np.where(J01_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlim(0,10000)
plt.xlabel('J0')
plt.show()

normJpeak_array = np.sort([Jpeak_array[i]/Ja_array[i] for i in range(csnum)])
plt.subplot(2,1,1)
plt.title('Historgram/CDF of Normalized Jpeak'+plot_title)
plt.hist(normJpeak_array,bins = 500)
plt.xlim(0,1)
x = np.linspace(normJpeak_array[0],normJpeak_array[-1],1000)
y = [np.where(normJpeak_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlim(0,1.5)
plt.xlabel('Jpeak/Ja')
plt.show()

normJ0_array = np.sort([J0_array[i]/Ja_array[i] for i in range(csnum)])
plt.subplot(2,1,1)
plt.title('Historgram/CDF of Normalized J0'+plot_title)
plt.hist(normJ0_array,bins = 500)
plt.xlim(0,1)
x = np.linspace(normJ0_array[0],normJ0_array[-1],1000)
y = [np.where(normJ0_array<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlim(0,1)
plt.xlabel('J0/Ja')
plt.show()

plt.scatter(Ja_array,Jpeak_array,s=3)
plt.title('Jpeak vs Ja'+plot_title)
plt.xlabel('Ja')
plt.ylabel('Jpeak')
plt.yscale('log')
plt.show()

plt.scatter(Ja_array,J0_array,s=3)
plt.title('J0 vs Ja'+plot_title)
plt.xlabel('Ja')
plt.ylabel('J0')
plt.yscale('log')
plt.show()

plt.scatter(Jth_array,Jpeak_array,s=3)
plt.title('Jpeak vs Jth'+plot_title)
plt.xlabel('Jth')
plt.ylabel('Jpeak')
plt.yscale('log')
plt.xscale('log')
plt.xlim(500)
plt.show()

plt.scatter(Jth_array,J0_array,s=3)
plt.title('J0 vs Jth'+plot_title)
plt.xlim(500)
plt.xlabel('Jth')
plt.ylabel('J0')
plt.yscale('log')
plt.xscale('log')
plt.show()

plt.scatter(beta_array,np.abs(assymx_array),s=2)
plt.title('Current Sheet Assymettry vs Proton Beta')
plt.xlabel('Beta')
plt.ylabel('CS Assymettry')
plt.xscale('log')
#%%
logbeta = np.log10(beta_array)
plt.subplot(2,1,1)

#plt.xticks(ticks =[0.001,0.05,0.01,0.05,0.1,0.5])
q  = np.where(np.array(beta_array)>0.002)
plt.hist(logbeta[q],bins=20)
logbeta1 = np.sort(np.array(logbeta)[q])
plt.title('Beta Histogram'+plot_title)
x = np.linspace(logbeta1[0],logbeta1[-1],1000)
y = [np.where(logbeta1<=xx)[0][-1]/csnum for xx in x]
plt.subplot(2,1,2)
plt.plot(x,y)
plt.xlabel('log(Beta)')
plt.show()
#%%
plt.hist(dist_array,bins=20)
plt.title('Histogram of CS Distance from Sun')
plt.xlabel('Distance (Solar Radii)')
plt.show()
#%%
#NORMALIZED THICKNESS VS PROTON BETA
norm_thickness = np.array([lamda_array[i]/inertial_length_array[i] for i in range(csnum)])
for u in [0,1]:
    plt.scatter(betas_array[u],norm_thickness,s=1)
    #plt.xlim(0,5)
    plt.xscale("log")
    plt.yscale('log')
    #plt.yscale("log")
    plt.title('Normalized Thickness vs. '+beta_labels[u])
    plt.xlabel('Beta')
    plt.ylabel('Lamda/Lamda p')
    
    
    x10 = np.linspace(np.log10(np.min(betas_array[u])),np.log10(np.max(betas_array[u])),10)
    xdiff = x10[1]-x10[0]
    xb = 10**x10
    xav = 10**(x10[1:]-xdiff)
    yav = []
    for i in range(len(xb)-1):
        q1 = np.where(betas_array[u]>xb[i])
        q2 = np.where(betas_array[u]<xb[i+1])
        q = np.intersect1d(q1,q2)
        yav.append(np.median(norm_thickness[q]))
        
    plt.scatter(xav,yav)
    
    
    def model(x,a,b):
        return a*x**b
    
    par0  = [1,1]
    
    def model1(x,a,b):
        return b*x + a
    
    #par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

    par,cov = fitter.curve_fit(model1,np.log10(betas_array[u]),np.log10(norm_thickness),par0)
    
    par[0] = 10**par[0]
    
    
    x = np.linspace(np.min(betas_array[u]),np.max(betas_array[u]),100)
    y = model(x,par[0],par[1])
    
    plt.plot(x,y,color  = 'orange')
    
    print(par)
    
    plt.show()
    

#%%
#NORMALIZED THICKNESS (LOWEST PERCENTILES) VS PROTON BETA

norm_thickness = np.array([lamda_array[i]/inertial_length_array[i] for i in range(csnum)])

#q8 = np.where(norm_thickness<np.percentile(norm_thickness,10))
q9 = []
for u in [0,1]:
    x10 = np.linspace(np.log10(np.min(betas_array[u])),np.log10(np.max(betas_array[u])),30)
    beta_bins = 10**x10
    for i in range(len(beta_bins)-1):
        q = np.where(np.logical_and(betas_array[u]>beta_bins[i],betas_array[u]<beta_bins[i+1]))
        if len(q[0]>0):
            q1 = np.where(norm_thickness<np.percentile(norm_thickness[q],10))
        q2 = np.intersect1d(q,q1)
        if len(q2>0):
            q9.append(q2)
    q8 = (np.concatenate(q9),)
    plt.scatter(betas_array[u][q8],norm_thickness[q8],s=1)
    #plt.xlim(0,5)
    plt.xscale("log")
    plt.yscale('log')
    #plt.yscale("log")
    plt.title('Normalized Thickness (Lowest 10 percentile) vs. '+beta_labels[u])
    plt.xlabel('Beta')
    plt.ylabel('Lamda/Lamda p')
    
    
    x10 = np.linspace(np.log10(np.min(betas_array[u])),np.log10(np.max(betas_array[u])),10)
    xdiff = x10[1]-x10[0]
    xb = 10**x10
    xav = 10**(x10[1:]-xdiff)
    yav = []
    for i in range(len(xb)-1):
        q1 = np.where(betas_array[u][q8]>xb[i])
        q2 = np.where(betas_array[u][q8]<xb[i+1])
        q = np.intersect1d(q1,q2)
        yav.append(np.median(norm_thickness[q8][q]))
        
    plt.scatter(xav,yav)
    
    
    def model(x,a,b):
        return a*x**b
    
    par0  = [1,1]
    
    def model1(x,a,b):
        return b*x + a
    
    #par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

    par,cov = fitter.curve_fit(model1,np.log10(betas_array[u][q8]),np.log10(norm_thickness[q8]),par0)
    
    par[0] = 10**par[0]
    
    
    x = np.linspace(np.min(betas_array[u]),np.max(betas_array[u]),100)
    y = model(x,par[0],par[1])
    
    plt.plot(x,y,color  = 'orange',label = 'Curve Fit')
    
    print(par)
    
    if u == 1:
    
        y = x**(5/18)
        
        plt.plot(x,y,color = 'red', label = 'Theoretical Prediction')
        
        plt.legend()

    
    
    plt.plot()
    
    plt.show()

#%%
#SPECIFIC REALATIONSHIP VS ELECTRON BETA
norm_thickness = np.array([lamda_array[i]/inertial_length_array[i] for i in range(csnum)])
beta_altered = [beta_array[i]*eTemp_array[i]/temp_array[i] for i in range(csnum)]
plt.scatter(beta_altered,norm_thickness,s=1)
#plt.xlim(0,5)
plt.xscale("log")
plt.yscale('log')
plt.xlim(0.001)
#plt.yscale("log")
plt.title('Normalized Thickness vs. Proton Beta (with Temperature Adjustment)')
plt.xlabel('Beta')
plt.ylabel('Lamda/Lamda p')


x10 = np.linspace(np.log10(np.min(beta_altered)),np.log10(np.max(beta_altered)),10)
xdiff = x10[1]-x10[0]
xb = 10**x10
xav = 10**(x10[1:]-xdiff)
yav = []
for i in range(len(xb)-1):
    q1 = np.where(beta_altered>xb[i])
    q2 = np.where(beta_altered<xb[i+1])
    q = np.intersect1d(q1,q2)
    yav.append(np.median(norm_thickness[q]))
    
plt.scatter(xav,yav)



def model(x,a,b):
    return a*x**b

par0  = [1,1]

def model1(x,a,b):
    return b*x + a

#par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

par,cov = fitter.curve_fit(model1,np.log10(beta_altered),np.log10(norm_thickness),par0)

par[0] = 10**par[0]

print(par)
    
    
x = np.linspace(np.min(beta_altered),np.max(beta_altered),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

plt.show()

#%%
#NORMALIZED CURRENT (Jpeak) VS PROTON BETA
Jpeak_array = np.array(Jpeak_array)
qq = np.where(np.logical_not(Jpeak_array<0))
norm_current = np.array([Jpeak_array[i]/Ja_array[i] for i in range(csnum)])[qq]
for u in [0,1]:
    plt.scatter(betas_array[u][qq],norm_current,s=1)
    #plt.xlim(0,5)
    plt.xscale("log")
    plt.yscale('log')
    #plt.yscale("log")
    plt.title('Normalized Jpeak vs. '+beta_labels[u])
    plt.xlabel('Beta')
    plt.ylabel('Current')
    
    
    x10 = np.linspace(np.log10(np.min(betas_array[u][qq])),np.log10(np.max(betas_array[u][qq])),10)
    xdiff = x10[1]-x10[0]
    xb = 10**x10
    xav = 10**(x10[1:]-xdiff)
    yav = []
    for i in range(len(xb)-1):
        q1 = np.where(betas_array[u][qq]>xb[i])
        q2 = np.where(betas_array[u][qq]<xb[i+1])
        q = np.intersect1d(q1,q2)
        yav.append(np.median(norm_current[q]))
        
    plt.scatter(xav,yav)
    
    
    def model(x,a,b):
        return a*x**b
    
    par0  = [1,1]
    
    def model1(x,a,b):
        return b*x + a
    
    #par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

    par,cov = fitter.curve_fit(model1,np.log10(betas_array[u][qq]),np.log10(norm_current),par0)
    
    par[0] = 10**par[0]
    
    
    
    x = np.linspace(np.min(betas_array[u][qq]),np.max(betas_array[u][qq]),100)
    y = model(x,par[0],par[1])
    
    plt.plot(x,y,color  = 'orange')
    
    plt.show()
    
    print(par)
#%%
#NORMALIZED CURRENT (J0) VS PROTON BETA
norm_current = np.array([J0_array[i]/Ja_array[i] for i in range(csnum)])
for u in [0,1]:
    plt.scatter(betas_array[u],norm_current,s=1)
    #plt.xlim(0,5)
    plt.xscale("log")
    plt.yscale('log')
    plt.title('Normalized J0 vs. '+beta_labels[u])
    plt.xlabel('Beta')
    plt.ylabel('Current')
    
    
    x10 = np.linspace(np.log10(np.min(betas_array[u])),np.log10(np.max(betas_array[u])),10)
    xdiff = x10[1]-x10[0]
    xb = 10**x10
    xav = 10**(x10[1:]-xdiff)
    yav = []
    for i in range(len(xb)-1):
        q1 = np.where(betas_array[u]>xb[i])
        q2 = np.where(betas_array[u]<xb[i+1])
        q = np.intersect1d(q1,q2)
        yav.append(np.median(norm_current[q]))
        
    plt.scatter(xav,yav)
    
    def model(x,a,b):
        return a*x**b
    
    par0  = [1,1]
    
    def model1(x,a,b):
        return b*x + a
    
    #par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

    par,cov = fitter.curve_fit(model1,np.log10(betas_array[u]),np.log10(norm_current),par0)
    
    par[0] = 10**par[0]
    
    print(par)
    
    x = np.linspace(np.min(betas_array[u]),np.max(betas_array[u]),100)
    y = model(x,par[0],par[1])
    
    plt.plot(x,y,color  = 'orange')
    
    plt.show()

#%%
beta_array = np.array(beta_array)
qlow = np.where(beta_array<0.01)
qmed1 = np.where(beta_array>0.01)
qmed2 = np.where(beta_array<0.1)
qmed = np.intersect1d(qmed1, qmed2)
qhigh = np.where(beta_array>0.1)

assymx_array = np.abs([Bxav_array[i]/deltaBx_array[i] for i in range(csnum)])
assym_array = np.array([deltaB_array[i]/Bav_array[i] for i in range(csnum)])
assymax_array = np.array([deltaBmax_array[i]/Bav_array[i] for i in range(csnum)])

assymx_split = [assymx_array[qq] for qq in [qlow,qmed,qhigh]]
assym_split = [assym_array[qq] for qq in [qlow,qmed,qhigh]]
assymax_split = [assymax_array[qq] for qq in [qlow,qmed,qhigh]]
labels = ['low','med','high']
normalization = [np.ones(len(assymx_split[i]))/len(assymx_split[i]) for i in range(3)]
plt.hist(assymx_split,weights = normalization,label = labels)
plt.title('Histogram of Current Sheet Assymettry in different Beta Regimes')
plt.xlabel('Bxav/deltaBx')
plt.legend()
plt.show()
plt.hist(assym_split,weights = normalization,label = labels,range = [0,0.1])
plt.title('Histogram of Current Sheet Compressability in different Beta Regimes')
plt.xlabel('Bav/deltaB')
plt.legend()
plt.show()
plt.hist(assymax_split,weights = normalization,label = labels,range = [0,0.1])
plt.title('Histogram of Current Sheet Compressibility (MAX) in different Beta Regimes')
plt.xlabel('Bav/deltaBmax')
plt.legend()
plt.show()

#%%

shear_angle_array = np.array(shear_angle_array)
shear_angle_split = [shear_angle_array[qq] for qq in [qlow,qmed,qhigh]]  
plt.hist(shear_angle_split,weights = normalization,label = labels)
plt.title('Histogram of Shear Angle in Different Beta')
plt.xlabel('Theta')
plt.show()

for u in [0,1]:
    plt.scatter(betas_array[u],shear_angle_array,s=1)
    #plt.xlim(0,5)
    plt.xscale("log")
    plt.yscale('log')
    #plt.yscale("log")
    plt.title('Shear Angle vs. '+beta_labels[u])
    plt.xlabel('Beta')
    plt.ylabel('Current')
    
    
    x10 = np.linspace(np.log10(np.min(betas_array[u])),np.log10(np.max(betas_array[u])),10)
    xdiff = x10[1]-x10[0]
    xb = 10**x10
    xav = 10**(x10[1:]-xdiff)
    yav = []
    for i in range(len(xb)-1):
        q1 = np.where(betas_array[u]>xb[i])
        q2 = np.where(betas_array[u]<xb[i+1])
        q = np.intersect1d(q1,q2)
        yav.append(np.median(shear_angle_array[q]))
        
    plt.scatter(xav,yav)
    
    
    def model(x,a,b):
        return a*x**b
    
    par0  = [1,1]
    
    def model1(x,a,b):
        return b*x + a
    
    #par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

    par,cov = fitter.curve_fit(model1,np.log10(betas_array[u]),np.log10(shear_angle_array),par0)
    
    par[0] = 10**par[0]
    
    x = np.linspace(np.min(betas_array[u]),np.max(betas_array[u]),100)
    y = model(x,par[0],par[1])
    
    plt.plot(x,y,color='orange')
    
    print(beta_labels[u],par)
    
    plt.show()

#%%

plt.scatter(eTemp_array1,eTemp_array2,s=3)
plt.ylabel('T (Jasper)')
plt.xlabel('T (QTN)')
plt.yscale('log')
plt.xscale('log')
plt.title('Comparison of Electron Temperature Measurments')
#%%
#PROTON BETA VS ELECTRON BETA PLOTS
#dates = np.array(dates)
#q = np.where(np.logical_not(np.logical_or(dates == '2021-08-10',dates=='2021-04-29')))
ebeta_array = np.array(ebeta_array)
plt.scatter(ebeta_array,beta_array,s=1)
plt.ylabel('Proton Beta')
plt.xlabel('Electron Beta')
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.001)
plt.ylim(0.001)
plt.title('Comparison of Proton and Electron Beta')

def model(x,a,b):
    return a*x**b

par0  = [1,1]
par, cov = fitter.curve_fit(model, ebeta_array, beta_array, par0)

x = np.linspace(np.min(ebeta_array),np.max(ebeta_array),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color  = 'orange')

plt.show()


plt.show()
#%%
#TESTING THEORETICAL RELATIONS

alpha = 0.26170057
gamma = 0.36645062
psi = 0.1014982



lamdanorm = np.array([lamda_array[i]/(inertial_length_array[i]*beta_array[i]**alpha) for i in range(csnum)])
shearanglenorm = np.array([shear_angle_array[i]/(beta_array[i]**gamma) for i in range(csnum)])
currentnorm = [J0_array[i]/(Ja_array[i]*beta_array[i]**psi) for i in range(csnum)]

q = np.where(lamdanorm>0.1)

plt.scatter(lamdanorm,shearanglenorm,s=1)

def model(x,a,b):
    return a*x**b

par0  = [1,1]

def model1(x,a,b):
    return b*x + a

#par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

#par,cov = fitter.curve_fit(model1,np.log10(lamdanorm[q]),np.log10(shearanglenorm[q]),par0)

#par[0] = 10**par[0]


par,cov = fitter.curve_fit(model,lamdanorm,shearanglenorm,par0)


print(par)

x = np.linspace(np.min(lamdanorm),np.max(lamdanorm),100)
y = model(x,par[0],par[1])

#y = model(x,par[0]/1.5,par[1]*2)

plt.plot(x,y,color='orange')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Normalized Thickness')
plt.ylabel('Nomralized Shear Angle')
plt.title('Shear Angle vs. CS Thickness, Normalized with respect to Beta Variance')
plt.show()

plt.scatter(lamdanorm,currentnorm,s = 1)

def model(x,a,b):
    return a*x**b

par0  = [1,1]

def model1(x,a,b):
    return b*x + a

#par,cov = fitter.curve_fit(model1,np.log10(xav),np.log10(yav),par0)

par,cov = fitter.curve_fit(model,lamdanorm,currentnorm,par0)
print(par)

x = np.linspace(np.min(lamdanorm),np.max(lamdanorm),100)
y = model(x,par[0],par[1])

plt.plot(x,y,color='orange')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Normalized Thickness')
plt.ylabel('Nomralized Current')
plt.title('Current vs. CS Thickness, Normalized with respect to Beta Variance')
plt.show()

#%%
norm_thickness = np.array([lamda_array[i]/inertial_length_array[i] for i in range(csnum)])
n = np.percentile(norm_thickness,1)
qnique1 = np.where(norm_thickness<n)
ex = np.array(exempt_array)[qnique1]
qnique = (np.array([qnique1[0][i]+ex[i] for i in range(len(ex))]),)
#%%
plt.hist(Bzdiff_array)
plt.title('Max Normalized Bz Magnitude Histogram')
plt.xlabel('Bz/Bx')
plt.show()
plt.hist(Bzav_array)
plt.title('Average Normalized Bz Magnitude Histogram')
plt.xlabel('Bz/Bx')
plt.show()
#%%
plt.scatter(cs_array,vd_array,s=1)
plt.title('Electron-ion drift versus the Ion-acoustic speed')
plt.xlabel('Acoustic Speed')
plt.ylabel('Drift Velocity')
plt.show()

ratiospeed = [vd_array[i]/cs_array[i] for i in range(csnum)]
ratiotemp = [eTemp_array[i]/temp_array[i] for i in range(csnum)]

plt.scatter(ratiotemp,ratiospeed,s=1)
plt.title('')
plt.xlabel('Te/Ti')
plt.ylabel('Vd/Cs')

mi = 1
me = 1

x = np.linspace(np.min(ratiotemp),np.max(ratiotemp),100)
y1 = 1 + np.sqrt((mi/me)*x**3)*np.exp(-(3/2)-x/2)  

plt.show()