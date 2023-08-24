#!/usr/bin/env python
# coding: utf-8
# %%
import numpy as np
import matplotlib.pyplot as plt

def b_init(zrange,zprecision):
    zlist = np.arange(0,zrange,zprecision)
    b = []
    for z in zlist:
        b.append(np.tanh(z))
    return b

def func(z1,z,vy,vz,step):
    indexz1 = int(np.abs(z1/step))
    indexz = int(np.abs(z/step))
    return vy**2+vz**2 - (vy - (Ub[indexz]-Ub[indexz1]))**2

def func1(z1,z,vy,vz,step):
    indexz1 = int(np.abs(z1/step))
    indexz = int(np.abs(z/step))
    a = vy**2+vz**2 - (vy - (Ub[indexz]-Ub[indexz1]))**2
    if a<0:
        return 0
    return np.sqrt(a)

def current_density(b,zrange,ed,n0,vrange,vinterval,zinterval,zprecision):
    global Ub
    Ub = []    
    index = range(0,50*int(1/zprecision))
    for i in index:
        area = 0
        for j in range(i):
            if j<len(b):
                area = area + b[j]*zprecision 
            else:
                area = area + zprecision
        Ub.append(area)
    if zrange<=2:
        z = np.arange(0,zrange,zinterval)
    else:
        z = np.concatenate((np.arange(0,2,zinterval),np.arange(2,zrange,zinterval*4)))
    J = []
    N = []
    vx = np.arange(0,vrange,vinterval)
    vy = np.arange(-vrange,vrange,vinterval)
    vz = np.arange(0,vrange,vinterval)
    for zz in z:
        #Iz = calc_Iz(zz,vy,vz,vinterval,zprecision)
        sum = 0
        sum1 = 0
        for j in range(len(vy)):
            for k in range(len(vz)):
                zbounds = []
                for n in [-1,1]:
                    zi = zz
                    while(True):
                        zi = zi + zprecision*n
                        if func(zi,zz,vy[j],vz[k],zprecision)<0:
                            z1 = zi - zprecision
                            break
                    zbounds.append(z1)
                if zbounds[0]*zbounds[1]<0:
                    zbounds[0] = 0
                if (zbounds[1]-zbounds[0])/100>zprecision:  
                    area = 0
                    zlist = np.arange(zbounds[0],zbounds[1]-zprecision,zprecision)
                    for zed in zlist:
                        area = area + func1(zed,zz,vy[j],vz[k],zprecision)*zprecision 
                elif zbounds[1]-zbounds[0]<zprecision/10:
                    area = 0
                else:
                    a = (zbounds[1]-zbounds[0])/100
                    fp = []
                    zp = []
                    zlist = np.arange(zbounds[0],zbounds[1],zprecision)
                    for zed in zlist:
                        fp.append(func1(zed,zz,vy[j],vz[k],zprecision))
                        zp.append(zed)
                    zp.append(zed+zprecision)
                    fp.append(-np.sqrt(-func(zbounds[1]+zprecision,zz,vy[j],vz[k],zprecision)))
                    area = 0
                    zlist = np.arange(zbounds[0],zbounds[1]+(zprecision+a),a)
                    f = np.interp(zlist,zp,fp)
                    count = 0
                    while zlist[count]<zbounds[1]:
                        area = area + f[count]*zprecision
                        count+=1
                    while(f[count]>0):
                        area = area + f[count]*zprecision
                        count +=1        
                iz = area/np.pi 
                    
                for i in range(len(vx)):
                    vsquare = vx[i]**2+vy[j]**2+vz[k]**2
                    a = vsquare-2*iz
                    if a>0:
                        sum += vy[j]*np.exp(-(np.sqrt(a)-ed)**2-2*iz)
                        sum1 += np.exp(-(np.sqrt(a)-ed)**2-2*iz)
        if zz<zinterval/10:
            c = n0/(sum1*4*vinterval**3)
        J.append(c*4*sum*vinterval**3)
        N.append(c*4*sum1*vinterval**3)
    zminus = []
    for i in range(len(z)-1):
        zminus.append(-z[len(z)-1-i])
    Jminus = J[::-1][:-1]
    zfull = np.concatenate((zminus,z))
    Jfull = np.concatenate((Jminus,J))
    plt.plot(zfull,Jfull)
    plt.xlabel('z')
    plt.ylabel('Current Density')
    plt.title('Current Density vs z $vrange$={}$ \epsilon_d$={}'.format(vrange,ed))
    plt.show()
    return J,N

def integrateJ(J,endindex,step):
    area = 0
    for i in range(endindex):
        area = area + J[i]*step 
    return area


b = b_init(5,0.02)
zprecision = 0.02
zrange = 5
vrange = 5
vinterval = 0.1
zinterval = 0.05
ed = 1
n0 = 1

allb = []
allJ = []
allN = []
iterations = 6
allb.append(b)
for n in range(iterations):

    J,N = current_density(b,zrange,ed,n0,vrange,vinterval,zinterval,zprecision)
    allJ.append(J)
    allN.append(N)
    zp = np.concatenate((np.arange(0,2,zinterval),np.arange(2,zrange,4*zinterval)))
    z = np.arange(0,zrange,zprecision)
    Jinterp = np.interp(z,zp,J)
    beta = 1/integrateJ(Jinterp,len(Jinterp),zprecision)
    b = []
    for i in range(len(Jinterp)):
        if n==0:
            b.append(beta*integrateJ(Jinterp,i,zprecision)) 
        else:
            b.append((beta*integrateJ(Jinterp,i,zprecision)+allb[n-1][i])/2) 
            
    plt.plot(z,b)
    allb.append(b)
    plt.xlabel('z')
    plt.ylabel('b(z)')
    plt.title('Calcultaed b for iteration {}'.format(n+1))
    plt.show()

for i in range(len(allb)):
    plt.plot(z,allb[i],label = '{}'.format(i+1))
    plt.xlabel('z')
    plt.ylabel('b(z)')
    plt.title('Calcultaed b for different iterations')
    plt.legend()
    plt.xlim([1.5,5])
    plt.ylim([0.9,1.05])
plt.show()
    
for i in range(len(allJ)):
    plt.plot(zp,allJ[i],label = '{}'.format(i))
    plt.xlabel('z')
    plt.ylabel('J (in terms of n0)')
    plt.title('Current Density Profiles for different iterations')
    plt.legend()
plt.show()
    
for i in range(len(allN)):
    plt.plot(zp,allN[i],label = '{}'.format(i))
    plt.xlabel('z')
    plt.ylabel('Density (in terms of n0)')
    plt.title('Density Profiles for different iterations')
    plt.legend()
plt.show()

# %%
for i in range(len(allb)):
    plt.plot(z,allb[i],label = '{}'.format(i))
    plt.xlabel('z')
    plt.ylabel('b(z)')
    plt.title('Calcultaed b for different iterations')
    plt.legend()
    #plt.xlim([1.5,5])
    #plt.ylim([0.8,1.08])
plt.show()


