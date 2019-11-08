#! /usr/bin/python3
import numpy as np
from sm_func import *
from glob import glob
import matplotlib.pyplot as plt

path = '../data/ArcticArc2007/'
out_path = '../plots/'

#initialize plots
fig1 = plt.figure(figsize=(9,7))

ax = fig1.add_subplot(211)
name = 'ArcticArc2007 vs. SnowModel'
ax.set_title(name,fontsize=25, loc='left')
ax.set_xlabel(r"Latitude",fontsize=20)
ax.set_ylabel(r"Snow depth (cm)",fontsize=20)

bx = fig1.add_subplot(212)
bx.set_xlabel(r"Latitude",fontsize=20)
bx.set_ylabel(r"Spatial derivative (cm)",fontsize=20)


#At each of 20 sites, 100 snow-thickness measurements were collected over a length of 400 m with an average spacing of 4 m (Gerland and Haas, 2011)
#open ski ArcArc2007 file
fname = 'ArcticArc2007_snow.csv'
year = getColumn(path+fname,0, delimiter=',')
mon = getColumn(path+fname,1, delimiter=',')
day = getColumn(path+fname,2, delimiter=',')
pn = getColumn(path+fname,3, delimiter=',')
snod = getColumn(path+fname,4, delimiter=',')
snod_sd = getColumn(path+fname,5, delimiter=',')

snod = np.array(snod,dtype=np.float)#/10         #change from mm to cm
snod_sd = np.array(snod_sd,dtype=np.float)#/10   #change from mm to cm
pn = np.array(pn,dtype=np.int)
print(pn)

#plot observations
ax.errorbar(pn,snod,snod_sd,label='obs (local floe)', capsize=5, capthick=3)

print('Observations')
#print('Mean:')
#print(snod)
#print(np.mean(snod))
#print('Rate of change:')
#print(np.gradient(snod))
print(np.mean(np.abs(np.gradient(snod))))


grad = np.gradient(snod)
grad_m = np.mean(np.abs(grad))

bx.plot(pn,grad,label='obs')


#open SnowModel file for same traverse
#k,iyr,imo,idy,ylat,xlon,snod_obs/10.0,sdev/10.0,(100.0 * snod_m(kk,k),kk=1,nnearest)
fname_part = 'icearc_nearest_25stns_snod_*.dat'
forcing_files = glob(path+fname_part)
for ff in forcing_files:

    pn = getColumn(ff,0, delimiter=' ', skipinitialspace=True, skipheader=False)
    lat = getColumn(ff,4, delimiter=' ', skipinitialspace=True, skipheader=False)
    
    pn = np.array(pn,dtype=np.int)
    lat = np.array(lat,dtype=np.float)
    lat = np.round(lat,1)

    snod_m = np.empty((20,25),dtype=np.float)

    for i in range(8,33):
        st = getColumn(ff,i, delimiter=' ', skipinitialspace=True, skipheader=False)
        st = np.array(st,dtype=np.float)
        j = i-8
        snod_m[:,j]=st

    #calculate mean and standard deviation for the SnowModel
    snod_m_m = np.mean(snod_m,axis=1)
    snod_m_sd = np.std(snod_m,axis=1)
    #print(snod_m_m.shape)

    #just closest station
    snod_m_m1 = np.mean(snod_m[:,:1],axis=1)
    snod_m_sd1 = np.std(snod_m[:,:1],axis=1)

    #just closest 5 stations
    snod_m_m5 = np.mean(snod_m[:,:6],axis=1)
    snod_m_sd5 = np.std(snod_m[:,:6],axis=1)

    #plot model with this forcing
    forcing = ff.split('_')[-1].split('.')[0]
    ax.errorbar(pn,snod_m_m1,snod_m_sd1,label=forcing+' (25 nearest)', capsize=3, capthick=2)
    
    print('Snow Model')
    print(forcing)
    #print('Mean:')
    #print(snod_m_m)
    #print(np.mean(snod_m_m))
    #print('Rate of change:')
    
    grad = np.gradient(snod_m_m1)
    grad_m = np.mean(np.abs(grad))
    print(grad_m)
    
    bx.plot(pn,grad,label=forcing+' (6 nearest)')

    


    
    

ax.set_xticks(pn)
ax.set_xticklabels(lat)

bx.set_xticks(pn)
bx.set_xticklabels(lat)


ax.legend()
bx.legend()
fig1.tight_layout()
fig1.savefig(out_path+'traverse_2007_all')


##scatter plot
#fig2 = plt.figure(figsize=(7,6))
#ax = fig2.add_subplot(111)
#name = 'Variance scatter plot'
#ax.set_title(name,fontsize=29, loc='left')
#ax.set_xlabel(r"ArcticArc2007",fontsize=20)
#ax.set_ylabel(r"SnowModel",fontsize=20)

#sc = ax.scatter(snod_sd**2, snod_m_sd**2, c=pn, cmap=plt.cm.Reds)
#plt.colorbar(sc)

##plt.show()
##exit()

#fig2.tight_layout()
#fig2.savefig(out_path+'traverse_2007_scatter_merra2')

