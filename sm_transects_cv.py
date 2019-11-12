#! /usr/bin/python3
import numpy as np
from sm_func import *
from glob import glob
import matplotlib.pyplot as plt

path = '../data/CryoVEX_2017/'
out_path = '../plots/'

#initialize plots
fig1 = plt.figure(figsize=(12,9))

ax = fig1.add_subplot(111)
name = 'CryoVEX2017 vs. SnowModel'
#ax.set_title(name,fontsize=25, loc='left')
ax.set_xlabel(r"Latitude",fontsize=20)
ax.set_ylabel(r"Snow depth (cm)",fontsize=20)

#bx = fig1.add_subplot(212)
#bx.set_xlabel(r"Latitude",fontsize=20)
#bx.set_ylabel(r"Spatial derivative (cm)",fontsize=20)

#open obs file
fname = 'CryoVex2017_snow.csv'
fname = 'CryoVex2017_snow_all.csv'
year = getColumn(path+fname,0, delimiter=',')
mon = getColumn(path+fname,1, delimiter=',')
day = getColumn(path+fname,2, delimiter=',')
pn = getColumn(path+fname,3, delimiter=',')
snod = getColumn(path+fname,6, delimiter=',')
snod_sd = getColumn(path+fname,7, delimiter=',')
snod_mo = getColumn(path+fname,8, delimiter=',')
snod_mo2 = getColumn(path+fname,9, delimiter=',')
lat = getColumn(path+fname,4, delimiter=',')

snod = np.array(snod,dtype=np.float)#*100         #change from m to cm
snod_sd = np.array(snod_sd,dtype=np.float)#*100   #change from m to cm
snod_mo = np.array(snod_mo,dtype=np.float)
snod_mo2 = np.array(snod_mo2,dtype=np.float)
pn = np.array(pn,dtype=np.int)
print(pn)
print(snod)
lat = np.array(lat,dtype=np.float)
lat = np.round(lat,1)

#sort modes into high and low modes
snod_mo_h = np.where(snod_mo2>22,snod_mo2,snod_mo)
snod_mo_l = np.where(snod_mo2<22,snod_mo2,snod_mo)

#OIB data
fname = 'OIB_CryoVex2017_snow.csv'
snod_oib_m = getColumn(path+fname,3, delimiter=',')
snod_oib_sd = getColumn(path+fname,4, delimiter=',')
snod_oib_mo = getColumn(path+fname,6, delimiter=',')
snod_oib_mo2 = getColumn(path+fname,7, delimiter=',')

snod_oib_m = np.array(snod_oib_m,dtype=np.float)*100         #change from m to cm
snod_oib_sd = np.array(snod_oib_sd,dtype=np.float)*100         #change from m to cm
snod_oib_mo = np.array(snod_oib_mo,dtype=np.float)*100
snod_oib_mo2 = np.array(snod_oib_mo2,dtype=np.float)*100

snod_oib_mo_h = np.where(snod_oib_mo2>22,snod_oib_mo2,snod_oib_mo)
snod_oib_mo_l = np.where(snod_oib_mo2<22,snod_oib_mo2,snod_oib_mo)


#plot observations
ax.errorbar(pn,snod,snod_sd,label='CryoVEx 2017 mean and standard deviation', capsize=5, capthick=3)
ax.errorbar(pn,snod_oib_m,snod_oib_sd,label='OIB 2017 mean and standard deviation', capsize=5, capthick=3,c='k')
ax.plot(pn,snod_mo_h,label='CryoVEx 2017 high mode')
ax.plot(pn,snod_mo_l,label='CryoVEx 2017 low mode')
ax.plot(pn,snod_oib_mo_h,label='OIB 2017 high mode',c='.3',ls='--')
ax.plot(pn,snod_oib_mo_l,label='OIB 2017 low mode',c='.3',ls=':')

#grad = np.gradient(snod)
#grad_m = np.mean(np.abs(grad))

#grad_mo = np.gradient(snod_mo)
#grad_mo_m = np.mean(np.abs(grad_mo))

#bx.plot(pn,grad,label='obs mean')
#bx.plot(pn,grad_mo,label='obs mode')

#print('Observations')
#print(grad_m)
#print(grad_mo_m)


colors = ['r','m']
ci=0

#open SnowModel file for same traverse
#k,iyr,imo,idy,ylat,xlon,snod_obs/10.0,sdev/10.0,(100.0 * snod_m(kk,k),kk=1,nnearest)
fname_part = 'cryovex_nearest_25stns_snod_*.dat'
forcing_files = glob(path+fname_part)
for ff in forcing_files:

    print(ff)
    pn = getColumn(ff,0, delimiter=' ', skipinitialspace=True, skipheader=False)
    lat = getColumn(ff,4, delimiter=' ', skipinitialspace=True, skipheader=False)
    
    pn = np.array(pn,dtype=np.int)
    lat = np.array(lat,dtype=np.float)
    lat = np.round(lat,1)

    snod_m = np.empty((10,25),dtype=np.float)

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

    #just closest 6 stations
    snod_m_m3 = np.mean(snod_m[:,:3],axis=1)
    snod_m_sd3 = np.std(snod_m[:,:3],axis=1)
    
    #closest 3 stations if snow depth is low
    snod_m_smooth = np.where(snod_m_m1<15,snod_m_m3,snod_m_m1)

    #plot model with this forcing
    forcing = ff.split('_')[-1].split('.')[0]
    ax.plot(pn,snod_m_m1,label='SnowModel '+forcing+' (1 nearest)',c=colors[ci])
    ax.plot(pn,snod_m_smooth,label='SnowModel '+forcing+' (adjusted low - 3 nearest)',ls='--',c=colors[ci])
    
    #print('Snow Model')
    #print(forcing)
    ##print('Mean:')
    ##print(snod_m_m)
    ##print(np.mean(snod_m_m))
    ##print('Rate of change:')
    
    #grad = np.gradient(snod_m_m1)
    #grad_m = np.mean(np.abs(grad))
    #print(grad_m)
    
    #bx.plot(pn,grad,label=forcing+' (1 nearest)')
    
    ci = ci+1

ax.set_xticks(pn)
ax.set_xticklabels(lat)

#bx.set_xticks(pn)
#bx.set_xticklabels(lat)


ax.legend()
#bx.legend()
fig1.tight_layout()
fig1.savefig(out_path+'traverse_2017_all')


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

