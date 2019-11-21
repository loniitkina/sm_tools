#! /usr/bin/python3
import numpy as np
from glob import glob
from scipy.stats import mode
from sm_func import *
import matplotlib.pyplot as plt

out_path = '../plots/'

#CryoVEX_2017
path = '../data/CryoVEX_2017/CryoVEX2017-OpenData-master/'

#open data file
#fname = 'Alert88N_Snow.csv'            #this files also includes sites like NE3 and nost just the core transect
fname = 'Alert88N_Snow_core.csv'

time = getColumn(path+fname,0, delimiter=',')
lat = getColumn(path+fname,1, delimiter=',')
lon = getColumn(path+fname,2, delimiter=',')
snod = getColumn(path+fname,3, delimiter=',')
pn = getColumn(path+fname,5, delimiter=',')

snod = np.array(snod,dtype=np.float)
lat = np.array(lat,dtype=np.float)
lon = np.array(lon,dtype=np.float)
pn = np.array(pn,dtype=np.int)
time = np.array(time,dtype=np.str)

#get the OIB data
path = '../data/OIB_2017/'
fname = 'IDCSI4_20170411.txt'
lat_oib = getColumn(path+fname,0, delimiter=',')
lon_oib = getColumn(path+fname,1, delimiter=',')
snod_oib = getColumn(path+fname,7, delimiter=',')
snod_unc = getColumn(path+fname,8, delimiter=',')
elapsed = getColumn(path+fname,16, delimiter=',')

lat_oib = np.array(lat_oib,dtype=np.float)
lon_oib = np.array(lon_oib,dtype=np.float)
snod_oib = np.array(snod_oib,dtype=np.float)
snod_unc = np.array(snod_unc,dtype=np.float)
elapsed = np.array(elapsed,dtype=np.float)

#get only valid data from the flight northwards
mask = (snod_oib != -99999.) & (elapsed>56000)
snod_oib = snod_oib[mask] *100                  #change from m to cm
lat_oib = lat_oib[mask]
lon_oib = lon_oib[mask]


#SnowModel MERRA-2
path = '../data/CryoVEX_2017/'
ff = path+'cryovex_nearest_25stns_snod_merra2.dat'

snod_merra = np.empty((10,25),dtype=np.float)
for i in range(8,33):
    st = getColumn(ff,i, delimiter=' ', skipinitialspace=True, skipheader=False)
    st = np.array(st,dtype=np.float)
    j = i-8
    snod_merra[:,j]=st

#SnowModel ERA-5
ff = path+'cryovex_nearest_25stns_snod_era5.dat'
snod_era = np.empty((10,25),dtype=np.float)
for i in range(8,33):
    st = getColumn(ff,i, delimiter=' ', skipinitialspace=True, skipheader=False)
    st = np.array(st,dtype=np.float)
    j = i-8
    snod_era[:,j]=st


#start the histogram spagetti plot
fig1 = plt.figure(figsize=(8,8))

ax = fig1.add_subplot(111)
name = 'PDFs of all CryoVEX2017'
ax.set_title(name,fontsize=25, loc='left')
ax.set_xlabel(r"Snow depth (cm)",fontsize=20)
ax.set_ylabel(r"Probability",fontsize=20)

years = []
months = []
days = []
profile = []
lats = []
lons = []
means = []
sdevs = []
modes = []
modes2 = []
ns = []

#group in 3 regions + shear zone (no data from CryoVEx)
#calculate statistics for each site
for i in range(4,5):                   #10 sites
    print(i)
    if i == 0:                      #stations 1-4 on realtivelly immobile ice
        region = 'ice bridge'
        print(region)
        snodl = snod[(pn<5)&(snod<120)]    #select site and remove all values of 120cm - that is the instrument top limit it accumulates high values, better to simply remote that and say we have excluded ridges
        
        snodl_oib = snod_oib[(lat_oib>83.5)&(lat_oib<84)]    #south of the shear zone
        
        model_merra = np.hstack(snod_merra[:4,:8])
        #ax.scatter(model,np.zeros_like(model),c='m')
        
        model_era = np.hstack(snod_era[:4,:8])
        #ax.scatter(model,np.zeros_like(model),c='c')
        
        
    if i == 1:                      #stations 1-4 on realtivelly immobile ice
        region = 'shear zone'
        print(region)
        snodl = snod[(pn<5)&(pn>5) &(snod<120)]    #select site and remove all values of 120cm - that is the instrument top limit it accumulates high values, better to simply remote that and say we have excluded ridges
        
        snodl_oib = snod_oib[(lat_oib>84)&(lat_oib<85)]    #south of the shear zone
               
    if i == 2:                      #stations 1-4 on realtivelly immobile ice
        region = 'SYI zone'
        print(region)
        snodl = snod[(pn>5)&(pn<10) &(snod<120)]    #select site and remove all values of 120cm - that is the instrument top limit it accumulates high values, better to simply remote that and say we have excluded ridges
        
        snodl_oib = snod_oib[(lat_oib>85)&(lat_oib<86.9)]    #south of the shear zone
        
        model_merra = np.hstack(snod_merra[4:9,:8])
        #ax.scatter(model,np.zeros_like(model),c='m')
        
        model_era = np.hstack(snod_era[4:9,:8])
        #ax.scatter(model,np.zeros_like(model),c='c')
       
        
    if i == 3:                      #stations 1-4 on realtivelly immobile ice
        region = 'FYI zone'
        print(region)
        snodl = snod[(pn==10)&(snod<120)]    #select site and remove all values of 120cm - that is the instrument top limit it accumulates high values, better to simply remote that and say we have excluded ridges
        
        snodl_oib = snod_oib[(lat_oib>83.5)&(lat_oib<87.1)]    #south of the shear zone
        
        model_merra = np.hstack(snod_merra[9,:8])
        #ax.scatter(model,np.zeros_like(model),c='m')
        
        model_era = np.hstack(snod_era[9,:8])
        #ax.scatter(model,np.zeros_like(model),c='c')   
        
    if i ==4:
        region = 'whole transect'
        print(region)
        snodl = snod[(snod<120)]    #select site and remove all values of 120cm - that is the instrument top limit it accumulates high values, better to simply remote that and say we have excluded ridges
        
        snodl_oib = snod_oib[(lat_oib>86.9)&(lat_oib<87.1)]    #south of the shear zone
        
        model_merra = np.hstack(snod_merra[:,:8])
        model_merra = model_merra[model_merra!=0]
        #ax.scatter(model,np.zeros_like(model),c='m')
        
        model_era = np.hstack(snod_era[:,:8])
        model_era = model_era[model_era!=0]
        #ax.scatter(model,np.zeros_like(model),c='c')         
        
        
    number = snodl.shape[0]
    snod_m = np.mean(snodl)
    snod_sd = np.std(snodl)
    print(snod_m)
    
    #find mode
    hist = np.histogram(snodl,bins=25)
    #print(hist)
    srt = np.argsort(hist[0])                           #indexes that would sort the array
    mm = srt[-1]                                        #same as: np.argmax(hist[0])
    #print(srt)
    #print(hist[1][mm])
    mm1 = np.argmax(hist[0])
    #print(mm,mm1)
    #print(hist[1][mm1])
    snod_mo = (hist[1][mm] + hist[1][mm+1])/2           #take mean of the bin for the mode value
    print(snod_mo)
    
    
    
    #find secondary mode
    #make a gradient of sorting indexes, first sharp gradient (from the back) is the second mode. 
    #print(srt)
    #grad = np.gradient(srt)
    #print(grad)
    #print(srt[1:]-srt[:-1])
    grad = abs(srt[1:]-srt[:-1])                          #treat negative and positive order changes equally
    #print(grad)
    srt2 = np.argsort(grad[:4])                           #Only search the first 1/4 of values from the back
    #print(srt2)
    mmt = srt2[-1]+2                                      #add 2: one for lost boundary in gradient and one for indexing from the back
    #print(mmt)
    mm2 = srt[-mmt]
    #print(mm2)
    #mm2 = np.argmax(grad[:6])
    snod_mo2 = (hist[1][mm2] + hist[1][mm2+1])/2
    print(snod_mo2)
    
    ##exit()
    ##some secondary means seem to be wrong, here manual fix:
    #if i==1: snod_mo2=50
    #if i==4: snod_mo2=41
    #if i==6: snod_mo2=8  
    #if i==9: snod_mo2=28 
    
    
    #plot
    num_bins = np.arange(0,125,5)
    ncv, bins, patches = ax.hist(snodl, num_bins, histtype='step', density=True, label = region+' CV - n: '+str(number),lw=3)
    ##print mean and mode
    #ax.plot(snod_mo,hist[0][mm],'o',c='k')
    ##plot second index if it is very different from the first one
    #if np.abs(snod_mo-snod_mo2)>10:
        #ax.plot(snod_mo2,hist[0][mm2],'o',c='0.5')


    noib, bins, patches = ax.hist(snodl_oib, num_bins, histtype='step', density=True, label = region+' OIB - n: '+str(snodl_oib.shape[0]),lw=3)

    #save histogram data for Glen
    tt = [bins,ncv,noib]
    table = list(zip(*tt))
    outname = '../data/CryoVEX_2017/CryoVex2017_hist.csv'
    print(outname)
    with open(outname, 'wb') as f:
        #header
        f.write(b'bin min,probability CV, probability OIB\n')
        np.savetxt(f, table, fmt="%s", delimiter=",")
    
    
    num_bins = np.arange(0,127,7)
    nm, bins, patches = ax.hist(model_merra, num_bins, histtype='step', density=True, label = region+' MERRA-2 - n: '+str(model_merra.shape[0]),lw=3)
    ne, bins, patches = ax.hist(model_era, num_bins, histtype='step', density=True, label = region+' ERA-5 - n: '+str(model_era.shape[0]),lw=3)


    #save histogram data for Glen
    tt = [bins,nm,ne]
    table = list(zip(*tt))
    outname = '../data/CryoVEX_2017/CryoVex2017_hist_snowmodel.csv'
    print(outname)
    with open(outname, 'wb') as f:
        #header
        f.write(b'bin min, probability merra, probability era\n')
        np.savetxt(f, table, fmt="%s", delimiter=",")
    

    #estimate date for each point
    timel = time[(pn<5)&(snod<120)][0]
    #print(timel)
    year = timel.split('-')[0]
    mon = timel.split('-')[1]
    day = timel.split('-')[2].split(' ')[0]
    print(year,mon,day)
    
    #get the point lat,lon
    latl = lat[(pn<5)&(snod<120)][0]
    lonl = lon[(pn<5)&(snod<120)][0]
    print(latl,lonl)
    
    print('-------')
    
    #store all the data in the lists
    years.append(year)
    months.append(mon)
    days.append(day)
    profile.append(i)
    lats.append(latl)
    lons.append(lonl)
    means.append(snod_m)
    sdevs.append(snod_sd)
    modes.append(snod_mo)
    if np.abs(snod_mo-snod_mo2)>10:
        modes2.append(snod_mo2)
    else:
        modes2.append(snod_mo)  #if no second mode, then same as primary mode
    ns.append(number)
    
    break

ax.legend()
fig1.tight_layout()
fig1.savefig(out_path+'hist_group5')

#write out derived data
tt = [years,months,days,profile,lats,lons,means,sdevs,modes,modes2,ns]
table = list(zip(*tt))


outname = '../data/CryoVEX_2017/CryoVex2017_snow_regions.csv'
print(outname)
with open(outname, 'wb') as f:
  #header
  f.write(b'Year,Month,Day,ProfileNo,Latitude,Longitude,Zs_mean,Zs_sdev,Zs_mode,Zs_mode2,N\n')
  np.savetxt(f, table, fmt="%s", delimiter=",")
