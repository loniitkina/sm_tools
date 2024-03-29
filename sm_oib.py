#! /usr/bin/python3
import numpy as np
from sm_func import *
from glob import glob
from netCDF4 import Dataset
import pyresample as pr
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

path = '../data/OIB_2017/'
path_cv = '../data/CryoVEX_2017/'
out_path = '../plots/'

#get the OIB data
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

#get only valid data from the flight southwards (NP-Alert)
mask = (snod_oib != -99999.) & (elapsed>56000)
snod_oib = snod_oib[mask] *100
lat_oib = lat_oib[mask]
lon_oib = lon_oib[mask]

##NSIDC grid
#fn = path+'asicd25e2_20121231_v01r02.nc'      #sample file with EASE 25-km grid
#f = Dataset(fn)
#lone = f.variables['longitude'][:]
#late = f.variables['latitude'][:]

#get CryoVEx coordinates
fname = 'CryoVex2017_snow_all.csv'
pn_cv = getColumn(path_cv+fname,3, delimiter=',')
lat_cv = getColumn(path_cv+fname,4, delimiter=',')
lon_cv = getColumn(path_cv+fname,5, delimiter=',')

pn_cv = np.array(pn_cv,dtype=np.int)
lat_cv = np.array(lat_cv,dtype=np.float)
lon_cv = np.array(lon_cv,dtype=np.float)

#project all coordinates
import pyproj
#use OSI-SAF projection: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
wgs84=pyproj.Proj("+init=EPSG:4326") 
nh_stere=pyproj.Proj("+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45")
##transform EASE grid locations (center ponts)
#xe,ye = pyproj.transform(wgs84, nh_stere,lone,late)
#transform CryoVEx points
xcv,ycv = pyproj.transform(wgs84, nh_stere,lon_cv,lat_cv)
#transform OIB points
xoib,yoib = pyproj.transform(wgs84, nh_stere,lon_oib,lat_oib)

##search for the EASE grid points that correrpond to the CryoVEx transect
#import scipy.spatial as spatial
#xye = np.dstack([xe.ravel(),ye.ravel()])[0]
#tree = spatial.KDTree(xye)

hd = 12500  #half-distance of EASE grid point, radius, 12.5km
#hd = 25000

#start the histogram spagetti plot
fig1 = plt.figure(figsize=(8,8))

ax = fig1.add_subplot(111)
name = 'PDFs of OIB'
ax.set_title(name,fontsize=25, loc='left')
ax.set_xlabel(r"Snow depth (cm)",fontsize=20)
ax.set_ylabel(r"Frequency",fontsize=20)

#prepare storing lists
snod_oib_m = []
snod_oib_sd = []
snod_oib_n = []
snod_oib_mo = []
snod_oib_mo2 = []

#loop though all the CryoVEx observation points
for i in range(0,len(lat_cv)):
    print(lat_cv[i],lon_cv[i])

    ##get closest center point of the 25-km EASE grid
    #dist, index = tree.query([xcv[i],ycv[i]])
    #print(dist)

    #center = xye[index]
    ##print(center)
    
    ##double check coordinates
    #lon_center,lat_center = pyproj.transform(nh_stere, wgs84, center[0],center[1])
    #print(lat_center,lon_center)
    
    #replace that EASE center with the CryoVEx observation points
    center = [xcv[i],ycv[i]]
    
    
    #get all OIB coordinates inside 12.5km radius of that center point
    mask = (xoib > center[0]-hd) & (xoib < center[0]+hd) & (yoib > center[1]-hd) & (yoib < center[1]+hd)
    snod_in_grid  = snod_oib[mask]
    #print(snod_in_grid)
    #print(snod_in_grid.shape)
    
    #calculate mean, sd
    snod_m = np.mean(snod_in_grid)
    snod_sd = np.std(snod_in_grid)
    snod_n = snod_in_grid.shape[0]
    
    
    #find mode
    hist = np.histogram(snod_in_grid,bins=25)
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

    #only store secondary mondes that are quite different from primary
    if np.abs(snod_mo-snod_mo2)<15:
        snod_mo2 = snod_mo

    #some secondary means seem to be wrong, here manual fix:
    #if i==4: snod_mo2=35                                        #for 50km box, 25km radius
    if i==7: snod_mo2=12                                        #for 25km box, 12.5km radius
    if i==8: snod_mo2=20
    
    #plot
    num_bins = 25
    n, bins, patches = ax.hist(snod_in_grid, num_bins, histtype='step', label = str(i+1)+'  - measurements: '+str(snod_n))
    #print mean and mode
    ax.plot(snod_mo,hist[0][mm],'o',c='k')
    #plot second index if it is very different from the first one
    if np.abs(snod_mo-snod_mo2)>.1:
        ax.plot(snod_mo2,hist[0][mm2],'o',c='0.5')
    
    
    
    
    
    #store data
    snod_oib_m.append(snod_m)
    snod_oib_sd.append(snod_sd)
    snod_oib_n.append(snod_n)
    snod_oib_mo.append(snod_mo)
    snod_oib_mo2.append(snod_mo2)
    
    print('-----------------------------------------------------')
    #break

ax.legend()
fig1.tight_layout()
fig1.savefig(out_path+'hist_oib_12km')


    
#write out derived data
tt = [pn_cv,lat_cv,lon_cv,snod_oib_m,snod_oib_sd,snod_oib_n,snod_oib_mo,snod_oib_mo2]
table = list(zip(*tt))


outname = path_cv+'OIB_CryoVex2017_snow_12km.csv'
print(outname)
with open(outname, 'wb') as f:
  #header
  f.write(b'ProfileNo,Latitude,Longitude,Zs_mean,Zs_sdev,N,Zs_mode,Zs_mode2\n')
  np.savetxt(f, table, fmt="%s", delimiter=",")
    
    
    
    
    
    
    

