#! /usr/bin/python3
import numpy as np
from glob import glob
from sm_func import *
import matplotlib.pyplot as plt

out_path = '../plots/'

#CryoVEX_2017
path = '../data/CryoVEX_2017/CryoVEX2017-OpenData-master/'

#open data file
fname = 'Alert88N_Snow.csv'
#fname = 'Alert88N_Snow_core.csv'

time = getColumn(path+fname,0, delimiter=',')
lat = getColumn(path+fname,1, delimiter=',')
lon = getColumn(path+fname,2, delimiter=',')
snod = getColumn(path+fname,3, delimiter=',')
pn = getColumn(path+fname,5, delimiter=',')

snod = np.array(snod,dtype=np.float)
lat = np.array(lat,dtype=np.float)
#pn = np.array(pn,dtype=np.int)

#just pack ice N of the shear zone (no ice bridge, landfast ice etc)
snod_pack = snod[(lat>84)]
print(snod.shape)

#plot
fig1 = plt.figure(figsize=(14,5))
ax = fig1.add_subplot(131)
ax.set_title('CryoVex 2017 N of the shear zone')


num_bins = 30
n, bins, patches = ax.hist(snod_pack, num_bins, lw=3, color='r', alpha=1, density=True, label = 'all pack ice', histtype='step')


# make a fit to the samples
from scipy import stats 
samples = snod
shape, loc, scale = stats.lognorm.fit(samples, floc=0)
x_fit       = np.linspace(samples.min(), samples.max(), 100)
samples_fit = stats.lognorm.pdf(x_fit, shape, loc=loc, scale=scale)

loc, scale = stats.norm.fit(samples)
samples_fit_n = stats.norm.pdf(x_fit, loc=loc, scale=scale)

# plot fit into histogram
ax.plot(x_fit, samples_fit, label='log-normal PDF', linewidth=2)
ax.plot(x_fit, samples_fit_n, label='normal PDF', linewidth=2)


#just pack ice N of the shear zone (no ice bridge, landfast ice etc)
#MYI
snod_myi = snod[(lat>84)&(lat<87.)]
#FYI
snod_fyi = snod[lat>87.]

n, bins, patches = ax.hist(snod_fyi, num_bins, facecolor='blue', alpha=0.3, density=True, label = 'FY pack ice')
n, bins, patches = ax.hist(snod_myi, num_bins, facecolor='green', alpha=0.3, density=True, label = 'MY pack ice')

ax.legend()


#N-ICE-2015 April/May data
path = '../data/N-ICE2015_MP_v1/'

#open data file
fname1 = '2015_01_*_MP_i.txt'
fname2 = '2015_02_*_MP_i.txt'
fname3 = '2015_03_*_MP_i.txt'
fname4 = '2015_04_*_MP_i.txt'
fname5 = '2015_05_*_MP_i.txt'
fname6 = '2015_06_*_MP_i.txt'

#all measurements
fl = sorted(glob(path+fname1)+glob(path+fname2)+glob(path+fname3)+glob(path+fname4)+glob(path+fname5)+glob(path+fname6))

#spring only
fl = sorted(glob(path+fname4)+glob(path+fname5)+glob(path+fname6))

print(fl)

bx = fig1.add_subplot(132)
bx.set_title('N-ICE 2015 spring transects')

sl = []
for i in fl:
    print(i)

    time = getColumn(i,0, delimiter=' ')
    idx = getColumn(i,1, delimiter=' ')
    s = getColumn(i,2, delimiter=' ')
    lat = getColumn(i,3, delimiter=' ')
    lon = getColumn(i,4, delimiter=' ')
    
    s = np.array(s,dtype=np.float)
    s = s[~np.isnan(s)]            #remove nan values
    
    #check that no values are higher than 120 (max depth possible with MagnaProbe)
    s = np.where(s>120,120,s)
    
    #exclude short transects (limited to one floe or special areas e.g. refrozen leads)
    print(s.shape)
    #if s.shape[0] < 2500: continue

    ##exclude the one with a long lead
    #if i == '../data/N-ICE2015_MP_v1/2015_05_08_MP_i.txt': print('gone!');continue
    #if i == '../data/N-ICE2015_MP_v1/2015_01_31_MP_i.txt': print('gone!');continue
    ##ALS overflight (very limited field with level ice)
    #if i == '../data/N-ICE2015_MP_v1/2015_04_18_MP_i.txt': print('gone!');continue

    
    num_bins = 30
    n, bins, patches = bx.hist(s, num_bins, facecolor='blue', alpha=0.3, density=True)
    #plt.show()

    sl.extend(s)

snod = np.array(sl)
print(snod.shape)
   
num_bins = 30
n, bins, patches = bx.hist(snod, num_bins, lw=3, color='r', alpha=1, histtype='step', density=True)


# make a fit to the samples
from scipy import stats 
samples = snod
shape, loc, scale = stats.lognorm.fit(samples, floc=0)
x_fit       = np.linspace(samples.min(), samples.max(), 100)
samples_fit = stats.lognorm.pdf(x_fit, shape, loc=loc, scale=scale)

loc, scale = stats.norm.fit(samples)
samples_fit_n = stats.norm.pdf(x_fit, loc=loc, scale=scale)

# plot fit into histogram
bx.plot(x_fit, samples_fit, label='log-normal PDF', linewidth=2)
bx.plot(x_fit, samples_fit_n, label='normal PDF', linewidth=2)
bx.set_xlabel('Snow Depth (cm)')
bx.legend()

#Barneo 2016
path = '../data/Barneo/data/data/MagnaProbe/'

#open data file
fname = 'NPI*.dat'

fl = sorted(glob(path+fname))
print(fl)

cx = fig1.add_subplot(133)
cx.set_title('Barneo 2016 transects')

sl = []
for i in fl:
    print(i)

    s = getColumn(i,3, delimiter=',', magnaprobe=True)
    
    s = np.array(s,dtype=np.float)
    s = s[~np.isnan(s)]            #remove nan values
    
    ##exclude short transects (limited to one floe or special areas e.g. refrozen leads)
    #print(s.shape)
    #if s.shape[0] < 2000: continue
    
    num_bins = 30
    n, bins, patches = cx.hist(s, num_bins, facecolor='blue', alpha=0.3, density=True)
    #plt.show()

    sl.extend(s)

snod = np.array(sl)  
   
num_bins = 30
n, bins, patches = cx.hist(snod, num_bins, lw=3, color='r', alpha=1, histtype='step', density=True)
  
# make a fit to the samples
from scipy import stats 
samples = snod
shape, loc, scale = stats.lognorm.fit(samples, floc=0)
x_fit       = np.linspace(samples.min(), samples.max(), 100)
samples_fit = stats.lognorm.pdf(x_fit, shape, loc=loc, scale=scale)

loc, scale = stats.norm.fit(samples)
samples_fit_n = stats.norm.pdf(x_fit, loc=loc, scale=scale)

# plot fit into histogram
cx.plot(x_fit, samples_fit, label='log-normal PDF', linewidth=2)
cx.plot(x_fit, samples_fit_n, label='normal PDF', linewidth=2)
cx.set_xlabel('Snow Depth (cm)')
cx.legend()


fig1.tight_layout()
fig1.savefig(out_path+'hist')
    
