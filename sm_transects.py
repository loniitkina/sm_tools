#! /usr/bin/python3
import numpy as np
from sm_func import *
import matplotlib.pyplot as plt

path = 'data/'
out_path = 'plots/'

#open ski ArcArc2007 file
fname = 'ArcticArc2007_snow.csv'
year = getColumn(path+fname,0, delimiter=',')
mon = getColumn(path+fname,1, delimiter=',')
day = getColumn(path+fname,2, delimiter=',')
pn = getColumn(path+fname,3, delimiter=',')
snod = getColumn(path+fname,4, delimiter=',')
snod_sd = getColumn(path+fname,5, delimiter=',')

snod = np.array(snod,dtype=np.float)/10         #change from mm to cm
snod_sd = np.array(snod_sd,dtype=np.float)/10   #change from mm to cm
pn = np.array(pn,dtype=np.int)
print(pn)


#open SnowModel file for same traverse
fname = 'nearest_25stns_snod.dat'
pn = getColumn(path+fname,0, delimiter=' ', skipinitialspace=True, skipheader=False)

snod_m = np.empty((20,25),dtype=np.float)

for i in range(1,26):
    st = getColumn(path+fname,i, delimiter=' ', skipinitialspace=True, skipheader=False)
    st = np.array(st,dtype=np.float)
    j = i-1
    snod_m[:,j]=st
    #print(st)
    

pn = np.array(pn,dtype=np.int)
print(pn)

#calculate mean and standard deviation for the SnowModel
snod_m_m = np.mean(snod_m,axis=1)
snod_m_sd = np.std(snod_m,axis=1)
print(snod_m_m.shape)

#just closest 15 stations
snod_m_m15 = np.mean(snod_m[:,:15],axis=1)
snod_m_sd15 = np.std(snod_m[:,:15],axis=1)

#just closest 5 stations
snod_m_m5 = np.mean(snod_m[:,:5],axis=1)
snod_m_sd5 = np.std(snod_m[:,:5],axis=1)

#plot
fig1 = plt.figure(figsize=(8,6))
ax = fig1.add_subplot(211)
name = 'ArcticArc2007 vs. SnowModel'
ax.set_title(name,fontsize=25, loc='left')
ax.set_xlabel(r"Station number",fontsize=20)
ax.set_ylabel(r"Snow depth (cm)",fontsize=20)

ax.errorbar(pn,snod,snod_sd,label='obs', capsize=5, capthick=3)
ax.errorbar(pn,snod_m_m,snod_m_sd5,label='obs + model SD 5 nearest', capsize=3, capthick=2)
ax.legend()

bx = fig1.add_subplot(212)
bx.set_xlabel(r"Station number",fontsize=20)
bx.set_ylabel(r"Snow depth (cm)",fontsize=20)

bx.errorbar(pn,snod_m_m,snod_m_sd,label='model 25 nearest', capsize=3, capthick=2)
bx.errorbar(pn,snod_m_m15,snod_m_sd15,label='model 15 nearest', capsize=3, capthick=2)
bx.errorbar(pn,snod_m_m5,snod_m_sd5,label='model 5 nearest', capsize=3, capthick=2)
bx.legend()

fig1.tight_layout()
fig1.savefig(out_path+'traverse_2007')

#scatter plot
fig2 = plt.figure(figsize=(7,6))
ax = fig2.add_subplot(111)
name = 'Variance scatter plot'
ax.set_title(name,fontsize=29, loc='left')
ax.set_xlabel(r"ArcticArc2007",fontsize=20)
ax.set_ylabel(r"SnowModel",fontsize=20)

sc = ax.scatter(snod_sd**2, snod_m_sd**2, c=pn, cmap=plt.cm.Reds)
plt.colorbar(sc)

#plt.show()
#exit()

fig2.tight_layout()
fig2.savefig(out_path+'traverse_2007_scatter')

