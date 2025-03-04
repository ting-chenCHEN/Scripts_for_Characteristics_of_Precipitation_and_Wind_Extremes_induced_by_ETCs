#
# Before opening Python execute:
#   module load python2
#   module load development/python-rpn
#

import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np                         # Good module for matrix and matrix operation
import matplotlib.pyplot as plt            # Module to produce figure
import matplotlib.colors as colors
import os                                  # Used to convert png to other format

import funcs_map
from funcs_map import lambert_map
from matplotlib.path import Path
import matplotlib.patches as patches
import netCDF4
from scipy.ndimage.filters import gaussian_filter
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
from matplotlib import gridspec,rc

sans_font = {'fontname':'sans-serif'}
mono_font = {'fontname':'monospace'}

fs=21
fs0=16
#Neu et al.'s(2013) colorbar
c1='#F8EC5A'
c2='#ECB249'
c3='#21A758'
c4='#40B5D9'
c5='#1E90CA'
c6='#355B9F'
c7='#323374'
c8='#9C3D90'
c9='#D94039'
c10='#A4AAB1'

def get_verts_new(categ_box):
    verts = [(categ_box['west'], categ_box['south']), 
     (categ_box['west'], categ_box['north']), 
     (categ_box['east'], categ_box['north']), 
     (categ_box['east'], categ_box['south']), 
     (0., 0.)]
    return verts

box = {'west':-180, 'east':0, 'north':80, 'south':20}

ITrack=netCDF4.Dataset('TPExtremes_p99p0_ETC_association_2001_2020_IMERG.nc','r')
Ilat0 = ITrack.variables['latitude'][:]
Ilon0 = ITrack.variables['longitude'][:]

IEXcase1 = ITrack.variables['Extreme_cases_sea'][0,:,:]
IEXcase2 = ITrack.variables['Extreme_cases_sea'][1,:,:]
IEXcase3 = ITrack.variables['Extreme_cases_sea'][2,:,:]
IEXcase4 = ITrack.variables['Extreme_cases_sea'][3,:,:]

IEXcaseall=IEXcase1+IEXcase2+IEXcase3+IEXcase4

IETCprob1 = ITrack.variables['Prob_Extreme_ass_ETC_sea'][0,:,:]
IETCprob2 = ITrack.variables['Prob_Extreme_ass_ETC_sea'][1,:,:]
IETCprob3 = ITrack.variables['Prob_Extreme_ass_ETC_sea'][2,:,:]
IETCprob4 = ITrack.variables['Prob_Extreme_ass_ETC_sea'][3,:,:]

Ilon, Ilat = np.meshgrid(Ilon0, Ilat0)
print (Ilat.shape[0])
print (Ilon.shape[1])
Inx= Ilon.shape[0]
Iny= Ilon.shape[1]

Imask=np.ones((Inx,Iny),dtype=bool)
for i in np.arange(Inx):
    for j in np.arange (Iny):
        if (Ilat[i,j]<=60 and Ilat[i,j]>=40 and Ilon[i,j]>=255 and Ilon[i,j]<=300):
             Imask[i,j] = False
#---------------------------------
Track=netCDF4.Dataset('TPExtremes_p99p0_ETC_association_2001_2020_ERA5.nc','r')

lat0 = Track.variables['latitude'][:]
lon0 = Track.variables['longitude'][:]

EXcase1 = Track.variables['Extreme_cases_sea'][0,:,:]
EXcase2 = Track.variables['Extreme_cases_sea'][1,:,:]
EXcase3 = Track.variables['Extreme_cases_sea'][2,:,:]
EXcase4 = Track.variables['Extreme_cases_sea'][3,:,:]

EXcaseall=EXcase1+EXcase2+EXcase3+EXcase4

ETCprob1 = Track.variables['Prob_Extreme_ass_ETC_sea'][0,:,:]
ETCprob2 = Track.variables['Prob_Extreme_ass_ETC_sea'][1,:,:]
ETCprob3 = Track.variables['Prob_Extreme_ass_ETC_sea'][2,:,:]
ETCprob4 = Track.variables['Prob_Extreme_ass_ETC_sea'][3,:,:]

lon, lat = np.meshgrid(lon0, lat0)
print (lat.shape[0])
print (lon.shape[1])
nx= lon.shape[0]
ny= lon.shape[1]

mask=np.ones((nx,ny),dtype=bool)
mask60d=np.ones((nx,ny),dtype=bool)
for i in np.arange(nx):
    for j in np.arange (ny):
        if (lat[i,j]<=60) :
            mask60d[i,j]= False
        if (lat[i,j]<=60 and lat[i,j]>=40 and lon[i,j]>=255 and lon[i,j]<=300):
             mask[i,j] = False

EXfreq1= EXcase1/EXcaseall
ETCcont1=ETCprob1
sea1='JJA'
IEXfreq1= IEXcase1/IEXcaseall
IETCcont1=IETCprob1

EXfreq2= EXcase2/EXcaseall
ETCcont2=ETCprob2
sea2='SON'
IEXfreq2= IEXcase2/IEXcaseall
IETCcont2=IETCprob2

EXfreq3= EXcase3/EXcaseall
ETCcont3=ETCprob3
sea3='DJF'
IEXfreq3= IEXcase3/IEXcaseall
IETCcont3=IETCprob3

EXfreq4= EXcase4/EXcaseall
ETCcont4=ETCprob4
sea4='MAM'
IEXfreq4= IEXcase4/IEXcaseall
IETCcont4=IETCprob4

#-----------------
# For Subdomain over NNA
ETCcont1=np.ma.masked_where(EXfreq1*100<0.5,EXfreq1)
ETCcont2=np.ma.masked_where(EXfreq2*100<0.5,EXfreq2)
ETCcont3=np.ma.masked_where(EXfreq3*100<0.5,EXfreq3)
ETCcont4=np.ma.masked_where(EXfreq4*100<0.5,EXfreq4)
ETCcont_masked1=np.ma.masked_where(mask>=1,ETCcont1)
ETCcont_masked2=np.ma.masked_where(mask>=1,ETCcont2)
ETCcont_masked3=np.ma.masked_where(mask>=1,ETCcont3)
ETCcont_masked4=np.ma.masked_where(mask>=1,ETCcont4)
  
IETCcont1=np.ma.masked_where(IEXfreq1*100<0.5,IEXfreq1)
IETCcont2=np.ma.masked_where(IEXfreq2*100<0.5,IEXfreq2)
IETCcont3=np.ma.masked_where(IEXfreq3*100<0.5,IEXfreq3)
IETCcont4=np.ma.masked_where(IEXfreq4*100<0.5,IEXfreq4)
IETCcont_masked1=np.ma.masked_where(Imask>=1,IETCcont1)
IETCcont_masked2=np.ma.masked_where(Imask>=1,IETCcont2)
IETCcont_masked3=np.ma.masked_where(Imask>=1,IETCcont3)
IETCcont_masked4=np.ma.masked_where(Imask>=1,IETCcont4)
  
ETCcont_subD1 = np.average(ETCcont_masked1)*100.
ETCcont_subD2 = np.average(ETCcont_masked2)*100.
ETCcont_subD3 = np.average(ETCcont_masked3)*100.
ETCcont_subD4 = np.average(ETCcont_masked4)*100.

IETCcont_subD1 = np.average(IETCcont_masked1)*100.
IETCcont_subD2 = np.average(IETCcont_masked2)*100.
IETCcont_subD3 = np.average(IETCcont_masked3)*100.
IETCcont_subD4 = np.average(IETCcont_masked4)*100.

#----------------------

fig = plt.figure(figsize=(6,3))
gs = gridspec.GridSpec(1,1)

ax0=plt.subplot(gs[0])
c='crimson'
plt.scatter(0, ETCcont_subD1, s=350, c='white', marker="D",alpha=0.6)
plt.scatter(0, ETCcont_subD1, s=350, c=c, marker="D")
plt.scatter(0.2, IETCcont_subD1, s=350, c='white', marker="o",alpha=0.6)
plt.scatter(0.2, IETCcont_subD1, s=350, c=c, marker="o",alpha=0.5)


c='orange'
plt.scatter(0.6, ETCcont_subD2, s=350, c='white', marker="D",alpha=0.6)
plt.scatter(0.6, ETCcont_subD2, s=350, c=c, marker="D")
plt.scatter(0.8, IETCcont_subD2, s=350, c='white', marker="o",alpha=0.6)
plt.scatter(0.8, IETCcont_subD2, s=350, c=c, marker="o",alpha=0.5)

c='royalblue'
plt.scatter(1.2, ETCcont_subD3, s=350, c='white', marker="D",alpha=0.6)
plt.scatter(1.2, ETCcont_subD3, s=350, c=c, marker="D")
plt.scatter(1.4, IETCcont_subD3, s=350, c='white', marker="o",alpha=0.6)
plt.scatter(1.4, IETCcont_subD3, s=350, c=c, marker="o",alpha=0.5)

c='green'
plt.scatter(1.8, ETCcont_subD4, s=350, c='white', marker="D",alpha=0.6)
plt.scatter(1.8, ETCcont_subD4, s=350, c=c, marker="D")
plt.scatter(2.0, IETCcont_subD4, s=350, c='white', marker="o",alpha=0.6)
plt.scatter(2.0, IETCcont_subD4, s=350, c=c, marker="o",alpha=0.5)


for marker in ["D","o"]:
    if marker=="D":
        label='ERA5'
    else:
        label='IMERG'
    plt.scatter([], [], c='k', alpha=0.5, s=150, marker=marker,
                label=label)
plt.legend(ncol=2, fontsize="16",columnspacing =0.1,handletextpad=0)


ax0.grid(axis='y',linestyle=':',zorder=0)

plt.xlim(-0.3,2.2)
plt.ylim(0,58)
plt.yticks(fontsize=fs0)
ax0.set_xticks([0.1,0.7,1.3,1.9])
ax0.set_yticks([10,20,30,40,50])

ax0.set_xticklabels(['JJA','SON','DJF','MAM'],rotation=0,ha='center',fontsize=fs,**mono_font,weight='bold')
colors = ['crimson', 'orange', 'royalblue','green']
for xtick, color in zip(ax0.get_xticklabels(), colors):
    xtick.set_color(color)

ax0.text(0.01,0.985,'(m)',bbox={'facecolor':'white','edgecolor':'None','alpha':0.9,'pad':2.0},
              horizontalalignment='left', verticalalignment='top',
              transform=ax0.transAxes,fontsize=18,zorder=28,**sans_font)
plt.savefig('Fig4_m', bbox_inches='tight', dpi=200)
plt.show()
plt.close()
