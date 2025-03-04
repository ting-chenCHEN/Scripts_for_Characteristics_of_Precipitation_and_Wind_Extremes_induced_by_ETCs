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
from funcs_map import  lambert_map
from matplotlib.path import Path
import matplotlib.patches as patches
import netCDF4
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
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

mono_font = {'fontname':'monospace'}
sans_font = {'fontname':'sans-serif'}

def get_verts_new(categ_box):
    verts = [(categ_box['west'], categ_box['south']), 
     (categ_box['west'], categ_box['north']), 
     (categ_box['east'], categ_box['north']), 
     (categ_box['east'], categ_box['south']), 
     (0., 0.)]
    return verts

box = {'west':-180, 'east':0, 'north':80, 'south':20}

Topo  = netCDF4.Dataset('Geopotential_orography.nc','r')
elevation = Topo.variables['z'][0,:,:]
print(elevation.shape)

Track=netCDF4.Dataset('TPExtremes_p99p0_ETC_association_2001_2020_ERA5.nc','r')

lat0 = Track.variables['latitude'][:]
lon0 = Track.variables['longitude'][:]

EXcase1 = Track.variables['Extreme_cases_sea'][0,:,:]
EXcase2 = Track.variables['Extreme_cases_sea'][1,:,:]
EXcase3 = Track.variables['Extreme_cases_sea'][2,:,:]
EXcase4 = Track.variables['Extreme_cases_sea'][3,:,:]


EXcaseall=EXcase1+EXcase2+EXcase3+EXcase4

lon, lat = np.meshgrid(lon0, lat0)
print (lat.shape[0])
print (lon.shape[1])

nx= lon.shape[0]
ny= lon.shape[1]

mask=np.ones((nx,ny),dtype=bool)
maskzero=np.ones((nx,ny))
for i in np.arange(nx):
    for j in np.arange (ny):
        if (lat[i,j]<=60 and lat[i,j]>=40 and lon[i,j]>=255 and lon[i,j]<=300):
             mask[i,j] = False
             maskzero[i,j] = 1.


fig = plt.figure(figsize=(5,12))
proj = ccrs.LambertConformal(central_longitude=-87, central_latitude=35,
                                 standard_parallels=[35])

seasons=[0,1,2,3]
for s in seasons: 
    if s == 0:
       EXfreq= EXcase1/EXcaseall
       sea='JJA'
       label='(a)'
    elif s == 1:
       EXfreq= EXcase2/EXcaseall
       sea='SON'
       label='(d)'
    elif s == 2:
       EXfreq= EXcase3/EXcaseall
       sea='DJF'
       label='(g)'
    elif s == 3:
       EXfreq= EXcase4/EXcaseall
       sea='MAM'
       label='(j)'

    ax = fig.add_subplot(4, 1, s+1, projection=proj)

    EXfreq_masked=np.ma.masked_where(mask>=1,EXfreq)
  
    EXfreq_subD = np.average(EXfreq_masked)

    colors=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10]

    fig, ax = lambert_map(extent=(-133, -37, 17, 67), fig=fig,
            cent_lon=-87, ax=ax)
    
    cf = ax.contourf(lon, lat, EXfreq*100., transform=ccrs.PlateCarree(), zorder=6,
                 cmap=plt.cm.get_cmap('Spectral_r'),
                 levels=[0.5,2,5,10,20,30,40,50,60,70],extend='max')

    cf0 = ax.contourf(lon0,lat0,elevation/9.80665, transform=ccrs.PlateCarree(), levels=[1000,5000],
               colors='none',hatches=['////'],zorder=7)
    cf2 = ax.contour(lon,lat,mask, transform=ccrs.PlateCarree(), levels=[1],
               colors='k',linewidths=2.5,zorder=7,alpha=0.7)

    gl= ax.gridlines(draw_labels=True, x_inline=False, y_inline=False,
              color='gray', linestyle='dashed', linewidth=0.5, zorder=21)

    gl.top_labels=False
    gl.bottom_labels=False
    gl.left_labels=False
    gl.right_labels=False

    gl.ylocator = LatitudeLocator()

    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

    gl.xlocator = mticker.FixedLocator([-150,-120,-90,-60,-30])
    gl.ylocator = mticker.FixedLocator([25,45,65,85])


    verts = get_verts_new(box)
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY,]
    path = Path(verts, codes)
    import matplotlib.patheffects as path_effects
    patch = patches.PathPatch(path, facecolor='none', ec='crimson', lw=2,
                            transform=ccrs.PlateCarree(), zorder=30)
    ax.text(0.01,0.03,''+sea+'',horizontalalignment='left',transform=ax.transAxes,fontsize=18,zorder=28,weight='bold',**mono_font)
    ax.text(0.96,0.1,''+str('{:2.0f}'.format(EXfreq_subD*100.))+'%',bbox={'facecolor':'white','alpha':0.7,'boxstyle':'round'},
              horizontalalignment='right', transform=ax.transAxes,fontsize=16,zorder=28)
    
    ax.text(0.01,0.985,''+label+'',bbox={'facecolor':'white','edgecolor':'None','alpha':0.9,'pad':2.0},
              horizontalalignment='left', verticalalignment='top',
              transform=ax.transAxes,fontsize=15,zorder=28,**sans_font)

fig.subplots_adjust(hspace=0.04)
plt.savefig('Fig4_leftcolumn', bbox_inches='tight', dpi=200)
plt.show()
plt.close()
