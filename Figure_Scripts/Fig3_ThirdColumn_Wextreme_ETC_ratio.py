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
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)

sans_font = {'fontname':'sans-serif'}

c1='#edf8fb'
c2='#bfd3e6'
c3='#9ebcda'
c4='#8c96c6'
c5='#8856a7'
c6='#810f7c'

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

Track=netCDF4.Dataset('WSExtremes_p99p0_ETC_association_2001_2020_ERA5.nc','r')

lat0 = Track.variables['latitude'][:]
lon0 = Track.variables['longitude'][:]

EXprecp1 = Track.variables['Extremevalues_total_sea'][0,:,:]
EXprecp2 = Track.variables['Extremevalues_total_sea'][1,:,:]
EXprecp3 = Track.variables['Extremevalues_total_sea'][2,:,:]
EXprecp4 = Track.variables['Extremevalues_total_sea'][3,:,:]

EXprecpall=EXprecp1+EXprecp2+EXprecp3+EXprecp4

EXcase1 = Track.variables['Extreme_cases_sea'][0,:,:]
EXcase2 = Track.variables['Extreme_cases_sea'][1,:,:]
EXcase3 = Track.variables['Extreme_cases_sea'][2,:,:]
EXcase4 = Track.variables['Extreme_cases_sea'][3,:,:]

EXcaseall= EXcase1+EXcase2+EXcase3+EXcase4

ETCprecp1 = Track.variables['Extremevalues_ETC_sea'][0,:,:]
ETCprecp2 = Track.variables['Extremevalues_ETC_sea'][1,:,:]
ETCprecp3 = Track.variables['Extremevalues_ETC_sea'][2,:,:]
ETCprecp4 = Track.variables['Extremevalues_ETC_sea'][3,:,:]

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
       ETCcase= ETCprob1*EXcase1
       ETCcont= ETCprob1
       EXcase=EXcase1
       EXfreq= EXcase1/EXcaseall
       sea='JJA'
       label='(c)'
    elif s == 1:
       ETCcase= ETCprob2*EXcase2
       ETCcont= ETCprob2
       EXcase=EXcase2
       EXfreq= EXcase2/EXcaseall
       sea='SON'
       label='(f)'
    elif s == 2:
       ETCcase= ETCprob3*EXcase3
       ETCcont= ETCprob3
       EXcase=EXcase3
       EXfreq= EXcase3/EXcaseall
       sea='DJF'
       label='(i)'
    elif s == 3:
       ETCcase= ETCprob4*EXcase4
       ETCcont= ETCprob4
       EXcase=EXcase4
       EXfreq= EXcase4/EXcaseall
       sea='MAM'
       label='(l)'

    ax = fig.add_subplot(4, 1, s+1, projection=proj)

    ETCcont_mask1=np.ma.masked_where(EXfreq*100<0.5,ETCcont*100)
    ETCcont_masked=np.ma.masked_where(mask,ETCcont_mask1)

  
    ETCcont_subD = np.average(ETCcont_masked)
    print('averaged over subDomain:',ETCcont_subD)

    colors=[c1,c2,c3,c4,c5,c6]
    fig, ax = lambert_map(extent=(-133, -37, 17, 67), fig=fig,
            cent_lon=-87, ax = ax)



    if s >= 0:
       cf = ax.contourf(lon, lat, ETCcont_mask1, transform=ccrs.PlateCarree(), zorder=6,
                 colors=colors,
                 levels=[30,50,60,70,80,90,95],extend='max')
       cf2 = ax.contour(lon,lat,mask, transform=ccrs.PlateCarree(), levels=[1],
               colors='k',linewidths=2.5,zorder=7,alpha=0.7)

    cf0 = ax.contourf(lon0,lat0,elevation/9.80665, transform=ccrs.PlateCarree(), levels=[1000,5000],
               colors='none',hatches=['////'],zorder=7)

    ax.text(0.96,0.1,''+str('{:2.0f}'.format(ETCcont_subD))+'%',bbox={'facecolor':'white','alpha':0.7,'boxstyle':'round'},
              horizontalalignment='right', transform=ax.transAxes,fontsize=16,zorder=28)
    ax.text(0.01,0.985,''+label+'',bbox={'facecolor':'white','edgecolor':'None','alpha':0.9,'pad':2.0},
              horizontalalignment='left', verticalalignment='top',
              transform=ax.transAxes,fontsize=15,zorder=28,**sans_font)
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
fig.subplots_adjust(hspace=0.04)
plt.savefig('Fig3_rightcolumn', bbox_inches='tight', dpi=200)
plt.show()
plt.close()
