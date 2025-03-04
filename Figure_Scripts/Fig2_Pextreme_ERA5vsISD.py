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
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap

import funcs_map
from funcs_map import lambert_map
from matplotlib.path import Path
import matplotlib.patches as patches
import netCDF4
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
import datetime as dt

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

lat0 = Topo.variables['latitude'][:]
lon0 = Topo.variables['longitude'][:]

# ISD 

ISD = 'ISD_extremes_p99p0.nc4'
ds=xr.open_dataset(ISD)
latsI=ds.variables['lat']
lonsI=ds.variables['lon']

exces=ds.groupby('time.season').count(dim='time')
tot=ds.count(dim='time')

# ERA5 (regridded to ISD)
reference_date = dt.datetime(2001, 1, 1, 0, 0, 0)  # YYYY, MM, DD, HH, MM, SS
filet='ERA5_regridded_ISDstation_TP_p99p0.nc4'
da=xr.open_dataset(filet)
lats=da.variables['latitude']
lons=da.variables['longitude']
hours_since=da.variables['Extreme_time'].values
vals=da.variables['Extreme_time'].values

sea_months={'DJF':[12,1,2],'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11]}
#------------------------------------
fig = plt.figure(figsize=(6.5,8.4))
proj = ccrs.LambertConformal(central_longitude=-100, central_latitude=35,
                                 standard_parallels=[35])

seasons=[0,1,2,3]
for s in seasons: 
    if s == 0:
       sea='JJA'
       label='(c)'
       label2='(d)'
    elif s == 1:
       sea='SON'
       label='(g)'
       label2='(h)'
    elif s == 2:
       sea='DJF'
       label='(k)'
       label2='(l)'
    elif s == 3:
       sea='MAM'
       label='(o)'
       label2='(p)'

    # counting the fraction 
    zzera5=np.zeros((hours_since.shape[1]))
    for jj in range(hours_since.shape[1]):
        temp=[]
        for ii in range(hours_since.shape[0]):
            date=reference_date + dt.timedelta(hours=int(hours_since[ii,jj]))
            if date.month in sea_months[sea]:
                temp.append(1)
        zzera5[jj]=100.*np.sum(temp)/hours_since.shape[0]

    ax = fig.add_subplot(4, 2, s*2+1, projection=proj)


    colors=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10]

    fig, ax = lambert_map(extent=(-124, -74, 23, 51), fig=fig,
            cent_lon=-150, ax=ax)
    levels=[0.5,2,5,10,15,20,30,40,50,60]
    cmap = plt.cm.get_cmap("Spectral_r")
    norm = BoundaryNorm(levels, cmap.N)
    cf = plt.scatter(lons, lats, c=zzera5, s=20, cmap=cmap,norm=norm, edgecolor='gray',linewidth=0.5,
            transform=ccrs.PlateCarree(),zorder=11,alpha=0.9)
    cf0 = ax.contourf(lon0,lat0,elevation/9.80665, transform=ccrs.PlateCarree(), levels=[1000,6000],
               colors='silver', zorder=7)

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
    if s==0:
       ax.text(0.5,1.03,'ERA5',horizontalalignment='center',transform=ax.transAxes,fontsize=13,zorder=24,weight='bold',**mono_font)
    ax.text(0.01,0.01,sea,horizontalalignment='left',transform=ax.transAxes,fontsize=12,zorder=24,weight='bold',**mono_font)
    
    ax.text(0.013,0.981,''+label+'',bbox={'facecolor':'white','edgecolor':'None','alpha':0.9,'pad':2.0},
              horizontalalignment='left', verticalalignment='top',
              transform=ax.transAxes,fontsize=11,zorder=24,**sans_font)

#--------------------------------------------------------
    zz=100.*exces.sel(season=sea).variables['pr'].values/tot.variables['pr'].values

    ax = fig.add_subplot(4, 2, s*2+1+1, projection=proj)
    fig, ax = lambert_map(extent=(-124, -74, 23, 51), fig=fig,
            cent_lon=-150, ax=ax)
    levels=[0.5,2,5,10,15,20,30,40,50,60]
    cmap = plt.cm.get_cmap("Spectral_r")
    norm = BoundaryNorm(levels, cmap.N)
    cf = plt.scatter(lonsI, latsI, c=zz, s=20, cmap=cmap,norm=norm, edgecolor='gray',linewidth=0.5,
            transform=ccrs.PlateCarree(),zorder=11, alpha=0.9)
    cf0 = ax.contourf(lon0,lat0,elevation/9.80665, transform=ccrs.PlateCarree(), levels=[1000,6000],
               colors='silver', zorder=7)

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
 
    if s==0:
       ax.text(0.5,1.03,'ISD',horizontalalignment='center',transform=ax.transAxes,fontsize=13,zorder=24,weight='bold',**mono_font)
    ax.text(0.01,0.01,sea,horizontalalignment='left',transform=ax.transAxes,fontsize=12,zorder=24,weight='bold',**mono_font)
    ax.text(0.013,0.981,''+label2+'',bbox={'facecolor':'white','edgecolor':'None','alpha':0.9,'pad':2.0},
              horizontalalignment='left', verticalalignment='top',
              transform=ax.transAxes,fontsize=11,zorder=24,**sans_font)

    
fig.subplots_adjust(hspace=0.01,wspace=0.)

cax=plt.axes([0.16,0.085,0.65,0.012])
cb = plt.colorbar(cf, cax=cax, ticks=levels, orientation='horizontal') #, ticks=np.arange(0,3,1))
cb.ax.tick_params(labelsize=12)
cb.ax.set_xticklabels(levels[:])
cb.ax.set_title('(%)', x=1.05, y=-1, fontsize=12)
plt.savefig('Fig2_Precip_ex990values_freq_ISD_vs_ERA5', bbox_inches='tight', dpi=200)
plt.show()
plt.close()
da.close()
ds.close()
