#
# Before opening Python execute:
#   module load python3
#   development/python3-rpn-code
#   source activate base_plus
#

import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np                         # Good module for matrix and matrix operation
import matplotlib as mpl
import matplotlib.pyplot as plt            # Module to produce figure
import matplotlib.colors as mcolors
import os                                  # Used to convert png to other format

import funcs_map
from funcs_map import lambert_map
from matplotlib.path import Path
import matplotlib.patches as patches
import netCDF4
from datetime import datetime


YUL_lat= 45.25
YUL_lon= 286.25


cmap_data = ["#98CBFE", "#0498FE", "#0F2CFF",
             "#00FF65", "#00CC01", "#009801", "#006500",
             "#FFFE31", "#FECC00", "#FF9800", "#FF6500",
             "#FF0298", "#9833CB", "#650098",
             ]

clevs = [ 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0,
          4.5, 6.0, 8.0, 10.0, 13.0, 16.0  ]

clevs_ano = np.arange(-15,17.5,2.5)
clevs_ano2 = np.arange(-9,10.5,1.5)

cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
norm = mcolors.BoundaryNorm(clevs, cmap.N)


data = open('HalloweenStorm_track.txt','r')
snum=6265

P990data=netCDF4.Dataset('TPp99p0_2001_2020_ERA5.nc','r')
Precpdata=netCDF4.Dataset('ERA5_2019OctNov_TP_10mUV.nc','r')
MSLdata=netCDF4.Dataset('ERA5_2019OctNov_mslp.nc','r')

lonp = Precpdata.variables['longitude'][:]  
latp = Precpdata.variables['latitude'][:] 

nctime = Precpdata.variables['time'][:] # get values
t_unit = Precpdata.variables['time'].units
t_cal = Precpdata.variables['time'].calendar

num=[]
point=[]
date=[]
lat=[]
lon=[]
r_pres=[]
avgpr=[]

for i, line in enumerate(data):
    words=line.split()
    if i >=2 :
          num0   =float(words[0])
          point0 =int(words[2])
          date0  =words[3]
          lat0   =float(words[4])
          lon0   =float(words[5])
          if (num0 == snum):
            num.append(num0)
            point.append(point0)
            date.append(date0)
            lat.append(lat0)
            lon.append(lon0)

tvalue = netCDF4.num2date(nctime, units = t_unit, calendar = t_cal)
#Find the correspondint time to the ERA5 nc file
def get_verts_new(categ_box):
    verts = [(categ_box['west'], categ_box['south']), 
     (categ_box['west'], categ_box['north']), 
     (categ_box['east'], categ_box['north']), 
     (categ_box['east'], categ_box['south']), 
     (0., 0.)]
    return verts

box = {'west':-180, 'east':0, 'north':90, 'south':20}

for pointnow in [25,36]:
    datetime = datetime.strptime(date[pointnow], '%Y%m%d%H')
    index= int(np.where(tvalue==datetime)[0])

    tpe  = Precpdata.variables['tp'][index,:,:]
    msl  = MSLdata.variables['msl'][index,:,:]
    p990 = P990data.variables['p99p0'][:,:]

    actp = tpe*1000.-p990


    fig, ax = lambert_map(extent=(-90.6, -59, 31, 57), figsize=(7,7))
    cf0 = ax.contourf(lonp, latp, actp[:,:], 
                      clevs_ano2,
                      cmap=plt.cm.get_cmap('RdBu_r'),
                      extend='both',
                      transform=ccrs.PlateCarree(), zorder=4)#, 
    cf1 = ax.contour(lonp, latp, actp[:,:], 
                    [0],colors='gold',linewidths=1.2,
                    transform=ccrs.PlateCarree(), zorder=4)#, 
    cf2= ax.contour(lonp, latp, msl*1./100.,colors='k',transform=ccrs.PlateCarree(),zorder=18,linewidths=1.0)
    ax.clabel(cf2,fmt='%4.0f',fontsize=14)
    cf=plt.scatter(lon[pointnow],lat[pointnow],c='k',s=100,transform=ccrs.PlateCarree(),zorder=19)
    cf=plt.scatter(lon[pointnow],lat[pointnow],c='r',s=50,transform=ccrs.PlateCarree(),zorder=19)
    cf1=plt.scatter(YUL_lon,YUL_lat,c='k',s=320,marker='*',transform=ccrs.PlateCarree(),zorder=19)
    cf1=plt.scatter(YUL_lon,YUL_lat,c='w',s=150,marker='*',transform=ccrs.PlateCarree(),zorder=19)
    ax.tissot(rad_km=1000, lons=lon[pointnow],lats=lat[pointnow],n_samples=36, edgecolor='crimson',linewidths=5,facecolor='none',zorder=20,alpha=0.5)

    cbar= fig.colorbar(cf0, orientation="horizontal", pad=0.05, shrink=0.8)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig('Fig7_HalloweenStorm_TPexceedance_point'+str(pointnow)+'', bbox_inches='tight', dpi=200)
    plt.show()
