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
from matplotlib import gridspec,rc
import funcs_map
from funcs_map import open_var_2d, lambert_map
from matplotlib.path import Path
import matplotlib.patches as patches
import netCDF4

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

fs=10

def get_verts_new(categ_box):
    verts = [(categ_box['west'], categ_box['south']), 
     (categ_box['west'], categ_box['north']), 
     (categ_box['east'], categ_box['north']), 
     (categ_box['east'], categ_box['south']), 
     (0., 0.)]
    return verts

box = {'west':-180, 'east':0, 'north':80, 'south':20}

WSc=netCDF4.Dataset('WSp99p0_2001_2020_ERA5.nc','r')
TPc=netCDF4.Dataset('TPp99p0_2001_2020_ERA5.nc','r')


ig=305 #from 210, 0.25 interval
jg=119 #from 75N, 0.25 innterval
WS99=WSc.variables['p99p0'][jg,ig]
TP99=TPc.variables['p99p0'][jg,ig]
lat=TPc.variables['latitude'][jg]
lon=TPc.variables['longitude'][ig]


#Halloween storm, original number: 11649, new number for newNNA: 2312
Track=netCDF4.Dataset('Timeseries_WSnTP_local_record_Montreal_2019HalloweenStorm.nc','r')
lat2=Track.variables['latitude'][jg]  #46.75
lon2=Track.variables['longitude'][ig] #288.75
print(lat2)
print(lon2)

liststorm=[1]
extremeTPh_list=[]
stormnum=0
for n in liststorm:

    TP0 = Track.variables['TP_series'][:,0,jg,ig]
    px= TP0.shape[0] #include missing vlaues
    pvalid=[]
    for p in np.arange(px):
        if (TP0[p]):
            pvalid.append(p)

    extremeTPh= 0 
    if (bool(pvalid)): #FALSE if the list is empty
      pmin= min(pvalid)
      pmax= max(pvalid)+1
            
      print(pmin)
      print(pmax)
      Dist = Track.variables['Distance_to_storm'][pmin:pmax,0,jg,ig]
      Dirc = Track.variables['Re_direction_of_storm'][pmin:pmax,0,jg,ig]

      TP = Track.variables['TP_series'][pmin:pmax,0,jg,ig]
      WS = Track.variables['WS_series'][pmin:pmax,0,jg,ig]
      lifetime = Track.variables['lifetime'][pmin:pmax,0,jg,ig]
      print (lifetime)
      extremeTPh= 0 
      for p in np.arange(pmax-pmin):
         if TP[p]>=TP99:
            if extremeTPh==0: 
               extremeTPh=extremeTPh+1
            else:
               if TP[p-1]>=TP99: # continuous   
                  extremeTPh=extremeTPh+1
               elif TP[p-2]>=TP99: # only skip one hour
                  extremeTPh=extremeTPh+1
               else:
                  extremeTPh_list.append(extremeTPh)  # record the previous one first
                  extremeTPh =0
                  extremeTPh=extremeTPh+1

               
      #break
    if extremeTPh >= 0:
        extremeTPh_list.append(extremeTPh)  # record the previous one first
        stormnum+=1
        print ('extremeTPh list=',extremeTPh_list)
        hours=pmax-pmin
        point=np.arange(pmin,pmax)
        fig = plt.figure(figsize=(4,7))
        gs = gridspec.GridSpec(3,1)
        ax0=plt.subplot(gs[0])
        ax0.plot(lifetime,Dist*1./1000.,color='k',alpha=0.8)   
        ax0.set_title('Distance and direction of ETC',loc='left',fontdict={'fontsize': fs})
        ax0.set_ylabel('(1000 km)',fontsize=fs)
        ax0.set_yticks([0.5,1,1.5,2])
        ax0.set_xlim(0,58)
        ax0.set_xticks(np.arange(6,58,6))
        plt.grid(axis='x',linestyle=':')
        ax02=ax0.twinx()
        for xx in range( len(point) ):
            if Dirc[xx]==8:
                Dirc[xx]=nan
            elif Dirc[xx]>=4:
                Dirc[xx]=Dirc[xx]-4
            else:
                Dirc[xx]=Dirc[xx]+4

        ax02.scatter(lifetime,Dirc,s=10,c='k',alpha=0.8)   
        ax02.set_yticks([0,1,2,3,4,5,6,7])
        ax02.set_yticklabels(['S','SW','W','NW','N','NE','E','SE'])
        ax02.yaxis.label.set_color('k')
        ax02.tick_params(axis='y',colors='k')
        ax02.set_xlim(0,58)
        
        for xx in range( len(point) ):
            if Dist[xx] <=1000:
               ax0.axvspan(xx-0.5, xx+0.5, alpha=0.25,color='gray',lw=0, zorder=0)

        ax1=plt.subplot(gs[1])
        ax1.bar(lifetime,TP,color='cornflowerblue',alpha=0.8)  
        print('TP=',TP-TP99)

        ax1.axhline(y=TP99,color='blue')
        ax1.text(-1.,TP99,"{:.1f}".format(TP99),color='blue',ha='right',va='center')
        ax1.set_title('Precipitation rate',loc='left',fontdict={'fontsize': fs})
        for xx in range( len(point) ):
            if Dist[xx] <=1000:
               ax1.axvspan(xx-0.5, xx+0.5, alpha=0.25,color='gray',lw=0,zorder=0)
        ax1.set_xlim(0,58)
        ax1.set_yticks([3,6,9])
        ax1.set_ylabel('(mm/hr)',fontsize=fs)
        ax1.set_xticks(np.arange(6,58,6))
        plt.grid(axis='x',linestyle=':')

        ax2=plt.subplot(gs[2])
        ax2.bar(lifetime,WS,color='crimson',alpha=0.8)   
        ax2.axhline(y=WS99,color='red')
        print('WS=',WS-WS99)
        ax2.text(-1.,WS99,"{:.1f}".format(WS99),color='red',ha='right',va='center')
        ax2.set_title('10-m wind speed',loc='left',fontdict={'fontsize': fs})
        ax2.set_xlabel('Cyclone lifetime (hour)',fontsize=fs)
        for xx in range( len(point) ):
            if Dist[xx] <=1000:
               ax2.axvspan(xx-0.5, xx+0.5, alpha=0.25,color='gray',lw=0,zorder=0)
        ax2.set_xlim(0,58)
        ax2.set_ylabel('(m/s)',fontsize=fs)
        ax2.set_xticks(np.arange(6,58,6))
        plt.grid(axis='x',linestyle=':')

        plt.subplots_adjust(hspace=0.4)
    
        plt.savefig('./Fig7_efg_time_series_TP_WS_duringETC_2019HalloweenStorm_Montreal', bbox_inches='tight', dpi=150)
        plt.show()
        plt.close()
