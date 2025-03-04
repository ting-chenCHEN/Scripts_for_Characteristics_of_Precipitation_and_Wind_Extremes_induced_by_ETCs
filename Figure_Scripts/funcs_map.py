#

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.path import Path
import matplotlib.patches as patches


def lambert_map(extent=(-82, -75, 41, 46), cent_lon =-80,figsize=(14, 12), fig = None, ax = None):

    proj = ccrs.LambertConformal(central_longitude=cent_lon, central_latitude=35,
                                 standard_parallels=[35])
    if ax == None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=proj)

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',  name='admin_1_states_provinces_lines',
        scale='50m', facecolor='none')
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='None',
                                            #edgecolor='face',
                                            #facecolor='0.9')
                                            facecolor='0.98')
    lakes_50m = cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                            edgecolor='None',
                                            facecolor=[(0.59375 , 0.71484375, 0.8828125)])
    lakes_50m_edge= cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                            edgecolor='k',
                                            facecolor='None',lw=0.8)
    #(0.59375 , 0.71484375, 0.8828125)
    ax.add_feature(land_50m); 
    #ax.add_feature(lakes_50m, zorder=3); 
    ax.add_feature(lakes_50m_edge, zorder=10,lw=0.8); 
    #ax.add_feature(cfeature.LAKES, edgecolor='white', zorder=10);
    ax.add_feature(cfeature.BORDERS, zorder=10,edgecolor='dimgray'); 
    #ax.add_feature(states_provinces, edgecolor='gray', zorder=10)
    #ax.coastlines('50m', zorder=15, color='k',lw=0.5)
    ax.coastlines('50m', zorder=15, color='k',linewidths=0.8)
    # Set plot bounds
    ax.set_extent(extent)
    return fig, ax
