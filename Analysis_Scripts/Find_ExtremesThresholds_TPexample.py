#
# Before opening Python execute:
#   module load python2
#   module load development/python-rpn
#

import xarray as xr
import numpy as np                         # Good module for matrix and matrix operation
import os                                  # Used to convert png to other format
import netCDF4

path='/home/archive/REANALYSES/ERA5/1h/tp/ll/nc4'

nc_files = ['%s/%s/%s/era5_tp_ll_%s%s_1h.nc4' % (path,str(yyyy),str(mm).rjust(2,"0"),str(yyyy),str(mm).rjust(2,"0")) 
            for yyyy in range(2001,2020) for mm in range(1,13)]

#print (nc_files)

ds = xr.open_mfdataset(nc_files, concat_dim="time", combine='nested',
                 data_vars='minimal',coords='minimal',compat='override',
                 parallel=True, autoclose=True)

# Cut into 13 patches of sub-domains:
#for i in range(0, 13):
for i in range(0, 26):
#for i in range(15, 26):
    print('i=',i)
    if (i == 0):
       #lwest=210+i*10
       lwest=210+i*5
    else:
       lwest=210+i*5+0.1
       #lwest=210+i*10+0.1
    least=lwest+5
    subtp=ds.tp.sel(latitude=slice(75,25), longitude=slice(lwest,least))
#    subtp=ds.tp.sel(latitude=slice(75,60), longitude=slice(lwest,least))

    subtp.load()
    subtp=subtp*1000.

    print('calculating p98.0')
    p980=subtp.quantile(0.98, dim="time",interpolation="lower")
    print('calculating p99.0')
    p990=subtp.quantile(0.99, dim="time",interpolation="lower")
    print('calculating p99.9')
    p999=subtp.quantile(0.999, dim="time",interpolation="lower")

    newp = xr.Dataset({})
    newp['p98p0'] = p980
    newp['p99p0'] = p990
    newp['p99p9'] = p999
    print('loading to netcdf...')
    newp.to_netcdf('TPperc_2001_2020_'+str(i)+'.nc',format="NETCDF3_64BIT")

    subtp.close()

print('merging...')

files = ['TPperc_2001_2020_%s.nc' % (str(i))
         for i in range(0,26)]

ds = [xr.open_dataset(f) for f in files]

merge980=xr.concat( [ds[i].p98p0.to_dataset() for i in range(0, 26)], dim='longitude')
merge980.to_netcdf('TPp98p0_2001_2020_ERA5.nc', format="NETCDF3_64BIT")
merge990=xr.concat( [ds[i].p99p0.to_dataset() for i in range(0, 26)], dim='longitude')
merge990.to_netcdf('TPp99p0_2001_2020_ERA5.nc', format="NETCDF3_64BIT")
merge999=xr.concat( [ds[i].p99p9.to_dataset() for i in range(0, 26)], dim='longitude')
merge999.to_netcdf('TPp99p9_2001_2020_ERA5.nc', format="NETCDF3_64BIT")

