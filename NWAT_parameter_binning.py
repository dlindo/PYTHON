# Code to take ROMS modeled data in the NWAT region and bin it by
# spatial location and mm/yy. Parameters include eddy kinetic energy
# and vorticity. 

from netCDF4 import Dataset
import pandas as pd
import numpy as np
import datetime
import time
from collections import defaultdict
from collections import OrderedDict
import warnings

## Convert VRT from psi to rho grid
## Two approaches: (1) calculate lat/lon grid for psi
## (2) Average/calculate vrt values for rho grid

def rho2psi(array, dim = '2D'):
    if dim == '2D':
        array_psi = .25*(array[:-1, :-1] + array[:-1, 1:] + array[1:,:-1] + array[1:,1:])
    elif dim == '3D':
        array_psi = .25*(array[:,:-1, :-1] + array[:,:-1, 1:] + array[:,1:,:-1] + array[:,1:,1:])
    return array_psi
    
## Also define rho to u and v grid conversions!

def rho2v(array, dim = '2D'):
    if dim == '2D':
        array_v = .5*(array[:-1,:] + array[1:,:])
    elif dim == '3D':
        array_v = .5*(array[:,:-1,:] + array[:,1:,:])
    return array_v

def rho2u(array, dim = '2D'):
    if dim == '2D':
        array_u = .5*(array[:,:-1] + array[:,1:])
    elif dim == '3D':
        array_u = .5*(array[:,:,:-1] + array[:, :,1:])
    return array_u
    
## Note: Because both grids are curvillinear, distances between points vary. This, coupled 
## with the fact that the rho grid contains the psi grid, means creating a "psi2rho" function would either 
## shrink the rho grid (by removing the edge points) or require inexact extrapolation to the edges of the grid

nc_file_format = '../vrt/nwat_sfc_vrt.{:04d}.nc'
grid = '../../../NWAT/nwat_grd_latlon.nc'
date_ref = datetime.date(year = 1, month = 9, day = 21) # Sept. 21, 0001 ("first" day of model)
grid_ds = Dataset(grid)
lat_rho = grid_ds.variables['lat_rho'][:].data
lon_rho = grid_ds.variables['lon_rho'][:].data
lat_psi = rho2psi(lat_rho)  
lon_psi = rho2psi(lon_rho)
parameter = 'ke'

#if parameter is vrt, use psi grid, if ke, use rho grid
if parameter == 'ke':
    lat, lon = lat_rho, lon_rho
elif parameter == 'vrt':
    lat, lon = lat_psi, lon_psi

## Create dictionary of indices for each spatial bin
minlat, maxlat = int(np.min((lat))),int(np.max(lat))
minlon, maxlon = int(np.min((lon))),int(np.max(lon))


d_spatial = {}

for i_lat in range(minlat, maxlat):
    for i_lon in range(minlon, maxlon):
        latlat = (lat >= i_lat) & (lat < i_lat + 1)
        lonlon = (lon >= i_lon) & (lon < i_lon + 1)
        intrsct = latlat & lonlon
        if np.sum(intrsct) > 0:
            d_spatial[(i_lat, i_lon)] = intrsct
            
od_spatial = OrderedDict(d_spatial)
indices = list(od_spatial.values())
## Loop through each vrt file and set of observations, take mean of vrt for each
## spatial bin, and add to a dataframe
def param_binning(parameter):
    df = pd.DataFrame({'latlon_bins': list(od_spatial.keys())})
    df.set_index('latlon_bins', inplace = True)
    for obs in range(0, 1230, 5):  # number of unique vrt files
        ds = Dataset(nc_file_format.format(obs))
        param = ds.variables[parameter][:].data
        param[param == 0] = np.nan  #If 0, actually coastline i.e. masked
        ds.close()
        param_ls = []
        for val in indices:
            ind = param[:,val]
            mean_array = np.nanmean(ind, axis = 1)
            param_ls.append(mean_array)
        if obs == 1225:  #only four observations in the last data file?
            file_obs = 4
        else:
            file_obs = 5
        for i in range(file_obs):
            days_passed = 5*(obs + i)  #i.e. in file 5, obs1, 30 days have passed
            mmddyy = date_ref + datetime.timedelta(days = days_passed)
            df[mmddyy] = [val[i] for val in param_ls]
    return df


start_time = time.time()
df = param_binning(parameter)
print(time.time() - start_time)

df = df.dropna() #If NA, it's masked coastline/land

## Transpose for grouping, and add month-year info
dfT = df.T
dfT['month-year'] = dfT.index.map(lambda x: (x.month, x.year))
dfT.reset_index(drop = True, inplace = True)
df_agg = dfT.groupby('month-year').agg('mean')
df_agg = df_agg.T

df_agg.to_csv('../../data_matrices/{}_binned.csv'.format(parameter))
