# Code to take identified eddies and tracks
# and calculate statistics based on 1-degree by 1-degree bins

from netCDF4 import Dataset
import pandas as pd
import numpy as np
import datetime
import time

start_time = time.time()

## User inputs
nc_file = 'cosmo-output/tracks/Anticyclonic.nc' #User enters in file needed
cyc = False  #False if anticyclonic, True if cyclonic

## Constants
day_ref = 1721424 #in days, corresponds to year 1, month 1, day 1 (01-01-01)
date_ref = datetime.datetime(year=1, month=1, day=1)

## Functions
def jdate_conversion(jdate):
    delta_t = datetime.timedelta(days = jdate - day_ref)
    delta_t
    date = date_ref + delta_t
    date_str = datetime.datetime.strftime(date, '%y-%m-%d')
    return date  #convert from a Julian date to datetime object

def df_from_nc(nc_file, date_fn = None):
    ds = Dataset(nc_file)
    lat = list(ds['lat'][:])  #latitude in degrees
    lon = list(ds['lon'][:])  #longitude in degrees
    radius = list(ds['radius_e'][:]) #radius in km
    amp = list(ds['A'][:])  #amplitude in m (not cm)
    date = list(ds['j1'][:]) #Julian Date, with first date (1721687) corresponding to 01-09-21 (year 1 of model, Sept. 21)
    eddy = list(ds['track'][:])  #eddy ID
    ext_contour_height = list(ds['height_external_contour'][:]) #height at outer contour
    inn_contour_height = list(ds['height_inner_contour'][:]) #height at innermost contour of eddy
    obs_number = list(ds['n'][:])
    df = pd.DataFrame([eddy, date, lat, lon, radius, amp, ext_contour_height, inn_contour_height,obs_number]).T
    df.columns = ['eddy_ID', 'date', 'lat', 'lon', 'radius', 'amp', 'outer_height', 'inner_height','obs_number']
    if date_fn:
        df['date'] = df['date'].apply(date_fn)
    df['eddy_ID'] = df['eddy_ID'].astype(int)
    ds.close()
    return df  #convert netcdf into pandas dataframe
    
def prev_month(tuple):
    if tuple[0] == 1: #i.e. January
        return (12, tuple[1] - 1)
    else:
        return (tuple[0] - 1, tuple[1])
        
def calc_duration(df):
    """
    This function calculates the number of days that a unique eddy was present in a lat/lon bin in a given month_year
    Note: if/else statement addresses scenario in which eddy is seen in one spatial bin for two or months in a row
    i.e., if it was obs. Jan 29 and Feb. 2 in the same spatial bin, two days of observation must be added to the dataframe 
    for January and 3 to the dataframe for February (algorithm otherwise doesn't account for 4 days between (not inclusive) Jan. 29 and Feb. 3,
    even though we know what bin the eddy was located in)
    """
    agg = df.groupby(['month-year', 'lat_lon_bins','eddy_ID']).agg({'eddy_ID':'count'}) #aggregate counts of unique eddies
    agg = agg.rename(columns = {'eddy_ID':'eddy_obs'}) 
    agg_copy = agg.copy()
    agg['duration'] = 0
    df1 = df.set_index(['month-year','lat_lon_bins','eddy_ID'])  #allows for easier indexing
    agg['duration'] = agg.apply(lambda row: (row - 1)*5 + 1) #1 obs/month -> 1 day, 2 obs -> 6 days duration, etc.
    for index, obs in agg_copy.itertuples():  #iterate over a copy
        my, latlon, ID = index
        try:
            date_prev = np.max([df1.at[(prev_month(my),latlon, ID),'date']]) #was eddy present in latlon bin in previous month?
        except KeyError:
            date_prev = None        
        if date_prev is not None: 
            date = np.min([df1.at[(my,latlon,ID),'date']])
            if (date - date_prev) == 5:  ##i.e. two consecutive obs
                date_delta = date - datetime.datetime(date.year, date.month, day = 1)
                agg.at[(my, latlon, ID),'duration'] += date_delta.days
                agg.at[(prev_month(my), latlon, ID),'duration'] += (4 - date_delta.days)
    return agg

## Main routine
df = df_from_nc(nc_file, date_fn = jdate_conversion)  #Load in a pandas dataframe]
df_binned = df.copy()
df_binned['month-year'] = df_binned.apply(lambda row: (row['date'].month, row['date'].year), axis = 1)
df_binned['lat_lon_bins'] = df_binned.apply(lambda row: (int(row['lat']), int(row['lon'])), axis = 1)

##Calculate duration statistics
duration = df_binned[['month-year','lat_lon_bins','eddy_ID','date']]
agg_dur = calc_duration(duration)
agg_dur.reset_index(inplace = True)  #ungroup and then regroup on appropriate columns
duration_stats = agg_dur.groupby(['lat_lon_bins', 'month-year']).agg({'duration': 'sum'})

#Calculate spatial/temporal statistics for all params except duration
bins_groupby = df_binned.groupby(['lat_lon_bins', 'month-year']) #group by lat/lon and months
stats = bins_groupby.agg({'eddy_ID': lambda e:e.nunique(), 'month-year': 'first', #get unique nObs, count of dates, and average everything else
    'lat_lon_bins': 'first', 'radius': np.mean, 'amp': np.mean, 
    'outer_height': np.mean, 'inner_height': np.mean})
stats.rename(columns = {'eddy_ID': 'obs'}, inplace = True)

#Join duration stats with all other params
stats = stats.join(duration_stats)

#Create output dataframe axes
params = ['obs','radius', 'amp', 'outer_height','inner_height','duration']

for param in params:
    matrix = stats[param].unstack(-1)
    matrix = matrix.reindex(sorted(matrix.columns, key=lambda label: (label[1],label[0])), axis = 1)
    matrix.to_csv("NWAT_eddy_stats_{0}_{1}.csv".format(param, (lambda x: 'cc' if x else 'ac')(cyc)))


