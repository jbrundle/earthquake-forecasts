#!/opt/local/bin python

#Note about the program:

#The larger the range of lat/lons entered below, the more probabilities will need to
#be returned by the web API and the slower the program will be. 

#It takes about 1 minute per 100 lat/lon points. The number of lat/lon points is 'h' below.

import sys
import math
import matplotlib as mpl
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from numpy import *
import scipy.interpolate
import matplotlib.mlab as ml
import requests
from urllib2 import Request, urlopen, URLError
import time
from itertools import product
import datetime

Location = 'California'

NELat       = float(sys.argv[1])
NELng       = float(sys.argv[2])
SWLat       = float(sys.argv[3])
SWLng       = float(sys.argv[4])
MagLo       = float(sys.argv[5])
Location    = str(sys.argv[6])  

lat_0_center = 0.5*(NELat + SWLat)
lon_0_center = 0.5*(NELng + SWLng)

MagLo = 5.5
Forecast_Mag = '6.0'
Forecast_Rad = '100'
Forecast_Wdw = '12'
Months       = '12'
Days         = '365'

    #......................................

print ' '
print '     Forecast Magnitude is: ', Forecast_Mag
print '     Change it? (y/n)'
resp    =   raw_input()
if resp == 'y':
    print ' '
    print '     Enter Forecast Magnitude: '
    Forecast_Mag = raw_input()

print '     Forecast Radius is: ', Forecast_Rad, 'km'
print '     Change it? (y/n)'
resp    =   raw_input()
if resp == 'y':
    print ' '
    print '     Enter Forecast Radius (km): '
    Forecast_Rad = raw_input()

print '     Forecast Time Window is: ', Forecast_Wdw, 'months'
print '     Change it? (y/n)'
resp    =   raw_input()
if resp == 'y':
    print ' '
    print '     Enter Forecast Time Window (months): '
    Forecast_Wdw        = raw_input()
    Months              = Forecast_Wdw
    Forecast_Months     = float(Forecast_Wdw)
    Forecast_Days       = Forecast_Months * 365.0/12.0
    Days                = str(Forecast_Days)

resp_map = ''
print ' '
print '     Default map is eTopo.  Change to Other map? (y/n)'
resp_map    =   raw_input()
if resp_map == 'y':
    print ' '
    print '     Shaded Relief (s) or Flat (f = No Relief)...?'
    resp_map_type   = raw_input()
 
    #......................................

lower_limit_flag = 'n'    # Flag to allow lower bound on forecast probability

    #......................................



todays_time_and_date = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
todays_date          =  datetime.date.today().strftime("%B %d, %Y")

print ' '
print ' Todays Time and Date: ', todays_time_and_date


latres = 0.5
lonres = 0.5

if Location == 'World':
    latres=2.0
    lonres=2.0

llcrnrlat = SWLat
llcrnrlon  = SWLng
urcrnrlat  = NELat
urcrnrlon  = NELng

    #   Variable y is coordinates of latitude, (min_lat, max_lat, increments)
    #   Variable x is coordinates of longitude, (min_lon, max_lon, increments)
    #   HTTP error 500 will occur if values are entered that include a lon or lat of 0


y =  np.arange(SWLat-2,NELat+2, latres)
x = np.arange(SWLng-2,NELng+2, lonres)

    #   Variable c joins x and y variables into latitude/longitude pairs

c = list(product(x, y))

nlat = int((NELat-SWLat+4)/latres)
nlon = int((NELng -SWLng+4)/lonres)

    #   Variable d makes c into a list

d = list(c)

    #   Variable h is the length of c

h = len(c)

    #   Calls the web API and causes it to return only the probability for each
    #   latitude/longitude pair, then assign the value to the probability array:
    #   - probability -

probability = np.zeros((nlat,nlon))

ix=0
iy=0
          
calls = 0
h=h-1

print ' '
print ' Total Number of API Downloads Will Be: ', h

#while calls < h:
for calls in range(h):    
#   calls = calls + 1
#   print ' Completed API Download Number: ', calls+1, "                                 \r",
    percent_complete = 100.0*float(calls+1)/float(h) 
    print ' API Downloads Completed: ', "{0:.2f}".format(percent_complete),'%', "                                 \r",   
    sys.stdout.flush()
    C = calls
    f = str(d[C-1])
    K = f.strip('() ')
    J = K.replace(" ", "")
    url = 'http://api.openhazards.com/GetEarthquakeProbability?q=' + J +'&m=' + Forecast_Mag + '&r=' + Forecast_Rad + '&w=' + Days
    request = Request(url)
    try:
        response = urlopen(request)
        probs = response.read()
        probs = probs.translate(None, '%')
        probs = probs[probs.find('<prob>')+6:probs.find('</prob>')]
        probability[ix][iy] = probs
        ix += 1
        ix = ix%(nlat)
        if ix == 0:
            iy += 1    # increment iy if ix is reduced back to 0 by the mod operator %
    
    except URLError, e:
        print 'error', e
    
    if calls >= h:
        break

print ' '
print ' Finished API Downloads, Now Generating Contour Map'
print ' '

    #   Create a basemap on which to plot probability data

map = Basemap(projection='cyl',llcrnrlat=SWLat, urcrnrlat=NELat,
            llcrnrlon=SWLng, urcrnrlon=NELng, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')
#map.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')   # Contours will not display on land if this is used

map_params_meridians_min    = 0
map_params_meridians_max    = 360
map_params_meridians_delta  = 4

map_params_parallels_min    = -90
map_params_parallels_max    =  90
map_params_parallels_delta  =  2

if Location == 'World':
    map_params_meridians_min    = 0
    map_params_meridians_max    = 360
    map_params_meridians_delta  = 60

    map_params_parallels_min    = -90
    map_params_parallels_max    =  90
    map_params_parallels_delta  =  30

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
map.drawmeridians(np.arange(map_params_meridians_min, map_params_meridians_max, map_params_meridians_delta),labels=[0,0,0,1])
map.drawparallels(np.arange(map_params_parallels_min, map_params_parallels_max, map_params_parallels_delta),labels=[1,0,0,0])
map.drawstates()
#map.drawrivers()

map.etopo()
if resp_map == 'y':
    if resp_map_type == 's':
        map.shadedrelief()
    if resp_map_type == 'f':
        map.drawlsmask(land_color='#ffe8a0', ocean_color='#e6e6ff')  # Use this if map.etopo() is  not used
        map.drawmapboundary(fill_color='#e6e6ff')

    #   Create Lat and Lon 1-d arrays for plotting

y =  np.arange(SWLat-2,NELat+2, latres)
x = np.arange(SWLng-2,NELng+2, lonres)


print ' '

    #   Some data for configuring the plot

probmax = int(np.amax(probability) + 2)

if probmax > 90:
    probmax = 101

if probmax < 15:
    probmax = 16

levels = np.arange(1,probmax+1,1)
ticklevels = [1,5,10, 20, 30, 40, 50, 60, 70, 80, 90]

if lower_limit_flag == 'y':
    levels = np.arange(25,probmax+1,1)              # Lower limit currently set to 25%
    ticklevels = [25, 30, 40, 50, 60, 70, 80, 90]   # Change tick levels for color contours

    #   Contour lat, lon, and probability data so that it appears on basemap. cmap=colormap.

cs = plt.contourf(x, y, probability, levels, cmap = 'jet', latlon='True', alpha=0.40)

plt.colorbar(ticks=ticklevels, shrink=0.55)

Sup_title_text = 'Earthquake Probability (%) for M>' + Forecast_Mag + '\nWithin ' + Months + ' months and ' + Forecast_Rad + ' km Radial Distance in ' + Location
plt.suptitle(Sup_title_text, fontsize = 14)

Title_text = todays_date

plt.title(Title_text, fontsize=12)

plt.savefig('Forecast-Contours.png')

print ' '
print '     Close plot window to continue...'
print ' '
print '     .......................................'

                 
plt.show()







