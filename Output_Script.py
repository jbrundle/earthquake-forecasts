##!/opt/local/bin python
#import sys
#sys.path.reverse()
##############################################################################

import sys
import os
import urllib
import datetime
import dateutil.parser# Open the data file

data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")

print 'Data from the ANSS Catalog for Large Earthquakes in the Selected Region'
print ' '

print '#  DateString    Time     DecimalDate    Long.   Lat.   Mag  Depth'
print ' '

i=0
for line in data_file:
    i+=1
    items = line.strip().split()

    date_string         = items[0]
    time_string         = items[1]
    ts                  = items[2]
    lon                 = items[3]
    lat                 = items[4]
    mag                 = items[5]
    dep                 = items[6]

    event = str(i)
    if i <10:
        event = '0'+str(i)
    print event, date_string, time_string, ts, lon, lat, mag, dep

# Finalize the output file

print ' '
print '#  DateString    Time     DecimalDate    Long.   Lat.   Mag  Depth'

data_file.close()

