#!/opt/local/bin python

    #   Python Menu Code to call the various procedures to plot EQ data and forecasts
    #
    #   This code downloads data from various web sites and uses matplotlib to construct plots

import sys
import os
import numpy as np
from array import array

from EQUtilities import *

    #.....................................
    #
    # Read pre-defined locations file
    #
    #.....................................


input_file = open("locations.txt", "r")
i=0
for line in input_file:
    i +=  1
input_file.close()  # Put the file pointer back at top

number_locations = i

    # Create arrays of length i filled with zeros

Location_file   = ["" for x in range(i)]
MinLat_file     = np.zeros(i)
MaxLat_file     = np.zeros(i)
MinLng_file     = np.zeros(i)
MaxLng_file     = np.zeros(i)

input_file = open("locations.txt", "r")

i=-1
for line in input_file:
    i+=1
    line    = line.strip()
    items   = line.split(',')

    items_array = np.asarray(items)

    Location_file[i]   = items_array[0]
    MinLat_file[i]     = items_array[1]
    MaxLat_file[i]     = items_array[2]
    MinLng_file[i]     = items_array[3]
    MaxLng_file[i]     = items_array[4]

Location    = Location_file[0]
SWLat       = MinLat_file[0]
NELat       = MaxLat_file[0]
SWLng       = MinLng_file[0]
NELng       = MaxLng_file[0]

Last_location = Location

input_file.close()  # Put the file pointer back at top

    #.....................................
    
MagLo = 5.0 # Initial default minimum magnitude

    #.....................................

#   Assume initially that the local region is a circle

completeness_mag    =   2.99
Circle_Location     =   'None'
Circle_Lat          =   0.0
Circle_Lng          =   0.0
Radius_float        =   0.0
earthquake_depth    =   1000.0
region_type         =   'Circle'

settings_params     =   []
settings_params.append(region_type)
settings_params.append(completeness_mag)
settings_params.append(earthquake_depth)
settings_params.append(Circle_Location)
settings_params.append(Circle_Lat)
settings_params.append(Circle_Lng)
settings_params.append(Radius_float)

save_settings(settings_params)

settings_params = get_settings()

    #.....................................

#   Set the intial values of completeness_mag and Circle_Location here
#   Write these to a file "current_settings.txt"

    #.....................................
print ' '
print 'Downloading default data set'
os.system("python generate_catalog_file.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    #.....................................

rmenu =  ''
print ' '
print ' This is the basic code to set up the menu of choices'


while rmenu != 'Z':
    print ' '
    print '     ...........................................'
    print ' '
    print '     Current location is: ', Location
    print ' '
    print '     A: Input Location and Catalog Parameters'
    print '     B: Generate or Re-Generate Working Catalog File'
    print '     C: Edit Working Catalog File (for Macs Only)'
    print '     D: Plot EQ Magnitude vs. Time'
    print '     E: Plot EQ Epicenters on a Map'
    print '     F: Plot NTW Earthquake Forecast on a Map'
    print '     G: Plot Gutenberg-Richter Relation for an Inter-Earthquake Sequence'
    print '     H: Plot Earthquake NowCast for Regional Area'
    print '     I: Plot Earthquake NowCast for Circle within Regional Area'
    print '     J: Plot Earthquake NowCast for Polygon within Regional Area'
    print '     K: Plot Proxy Strain vs. Time in the Region'
    print '     L: Plot Proxy Strain vs. Time in the Circle'
#    print '     M: Plot Earthquake Forecast vs. Time in the Circle'
#    print '     N: Plot Filtered Earthquake NowCast for Regional Area (Same as H but Log10-Linear)'
    print ' '
    print '     Z: Get Me Out of Here!'
    print ' '
    print '     (You must choose ... but choose wisely!'
    print '     For a true choice will bring you joy,'
    print '     but a poor choice will take it from you...)'
    print ' '
    rmenu = ''
    rmenu = raw_input("Enter a Choice: \n") 
    print ' '

    #.....................................
    
    if rmenu == 'B':
        print 'Choice B: Generate or Re-Generate Working Catalog File'
        os.system("python generate_catalog_file.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'C':
        print 'Choice C: Open Working Catalog File in Text Editor'
#       os.system("open -t EQ_Working_Catalog.txt") # Only for Macs.  For other systems, substitute something like "gedit EQ_Working_Catalog.txt"
        os.system("open -e EQ_Working_Catalog.txt") # Forces file open with Mac app Textedit                      

    if rmenu == 'D':
        print 'Choice D: Plot EQ Magnitude vs. Time'
        os.system("python plot_ANSS_seismicity.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))
    
    if rmenu == 'E':
        print 'Choice E: Plot EQ Epicenters on a Map'
        os.system("python plot_ANSS_epicenters.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'F':
        print 'Choice F: Plot NTW EQ Forecast on a Map'
        os.system("python contour_eq_probs.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'G':
        print 'Choice G: Plot Gutenberg-Richter relation (Cumulative frequency-magnitude)'
        os.system("python plot_GR_relation.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'H':
        print 'Choice H: Plot Earthquake EPS in Region (Defined by histogram of small earthquake counts)'
        os.system("python plot_EQ_EPS_Region.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'I':
        print 'Choice I: Plot Earthquake EPS in Circle within Region (Defined by histogram of small earthquake counts)'
        os.system("python plot_EQ_EPS_Region_Circle.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'J':
        print 'Choice J: Plot Earthquake EPS in Polygon within Region (Defined by histogram of small earthquake counts)'
        os.system("python plot_EQ_EPS_Region_Polygon.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'K':
        print 'Choice K: Plot Proxy Strain vs. Time in the Region'
        os.system("python plot_proxy_strain_vs_time_region.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'L':
        print 'Choice L: Plot Proxy Strain vs. Time in the Circle'
        os.system("python plot_proxy_strain_vs_time_circle.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'M':
        print 'Choice M: Plot Forecast in Circle within Region (Defined by histogram of small earthquake counts)'
        os.system("python plot_Forecast_EPS_Region_Circle.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'N':
        print 'Choice N: Plot Filtered Earthquake EPS in Region (Defined by histogram of small earthquake counts)'
        os.system("python plot_Filtered_EQ_EPS_Region.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'A':
        print 'Choice A: Enter the map and plot parameters (Default is California)'
        print ' '
        print ' Enter a different predefined location or large earthquake magnitude? (y/n)'
        respl = raw_input()
        if respl == 'y':
    #......................................

            print ' '
            print '     Current pre-defined Locations are: '
            print ' '
            for j in range(number_locations):
                print '        ', Location_file[j]

    #......................................

            print ' '
            print ' Current location is: ', Location
            print ' Enter new location (Case sensitive: Overrides previous parameter set):'
            Location = raw_input()

            location_flag = 0
            for j in range(number_locations):
                if Location == Location_file[j]:
                    SWLat = MinLat_file[j]
                    NELat = MaxLat_file[j]
                    SWLng = MinLng_file[j]
                    NELng = MaxLng_file[j]
                    location_flag = 1
                    Last_location = Location

            if location_flag == 0:
                print ' '
                print '     Invalid location, try again...'
                Location = Last_location

        print ' '
        print ' Minimum magnitude (for large earthquakes) is currently set at: ', MagLo
        print ' Enter new minimum magnitude? (y/n)'
        respm = raw_input()
        if respm == 'y':
            print ' '
            print ' Current minimum magnitude is: ', MagLo
            print ' Enter new minimum magnitude (must be M>5.0):'
            MagLo = raw_input()

        if respl !='y':
            print ' '
            print ' Enter parameters (Lats, Longs, Location Name?) (y/n)'
            respp = raw_input()
            if respp == 'y':
                print ' '
                print ' Enter Min Lat, Max Lat, Min Long, Max Long, Location Name'
                print '>>>>>>   (Requires 5 parameters separated by commas)'
                items = raw_input().split(',')
                print ' '
                print ' You Entered: ', items
                print ' '
                SWLat       = items[0]
                NELat       = items[1]
                SWLng       = items[2]
                NELng       = items[3]
                Location    = items[4]


        print 'Downloading catalog data for new location and/or with new minimum magnitude'
        os.system("python generate_catalog_file.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))


        

        



    
