#!/opt/local/bin python

    #   Python Menu Code to call the various procedures to plot EQ data and forecasts
    #
    #   This code downloads data from various web sites and uses matplotlib to construct plots

import sys
import os
import plot_ANSS_seismicity
import numpy as np
from array import array

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
#   print items  # For testing purposes
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
    
MagLo = 6.0 # Initial default minimum magnitude

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
    print ' '
    print '     Z: Get Me Out of Here!'
    print ' '
    print '     (You must choose ... but choose wisely!'
    print '     For a true choice will bring you joy,'
    print '     but a poor choice will take it from you...)'
    print ' '
    rmenu = raw_input( ) 
    print ' '

    #.....................................
    
    if rmenu == 'B':
        print 'Choice B: Generate or Re-Generate Working Catalog File'
        os.system("python generate_catalog_file.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'C':
        print 'Choice C: Open Working Catalog File in Text Editor'
        os.system("open -t EQ_Working_Catalog.txt") # Only for Macs.  For other systems, substitute something like "gedit EQ_Working_Catalog.txt"                                       

    if rmenu == 'D':
        print 'Choice D: Plot EQ Magnitude vs. Time'
        os.system("python plot_ANSS_seismicity.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))
    
    if rmenu == 'E':
        print 'Choice E: Plot EQ Epicenters on a Map'
        os.system("python plot_ANSS_epicenters.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))

    if rmenu == 'F':
        print 'Choice F: Plot NTW EQ Forecast on a Map'
        os.system("python contour_eq_probs.py {0} {1} {2} {3} {4} {5}".format(NELat, NELng, SWLat, SWLng, MagLo, Location))



    if rmenu == 'A':
        print 'Choice A: Enter the map and plot parameters (Default is California)'
        print ' '
        print ' Enter a different predefined location? (y/n)'
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
            print ' Enter new minimum magnitude (must be M>5.5):'
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

        

        



    
