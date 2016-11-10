#!/opt/local/bin python
#import sys
#sys.path.reverse()

    #   Earthquake Methods library of methods and functions
    #   
    #   This code base collects the methods and functions used to make
    #   plots and maps of earthquake data and activity
    #
    ######################################################################

import sys
import matplotlib
import matplotlib.mlab as mlab
import numpy as np
from numpy import *
from mpl_toolkits.basemap import Basemap
from array import array
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches

import datetime
import dateutil.parser

import urllib
import os

import math
from math import exp

from EQUtilities import *


    ######################################################################

    # Set the catalog parameters

def get_catalog(NELat, NELng, SWLat, SWLng, MagLo):

#    completeness_mag = 2.99

    settings_params = get_settings()
    region_type = settings_params[0]
    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])

#   Read completeness_mag from "current_settings.txt"

    print ' '
    print ' Current value of catalog completeness magnitude is: ', completeness_mag
    print ' Do you want to change it? (y/n)'
    resp_completeness = raw_input()

#   Write new value of completeness_mag to "Settngs_File.txt"

    if resp_completeness == 'y':
        print ' Enter new value of completeness magnitude'
        completeness_mag = float(raw_input())

#   Read earthquake depth from "current_settings.txt"

    print ' '
    print ' Current value of max earthquake depth is: ', earthquake_depth
    print ' Do you want to change it? (y/n)'
    resp_depth = raw_input()

#   Write new value of max earthquake depth to "current_settings.txt"

    if resp_depth == 'y':
        print ' Enter new value of max earthquake depth (km):'
        earthquake_depth = float(raw_input())

    if region_type == 'Circle':
        settings_params[1] = completeness_mag
        settings_params[2] = earthquake_depth
        save_settings(settings_params)
    data = {

    # Adjust these  -   CA-NV
        "minmag": completeness_mag,
        "minlat": SWLat,
        "maxlat": NELat,
        "minlon": SWLng,
        "maxlon": NELng,
#        "mindepth": 0,         # Leaving these depth params in leads to missing some earthquake events
#        "maxdepth": 1000,


    # Date parameters
        "mintime": "1900/01/01",
#        "maxtime": "1994/01/16",
        "maxtime": "2110/01/01",


    # Leave these
        "format": "cnss",
        "output": "readable",
        "etype": "E",
        "searchlimit": 99999999,
    }

    spacer = ' '

    # Fetch the data
    url = "http://www.ncedc.org/cgi-bin/catalog-search2.pl"

    params = urllib.urlencode(data)
    catalog = urllib.urlopen(url, params).readlines()

    #############################################################

    #   Write output record file

    # Open the output file
    output = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "w")

    # Format the data

    date_string = ' '
    time_string = ' '
    ts = 0.

    i=-1
    for line in catalog:
        i+=1
        items = line.strip().split()

        try:
            date = dateutil.parser.parse(items[0] + ' ' + items[1])     # gives year-month-day (space) hours:minutes:seconds
    #       date = dateutil.parser.parse(items[0])    <-- gives year-month-day only
            ts = date.year + float(date.strftime("%j"))/366

            lat                 = items[2]
            lon                 = items[3]
            dep                 = items[4]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]

    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

            if dep == 'Unk':
                dep = '0.0'

            if mag == 'Unk':
                mag = '0.0'

            if mag == 'ML':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mb':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mw':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mc':
                mag = items[4]
                dep = '0.0'

            if mag == 'Md':
                mag = items[4]
                dep = '0.0'

            if mag == 'Ms':
                mag = items[4]
                dep = '0.0'

    #       print lat,lon,dep,mag,date_string,time_string

            if (float(dep) <= float(earthquake_depth)):
                output.write("%s\t%s\t%f\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

        except:
           pass

    # Finalize the output file
    output.close()

    print ' '
    print '     Read in Catalog Date, now Building Working File'

    #############################################################

    #   Now Write the Output Working File

    #############################################################

    #   First count number of lines in master data file that we just downloaded above   

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")

    #   Now open the working file where we will re-write the data

    j=0
    for line in data_file:
        j+=1
    h=j

    data_file.close()
    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r") # Put the file pointer back at the top


    working_file = open("EQ_Working_Catalog.txt", "w")

    i=0
    j=0
    for line in data_file:
        j+=1
        items = line.strip().split()
    #
        percent_complete = 100.0*float(j)/float(h) 
        print '     Percent File Completed: ', "{0:.2f}".format(percent_complete),'%', "                                 \r",   

    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #

        date_string         = items[0]
        time_string         = items[1]
        ts                  = items[2]
        lon                 = items[3]
        lat                 = items[4]
        mag                 = items[5]
        dep                 = items[6]

        lat     = 0.0001*float(int(10000.0*float(lat + '000000')))

        if float(mag) >= MagLo:
            i+=1
            event = str(i)
            if i <10:
                event = '0'+str(i)

            working_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()

    # Finalize the output file

    working_file.close()
    data_file.close()

    return None

    ######################################################################

    # Set the catalog parameters for circular region

def get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, MagLo):


    settings_params = get_settings()
    region_type = settings_params[0]
    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])


#   Read completeness_mag and max earthquake depth from "current_settings.txt"

    data = {

    # Adjust these
        "minmag": completeness_mag,
#       "mindepth": 0,
#       "maxdepth": 1000,

    # Date parameters
        "mintime": "1900/01/01",
#        "maxtime": "1994/01/16",
        "maxtime": "2110/01/01",


    # Leave these
        "delta": "%s,%s,0,%s" % (Circle_Lat, Circle_Lng, Radius_float),
        "format": "cnss",
        "output": "readable",
        "etype": "E",
        "searchlimit": 999999,
}

    spacer = ' '

    # Fetch the data
    url = "http://www.ncedc.org/cgi-bin/catalog-search2.pl"

    params = urllib.urlencode(data)
    catalog = urllib.urlopen(url, params).readlines()

    #############################################################

    #   Write output record file


    # Open the output file
    output = open("ANSS_%s.circle.catalog" % datetime.date.today().strftime("%F"), "w")

    # Format the data

    date_string = ' '
    time_string = ' '
    ts = 0.

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []

    i=-1
    for line in catalog:
        i+=1
        items = line.strip().split()

        try:
            date = dateutil.parser.parse(items[0] + ' ' + items[1])     # gives year-month-day (space) hours:minutes:seconds
    #       date = dateutil.parser.parse(items[0])    <-- gives year-month-day only
            ts = date.year + float(date.strftime("%j"))/366

            lat                 = items[2]  #  These are the items as formatted in the online catalog file
            lon                 = items[3]
            dep                 = items[4]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]

    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

            if dep == 'Unk':
                dep = '0.0'

            if mag == 'Unk':
                mag = '0.0'

            if mag == 'ML':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mb':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mw':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mc':
                mag = items[4]
                dep = '0.0'

            if mag == 'Md':
                mag = items[4]
                dep = '0.0'

            if mag == 'Ms':
                mag = items[4]
                dep = '0.0'

            mag_array.append(mag)           #   List of magnitudes
            date_array.append(date_string)
            time_array.append(time_string)
            year_array.append(ts)

            if (float(dep) <= float(earthquake_depth)):
                output.write("%s\t%s\t%f\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

        except:
           pass

    # Finalize the output file
    output.close()

    #############################################################

    #   Now Write the Output Working File

    #############################################################

    #   First count number of lines in master data file that we just downloaded above

    data_file = open("ANSS_%s.circle.catalog" % datetime.date.today().strftime("%F"), "r")

    #   Now open the working file where we will re-write the data

    j = 0
    for line in data_file:
        j += 1
    h = j

    data_file.close()
    data_file = open("ANSS_%s.circle.catalog" % datetime.date.today().strftime("%F"),
                     "r")  # Put the file pointer back at the top

    os.system("rm EQ_Working_Circle_Catalog.txt")

    working_file = open("EQ_Working_Circle_Catalog.txt", "w")

    i = 0
    j = 0
    for line in data_file:
        j += 1
        items = line.strip().split()
        #
        percent_complete = 100.0 * float(j) / float(h)
        print '     Percent File Completed: ', "{0:.2f}".format(
            percent_complete), '%', "                                 \r",

        #
        #   Remember that at this point, all the below variables are string variables, not floats
        #

        date_string = items[0]
        time_string = items[1]
        ts = items[2]
        lon = items[3]
        lat = items[4]
        mag = items[5]
        dep = items[6]

        lat = 0.0001 * float(int(10000.0 * float(lat + '000000')))

        if float(mag) >= MagLo:
            i += 1
            event = str(i)
            if i < 10:
                event = '0' + str(i)

            working_file.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()

    # Finalize the output file

    working_file.close()
    data_file.close()

    return mag_array, date_array, time_array, year_array

    ######################################################################

def get_polygon_catalog(polygon_vertex_data,MagLo):


    settings_params = get_settings()
    region_type = settings_params[0]

    if region_type == 'Polygon':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])
        Polygon_Location = settings_params[3]

    number_polygon_vertices = (len(settings_params)-4)/2

    #   Construct the string of polygon vertices

    polygon_vertices = settings_params[4] + ',' + settings_params[5]    #   First vertex

    for i in range(0,number_polygon_vertices-1):
        polygon_vertices += ',' + settings_params[6+2*i] + ',' + settings_params[7+2*i]

    #   Close the polygon

    polygon_vertices += ',' + settings_params[4] + ',' + settings_params[5]

#   Read completeness_mag and max earthquake depth from "current_settings.txt"

    data = {

    # Adjust these
        "minmag": completeness_mag,
#       "mindepth": 0,
#       "maxdepth": 1000,

    # Date parameters
        "mintime": "1900/01/01",
#        "maxtime": "1994/01/16",
        "maxtime": "2110/01/01",


    # Leave these
        "polygon": "%s" % (polygon_vertices),
        "format": "cnss",
        "output": "readable",
        "etype": "E",
        "searchlimit": 999999,
}

    spacer = ' '

    # Fetch the data
    url = "http://www.ncedc.org/cgi-bin/catalog-search2.pl"

    params = urllib.urlencode(data)
    catalog = urllib.urlopen(url, params).readlines()

    #############################################################

    #   Write output record file


    # Open the output file
    output = open("ANSS_%s.polygon.catalog" % datetime.date.today().strftime("%F"), "w")

    # Format the data

    date_string = ' '
    time_string = ' '
    ts = 0.

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []

    i=-1
    for line in catalog:
        i+=1
        items = line.strip().split()

        try:
            date = dateutil.parser.parse(items[0] + ' ' + items[1])     # gives year-month-day (space) hours:minutes:seconds
    #       date = dateutil.parser.parse(items[0])    <-- gives year-month-day only
            ts = date.year + float(date.strftime("%j"))/366

            lat                 = items[2]  #  These are the items as formatted in the online catalog file
            lon                 = items[3]
            dep                 = items[4]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]

    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

            if dep == 'Unk':
                dep = '0.0'

            if mag == 'Unk':
                mag = '0.0'

            if mag == 'ML':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mb':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mw':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mc':
                mag = items[4]
                dep = '0.0'

            if mag == 'Md':
                mag = items[4]
                dep = '0.0'

            if mag == 'Ms':
                mag = items[4]
                dep = '0.0'

            mag_array.append(mag)           #   List of magnitudes
            date_array.append(date_string)
            time_array.append(time_string)
            year_array.append(ts)

            if (float(dep) <= float(earthquake_depth)):
                output.write("%s\t%s\t%f\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

        except:
           pass

    # Finalize the output file
    output.close()

    #############################################################

    #   Now Write the Output Working File

    #############################################################

    #   First count number of lines in master data file that we just downloaded above

    data_file = open("ANSS_%s.polygon.catalog" % datetime.date.today().strftime("%F"), "r")

    #   Now open the working file where we will re-write the data

    j = 0
    for line in data_file:
        j += 1
    h = j

    data_file.close()
    data_file = open("ANSS_%s.polygon.catalog" % datetime.date.today().strftime("%F"),
                     "r")  # Put the file pointer back at the top

    os.system("rm EQ_Working_Polygon_Catalog.txt")

    working_file = open("EQ_Working_Polygon_Catalog.txt", "w")

    i = 0
    j = 0
    for line in data_file:
        j += 1
        items = line.strip().split()
        #
        percent_complete = 100.0 * float(j) / float(h)
        print '     Percent File Completed: ', "{0:.2f}".format(
            percent_complete), '%', "                                 \r",

        #
        #   Remember that at this point, all the below variables are string variables, not floats
        #

        date_string = items[0]
        time_string = items[1]
        ts = items[2]
        lon = items[3]
        lat = items[4]
        mag = items[5]
        dep = items[6]

        lat = 0.0001 * float(int(10000.0 * float(lat + '000000')))

        if float(mag) >= MagLo:
            i += 1
            event = str(i)
            if i < 10:
                event = '0' + str(i)

            working_file.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()

    # Finalize the output file

    working_file.close()
    data_file.close()

    return mag_array, date_array, time_array, year_array

    ######################################################################


def magtime(NELat, NELng, SWLat, SWLng, MagLo, Location):


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(i)]
    yrs                 =   zeros(i)
    lng                 =   zeros(i)
    lat                 =   zeros(i)
    mag                 =   zeros(i)
    dep                 =   zeros(i)
    time_string         =   ["" for x in range(i)]
    date_string         =   ["" for x in range(i)]

    # Loop over lines and extract variables of interest

    i=-1
    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    range_limit = len(lat) - 1

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ' '

    #.................................................................

    plotmin = MagLo - 0.5
    plotmax = np.amax(mag) + 0.5

    plt.plot()
    plt.ylim(plotmin, plotmax)
    delta_years = 0.07*(yrs[range_limit]-yrs[0])
    plt.xlim(int(yrs[0]-delta_years), int(yrs[range_limit] +delta_years))
#    plt.xlim(1930.0, 2020.0)

    SupTitle_text = 'Earthquake Magnitude vs. Time for M>' + str(MagLo) + ' in ' + Location
    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'From: ' + date_string[0] + ' ' + time_string[0] + '     To: ' + date_string[range_limit] + ' ' + time_string[range_limit]
    plt.title(Title_text)

    plt.xlabel('Time (years)', fontsize=14)
    plt.ylabel('Magnitude', fontsize=14)

    plot_lines = 'TRUE'

    if plot_lines == 'TRUE':
        for i in range(0,range_limit+1):
            x1=yrs[i]
            x2=yrs[i]
            y1=0.0
            y2=mag[i]
            plt.plot([x1, x2], [y1, y2], 'r-', lw=1.15)


    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_limit):
            plt.plot(yrs, mag, 'bo', ms=6)


    plt.savefig('Magnitude-Time.pdf')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    return None

    ######################################################################

def freqmag(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'FALSE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    print ' '
    print '       Sequence', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude'
    print ' '

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1
            if eq_sequence_number <10:
                blank_space = '         '
            if eq_sequence_number >=10:
                blank_space = '        '
            if eq_sequence_number >=100:
                blank_space = '       '
            print blank_space, eq_sequence_number, '    ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i]

        last_eq = eq_sequence_number

    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ' '
 
    #.................................................................

    print ' '
    print ' Enter Earthquake Sequence Numbers for Plot (Initial, Final)'
    print ' (To plot data since Last large EQ, enter (Last EQ, Last EQ + 1):'
    print ''
    Sequence_Numbers_Input = raw_input()
    Seq = Sequence_Numbers_Input.split(',')

    #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(Seq[0]) - 1   
    second_eq   = int(Seq[1]) - 1

    if print_data == 'TRUE':
        print ' '
        print int(Seq[0]), int(Seq[1])

    number_of_eq_cycles = int(Seq[1]) - int(Seq[0])


    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array

    number_small_eq_array   =   zeros(number_eqs+10)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = 0             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query     = float(items[5])
        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1
        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]

    if print_data == 'TRUE':
        print 'total_number_small_eqs: ', total_number_small_eqs        
        print 'number_small_eq_array: ', number_small_eq_array[:]

    #.................................................................

    resp_scale_stack = 'n'
    normalize_stack = 'n'
    if number_of_eq_cycles > 1:
        print ''
        print ' Normalize the Stack of Multiple EQ Cycles by Number of Cycles? (y/n)'
        normalize_stack = raw_input()        

    #.................................................................

    #   Open the complete data file again and collect the small earthquakes for plotting

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    sequence =0                                 #   This will be the display large earthquake number, which is 1 more than the internal large earthquake number
                                                #   That is, internally, the first large earthquake is number 0.  For display, it is number 1.

    eq_id=0                                     #   This is the small earthquake ID number

    total_weight_factor = 0.0

    for line in data_file:
        items = line.strip().split()
    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #

        date_small         = items[0]
        time_small         = items[1]
        ts_small           = items[2]

        lon_small          = float(items[3])
        lat_small          = float(items[4])
        mag_small          = float(items[5])
        dep_small          = float(items[6])

        if mag_small >= MagLo:
            freq_mag_bins_pdf += freq_mag_pdf_working
            sequence += 1                                           # Increments the sequence number
            number_small_eqs = 0                                    # Reset number of small earthquakes = 0
            number_nonzero_bins=0
            freq_mag_pdf_working =  zeros(number_mag_bins)

        if first_eq <= sequence < second_eq:
            if mag_small < MagLo and mag_small >= min_mag:
                eq_id += 1                                          # Assigns an ID number to each small earthquake
                bin_number = int((mag_small - min_mag) * 10)
                freq_mag_pdf_working[bin_number] += 1                  # Increment the appropriate bin to compute the GR PDF
               

    data_file.close()                                               # Close the data file

    #.................................................................

    

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]

    if normalize_stack == 'y':
        freq_mag_bins_sdf = [i/float(number_of_eq_cycles) for i in freq_mag_bins_sdf]


    print ''

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff

    for i in range(0,number_mag_bins):
        mag_array[i]=min_mag + float(i)*bin_diff


    #.................................................................

    if print_data == 'TRUE':
        print ''
        print mag_bins
        print ''
        print freq_mag_bins_pdf
        print ''
        print freq_mag_bins_sdf
        print ''
        print log_freq_mag_bins
        print ''

    #.................................................................

#
#     Fit the line to find the b-value
#

    xr    =   zeros(number_mag_bins)
    yr    =   zeros(number_mag_bins)

    print  ' '
    print  ' Max number of magnitude bins is', number_mag_bins
    print  ' '
    print  ' Number nonzero bins is: ', number_nonzero_bins
    print  ' '
    print  ' Enter LOW and HIGH mag. values'
    print  ' to determine b-value (magnitudes)'
    print  ' (to bypass this put HIGH < LOW)'
    print  ' '

    mag1 = 0.0
    mag2 = 0.0

    Line_Magnitude_Limits = raw_input()
    LineMagLims = Line_Magnitude_Limits.split(',')

    #   The large earthquakes were renumbered to start from 1 so need to correct for that

    mag1 = float(LineMagLims[0])   
    mag2 = float(LineMagLims[1]) 

    if print_data == 'TRUE':
        print ' '
        print mag1, mag2
        print mag_bins

    #   Here us where we count the number of data to fit, then define xr[:] and yr[:]

    xr    =   zeros(1000)
    yr    =   zeros(1000)

    k = -1
    for i in range(0,number_nonzero_bins):
        if (mag_bins[i] >= mag1) and (mag_bins[i] <= mag2):
            k = k + 1
            xr[k] = mag_bins[i]
            yr[k] = log_freq_mag_bins[i]
            if print_data == 'TRUE':
                print ' k, xr, yr: ', k, xr[k], yr[k]

    kmax = k+1

    if kmax > 1:
        slope, cept, errs, errc, s2 = linfit(xr,yr, kmax)
        bval = - slope
        print ' '
        print ' b - value is:'
        print   bval,' +/- ',errs
        print ' '
        print ' Intercept is: ', cept
        print ' '
        print ' '


    #.................................................................

    resp_droplet = 'n'
    last_droplet_plot = 'n'

    if kmax > 1:
        print 'Plot the best fitting droplet curve? (y/n)'
        resp_droplet = raw_input()

        if resp_droplet == 'y' and number_of_eq_cycles == 1:
            print ''
            print 'Plot the last (possibly multicycle) droplet curve? (y/n)'
            last_droplet_plot = raw_input()

        if resp_droplet == 'y':
            x_droplet, y_droplet, cept, slope, droplet_bins, xi_droplet, sigma_droplet, sum_of_squares_low =                \
                    dropletFit(mag_bins, log_freq_mag_bins, mag1, mag2, number_nonzero_bins, MagLo)

            sdev_xi, sdev_sigma =                                                                                           \
                    deviationFit(mag_bins, log_freq_mag_bins, mag1, mag2, number_nonzero_bins, droplet_bins, cept, slope,   \
                    sum_of_squares_low, xi_droplet, sigma_droplet)        

    #.................................................................

    #   Store the last previously computed droplet

    droplet_store_flag = 'w'

    if last_droplet_plot == 'y':
        droplet_store_flag = 'r'
        x_droplet_last, y_droplet_last = storeDroplet(x_droplet,y_droplet,droplet_bins, droplet_store_flag)
 

    #.................................................................


    plotmin = int(min(log_freq_mag_bins) -1)
    if plotmin > 0: 
        plotmin = 0.0
    plotmax = int(log_freq_mag_bins[0] + 1.0)

    plt.plot()
    plt.ylim(plotmin, plotmax)
    plt.xlim(min_mag, MagLo)

    if normalize_stack != 'y':
        SupTitle_text = 'Number-Magnitude Relation for Small EQs in ' + Location
        plt.suptitle(SupTitle_text, fontsize=16)

    if normalize_stack == 'y':
        SupTitle_text = 'Normalized Number-Magnitude for Small EQs in ' + Location
        plt.suptitle(SupTitle_text, fontsize=16)

    if second_eq != range_limit:

        text_x      = (MagLo-min_mag)*0.025 + min_mag
        text_y      = (plotmax-plotmin) * 0.15 + plotmin
        text_string = 'Number of EQ Cycles: ' + str(int(Seq[1]) - int(Seq[0]))
        plt.text(text_x,text_y,text_string, fontsize=12)

        text_x      = (MagLo-min_mag)*0.025 + min_mag
        text_y      = (plotmax-plotmin) * 0.10 + plotmin
        text_string = 'Initial EQ: Date=' + date_string[first_eq] + ';  Time=' + time_string[first_eq] + ';   Mag='+ str(mag[first_eq])
        plt.text(text_x,text_y,text_string, fontsize=12)

        text_x      = (MagLo-min_mag)*0.025 + min_mag
        text_y      = (plotmax-plotmin) * 0.05 + plotmin
        text_string = 'Final EQ:  Date=' + date_string[second_eq] + ';  Time=' + time_string[second_eq] + ';   Mag='+ str(mag[second_eq])
        plt.text(text_x,text_y,text_string, fontsize=12)

    if second_eq == range_limit:

        text_x      = (MagLo-min_mag)*0.025 + min_mag
        text_y      = (plotmax-plotmin) * 0.10 + plotmin
        text_string = 'Data After Most Recent EQ'
        plt.text(text_x,text_y,text_string, fontsize=12)

        text_x      = (MagLo-min_mag)*0.025 + min_mag
        text_y      = (plotmax-plotmin) * 0.05 + plotmin
        text_string = 'Date=' + date_string[first_eq] + ';  Time=' + time_string[first_eq] + ';   Mag='+ str(mag[first_eq])
        plt.text(text_x,text_y,text_string, fontsize=12)

    if kmax > 1:

        projected_max_mag   = float(cept/bval)
        proj_max_max        = float((cept+errc)/(bval - errs)) - projected_max_mag
        proj_max_min        = projected_max_mag - float((cept-errc)/(bval + errs)) 
        error_max_eq        = 0.5*(proj_max_max + proj_max_min)

        text_x      = (MagLo-min_mag)*0.55 + min_mag
        text_y      = (plotmax-plotmin) * 0.95 + plotmin
        bval        = '%.2f'%(bval)
        error       = '%.2f'%(errs)
        text_string = 'b-value:  ' + bval + ' +/- ' + error
        plt.text(text_x,text_y,text_string, fontsize=12)    


        text_x      = (MagLo-min_mag)*0.55 + min_mag
        text_y      = (plotmax-plotmin) * 0.90 + plotmin
        proj_max    = '%.2f'%(projected_max_mag)
        error_max_eq = '%.2f'%(error_max_eq)
        text_string = 'Avg. Projected EQ Mag:  M' + proj_max + ' +/- ' + error_max_eq
        plt.text(text_x,text_y,text_string, fontsize=12)    

        text_x      = (MagLo-min_mag)*0.55 + min_mag
        text_y      = (plotmax-plotmin) * 0.85 + plotmin

        observed_mean_mag = mean_val(mag)
        observed_dev_mag, observed_variance_mag = std_var_val(mag)
        obs_mag    = '%.2f'%(observed_mean_mag)
        obs_dev = '%.2f'%(observed_dev_mag)

        text_string = 'Avg. Observed EQ Mag:  M' + obs_mag + ' +/- ' + obs_dev
        plt.text(text_x,text_y,text_string, fontsize=12)    


        if resp_droplet == 'y':
            text_x      = (MagLo-min_mag)*0.55 + min_mag
            text_y      = (plotmax-plotmin) * 0.80 + plotmin
            sigma        = '%.2f'%(sigma_droplet)
            error       = '%.2f'%(sdev_sigma)
            text_string = 'Sigma:  ' + sigma + ' +/- ' + error
            plt.text(text_x,text_y,text_string, fontsize=12)    

            text_x      = (MagLo-min_mag)*0.55 + min_mag
            text_y      = (plotmax-plotmin) * 0.75 + plotmin
            corner_mag  = '%.2f'%(xi_droplet)
            error       = '%.2f'%(sdev_xi)
            text_string = 'Corner Mag:  ' + corner_mag + ' +/- ' + error
            plt.text(text_x,text_y,text_string, fontsize=12)    

    plt.xlabel('Magnitude', fontsize=14)
    plt.ylabel('Log$_{10}$ ( Number )', fontsize=14)

    if plotmin < 0.0:               #   Draw a horizontal zero line
        zero_line_x = zeros(2)
        zero_line_x[0] = min_mag
        zero_line_x[1] = MagLo
        zero_line_y = zeros(2)

        for i in range(0,range_bins):
            plt.plot(zero_line_x, zero_line_y, 'k-', lw=0.75)

    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_bins):
            plt.plot(mag_bins, log_freq_mag_bins, 'bo', ms=6)

    if resp_droplet == 'y':
        plt.plot(x_droplet, y_droplet, 'g-', lw=1.15)   

    if last_droplet_plot == 'y':
        plt.plot(x_droplet_last, y_droplet_last, 'g--', lw=1.15)
    
    if (kmax > 1):
        x1 = xr[kmax-1]                 # Note that xr[kmax-1] is the largest magnitude that was used to fit the line
        y1 = slope*xr[kmax-1] + cept
        x2 = MagLo
        y2 = slope*MagLo + cept
        plt.plot([x1, x2], [y1, y2], 'r--', lw=1.15)

        x1 = mag1
        y1 = slope*mag1 + cept
        x2 = xr[kmax-1]
        y2 = slope*xr[kmax-1] + cept
        plt.plot([x1, x2], [y1, y2], 'r-', lw=1.15)



    plt.savefig('Number-Magnitude.pdf')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    return None

    ######################################################################

def histogram_eps_region(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number

    input_file.close()  # Put the file pointer back at top


    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq


    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    number_small_eq_array   =   np.zeros(number_eqs)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1
        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    print ' '
    print '       Large EQ', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude', '    Number Small EQs (After this large EQ)'
    print ' '

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1  
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]
        if i+1 <10:
            blank_space = '         '
        if i+1 >=10:
            blank_space = '        '
        if i+1 >=100:
            blank_space = '       '
        print blank_space, i+1, '     ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i],  \
                '          ',int(number_small_eq_array[i])

    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ''
    print '     Total Number of Small Earthquakes: ', int(total_number_small_eqs)        
    print ''


   #.................................................................

#   num_bins = 100
    num_bins = 500

    #   Histogram the data

    #   n[i] below is an array with the number of intervals in bin [i]
    #   bins[i] is an array with the location of the bin center
    #   patches[i] is an array with descriptions of each histogram rectangle
    #   set normed = 1 if you want the histogram to be a pdf

    #    n, bins, patches = plt.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    todays_count = number_small_eq_array[last_eq-1]

    number_small_eqs_excluding_last = np.zeros(number_eqs-1)

    for i in range(0,number_eqs-1):
        number_small_eqs_excluding_last[i] = number_small_eq_array[i]

    mean    =   mean_val(number_small_eqs_excluding_last)
    std_deviation, variance = std_var_val(number_small_eqs_excluding_last)

    mean_eqs = mean

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print ''

    #.................................................................

    #   Bin the data

    todays_probability = 0.0

    n, bins = histogram(number_small_eqs_excluding_last, num_bins)

    #   -----------------------------------------


    cum_poisson = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_small_EQs))

    cum_poisson[:] = (cum_poisson[:]) * 100.0

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:]/ float(last_eq)

    scale_factor = 100.0/cum_prob[len(bins)-1]

    cum_prob[:] = cum_prob[:]*scale_factor	#   Ensures that cum prob goes from 0%->100%

#    cum_prob[:] = cum_prob[:] *100.0
    cum_prob[:] = cum_prob[:]

    for i in range(1, len(bins)):
        if todays_count >= bins[i]:
            todays_probability = cum_prob[i]

#    mean_eqs = mean_eqs*1.45   #   Test statement



    tau_best, beta_best, sdev_tau, sdev_beta, sum_of_squares_low = weibullBetaFit(n, bins, cum_prob, mean_eqs)

    tau_string   = "{0:.2f}".format(tau_best) + ' +/- '  + "{0:.2f}".format(sdev_tau)
    beta_string  = "{0:.2f}".format(beta_best) + ' +/- ' + "{0:.2f}".format(sdev_beta)

    print ''
    print ' Weibull fit data: '
    print ''
    print ' Tau: ', tau_string
    print 'Beta: ', beta_string
    print ''

    cum_weibull =   np.zeros(len(bins)+1)

    for i in range(0,len(bins)):
        arg_weibull = ((bins[i]/tau_best)**beta_best)
        cum_weibull[i] = 1.0 - math.exp(- arg_weibull)

    cum_weibull[:] = cum_weibull[:] *100.0


    #.................................................................

    #   Define the Temperature

    temperature = "{0:.1f}".format(todays_probability)+'%'

    #.................................................................

    #  Plot the data

    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[14, 1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])

    #   Draw the bins

    ax0.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    #   Get the axis limits

    ymin, ymax = ax0.get_ylim()
    xmin, xmax = ax0.get_xlim()

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
    SupTitle_text = 'Potential for M>' + str(MagLo) + ' Earthquakes in ' + Location + ' on ' + todays_date
    SupTitle_text = 'Potential for M>' + str(MagLo) + ' Earthquakes in ' + Location
    Title_text    = 'After M' + '%.2f'%mag[last_eq-1] + ' on ' + date_string[last_eq-1] + ' at ' + time_string[last_eq-1]

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.90 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.85 + ymin

    count_string    = str(int(todays_count))
    text_string1 = 'Todays Small EQ Count:  ' + count_string 
    text_string2 = 'On:  ' + date_small_quake + '   at:  ' + time_small_quake
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.77 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string = 'Todays EPS Value:  ' + prob +'%'
    plt.text(text_x,text_y,text_string, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.70 + ymin
    text_string = 'Weibull Parameters: '
    plt.text(text_x,text_y,text_string, fontsize=12)    

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.65 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.60 + ymin

    count_string    = str(int(todays_count))
    text_string1 = '    Tau:  ' + tau_string 
    text_string2 = '   Beta:  ' + beta_string
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12) 

    bins        = np.append(bins,bins[num_bins])
    cum_prob    = np.append(cum_prob, 100.0)
    cum_poisson = np.append(cum_poisson, 100.0)

    ax1 = ax0.twinx()

    plt.ylim(ymax = 100, ymin = 0)
    ymin, ymax = ax1.get_ylim()

    ax1.plot(bins,cum_weibull, 'g-', lw=1.2)    #   Uncomment if you want to plot the Weibull curve with same mean
    ax1.plot(bins,cum_poisson, 'b--')    #   Uncomment if you want to plot the Poisson curve with same mean
    ax1.plot(bins,cum_prob, 'r-')

    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)
    x=[mean_small_EQs,mean_small_EQs]
    y=[0,100]
    ax1.plot(x,y, 'k', ls='dotted')


    #.................................................................

    #   Thermometer

    ax3 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([0,1])                         #   Show the y-axis labels
    plt.xlim([0,1])                         #   Set the x-axis limits  

    frame1.axes.get_xaxis().set_ticks([])   #   Hide the x-axis ticks and labels
    frame1.axes.get_yaxis().set_ticks([])

    ax4 = ax3.twinx()
    plt.ylim([0,100])                         #   Show the y-axis labels
    plt.xlim([0,1])

    ax4.get_yaxis().set_ticks([])
#    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted', lw=1)

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays EPS Value (Cumulative Probability) is: ', temperature
    print ''

    matplotlib.pyplot.savefig('EPS-Region.pdf')

    plt.show()

    return None

    ######################################################################

def histogram_eps_region_circle(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number


    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq

    input_file.close()  # Put the file pointer back at top

    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    number_small_eq_array   =   np.zeros(number_eqs)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1
        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    print ' '
    print '       Large EQ', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude', '    Number Small EQs'
    print ' '

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1  
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]
        if i+1 <10:
            blank_space = '         '
        if i+1 >=10:
            blank_space = '        '
        if i+1 >=100:
            blank_space = '       '
        print blank_space, i+1, '     ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i],  \
                '          ',int(number_small_eq_array[i])

    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ''
    print '     Total Number of Small Earthquakes: ', int(total_number_small_eqs)        
    print ''

   #.................................................................
   #
   #    Get the data on the circular region

    Circle_Location, Circle_Lat, Circle_Lng, Radius_float = get_circle_data(Location)

    mag_array, date_array, time_array, year_array = get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, MagLo)

   #    Reverse the magnitude array for the magnitudes of events in the circular region

    mag_array_reversed  =   list(reversed(mag_array))
    date_array_reversed =   list(reversed(date_array))
    time_array_reversed =   list(reversed(time_array))
    year_array_reversed =   list(reversed(year_array))

   #    Now find the number of small events within the circle since last large event

    todays_count = 0    

    mag_flag = 0
    for i in range(0,len(mag_array_reversed)):
        if float(mag_array_reversed[i]) >= MagLo and mag_flag == 0:
            last_eq_mag         = str(mag_array_reversed[i])
            last_eq_date        = date_array_reversed[i]
            last_eq_time        = time_array_reversed[i]  
            mag_flag = 1
        if float(mag_array_reversed[i]) < MagLo and mag_flag == 0:
            todays_count += 1

   #.................................................................

#   num_bins = 100
    num_bins = 200

    #   Histogram the data

    #   n[i] below is an array with the number of intervals in bin [i]
    #   bins[i] is an array with the location of the bin center
    #   patches[i] is an array with descriptions of each histogram rectangle
    #   set normed = 1 if you want the histogram to be a pdf

    #    n, bins, patches = plt.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    number_small_eqs_excluding_last = np.zeros(number_eqs-1)

    for i in range(0,number_eqs-1):
        number_small_eqs_excluding_last[i] = number_small_eq_array[i]

    mean    =   mean_val(number_small_eqs_excluding_last)
    std_deviation, variance = std_var_val(number_small_eqs_excluding_last)

    mean_eqs = mean

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print ''

    #.................................................................

    #   Bin the data

    todays_probability = 0.0

    n, bins = histogram(number_small_eqs_excluding_last, num_bins)

    cum_poisson = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_eqs))

    cum_poisson[:] = (cum_poisson[:]) * 100.0

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:]/ float(last_eq)
    cum_prob[:] = cum_prob[:] *100.0

    tau_best, beta_best, sdev_tau, sdev_beta, sum_of_squares_low = weibullBetaFit(n, bins, cum_prob, mean_eqs)

    tau_string  = str(tau_best) + ' +/- ' + str(sdev_tau)
    beta_string = str(beta_best) + ' +/- ' + str(sdev_beta)

    print ''
    print ' Weibull fit data: '
    print ''
    print ' Tau: ', tau_string
    print 'Beta: ', beta_string
    print ''

    cum_weibull =   np.zeros(len(bins)+1)

    for i in range(1,len(bins)):
        cum_weibull[i] = 1.0 - math.exp(-((bins[i]/tau_best)**beta_best))

    cum_weibull[:] = cum_weibull[:] *100.0
    

    #.................................................................

    for i in range(1, len(bins)):           #   Define today's temperature = Earthquake Potential Score (EPS)
        if todays_count >= bins[i]:
            todays_probability = cum_prob[i]

    #   Define the Temperature

    temperature = "{0:.1f}".format(todays_probability)+'%'

    #.................................................................

    #  Plot the data

    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[14, 1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])

    #   Draw the bins

    ax0.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    #   Get the axis limits

    ymin, ymax = ax0.get_ylim()
    xmin, xmax = ax0.get_xlim()

    Circle_Location_actual = Circle_Location.split('-')

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
#    SupTitle_text = 'EPS for M>' + str(MagLo) + ' Earthquakes within ' + str(Radius_float) + ' km of ' + Circle_Location_actual[0] + ' on ' + todays_date
    SupTitle_text = 'EPS for M>' + str(MagLo) + ' Earthquakes within ' + str(Radius_float) + ' km of ' + Circle_Location_actual[0]
    Title_text    = 'After M' + last_eq_mag + ' on ' + last_eq_date + ' at ' + last_eq_time

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.90 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.85 + ymin

    count_string    = str(int(todays_count))
    text_string1 = 'Todays Small EQ Count:  ' + count_string 
    text_string2 = 'On:  ' + date_small_quake + '   at:  ' + time_small_quake
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.77 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string = 'Todays EPS Value:  ' + prob +'%'
    plt.text(text_x,text_y,text_string, fontsize=12)    

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12) 

    bins        = np.append(bins,bins[num_bins])
    cum_prob    = np.append(cum_prob, 100.0)
    cum_poisson = np.append(cum_poisson, 100.0)

    ax1 = ax0.twinx()
    ax1.plot(bins,cum_prob, 'r-')
    plt.ylim(ymax = 100, ymin = 0)
    ymin, ymax = ax1.get_ylim()

#   ax1.plot(bins,cum_weibull, 'g-', lw=1.2)    #   Uncomment if you want to plot the Weibull curve with same mean
    ax1.plot(bins,cum_poisson, 'b--')    #   Uncomment if you want to plot the Poisson curve with same mean

    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)
    x=[mean_small_EQs,mean_small_EQs]
    y=[0,100]
    ax1.plot(x,y, 'k', ls='dotted')


    #.................................................................

    #   Thermometer

    ax3 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([0,1])                         #   Show the y-axis labels
    plt.xlim([0,1])                         #   Set the x-axis limits  

    frame1.axes.get_xaxis().set_ticks([])   #   Hide the x-axis ticks and labels
    frame1.axes.get_yaxis().set_ticks([])

    ax4 = ax3.twinx()
    plt.ylim([0,100])                         #   Show the y-axis labels
    plt.xlim([0,1])

    ax4.get_yaxis().set_ticks([])
#    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted', lw=1)

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays Earthquake Potential Score (Cumulative Probability) is: ', todays_probability
    print ''

    matplotlib.pyplot.savefig('Earthquake-Potential-Score.pdf')
    matplotlib.pyplot.savefig('Earthquake-Potential-Score.jpg')

    plt.show()

    return None


    ######################################################################

def histogram_eps_region_polygon(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number


    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq

    input_file.close()  # Put the file pointer back at top

    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    number_small_eq_array   =   np.zeros(number_eqs)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1
        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    print ' '
    print '       Large EQ', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude', '    Number Small EQs'
    print ' '

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1  
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]
        if i+1 <10:
            blank_space = '         '
        if i+1 >=10:
            blank_space = '        '
        if i+1 >=100:
            blank_space = '       '
        print blank_space, i+1, '     ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i],  \
                '          ',int(number_small_eq_array[i])

    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ''
    print '     Total Number of Small Earthquakes: ', int(total_number_small_eqs)        
    print ''

   #.................................................................
   #
   #    Get the data on the polygon region

    polygon_vertex_data = get_polygon_data(Location)

    mag_array, date_array, time_array, year_array = get_polygon_catalog(polygon_vertex_data,MagLo)

   #    Reverse the magnitude array for the magnitudes of events in the circular region

    mag_array_reversed  =   list(reversed(mag_array))
    date_array_reversed =   list(reversed(date_array))
    time_array_reversed =   list(reversed(time_array))
    year_array_reversed =   list(reversed(year_array))

   #    Now find the number of small events within the polygon since last large event

    todays_count = 0    

    mag_flag = 0
    for i in range(0,len(mag_array_reversed)):
        if float(mag_array_reversed[i]) >= MagLo and mag_flag == 0:
            last_eq_mag         = str(mag_array_reversed[i])
            last_eq_date        = date_array_reversed[i]
            last_eq_time        = time_array_reversed[i]  
            mag_flag = 1
        if float(mag_array_reversed[i]) < MagLo and mag_flag == 0:
            todays_count += 1

   #.................................................................

#   num_bins = 100
    num_bins = 200

    #   Histogram the data

    #   n[i] below is an array with the number of intervals in bin [i]
    #   bins[i] is an array with the location of the bin center
    #   patches[i] is an array with descriptions of each histogram rectangle
    #   set normed = 1 if you want the histogram to be a pdf

    #    n, bins, patches = plt.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    number_small_eqs_excluding_last = np.zeros(number_eqs-1)

    for i in range(0,number_eqs-1):
        number_small_eqs_excluding_last[i] = number_small_eq_array[i]

    mean    =   mean_val(number_small_eqs_excluding_last)
    std_deviation, variance = std_var_val(number_small_eqs_excluding_last)

    mean_eqs = mean

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print ''

    #.................................................................

    #   Bin the data

    todays_probability = 0.0

    n, bins = histogram(number_small_eqs_excluding_last, num_bins)

    cum_poisson = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_eqs))

    cum_poisson[:] = (cum_poisson[:]) * 100.0

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:]/ float(last_eq)
    cum_prob[:] = cum_prob[:] *100.0

    tau_best, beta_best, sdev_tau, sdev_beta, sum_of_squares_low = weibullBetaFit(n, bins, cum_prob, mean_eqs)

    tau_string  = str(tau_best) + ' +/- ' + str(sdev_tau)
    beta_string = str(beta_best) + ' +/- ' + str(sdev_beta)

    print ''
    print ' Weibull fit data: '
    print ''
    print ' Tau: ', tau_string
    print 'Beta: ', beta_string
    print ''

    cum_weibull =   np.zeros(len(bins)+1)

    for i in range(1,len(bins)):
        cum_weibull[i] = 1.0 - math.exp(-((bins[i]/tau_best)**beta_best))

    cum_weibull[:] = cum_weibull[:] *100.0
    

    #.................................................................

    for i in range(1, len(bins)):           #   Define today's temperature = Earthquake Potential Score (EPS)
        if todays_count >= bins[i]:
            todays_probability = cum_prob[i]

    #   Define the Temperature

    temperature = "{0:.1f}".format(todays_probability)+'%'

    #.................................................................

    #  Plot the data

    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[14, 1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])

    #   Draw the bins

    ax0.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    #   Get the axis limits

    ymin, ymax = ax0.get_ylim()
    xmin, xmax = ax0.get_xlim()

    Polygon_Location = polygon_vertex_data[0]
    Polygon_Location_actual = Polygon_Location.split('-')

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
#    SupTitle_text = 'EPS for M>' + str(MagLo) + ' Earthquakes within ' + Polygon_Location_actual[0] + ' on ' + todays_date
    SupTitle_text = 'EPS for M>' + str(MagLo) + ' Earthquakes within ' + Polygon_Location_actual[0]
    Title_text    = 'After M' + last_eq_mag + ' on ' + last_eq_date + ' at ' + last_eq_time

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.90 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.85 + ymin

    count_string    = str(int(todays_count))
    text_string1 = 'Todays Small EQ Count:  ' + count_string 
    text_string2 = 'On:  ' + date_small_quake + '   at:  ' + time_small_quake
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.77 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string = 'Todays EPS Value:  ' + prob +'%'
    plt.text(text_x,text_y,text_string, fontsize=12)    

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12) 

    bins        = np.append(bins,bins[num_bins])
    cum_prob    = np.append(cum_prob, 100.0)
    cum_poisson = np.append(cum_poisson, 100.0)

    ax1 = ax0.twinx()
    ax1.plot(bins,cum_prob, 'r-')
    plt.ylim(ymax = 100, ymin = 0)
    ymin, ymax = ax1.get_ylim()

#   ax1.plot(bins,cum_weibull, 'g-', lw=1.2)    #   Uncomment if you want to plot the Weibull curve with same mean
    ax1.plot(bins,cum_poisson, 'b--')    #   Uncomment if you want to plot the Poisson curve with same mean

    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)
    x=[mean_small_EQs,mean_small_EQs]
    y=[0,100]
    ax1.plot(x,y, 'k', ls='dotted')


    #.................................................................

    #   Thermometer

    ax3 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([0,1])                         #   Show the y-axis labels
    plt.xlim([0,1])                         #   Set the x-axis limits  

    frame1.axes.get_xaxis().set_ticks([])   #   Hide the x-axis ticks and labels
    frame1.axes.get_yaxis().set_ticks([])

    ax4 = ax3.twinx()
    plt.ylim([0,100])                         #   Show the y-axis labels
    plt.xlim([0,1])

    ax4.get_yaxis().set_ticks([])
#    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted', lw=1)

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays Earthquake Potential Score (Cumulative Probability) is: ', todays_probability
    print ''

    matplotlib.pyplot.savefig('Earthquake-Potential-Score.pdf')
    matplotlib.pyplot.savefig('Earthquake-Potential-Score.jpg')

    plt.show()

    return None


    ######################################################################

def forecast_eps_region_circle(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number


    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq

    input_file.close()  # Put the file pointer back at top

    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    number_small_eq_array   =   np.zeros(number_eqs)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1
        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    print ' '
    print '       Large EQ', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude', '    Number Small EQs'
    print ' '

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1  
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]
        if i+1 <10:
            blank_space = '         '
        if i+1 >=10:
            blank_space = '        '
        if i+1 >=100:
            blank_space = '       '
        print blank_space, i+1, '     ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i],  \
                '          ',int(number_small_eq_array[i])

    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ''
    print '     Total Number of Small Earthquakes: ', int(total_number_small_eqs)        
    print ''

    #.................................................................
    #
    #    Get the data on the circular region

    Circle_Location, Circle_Lat, Circle_Lng, Radius_float = get_circle_data(Location)

    mag_array, date_array, time_array, year_array = get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, MagLo)

    #    Reverse the magnitude array for the magnitudes of events in the circular region

    mag_array_reversed  =   list(reversed(mag_array))
    date_array_reversed =   list(reversed(date_array))
    time_array_reversed =   list(reversed(time_array))
    year_array_reversed =   list(reversed(year_array))

    #.................................................................

    #    Now count the number of small events within the circle since last large event
    #       and compute the number eq cycles, poisson rate, and other items

    todays_count = 0
    mag_flag = 0

    number_large_eqs = 0

    for i in range(0,len(mag_array_reversed)):
        if float(mag_array_reversed[i]) >= MagLo:
            number_large_eqs += 1

        if float(mag_array_reversed[i]) < MagLo and mag_flag == 0:
            todays_count += 1

        if float(mag_array_reversed[i]) >= MagLo and mag_flag == 0:
            interim_count_small = 0
            last_eq_mag         = str(mag_array_reversed[i])
            last_eq_date        = date_array_reversed[i]
            last_eq_time        = time_array_reversed[i]  
            last_eq_year        = year_array_reversed[i]
            mag_flag = 1

        if float(mag_array_reversed[i]) <= MagLo and mag_flag == 1:
            interim_count_small += 1

        if float(mag_array_reversed[i]) >= MagLo and mag_flag == 1:
            first_eq_year       = year_array_reversed[i]        #   Keep resetting this until the first large eq is found
            total_count_small = interim_count_small             #   Ditto

    number_large_eq_cycles = number_large_eqs - 1

    print last_eq_year, first_eq_year

    print last_eq_year, first_eq_year, number_large_eq_cycles

    if number_large_eq_cycles >= 1:
        avg_time_per_cycle      = (float(last_eq_year)-float(first_eq_year))/float(number_large_eq_cycles)
        avg_small_eqs_per_cycle = float(total_count_small) / float(number_large_eq_cycles)
        poisson_rate_small_eqs  =   avg_small_eqs_per_cycle / avg_time_per_cycle    #   Rate for small eqs in the circle

    if number_large_eq_cycles < 1:
        print 'No large earthquake cycles!  Return'
        return None

   #.................................................................

#   num_bins = 100
    num_bins = 200

    #   Histogram the data

    #   n[i] below is an array with the number of intervals in bin [i]
    #   bins[i] is an array with the location of the bin center
    #   patches[i] is an array with descriptions of each histogram rectangle
    #   set normed = 1 if you want the histogram to be a pdf

    #    n, bins, patches = plt.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    number_small_eqs_excluding_last = np.zeros(number_eqs-1)

    for i in range(0,number_eqs-1):
        number_small_eqs_excluding_last[i] = number_small_eq_array[i]

    mean    =   mean_val(number_small_eqs_excluding_last)
    std_deviation, variance = std_var_val(number_small_eqs_excluding_last)

    mean_eqs = mean

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print ''

    #.................................................................

    #   Bin the data

    todays_probability = 0.0

    n, bins = histogram(number_small_eqs_excluding_last, num_bins)  #   n is number of eq cycles in a bin

    cum_poisson         = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob            = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_eqs))

    cum_poisson[:] = (cum_poisson[:]) * 100.0   #   In percent

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:]/ float(last_eq)
    cum_prob[:] = cum_prob[:] *100.0            #   In percent

    tau_best, beta_best, sdev_tau, sdev_beta, sum_of_squares_low = weibullBetaFit(n, bins, cum_prob, mean_eqs)

    cond_weib_prob      = np.zeros(len(bins))   #   Conditional Weibull probability

    forecast_time = 3.0     #   In years

    arg_weib            = (poisson_rate_small_eqs*forecast_time/tau_best)**beta_best

    cond_weib_prob[0]   = 1.0 - math.exp(-arg_weib)

    for i in range(1,len(bins)):
        arg_weib_1  =   ((bins[i] + poisson_rate_small_eqs*forecast_time)/tau_best)**beta_best
        arg_weib_2  =   (bins[i]/tau_best)**beta_best
        arg_weib    = arg_weib_2 - arg_weib_1
        cond_weib_prob[i]  =  1.0 - math.exp(arg_weib)

    cond_weib_prob[:]   =   cond_weib_prob[:] * 100.0   #   In percent

    print cond_weib_prob[:]

    tau_string  = str(tau_best) + ' +/- ' + str(sdev_tau)
    beta_string = str(beta_best) + ' +/- ' + str(sdev_beta)

    print ''
    print ' Weibull fit data: '
    print ''
    print ' Tau: ', tau_string
    print 'Beta: ', beta_string
    print ''

    cum_weibull =   np.zeros(len(bins)+1)

    for i in range(1,len(bins)):
        cum_weibull[i] = 1.0 - math.exp(-((bins[i]/tau_best)**beta_best))

    cum_weibull[:] = cum_weibull[:] *100.0
    

    #.................................................................

    for i in range(1, len(bins)):           #   Define today's temperature = Earthquake Potential Score (EPS)
        if todays_count >= bins[i]:
            print i, cond_weib_prob[i]
            todays_probability = cond_weib_prob[i]

    print 'todays_probability: ', todays_probability

    #   Define the Temperature

    temperature = "{0:.1f}".format(todays_probability)+'%'

    #.................................................................

    #  Plot the data

    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[14, 1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])

    #   Draw the bins

    ax0.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    #   Get the axis limits

    ymin, ymax = ax0.get_ylim()
    xmin, xmax = ax0.get_xlim()

    Circle_Location_actual = Circle_Location.split('-')

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
#    SupTitle_text = 'EPS for M>' + str(MagLo) + ' Earthquakes within ' + str(Radius_float) + ' km of ' + Circle_Location_actual[0] + ' on ' + todays_date
    SupTitle_text = 'Forecast for M>' + str(MagLo) + ' EQs within ' + str(Radius_float) + ' km of ' + Circle_Location_actual[0] + ' and less than ' + str(forecast_time) + ' years'
    Title_text    = 'After M' + last_eq_mag + ' on ' + last_eq_date + ' at ' + last_eq_time

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.90 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.85 + ymin

    count_string    = str(int(todays_count))
    text_string1 = 'Todays Small EQ Count' + count_string 
    text_string2 = 'On:  ' + date_small_quake + '   at:  ' + time_small_quake
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.77 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string1 = 'Todays EPS Value:  ' + prob +'%'
    plt.text(text_x,text_y,text_string1, fontsize=12)   

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12) 

    bins                = np.append(bins,bins[num_bins])
    last_value          = len(cond_weib_prob)
    last_prob           = cond_weib_prob[last_value-1]
    cond_weib_prob      = np.append(cond_weib_prob, last_prob)
    cum_poisson         = np.append(cum_poisson, 100.0)

    ax1 = ax0.twinx()
    ax1.plot(bins,cond_weib_prob, 'r-')
    plt.ylim(ymax = 100, ymin = 0)
    ymin, ymax = ax1.get_ylim()

#   ax1.plot(bins,cum_weibull, 'g-', lw=1.2)    #   Uncomment if you want to plot the Weibull curve with same mean
#   ax1.plot(bins,cum_poisson, 'b--')    #   Uncomment if you want to plot the Poisson curve with same mean

    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)
    x=[mean_small_EQs,mean_small_EQs]
    y=[0,100]
    ax1.plot(x,y, 'k', ls='dotted')


    #.................................................................

    #   Thermometer

    ax3 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([0,1])                         #   Show the y-axis labels
    plt.xlim([0,1])                         #   Set the x-axis limits  

    frame1.axes.get_xaxis().set_ticks([])   #   Hide the x-axis ticks and labels
    frame1.axes.get_yaxis().set_ticks([])

    ax4 = ax3.twinx()
    plt.ylim([0,100])                         #   Show the y-axis labels
    plt.xlim([0,1])

    ax4.get_yaxis().set_ticks([])
#    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted')
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted')

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays Earthquake Potential Score (Cumulative Probability) is: ', todays_probability
    print ''

    matplotlib.pyplot.savefig('Earthquake-Forecast.pdf')
    matplotlib.pyplot.savefig('Earthquake-Forecast.pdf')

    plt.show()

    return None


    ######################################################################

def map_epicenters(NELat, NELng, SWLat, SWLng, MagLo, Location):

    resp_map = ''
    print ' '
    print '     Default map is eTopo.  Change to Other map? (y/n)'
    resp_map    =   raw_input()
    if resp_map == 'y':
        print ' '
        print '     Shaded Relief (s) or Flat (f = No Relief)...?'
        resp_map_type   = raw_input()

    lat_0_center = 0.5*(NELat + SWLat)
    lon_0_center = 0.5*(NELng + SWLng)


    llcrnrlat = SWLat
    llcrnrlon  = SWLng
    urcrnrlat  = NELat
    urcrnrlon  = NELng

    m = Basemap(projection='cyl',llcrnrlat=SWLat, urcrnrlat=NELat,
            llcrnrlon=SWLng, urcrnrlon=NELng, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,4),labels=[0,0,0,1])
    m.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0])

    #   m.shadedrelief()
    #   m.bluemarble()
    #   m.etopo()
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    #m.drawrivers()

    m.etopo()
    if resp_map == 'y':
        if resp_map_type == 's':
            m.shadedrelief()
        if resp_map_type == 'f':
            m.drawlsmask(land_color='#ffe8a0', ocean_color='#e6e6ff')  # Use this if map.etopo() is  not used
            m.drawmapboundary(fill_color='#e6e6ff')

    m.etopo()
    if resp_map == 'y':
        m.shadedrelief()

    settings_params = get_settings()

    region_type = settings_params[0]

    #   -------------------------------------------

    if region_type == 'Circle':
        Polygon_Location = 'None'

        Circle_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Circle_Location != 'None':
            earthquake_depth = float(settings_params[2])
            CircleCenterLat = float(settings_params[4])
            CircleCenterLng = float(settings_params[5])
            CircleRadius    = float(settings_params[6])

        print ' Set Circle Location to None? (y/n)'
        Location_response = raw_input()
        
        if Location_response == 'y':
            Circle_Location = 'None'
            region_type = 'Circle'
            settings_params[3]  =   Circle_Location
            settings_params[4]  =   0.0
            settings_params[5]  =   0.0
            settings_params[6]  =   0.0
            save_settings(settings_params)
  
         #   -------------------------------------------

    if region_type == 'Polygon':
        Circle_Location = 'None'

        number_polygon_vertices = (len(settings_params)-4)/2

        Polygon_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Polygon_Location != 'None':
            earthquake_depth = float(settings_params[2])

        x_poly = []
        y_poly = []

        for j in range(0,number_polygon_vertices):
            y_poly.append(settings_params[4+2*j])
            x_poly.append(settings_params[5+2*j])

    #   Close the polygon

        y_poly.append(settings_params[4])
        x_poly.append(settings_params[5])

        for i in range(0,len(x_poly)):
            x_poly[i] = float(x_poly[i])
            y_poly[i] = float(y_poly[i])

        print ' Set Polygon Location to None? (y/n)'
        Location_response = raw_input()
        
        if Location_response == 'y':
            settings_params[0] = 'None'
            settings_params[3] = 'None'
            for i in range(4,len(settings_params)):
                settings_params[i] = 0.0
  
    #   -------------------------------------------


    #   Open input file and find the number of lines in the file

    input_file = open("EQ_Working_Catalog.txt", "r")

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs=i

    #   Open input file again

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Loop over lines and extract variables of interest

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(i)]
    yrs                 =   zeros(i)
    lng                 =   zeros(i)
    lat                 =   zeros(i)
    mag                 =   zeros(i)
    dep                 =   zeros(i)
    time_string         =   ["" for x in range(i)]
    date_string         =   ["" for x in range(i)]

    # Loop over lines and extract variables of interest

    i=-1
    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    input_file.close()  # Put the file pointer back at top

    z_limit = len(lat)
    range_limit = z_limit - 1

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print_text_4 = '     at depths < ' + str(settings_params[2]) + ' km'
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ' '

    #.................................................................

    for z in range(z_limit):
        x,y = m(lng,lat)

    #   -----------------------
        if mag[z] >= 5.0:
            mark_size = 3
        if mag[z] >= 5.5:
            mark_size = 4
        if mag[z] >= 6:
            mark_size = 6
        if mag[z] >= 6.5:
            mark_size = 8
        if mag[z] >= 7:
            mark_size = 10
        if mag[z] >= 7.5:
            mark_size = 12
        if mag[z] >= 8.0:
            mark_size = 14
        if mag[z] >= 8.5:
            mark_size = 16
    #   ----------------------
            
    #   m.plot(x[z], y[z], "ro", mfc='none', mec='r', markeredgewidth=2.0, ms=mark_size[z])
        m.plot(x[z], y[z], "ro", ms=mark_size)

    if (Circle_Location != 'None'):
        x_circle_dg, y_circle_dg = draw_big_circle(CircleCenterLat, CircleCenterLng, CircleRadius)
        
        m.plot(x_circle_dg, y_circle_dg, "b-", lw=1.3)

    if (Polygon_Location != 'None'):
        m.plot(x_poly, y_poly, "b-", lw=1.6)


    #   ------------------------------------------------------
    #   Plot 150 km circle at latest event

    lat_latest = lat[i]
    lng_latest = lng[i]
    x,y = m(lng_latest,lat_latest)

    #   Draw a circle of radius 200 km

    Radius = 200.0  # km
    mark_size = int(155.0 * Radius/200.0)   #  approx at mid latitudes

    #   mark_size = 155
    #   m.plot(x,y, "o", mfc='none', ms=mark_size)

    #   ------------------------------------------------------

    sfv     =   (NELat-SWLat)/12.0
    sfh     =   (NELng-SWLng)/16.0

    sfh     =   sfh * 0.9
    sfht    =   sfh * 1.1  # move the text over a bit more

    lat=[llcrnrlat - 0.30*sfv, llcrnrlat + 0.25*sfv,llcrnrlat + 0.80*sfv, llcrnrlat + 1.30*sfv, llcrnrlat + 1.80*sfv, llcrnrlat + 2.40*sfv, llcrnrlat + 3.1*sfv, llcrnrlat + 3.95*sfv]
    lng=[urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht]
    x,y = m(lng,lat)

    latc=[llcrnrlat - 0.25*sfv, llcrnrlat + 0.30*sfv,llcrnrlat + 0.85*sfv, llcrnrlat + 1.375*sfv, llcrnrlat + 1.90*sfv, llcrnrlat + 2.50*sfv, llcrnrlat + 3.2*sfv, llcrnrlat + 4.0*sfv]
    lngc=[urcrnrlon + 0.72*sfh, urcrnrlon + 0.72*sfh, urcrnrlon + 0.73*sfh, urcrnrlon + 0.74*sfh, urcrnrlon + 0.76*sfh, urcrnrlon + 0.78*sfh, urcrnrlon + 0.80*sfh, urcrnrlon + 0.80*sfh]
    xc,yc = m(lngc,latc)

    mag_value=['5.0+','5.5+','6+','6.5+','7+','7.5+', '8+', '8.5+']

    mark_size = [3, 4, 6, 8, 10, 12, 14, 16]

    for z in range(8):
        plt.text(x[z],y[z],mag_value[z], fontsize=10,)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False)

    #   matplotlib.pyplot.savefig('map_bw.pdf')

    if region_type == 'Circle':

        if Circle_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' in ' + Location

        if Circle_Location != 'None':
            Circle_Location_actual = Circle_Location.split('-')
            SupTitle_text = 'Earthquakes M>' + str(MagLo) + ' in ' + Location + ' (Blue Circle = ' + str(CircleRadius) + ' km Radius)' 

    if region_type == 'Polygon':

        if Polygon_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) +  'km'

        if Polygon_Location != 'None':
            Polygon_Location_actual = Polygon_Location.split('-')
            SupTitle_text = 'Earthquakes M>' + str(MagLo) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) + ' km'

    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'From: ' + date_string[0] + '  ' + time_string[0] + '   To:  ' + date_string[range_limit] + '  ' + time_string[range_limit] 
    plt.title(Title_text, fontsize=12)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')

    matplotlib.pyplot.savefig('Seismicity-Map.pdf')
    matplotlib.pyplot.savefig('Seismicity-Map.jpg')

    print ' '
    print '     Close plot window to continue...'

    plt.show()

    return None

    #.................................................................

def get_circle_data(Location):

#   Read Circle_Location and data from "Settings_File.txt"

    settings_params = get_settings()
    default_settings_params = settings_params

    region_type = settings_params[0]
    completeness_mag    =   float(settings_params[1])
    earthquake_depth    =   float(settings_params[2])

    if region_type == 'Circle':     #   If still a circle, proceed
        Circle_Location     =   settings_params[3]

    if region_type != 'Circle':     #   If not a circle, reset settings_params to null
        Circle_Location     =   'None'
        Circle_Lat          =   0.0
        Circle_Lng          =   0.0
        Radius_float        =   0.0

        settings_params     =   []
        settings_params.append(region_type)
        settings_params.append(completeness_mag)
        settings_params.append(earthquake_depth)
        settings_params.append(Circle_Location)
        settings_params.append(Circle_Lat)
        settings_params.append(Circle_Lng)
        settings_params.append(Radius_float)

   #
   #    Get the data for the circle: Read pre-defined locations file
   #

    input_file = open("circlelocations.txt", "r")
    i=0
    for line in input_file:
        i +=  1
    input_file.close()  # Put the file pointer back at top

    number_circle_locations = i

    # Create arrays of length i filled with zeros

    Circle_Location_file   = ["" for x in range(i)]

    CircleLat       = np.zeros(i)
    CircleLng       = np.zeros(i)
    Radius          = np.zeros(i)

    input_file = open("circlelocations.txt", "r")

    i=-1
    for line in input_file:
        i+=1
        line    = line.strip()
        items   = line.split(',')

#   print items for testing purposes

        items_array = np.asarray(items)

        Circle_Location_file[i]   = items_array[0]
        CircleLat[i]       = float(items_array[1])
        CircleLng[i]       = float(items_array[2])
        Radius[i]          = float(items_array[3])

    input_file.close()  # Put the file pointer back at top


    print ' '
    print '     Current pre-defined Locations are: '
    print ' '
    for j in range(0,number_circle_locations):
        print '        ', Circle_Location_file[j]
    print ' '
    print '     Current region is: '
    print '        ', Location
    print ' '
    print '     Current circle location is: ', Circle_Location
    print ' '
    print '     Pick a new circle location? (y/n):'
    new_location_resp   =   raw_input()

    Current_Circle_Location = Circle_Location

    if new_location_resp == 'y':
        print ' Enter location (Case sensitive: Overrides previous parameter set):'
        Circle_Location = raw_input()
        region_type = 'Circle'

    if new_location_resp != 'y':
        Circle_Location = Current_Circle_Location
        region_type = 'Circle'

    circle_location_flag = 0
    counter = 0

    while (circle_location_flag == 0):
        for j in range(0,number_circle_locations):

            counter += 1

            if Circle_Location == Circle_Location_file[j]:

                CircleCenterLat = CircleLat[j]
                CircleCenterLng = CircleLng[j]
                CircleRadius    = Radius[j]

                Circle_Lat         = float(CircleCenterLat)
                Circle_Lng         = float(CircleCenterLng)
                Radius_float       = float(CircleRadius)
 
                settings_params[0]  =   region_type
                settings_params[1]  =   completeness_mag
                settings_params[2]  =   earthquake_depth
                settings_params[3]  =   Circle_Location
                settings_params[4]  =   Circle_Lat
                settings_params[5]  =   Circle_Lng
                settings_params[6]  =   Radius_float
          
                circle_location_flag = 1

        if circle_location_flag == 0 and counter == number_circle_locations:
            circle_location_flag = 1
            settings_params = default_settings_params
            Circle_Location = settings_params[4]
            Circle_Lat      = float(settings_params[5])
            Circle_Lng      = float(settings_params[6])
            Radius_float    = float(settings_params[7])
            print ' '
            print '     Invalid location, try again...'
            print ' '
            print '     (Press any key to continue)'
            resp = raw_input()
            print ' '



#   Write new value of Circle_Location to "Settings_File.txt"

    save_settings(settings_params)

    return (Circle_Location, Circle_Lat, Circle_Lng, Radius_float)

    #.................................................................

def get_polygon_data(Location):

#   Read Polygon_Location and data from "Settings_File.txt"

    settings_params = get_settings()
    default_settings_params = settings_params

    region_type = settings_params[0]
    completeness_mag    =   float(settings_params[1])
    earthquake_depth    =   float(settings_params[2])

    polygon_data = {}

    if region_type == 'Polygon':     #   If still a Polygon, proceed
        Polygon_Location     =   settings_params[3]

    if region_type != 'Polygon':     #   If not a Polygon, reset settings_params to null
        Polygon_Location    =   'None'

    polygon_data[0] = ['None', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0',]    # Dummy List = value of dict for index[0]
    default_vertex_points = polygon_data[0]


   #
   #    Get the data for the Polygon: Read pre-defined locations file
   #

    input_file = open("polygonlocations.txt", "r")
    i=0
    for line in input_file:
        i +=  1
    input_file.close()  # Put the file pointer back at top

    number_polygon_locations = i

    # Create arrays of length i filled with zeros

    Polygon_Location_List   = ["" for x in range(i)]

    input_file = open("polygonlocations.txt", "r")

    polygon_data = {}

    i=-1
    for line in input_file:
        i+=1
        line    = line.strip()
        items   = line.split(',')
        Polygon_Location_List[i] = items[0]

#   print items for testing purposes

        vertex_points = []

        for j in range(0,len(items)):
            vertex_points.append(items[j])

        polygon_data[i] = vertex_points

    input_file.close()  # Put the file pointer back at top

    print ' '
    print '     Current pre-defined Locations are: '
    print ' '

    for i in range(0,number_polygon_locations):
        print '        ', Polygon_Location_List[i]

    print ' '
    print '     Current region is: '
    print '        ', Location
    print ' '
    print '     Current polygon location is: ', Polygon_Location
    print ' '
    print '     Pick a new polygon location? (y/n):'
    new_location_resp   =   raw_input()

    Current_Polygon_Location = Polygon_Location

    if new_location_resp == 'y':
        print ' Enter location (Case sensitive: Overrides previous parameter set):'
        Polygon_Location = raw_input()
        region_type = 'Polygon'
        print 'Polygon_Location: ', Polygon_Location

    if new_location_resp != 'y':
        region_type = 'Polygon'
        print ' '
        print ' Then we use the current region (which could be the dummy region)'
        print ' '

    polygon_location_flag = 0
    counter = 0
    settings_params = []

    while (polygon_location_flag == 0):

        for i in range(0,number_polygon_locations):

            counter += 1

            if Polygon_Location == polygon_data[i][0]:
                polygon_location_flag = 1

                settings_params.append(region_type)
                settings_params.append(completeness_mag)
                settings_params.append(earthquake_depth)

                vertex_points = polygon_data[i]

                for j in range(0,len(vertex_points)):
                    settings_params.append(vertex_points[j])

        if polygon_location_flag == 0 and counter >= number_polygon_locations:
            polygon_location_flag = 1
            settings_params = default_settings_params
            vertex_points = default_vertex_points
            print ' '
            print '     Invalid location, try again...'
            print ' '
            print '     (Press any key to continue)'
            resp = raw_input()
            print ' '

#   Write new value of Polygon_Location to "Settings_File.txt"

    save_settings(settings_params)

    return (vertex_points)

    #.................................................................

def proxy_strain_region(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_large_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_large_eqs)]
    yrs                 =   zeros(number_large_eqs)
    lng                 =   zeros(number_large_eqs)
    lat                 =   zeros(number_large_eqs)
    mag                 =   zeros(number_large_eqs)
    dep                 =   zeros(number_large_eqs)
    max_number          =   zeros(number_large_eqs)

    time_string         =   ["" for x in range(number_large_eqs)]
    date_string         =   ["" for x in range(number_large_eqs)]

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number

    input_file.close()  # Put the file pointer back at top

    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq


    #.................................................................


    #   Now open the complete data file so that we can retrieve the small earthquakes

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    small_eq_count = 0
    mag_flag = 'FALSE'
    for line in data_file:
        items = line.strip().split()
        mag_query = float(items[5])
        if mag_query >= MagLo:
            mag_flag = 'TRUE'
        if mag_flag == 'TRUE':
            small_eq_count +=  1

    data_file.close()  # Put the file pointer back at top

    number_small_eq_array   =   np.ones(small_eq_count + number_large_eqs + 1)              #   Create small eq array and fill with 1's
    time_small_eq_array     =   np.zeros(small_eq_count+ number_large_eqs + 1)
    mag_small_eq_array      =   np.zeros(small_eq_count+ number_large_eqs + 1)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    index_small_eq = 0
    mag_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])
        time_small_eq       = float(items[2])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if (mag_query < MagLo) and (mag_flag == 'TRUE'):
            index_small_eq += 1
            number_small_eq_array[index_small_eq] = number_small_eq_array[index_small_eq - 1] + 1
            time_small_eq_array[index_small_eq] = time_small_eq
            mag_small_eq_array[index_small_eq]  =  mag_query

        if mag_query >= MagLo:
            mag_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1

            index_small_eq += 1
            number_small_eq_array[index_small_eq]   = number_small_eq_array[index_small_eq-1]
            time_small_eq_array[index_small_eq]     = time_small_eq
            mag_small_eq_array[index_small_eq]      = mag_query
            max_number[index_large_eq] = np.log10(number_small_eq_array[index_small_eq])

            index_small_eq += 1
            number_small_eq_array[index_small_eq] = 1   #   Re-initialize the small earthquake number
            time_small_eq_array[index_small_eq]     = time_small_eq
            mag_small_eq_array[index_small_eq] = mag_query

    data_file.close()               #   Close the data file

    max_number[0] = -100.   #   This one is zero, remove it

    number_small_eq_array = np.log10(number_small_eq_array)

#    fig, ax1 = plt.subplots()

    max_plot = max(number_small_eq_array) + 1

    range_limit = len(lat) - 1

    min_time = min(yrs) - 5.0
    max_time = max(time_small_eq_array) + 5.0


    plot_lines = 'TRUE'

    plt.subplot(211)

    SupTitle_text = 'Proxy Strain & Magnitude vs. Time for M>' + str(MagLo) + ' in ' + Location
    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'From: ' + date_string[0] + ' ' + time_string[0] + '     To: ' + date_string[range_limit] + ' ' + time_string[range_limit]
    plt.title(Title_text)


    if plot_lines == 'TRUE':
        for i in range(0,range_limit+1):
            x1=yrs[i]
            x2=yrs[i]
            y1=0.0
            y2=max_plot
    #        plt.plot([x1, x2], [y1, y2], 'b', ls='dotted', lw=0.9)

    #   Plot with multiple axes

    plt.plot(time_small_eq_array, number_small_eq_array, 'k')
#    ax=plt.gca()    #   1 plot with 1 row and 1 column

#    ax1 = ax.twinx()   #   If you want two axes with second one on right side

    plt.ylim([0.,max_plot])
    plt.xlim([min_time,max_time])
    plt.ylabel('Proxy Strain')

    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_limit):
            plt.plot(yrs, max_number, 'bo', ms=6)

    last_time_index     = len(time_small_eq_array) - 1.0
    last_time           = time_small_eq_array[last_time_index]
    last_proxy_strain   =  number_small_eq_array[last_time_index]

    if plot_lines != 'TRUE':
        plt.plot(last_time, last_proxy_strain, 'ro', lw=1.15)

#    ax.set_xlabel('Time (Years)')
#    ax1.set_ylim([5,8])
#    ax1.set_ylabel('Magnitude')
#    plt.plot(time_small_eq_array, number_small_eq_array, 'k')

    max_plot_mag = float(int(max(mag))) + 1.0
    min_plot_mag = MagLo - 0.5

    plt.subplot(212)
#   ax1=plt.gca()

    plt.ylim([min_plot_mag,max_plot_mag])
    plt.xlim([min_time,max_time])

    plt.ylabel('Magnitude')
    plt.xlabel('Time (Years)')

    plot_lines = 'TRUE'

    if plot_lines == 'TRUE':
        for i in range(0,range_limit+1):
            x1=yrs[i]
            x2=yrs[i]
            y1=min_plot_mag
            y2=max_plot_mag
    #        plt.plot([x1, x2], [y1, y2], 'b', ls='dotted', lw=0.9)


    if plot_lines == 'TRUE':
        for i in range(0,range_limit+1):
            x1=yrs[i]
            x2=yrs[i]
            y1=0.0
            y2=mag[i]
            plt.plot([x1, x2], [y1, y2], 'r-', lw=1.15)


    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_limit):
            plt.plot(yrs, mag, 'bo', ms=6, lw=1.0)

    plt.savefig('Proxy_Strain-Magnitude-Time.png')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    #.................................................................


    return


    #.................................................................


def proxy_strain_circle(NELat, NELng, SWLat, SWLng, MagLo, Location):

    # .................................................................
    #
    #    Get the data on the circular region

    Circle_Location, Circle_Lat, Circle_Lng, Radius_float = get_circle_data(Location)

    mag_array, date_array, time_array, year_array = get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, MagLo)

    #    Reverse the magnitude array for the magnitudes of events in the circular region

    print_data = 'TRUE'

    #   Open input file

    input_file = open("EQ_Working_Circle_Catalog.txt", "r")

    #   Find the number of lines in the file

    i = 0
    for line in input_file:
        i += 1

    input_file.close()  # Put the file pointer back at top

    number_large_eqs = i

    input_file = open("EQ_Working_Circle_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string = ["" for x in range(number_large_eqs)]
    yrs = zeros(number_large_eqs)
    lng = zeros(number_large_eqs)
    lat = zeros(number_large_eqs)
    mag = zeros(number_large_eqs)
    dep = zeros(number_large_eqs)
    max_number = zeros(number_large_eqs)

    time_string = ["" for x in range(number_large_eqs)]
    date_string = ["" for x in range(number_large_eqs)]

    # Loop over lines and extract variables of interest

    i = -1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i += 1

        indx_string[i] = data_array[0]
        date_string[i] = data_array[1]
        time_string[i] = data_array[2]

        yrs[i] = float(data_array[3])
        lng[i] = float(data_array[4])
        lat[i] = float(data_array[5])
        mag[i] = float(data_array[6])
        dep[i] = float(data_array[7])

        # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number

    input_file.close()  # Put the file pointer back at top

    # .................................................................

        #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq = int(0)
    second_eq = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq

    # .................................................................


    #   Now open the complete data file so that we can retrieve the small earthquakes

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    data_file = open("ANSS_%s.circle.catalog" % datetime.date.today().strftime("%F"),
                     "r")  # This is all the data and all the small events

    small_eq_count = 0
    mag_flag = 'FALSE'
    for line in data_file:
        items = line.strip().split()
        mag_query = float(items[5])
        if mag_query >= MagLo:
            mag_flag = 'TRUE'
        if mag_flag == 'TRUE':
            small_eq_count += 1

    data_file.close()  # Put the file pointer back at top

    number_small_eq_array = np.ones(small_eq_count + number_large_eqs + 1)  # Create small eq array and fill with 1's
    time_small_eq_array = np.zeros(small_eq_count + number_large_eqs + 1)
    mag_small_eq_array = np.zeros(small_eq_count + number_large_eqs + 1)

    data_file = open("ANSS_%s.circle.catalog" % datetime.date.today().strftime("%F"),
                     "r")  # This is all the data and all the small events

    index_large_eq = -1  # So the first large earthquake cycle will have index 0
    index_small_eq = 0
    mag_flag = 'FALSE'

    for line in data_file:
        items = line.strip().split()
        mag_query = float(items[5])
        time_small_eq = float(items[2])

        if mag_query < MagLo:
            date_small_quake = items[0]
            time_small_quake = items[1]

        if (mag_query < MagLo) and (mag_flag == 'TRUE'):
            index_small_eq += 1
            number_small_eq_array[index_small_eq] = number_small_eq_array[index_small_eq - 1] + 1
            time_small_eq_array[index_small_eq] = time_small_eq
            mag_small_eq_array[index_small_eq] = mag_query

        if mag_query >= MagLo:
            mag_flag = 'TRUE'  # Takes care of initial state where there is no initial large earthquake
            number_in_cycle = 0
            index_large_eq += 1

            index_small_eq += 1
            number_small_eq_array[index_small_eq] = number_small_eq_array[index_small_eq - 1]
            time_small_eq_array[index_small_eq] = time_small_eq
            mag_small_eq_array[index_small_eq] = mag_query
            max_number[index_large_eq] = np.log10(number_small_eq_array[index_small_eq])

            index_small_eq += 1
            number_small_eq_array[index_small_eq] = 1  # Re-initialize the small earthquake number
            time_small_eq_array[index_small_eq] = time_small_eq
            mag_small_eq_array[index_small_eq] = mag_query

    data_file.close()  # Close the data file

    max_number[0] = -100.  # This one is zero, remove it

    number_small_eq_array = np.log10(number_small_eq_array)

    #    fig, ax1 = plt.subplots()

    max_plot = max(number_small_eq_array) + 1

    range_limit = len(lat) - 1

    min_time = min(yrs) - 5.0
    max_time = max(time_small_eq_array) + 5.0

    plot_lines = 'TRUE'

    plt.subplot(211)

    Circle_Location_actual = Circle_Location.split('-')
    todays_date = datetime.date.today().strftime("%B %d, %Y")
    SupTitle_text = 'Proxy Strain & Magnitude vs. Time for M>' + str(MagLo)

    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'Within ' + str(Radius_float) + ' km of ' + Circle_Location_actual[0]
    plt.title(Title_text)

    if plot_lines == 'TRUE':
        for i in range(0, range_limit + 1):
            x1 = yrs[i]
            x2 = yrs[i]
            y1 = 0.0
            y2 = max_plot
    # plt.plot([x1, x2], [y1, y2], 'b', ls='dotted', lw=0.9)

    #   Plot with multiple axes

    plt.plot(time_small_eq_array, number_small_eq_array, 'k')
    ax = plt.gca()  # 1 plot with 1 row and 1 column

    #    ax1 = ax.twinx()   #   If you want two axes with second one on right side

    ax.set_ylim([0., max_plot])
    ax.set_xlim([min_time, max_time])
    ax.set_ylabel('Proxy Strain')

    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0, range_limit):
            plt.plot(yrs, max_number, 'bo', ms=6)

    last_time_index = len(time_small_eq_array) - 1
    last_time = time_small_eq_array[last_time_index]
    last_proxy_strain = number_small_eq_array[last_time_index]

    if plot_lines != 'TRUE':
        plt.plot(last_time, last_proxy_strain, 'r', marker='o', ms=6)


        #    ax.set_xlabel('Time (Years)')
        #    ax1.set_ylim([5,8])
        #    ax1.set_ylabel('Magnitude')
        #    plt.plot(time_small_eq_array, number_small_eq_array, 'k')

    max_plot_mag = float(int(max(mag))) + 1.0
    min_plot_mag = MagLo - 0.5

    plt.subplot(212)
    ax1 = plt.gca()

    ax1.set_ylim([min_plot_mag, max_plot_mag])
    ax1.set_xlim([min_time, max_time])

    ax1.set_ylabel('Magnitude')
    ax1.set_xlabel('Time (Years)')

    plot_lines = 'TRUE'

    if plot_lines == 'TRUE':
        for i in range(0, range_limit + 1):
            x1 = yrs[i]
            x2 = yrs[i]
            y1 = min_plot_mag
            y2 = max_plot_mag
    #       plt.plot([x1, x2], [y1, y2], 'b', ls='dotted', lw=0.9)

    if plot_lines == 'TRUE':
        for i in range(0, range_limit + 1):
            x1 = yrs[i]
            x2 = yrs[i]
            y1 = 0.0
            y2 = mag[i]
            plt.plot([x1, x2], [y1, y2], 'r-', lw=1.15)

    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0, range_limit):
            plt.plot(yrs, mag, 'bo', ms=6)

    plt.savefig('Proxy_Strain-Magnitude-Time.pdf')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    return

    ######################################################################

def filtered_histogram_eps_region(NELat, NELng, SWLat, SWLng, MagLo, Location):

    print_data = 'TRUE'

    last_droplet_plot   =   'n'


    #   Open input file

    input_file = open("EQ_Working_Catalog.txt", "r")

    #   Find the number of lines in the file

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs = i

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(number_eqs)]
    yrs                 =   zeros(number_eqs)
    lng                 =   zeros(number_eqs)
    lat                 =   zeros(number_eqs)
    mag                 =   zeros(number_eqs)
    dep                 =   zeros(number_eqs)

    time_string         =   ["" for x in range(number_eqs)]
    date_string         =   ["" for x in range(number_eqs)]

    # Bins for Number-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1

    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    # Loop over lines and extract variables of interest

    i=-1
    eq_sequence_number = 0

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    # Format and print the earthquake sequence numbers

        if mag[i] >= MagLo:
            eq_sequence_number += 1

        last_eq = eq_sequence_number

    mean_mainshock_interval =  (yrs[i] - yrs[0])/float(i)

    mean_mainshock_interval = mean_mainshock_interval * 12.0

    input_file.close()  # Put the file pointer back at top

    #.................................................................

     #   The large earthquakes were renumbered to start from 1 so need to correct for that

    first_eq    = int(0)
    second_eq   = int(last_eq) - 1

    number_of_eq_cycles = last_eq - first_eq

    earthquake_tag  =   np.zeros(number_eqs, dtype=int)

    #.................................................................

    #   Insert some code here to tag the cycles as "acceptable" or "not"
    #       depending on whether the calendar time constraint is satisfied.
    #   
    #   Constraint:  if the time interval between two large quakes is >=
    #       the specified interval
    #
    #   Acceptable: earthquake_tag[i] = 1    Unacceptable: earthquake_tag[i] = 0

    #   EQ_Condition: 
    #
    #                   = 1 if next large earthquake occurs after time_diff after current earthquake
    #                   = 2 if current large earthquake follows at least time_diff after previous earthquake
    #                   = 3 if next large earthquake occurs before time_diff after current earthquake
    #                   = 4 if current large earthquake follows no more than time_diff after previous


    print ''
    print 'Enter Time Delay Between Mainshocks (months): '
    time_fraction = raw_input()
    time_diff = float(time_fraction)/12.0    # Convert to fractions of a year

    print 'time_diff: ', time_diff

    EQ_Condition = 1

    if EQ_Condition == 1:
    
        for i in range(0,number_eqs-1):
            diff_yrs = yrs[i+1] - yrs[i]
            if diff_yrs >= time_diff:
                earthquake_tag[i] = 1

    if EQ_Condition == 2:
        for i in range(1,number_eqs-1):
            diff_yrs = yrs[i] - yrs[i-1]
            if diff_yrs >= time_diff:
                earthquake_tag[i] = 1

    if EQ_Condition == 3:
    
        for i in range(0,number_eqs-1):
            diff_yrs = yrs[i+1] - yrs[i]
            if diff_yrs <= time_diff:
                earthquake_tag[i] = 1

    if EQ_Condition == 4:
        for i in range(1,number_eqs-1):
            diff_yrs = yrs[i] - yrs[i-1]
            if diff_yrs <= time_diff:
                earthquake_tag[i] = 1


    number_acceptable_cycles = np.sum(earthquake_tag)

    #.................................................................
    

    #   Now open the complete data file so that we can retrieve the small earthquakes    

    #   Count the number of small earthquakes in each cycle and store in array for binning
    #       and plotting in histogram

    number_small_eq_array   =   np.zeros(number_eqs)

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    index_large_eq = -1             #   So the first large earthquake cycle will have index 0
    large_eq_flag = 'FALSE'

    for line in data_file:
        items   = line.strip().split()
        mag_query           = float(items[5])

        if mag_query < MagLo:
            date_small_quake    = items[0]
            time_small_quake    = items[1]

        if mag_query >= MagLo:
            large_eq_flag = 'TRUE'  #   Takes care of initial state where there is no initial large earthquake
            number_in_cycle =  0
            index_large_eq += 1

        if (mag_query < MagLo) and (large_eq_flag == 'TRUE'):
            number_small_eq_array[index_large_eq] += 1
    
    data_file.close()               #   Close the data file

    print ' '
    print '       Large EQ', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude', '    Number Small EQs (After this large EQ)'
    print ' '

    total_number_small_eqs = 0
    eq_limit = second_eq-first_eq +1  
    for i in range(first_eq,eq_limit):
        total_number_small_eqs += number_small_eq_array[i]
        if i+1 <10:
            blank_space = '         '
        if i+1 >=10:
            blank_space = '        '
        if i+1 >=100:
            blank_space = '       '
        print blank_space, i+1, '     ', date_string[i] + ' ' + time_string[i], '     ', '%.2f'%mag[i], '     ', '%.3f'%(lat[i]), '    ', '%.3f'%lng[i],  \
                '          ',int(number_small_eq_array[i])


    #.................................................................

    range_limit = len(lat)

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' in ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print ' '
    print print_text_1
    print ' '
    print print_text_2
    print print_text_3
    print ''
    print '     Total Number of Small Earthquakes: ', int(total_number_small_eqs)        
    print ''


   #.................................................................

    #   Histogram the data

    #   n[i] below is an array with the number of intervals in bin [i]
    #   bins[i] is an array with the location of the bin center
    #   patches[i] is an array with descriptions of each histogram rectangle
    #   set normed = 1 if you want the histogram to be a pdf

    #    n, bins, patches = plt.hist(number_small_eq_array, num_bins, facecolor='green', alpha=0.5)

    todays_count = number_small_eq_array[last_eq-1]

    number_small_eqs_excluding_last = np.zeros(number_acceptable_cycles)

    j = 0
    for i in range(0,number_eqs-1):
        if earthquake_tag[i] == 1:
            number_small_eqs_excluding_last[j] = number_small_eq_array[i]
            j += 1

    mean    =   mean_val(number_small_eqs_excluding_last)
    std_deviation, variance = std_var_val(number_small_eqs_excluding_last)

    mean_small_original = mean_val(number_small_eq_array)

    mean_eqs = mean

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  
    mean_small_EQs_original = '%.2f'%mean_small_original

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print '   Mean Number Small Earthquakes Original: ', mean_small_EQs_original
    print ''

    #.................................................................

    #   Bin the data

    num_bins = 500

    max_bins = max(number_small_eq_array)

    n, bins = histogram(number_small_eqs_excluding_last, range=(0,max_bins), bins=num_bins)

    #   -----------------------------------------

    todays_probability = 0.0

    cum_poisson = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_small_EQs))

    cum_poisson[:] = (cum_poisson[:]) * 100.0

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:]/ float(last_eq)

    scale_factor = 100.0/cum_prob[len(bins)-1]

    cum_prob[:] = cum_prob[:]*scale_factor	#   Ensures that cum prob goes from 0%->100%

#    cum_prob[:] = cum_prob[:] *100.0
    cum_prob[:] = cum_prob[:]

    for i in range(1, len(bins)):
        if todays_count >= bins[i]:
            todays_probability = cum_prob[i]

#    todays_probability = 100.0 - todays_probability


    tau_best, beta_best, sdev_tau, sdev_beta, sum_of_squares_low = weibullBetaFit(n, bins, cum_prob, mean_eqs)

    tau_string   = "{0:.2f}".format(tau_best) + ' +/- '  + "{0:.2f}".format(sdev_tau)
    beta_string  = "{0:.2f}".format(beta_best) + ' +/- ' + "{0:.2f}".format(sdev_beta)

    print ''
    print ' Weibull fit data: '
    print ''
    print ' Tau: ', tau_string
    print 'Beta: ', beta_string
    print ''

    cum_weibull =   np.zeros(len(bins))

    for i in range(0,len(bins)):
        arg_weibull = ((bins[i]/tau_best)**beta_best)
        cum_weibull[i] = 1.0 - math.exp( - arg_weibull)

#        	cum_weibull[i] = math.exp(- arg_weibull)

    for i in range(0,len(bins)):
#        cum_prob[i] = 100.0 - cum_prob[i]
        cum_prob[i] = cum_prob[i]        

    for i in range(0,len(bins)):
#        cum_poisson[i] = 100.0 - cum_poisson[i]
        cum_poisson[i] = cum_poisson[i]

    cum_weibull[:] = cum_weibull[:] *100.0

    #.................................................................

    explanation1 = '>> For EQ_Condition = 1, this chart shows the probability that, for a given '
    explanation2 = '   minimum current earthquake interoccurrance time of ' + time_fraction + ' month(s), '
    explanation3 = '   there is less than a ' + '%.1f'%(todays_probability) + '% probability that the '
    explanation4 = '   current cycle will contain only ' + str(todays_count) + ' small earthquakes'

    print ''
    print explanation1
    print explanation2
    print explanation3
    print explanation4
    print ''

    #.................................................................

    #   Define the Temperature

    temperature = "{0:.1f}".format(todays_probability)+'%'

    #.................................................................

    #  Plot the data

    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[14, 1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])

    #   Draw the bins

    ax0.hist(number_small_eqs_excluding_last, num_bins, facecolor='green', alpha=0.5)
    plt.ylim((0,100))

    #   Get the axis limits

    ymin, ymax = ax0.get_ylim()
    xmin, xmax = ax0.get_xlim()

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
    SupTitle_text = 'Potential for M>' + str(MagLo) + ' Earthquakes in ' + Location + ' on ' + todays_date
    SupTitle_text = 'Probability for M>' + str(MagLo) + ' Earthquakes with Mainshock Time Delay > ' + time_fraction + ' Month(s)'
    Title_text    = 'After M' + '%.2f'%mag[last_eq-1] + ' on ' + date_string[last_eq-1] + ' at ' + time_string[last_eq-1] + ' in ' + Location

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.90 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.85 + ymin

    text_x3          = (xmax-xmin)*0.40 + xmin
    text_y3          = (ymax-ymin) * 0.80 + ymin

    count_string    = str(int(todays_count))
    text_string1 = 'Todays Small EQ Count:  ' + count_string 
    text_string2 = 'On:  ' + date_small_quake + '   at:  ' + time_small_quake
    interval     = '%.3f'%(mean_mainshock_interval)
    text_string3 = 'Mean Mainshock Interval: ' + interval + ' Months'

    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    
    plt.text(text_x3,text_y3,text_string3, fontsize=12)    


    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.72 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string = 'Todays DEPS Value:  ' + prob +'%'
    plt.text(text_x,text_y,text_string, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.65 + ymin
    text_string = 'Weibull Parameters: '
    plt.text(text_x,text_y,text_string, fontsize=12)    

    text_x1          = (xmax-xmin)*0.40 + xmin
    text_y1          = (ymax-ymin) * 0.60 + ymin

    text_x2          = (xmax-xmin)*0.40 + xmin
    text_y2          = (ymax-ymin) * 0.55 + ymin

    count_string    = str(int(todays_count))
    text_string1 = '    Tau:  ' + tau_string 
    text_string2 = '   Beta:  ' + beta_string
    plt.text(text_x1,text_y1,text_string1, fontsize=12)    
    plt.text(text_x2,text_y2,text_string2, fontsize=12)    

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12) 

    bins        = np.append(bins,bins[num_bins])

    last_prob_value = cum_prob[len(bins)-2]
    cum_prob    = np.append(cum_prob, last_prob_value)

    last_poisson_value = cum_poisson[len(bins)-2]
    cum_poisson = np.append(cum_poisson, last_poisson_value)

    last_weibull_value = cum_weibull[len(bins)-2]
    cum_weibull = np.append(cum_weibull, last_weibull_value)

    ax1 = ax0.twinx()

    plt.ylim(ymax = 100, ymin = 0)
    ymin, ymax = ax1.get_ylim()

    ax1.plot(bins,cum_weibull, 'g-', lw=1.2)    #   Uncomment if you want to plot the Weibull curve with same mean
    ax1.plot(bins,cum_poisson, 'b--')    #   Uncomment if you want to plot the Poisson curve with same mean
    ax1.plot(bins,cum_prob, 'r-')

    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)
    x=[mean_small_EQs,mean_small_EQs]
    y=[0,100]
    ax1.plot(x,y, 'k', ls='dotted')


    #.................................................................

    #   Thermometer

    ax3 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([0,1])                         #   Show the y-axis labels
    plt.xlim([0,1])                         #   Set the x-axis limits  

    frame1.axes.get_xaxis().set_ticks([])   #   Hide the x-axis ticks and labels
    frame1.axes.get_yaxis().set_ticks([])

    ax4 = ax3.twinx()
    plt.ylim([0,100])                         #   Show the y-axis labels
    plt.xlim([0,1])

    ax4.get_yaxis().set_ticks([])
#    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted', lw=1)

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays EPS Value (Cumulative Probability) is: ', temperature
    print ''

    matplotlib.pyplot.savefig('EPS-Region-Filtered.pdf')

    plt.show()

    return None



# .................................................................




