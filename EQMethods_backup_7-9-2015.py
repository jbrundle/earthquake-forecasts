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
import numpy as np
from numpy import *
from mpl_toolkits.basemap import Basemap
from array import array
import matplotlib.pyplot as plt

import datetime
import dateutil.parser

import urllib
import os

import math


    ######################################################################



    # Set the catalog parameters

def get_catalog(NELat, NELng, SWLat, SWLng, MagLo):

    data = {

    # Adjust these  -   CA-NV
        "minmag": 3.0,
        "minlat": SWLat,
        "maxlat": NELat,
        "minlon": SWLng,
        "maxlon": NELng,
#        "mindepth": 0,         # Leaving these depth params in leads to missing some earthquake events
#        "maxdepth": 1000,


    # Date parameters
        "mintime": "1900/1/1",
        "maxtime": "2110/1/1",


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
#           print i
#           print i, date_string, time_string, ts, lon, lat, mag, dep
            working_file.write("%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()

    # Finalize the output file

    working_file.close()
    data_file.close()

    return None

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
    plt.xlim(yrs[0], yrs[range_limit])

    SupTitle_text = 'Earthquake Magnitude vs. Time for M>' + str(MagLo) + ' in ' + Location
    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'From: ' + date_string[0] + ' ' + time_string[0] + '     To: ' + date_string[range_limit] + ' ' + time_string[range_limit]
    plt.title(Title_text)

    plt.xlabel('Time (years)')
    plt.ylabel('Magnitude')

    plot_lines = 'TRUE'

    if plot_lines == 'TRUE':
        for i in range(0,range_limit):
            x1=yrs[i]
            x2=yrs[i]
            y1=0.0
            y2=mag[i]
            plt.plot([x1, x2], [y1, y2], 'r-', lw=1.15)


    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_limit):
            plt.plot(yrs, mag, 'bo', ms=4)


    plt.savefig('Magnitude-Time.png')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    return None

    ######################################################################

def freqmag(NELat, NELng, SWLat, SWLng, MagLo, Location):


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

    # Bins for frequency-magnitude plot

    min_mag = 3.0
    bin_diff= 0.1
    number_mag_bins     =   (MagLo - min_mag) / bin_diff + 1      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    freq_mag_bins_pdf   =   zeros(number_mag_bins)
    freq_mag_bins_sdf   =   zeros(number_mag_bins)

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

    #.................................................................

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

    print ' '
    print ' Enter Earthquake Sequence Numbers for Plot (Initial, Final):'
    Sequence_Numbers_Input = raw_input()
    Seq = Sequence_Numbers_Input.split(',')

    first_eq    = int(Seq[0]) - 1       # I renumbered the large earthquakes to start from 1 so need to correct for that
    second_eq   = int(Seq[1]) - 1
    print ' '
    print int(Seq[0]), int(Seq[1])

    number_of_eq_cycles = int(Seq[1]) - int(Seq[0])

    #   Now open the complete data file so that we can retrieve the small earthquakes        ----------------------------------------

    data_file = open("ANSS_%s.catalog" % datetime.date.today().strftime("%F"), "r")     #  This is all the data and all the small events

    sequence =0                                   # This will be the large earthquake number
    eq_id=0                                       # This is the small earthquake ID number

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
            sequence += 1                                           # Increments the sequence number

        if first_eq <= sequence < second_eq:
            if mag_small < MagLo:
                eq_id += 1                                          # Assigns an ID number to each small earthquake
                bin_number = int((mag_small - min_mag) * 10)
                freq_mag_bins_pdf[bin_number] += 1                  # Increment the appropriate bin to compute the GR PDF
    #           print sequence, eq_id, mag_small, dep_small, lat_small, lon_small, date_small, time_small

    data_file.close()                                               # Close the data file

    #.................................................................

    range_mag_bins = int(number_mag_bins - 1)

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]

    print ''

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                         # Find the number of nonzero bins
        if freq_mag_bins_pdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins-1)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins            =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff


    print ''
    print mag_bins
    print ''
    print freq_mag_bins_sdf
    print ''
    print log_freq_mag_bins
    print ''

    #.................................................................


    plotmin = 0.0
    plotmax = int(log_freq_mag_bins[0] + 1.0)

    plt.plot()
    plt.ylim(plotmin, plotmax)
    plt.xlim(min_mag, MagLo)

    SupTitle_text = 'Earthquake Frequency-Magnitude Relation in ' + Location
    plt.suptitle(SupTitle_text, fontsize=16)

#    Title_text =  +'     Final: ' + date_string[second_eq] + ' ' + time_string[second_eq] + ', Mag= '+ str(mag[second_eq])
#    plt.title(Title_text)

    text_x      = (MagLo-min_mag)*0.025 + min_mag
    text_y      = plotmax * 0.05
    text_string = 'Final EQ:  Date=' + date_string[second_eq] + ';  Time=' + time_string[second_eq] + ';   Mag='+ str(mag[second_eq])
    plt.text(text_x,text_y,text_string, fontsize=12)

    text_x      = (MagLo-min_mag)*0.025 + min_mag
    text_y      = plotmax * 0.10
    text_string = 'Initial EQ: Date=' + date_string[first_eq] + ';  Time=' + time_string[first_eq] + ';   Mag='+ str(mag[first_eq])
    plt.text(text_x,text_y,text_string, fontsize=12)

    text_x      = (MagLo-min_mag)*0.025 + min_mag
    text_y      = plotmax * 0.15
    text_string = 'Number of EQ Cycles: ' + str(int(Seq[1]) - int(Seq[0]))
    plt.text(text_x,text_y,text_string, fontsize=12)

    plt.xlabel('Magnitude')
    plt.ylabel('Frequency')

    plot_lines = ''
    if plot_lines != 'TRUE':
        for i in range(0,range_bins):
            plt.plot(mag_bins, log_freq_mag_bins, 'bo', ms=6)


    plt.savefig('Frequency-Magnitude.png')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

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

    z_limit = len(lat)
    range_limit = z_limit - 1

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

    for z in range(z_limit):
        x,y = m(lng,lat)

    #   -----------------------
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

    lat=[llcrnrlat + 0.25*sfv,llcrnrlat + 0.80*sfv, llcrnrlat + 1.30*sfv, llcrnrlat + 1.80*sfv, llcrnrlat + 2.40*sfv, llcrnrlat + 3.1*sfv, llcrnrlat + 3.95*sfv]
    lng=[urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht]
    x,y = m(lng,lat)

    latc=[llcrnrlat + 0.30*sfv,llcrnrlat + 0.85*sfv, llcrnrlat + 1.375*sfv, llcrnrlat + 1.90*sfv, llcrnrlat + 2.50*sfv, llcrnrlat + 3.2*sfv, llcrnrlat + 4.0*sfv]
    lngc=[urcrnrlon + 0.72*sfh, urcrnrlon + 0.73*sfh, urcrnrlon + 0.74*sfh, urcrnrlon + 0.76*sfh, urcrnrlon + 0.78*sfh, urcrnrlon + 0.80*sfh, urcrnrlon + 0.80*sfh]
    xc,yc = m(lngc,latc)

    mag_value=['5.5+','6+','6.5+','7+','7.5+', '8+', '8.5+']

    mark_size = [4, 6, 8, 10, 12, 14, 16]

    for z in range(7):
        plt.text(x[z],y[z],mag_value[z], fontsize=10,)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False)

    #   matplotlib.pyplot.savefig('map_bw.pdf')

    SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' in ' + Location
    plt.suptitle(SupTitle_text, fontsize=16)

    Title_text = 'From:  ' + date_string[0] + '  ' + time_string[0] + '   To:   ' + date_string[range_limit] + '  ' + time_string[range_limit]
    plt.title(Title_text, fontsize=14)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')

    matplotlib.pyplot.savefig('Seismicity-Map.png')

    print ' '
    print '     Close plot window to continue...'

    plt.show()

    return None

    #.................................................................


