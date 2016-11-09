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

    completeness_mag = 2.99

    print ' '
    print ' Current value of catalog completeness magnitude is: ', completeness_mag
    print ' Do you want to change it? (y/n)'
    resp_completeness = raw_input()

    if resp_completeness == 'y':
        print ' Enter new value of completeness magnitude'
        completeness_mag = float(raw_input())

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

    plt.xlabel('Time (years)', fontsize=14)
    plt.ylabel('Magnitude', fontsize=14)

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
            plt.plot(yrs, mag, 'bo', ms=6)


    plt.savefig('Magnitude-Time.png')

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



    plt.savefig('Number-Magnitude.png')

    print ' '
    print '     Close plot window to continue...'
    print ' '
    print '     .......................................'

    plt.show()

    return None

    ######################################################################

def histogram_clock(NELat, NELng, SWLat, SWLng, MagLo, Location):

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

#    print ' '
#    print '       Sequence', '      Date and Time', '       Magnitude', '   Latitude', '    Longitude'
#    print ' '

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
        mag_query     = float(items[5])
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

    num_bins = 100

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

    mean_small_EQs = '%.2f'%mean
    std_dev_small_EQs = '%.2f'%std_deviation  

    print '         Mean Number of Small Earthquakes: ', mean_small_EQs
    print '     Standard Deviation Small Earthquakes: ', std_dev_small_EQs
    print ''

    #.................................................................

    #   Bin the data

    n, bins = histogram(number_small_eqs_excluding_last, num_bins)

    cum_poisson = np.zeros(len(bins))   #   Cumulative Poisson distribution

    cum_prob = np.zeros(len(bins))

    for i in range(1,len(bins)):
        cum_prob[i] = cum_prob[i-1] +  n[i-1]
        cum_poisson[i] = 1.0 - math.exp(-bins[i]/float(mean_small_EQs))

    cum_poisson[:] = cum_poisson[:] * 100.0

    print cum_poisson[:]

    #   Calc cumulative probability

    cum_prob[:] = cum_prob[:] / float(last_eq)
    cum_prob[:] = cum_prob[:] *100.0

    for i in range(1, len(bins)):
        if todays_count > bins[i]:
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

    todays_date          =  datetime.date.today().strftime("%B %d, %Y")
    SupTitle_text = 'Risk of M>' + str(MagLo) + ' Earthquakes in ' + Location + ' on ' + todays_date
    Title_text    = 'After M' + '%.2f'%mag[last_eq-1] + ' on ' + date_string[last_eq-1] + ' at ' + time_string[last_eq-1]

    plt.xlabel('Number of Small Earthquakes Between Large Earthquakes')
    plt.ylabel('Number of Earthquake Intervals')

    #   Write the legends on the plot

    text_x          = (xmax-xmin)*0.40 + xmin
    text_y          = (ymax-ymin) * 0.90 + ymin
    count_string    = str(int(todays_count))
    text_string = 'Todays Small EQ Count:  ' + count_string
    plt.text(text_x,text_y,text_string, fontsize=12)    

    text_x      = (xmax-xmin)*0.40 + xmin
    text_y      = (ymax-ymin) * 0.85 + ymin
    prob        = '%.1f'%(todays_probability)
    text_string = 'Todays Probability:  ' + prob +'%'
    plt.text(text_x,text_y,text_string, fontsize=12)    

    #.................................................................

    # Main plot

    plt.suptitle(SupTitle_text, fontsize=14)
    plt.title(Title_text, fontsize=12)

    ax1 = ax0.twinx()
    ax1.plot(bins,cum_prob, 'r-')
#   ax1.plot(bins,cum_poisson, 'b-')    #   Uncomment if you want to plot the Poisson curve with same mean
    ax1.get_yaxis().set_ticks([0.,25, 50, 75, 100])
    ax1.plot([xmin,xmax],[50,50], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[25,25], 'b', ls='dotted')
    ax1.plot([xmin,xmax],[75,75], 'b', ls='dotted')
    ax1.plot([todays_count,xmax],[todays_probability,todays_probability], 'r', ls='dashed', lw=1)
    x=[todays_count, todays_count]
    y=[0,todays_probability]
    ax1.plot(x,y, 'r', ls='dashed', lw=1)
    ax1.plot([todays_count], [todays_probability], 'ro', ms=8)

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
    ax4.set_ylabel('Seismic "Temperature" = Current Cumulative Probability (%)', rotation=90)
    x=[0.0,1.0]
    y=[todays_probability,todays_probability]
    ax4.fill_between(x,0,y, alpha=0.75, facecolor='red')
    
    ax4.plot([0,1],[50,50], 'b', ls='dotted')
    ax4.plot([0,1],[25,25], 'b', ls='dotted', lw=1)
    ax4.plot([0,1],[75,75], 'b', ls='dotted')

    #.................................................................

    print ''
    print 'Todays small earthquake count is: ', int(todays_count)

    print 'Todays Temperature (Cumulative Probability) is: ', temperature
    print ''

    matplotlib.pyplot.savefig('Earthquake-Hazard.png')

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



