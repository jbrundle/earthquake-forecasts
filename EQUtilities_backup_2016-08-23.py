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
from array import array

import os

import math
import scipy.special

    ######################################################################

def get_settings():

    working_file = open("Settings_File.txt", "r")

    for line in working_file:
        settings_params = line.strip().split()

    working_file.close()

    return settings_params

def save_settings(completeness_mag, Circle_Location, Circle_Lat, Circle_Lng, Radius_float, earthquake_depth):

    #   Saving the settings in a way that they can be reset from time to time

    working_file = open("Settings_File.txt", "w")

    working_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (completeness_mag, Circle_Location, Circle_Lat, Circle_Lng, \
                                                 Radius_float, earthquake_depth))

    working_file.close()

    return None


def mean_val(number_array):

    #   The obvious: returns the mean of an array of numbers

    N_Dim_mean = len(number_array)
    mean = sum(number_array)/float(N_Dim_mean)

    return mean

    #   .................................................................

def std_var_val(number_array):

    #   The obvious: returns the standard deviation of an array of numbers

    N_Dim_var = len(number_array)
    mean_of_array = mean_val(number_array)

    adjusted_arguments = [((i-mean_of_array)**2)/float(N_Dim_var-1) for i in number_array]   # Sample variance

    variance = sum(adjusted_arguments)   

    standard_deviation = math.sqrt(variance) 

    return (standard_deviation, variance)


def linfit(x, y, n):

#
#     This program fits a straight line to N data.  The data are
#     assumed to be error-free.  
#
#     Definitions:
#
#     x[i]:  Abscissas of the data
#     y[i]:  Ordinates of the data
#     slope: The resulting best fit slope
#     cept:  The resulting best fit y-intercept
#     errs:  The standard error for the slope
#     errc:  The standard error for the intercept
#     n:     The exact number of data to fit
#
#
#   NOTE!!:     The *exact number* n of data to fit *must equal* = len(x) = len(y)
#               Otherwise, the code will blow up
#
#
#    n = int(len(x))

    n = int(n)

    ata = [zeros(2),zeros(2)]
    aty = [zeros(2),zeros(2)]
    atainv = [zeros(2),zeros(2)]
#    
    sumx2 = 0.
    xsum = 0.
    yxsum = 0.
    ysum = 0.

    for i in range(0,n):
        sumx2 = sumx2 + x[i] * x[i]
        xsum = xsum + x[i]
        yxsum = yxsum + y[i] * x[i]
        ysum = ysum + y[i]
#
#    ata[0][0] = sumx2
#    ata[0][1] = xsum
#    ata[1][0] = xsum
#    ata[1][1] = float(n)

    ata[0][0] = sumx2
    ata[1][0] = xsum
    ata[0][1] = xsum
    ata[1][1] = float(n)
#
    aty[0] = yxsum
    aty[1] = ysum

#
    det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0]
    atainv[0][0] = ata[1][1]/det
    atainv[0][1] = -ata[0][1]/det
    atainv[1][0] = -ata[1][0]/det
    atainv[1][1] = ata[0][0]/det
#
    slope = atainv[0][0] * aty[0] + atainv[0][1] * aty[1]
    cept = atainv[1][0] * aty[0] + atainv[1][1] * aty[1]

    s2 = 0
    for i in range(0,n):
        s2 = s2 + (y[i] - cept - slope * x[i])**2

    s2 = s2 / (float(n) - 2.)
      
    errs = math.sqrt( float(n) * s2 / det )
    errc = math.sqrt( s2 * sumx2 / det)

    print slope, cept, errs, errc, s2

    return (slope, cept, errs, errc, s2)

    #   Usage in calling program:  slope, cept, errs, errc, s2 = linfit(x,y)

    #.................................................................

def droplet(x, y, x1, x2, nonzero_bins, MagLo):

#
#     This program computes a Fisher droplet-type model to the data (power law + decaying
#       exponentialto N data.  The data are assumed to be error-free.  
#
#     Definitions:
#
#     x[i]:  Abscissas of all the nonzero data points, not just the scaling region
#     y[i]:  Ordinates of all the nonzero data points, not just the scaling region
#            These are the log10 of the bin frequencies
#
#     slope: The resulting best fit slope of the scaling part
#     cept:  The resulting best fit y-intercept of the scaling part
#     errs:  The standard error for the slope
#     errc:  The standard error for the intercept
#     x1:    Left end of the scaling range
#     x2:    Right end of the scaling range
#
#     n = int(len(x))

    #   Find slope and intercept of small magnitude scaling part

    scaling_bins = 0
    droplet_bins = 0

    for i in range(0,nonzero_bins):
        if (x[i] >= x1) and (x[i] <= x2):
            scaling_bins += 1
        if x[i] >= x1:
            droplet_bins += 1

    xscale    =   zeros(scaling_bins)
    yscale    =   zeros(scaling_bins)

    k = -1
    for i in range(0,nonzero_bins):
        if (x[i] >= x1) and (x[i] <= x2):
            k = k + 1
            xscale[k] = x[i]
            yscale[k] = y[i]
            print ' k, xr, yr: ', k, xscale[k], yscale[k]

    kmax = k+1

    if kmax > 1:
        slope, cept, errs, errc, s2 = linfit(xscale, yscale, scaling_bins)
        bval = - slope
        print  ' '
        print  ' b - value is:'
        print   bval,' +/- ',errs
        print ' '
        print ' '

    x_droplet   =   zeros(droplet_bins)
    y_droplet   =   zeros(droplet_bins)

    j=-1

    xi_droplet      = MagLo
    sigma_droplet   = 0.66666       #   MagLo = 6.0
    sigma_droplet   = 0.50000       #   MagLo = 6.5
    sigma_droplet   = 0.33333       #   MagLo = 7.0     Seems to be the best overall
    exp_constant    = 1.0                #   Can be <1 for a better fit to Japan and Chile, for example

    sigma_droplet   = 0.33333 * (8.0 - MagLo)

#   sigma_droplet   = 0.83333       #   MagLo = 5.5     Code blowing up on 5.5 probably due to too many bins on abscissa of f-m plot

    for i in range(0,nonzero_bins):
        if x[i] >= x1:
            j+=1
            x_droplet[j] = x[i]
#           exp_factor = 10.0**(1.5*sigma_droplet*(x[i] - xi_droplet**0.985))   #   With sigma_droplet = 0.50000 and MagLo = 7.0
            exp_factor = 10.0**(1.5*sigma_droplet*(x[i] - xi_droplet))
            y_droplet[j] = cept + slope*x[i] - exp_constant*exp_factor
            print x[i], exp_factor

    return (x_droplet, y_droplet)   # These are arrays:  y_droplet = a - b*x_droplet - (x_droplet/xi)**sigma

    #   Usage in calling program:  slope, cept, errs, errc, s2 = linfit(x,y)

    #.................................................................

def dropletFit(x, y, x1, x2, nonzero_bins, MagLo):

#
#     This program fits a Fisher droplet-type model to the data (power law + decaying
#       exponentialto N data.  The data are assumed to be error-free.  
#
#     Definitions:
#
#     x[i]:  Abscissas of all the nonzero data points, not just the scaling region
#     y[i]:  Ordinates of all the nonzero data points, not just the scaling region
#            These are the log10 of the bin frequencies
#
#     slope: The resulting best fit slope of the scaling part
#     cept:  The resulting best fit y-intercept of the scaling part
#     errs:  The standard error for the slope
#     errc:  The standard error for the intercept
#     x1:    Left end of the scaling range
#     x2:    Right end of the scaling range
#
#     n = int(len(x))

    #   Find slope and intercept of small magnitude scaling part

    scaling_bins = 0
    droplet_bins = 0

    for i in range(0,nonzero_bins):
        if (x[i] >= x1) and (x[i] <= x2):
            scaling_bins += 1
        if x[i] >= x1:
            droplet_bins += 1

    xscale    =   zeros(scaling_bins)
    yscale    =   zeros(scaling_bins)

    k = -1
    for i in range(0,nonzero_bins):
        if (x[i] >= x1) and (x[i] <= x2):
            k = k + 1
            xscale[k] = x[i]
            yscale[k] = y[i]
            print ' k, xr, yr: ', k, xscale[k], yscale[k]

    kmax = k+1

    if kmax > 1:
        slope, cept, errs, errc, s2 = linfit(xscale, yscale, scaling_bins)
        bval = - slope

    x_droplet   =   zeros(droplet_bins)
    y_droplet   =   zeros(droplet_bins)

    sum_of_squares_low  =   100000000.0     #   Some large number

    delta_xi                =   0.1
    delta_droplet           =   0.01

# 

    ixi_range       =   30
    droplet_range   =   200

    xi_droplet = MagLo
    xi_droplet = MagLo + 0.5*(ixi_range)*delta_xi

    #   Loop to find the best fitting values of xi and sigma

    for ixi in range(0,ixi_range):
        sigma_droplet  =   1.0 + delta_droplet
        xi_droplet -= delta_xi
        for idrop in range(0,droplet_range):
            sigma_droplet -= delta_droplet

            j=-1
            diff = zeros(droplet_bins)
            for i in range(0,nonzero_bins):     # Cycle through all the droplet points
                if x[i] >= x1:
                    j+=1
                    x_droplet[j] = x[i]
                    exp_factor = 10.0**(1.5*sigma_droplet*(x[i] - xi_droplet))
                    y_droplet[j] = cept + slope*x[i] - exp_factor
                    diff[j] = y_droplet[j] - y[i]

    #
            sum_of_squares = sum(diff[:]*diff[:])

            if (sum_of_squares < sum_of_squares_low):
    #           print 'ixi, idrop, sum_of_squares, sum_of_squares_low, xi_droplet_best, sigma_droplet_best: ', ixi, idrop, sum_of_squares, sum_of_squares_low, xi_droplet, sigma_droplet
                sum_of_squares_low  =   sum_of_squares
                xi_droplet_best     =   xi_droplet
                sigma_droplet_best  =   sigma_droplet

    print ''
    print ' Best Fitting Values of sum, xi, sigma: ', sum_of_squares_low, xi_droplet_best, sigma_droplet_best
    print ''

    #   Compute best fitting droplet curve

    xi_droplet     =   xi_droplet_best
    sigma_droplet  =   sigma_droplet_best

    j=-1

    for i in range(0,nonzero_bins):
        if x[i] >= x1:
            j+=1
            x_droplet[j] = x[i]
            exp_factor = 10.0**(1.5*sigma_droplet*(x[i] - xi_droplet))
            y_droplet[j] = cept + slope*x[i] - exp_factor

    # These are arrays:  y_droplet = a - b*x_droplet - (x_droplet/xi)**sigma

    return (x_droplet, y_droplet, cept, slope, droplet_bins, xi_droplet_best, sigma_droplet_best, sum_of_squares_low)   

    #   Usage in calling program:  slope, cept, errs, errc, s2 = linfit(x,y)

    #.................................................................

def deviationFit(x, y, x1, x2, nonzero_bins, droplet_bins, cept, slope, sum_of_squares_low, xi_droplet_best, sigma_droplet_best):

#
#     This program takes a fit of the data to the Fisher droplet-type model a (power law + decaying
#       exponentialto N data) and estimates the variance of the fit around the best fitting values
#       of parameters.  The data are assumed to be error-free.  
#
#     Definitions:
#
#     x[i]:  Abscissas of all the nonzero data points, not just the scaling region
#     y[i]:  Ordinates of all the nonzero data points, not just the scaling region
#            These are the log10 of the bin frequencies
#
#     slope: The resulting best fit slope of the scaling part
#     cept:  The resulting best fit y-intercept of the scaling part
#     errs:  The standard error for the slope
#     errc:  The standard error for the intercept
#     x1:    Left end of the scaling range
#     x2:    Right end of the scaling range
#

    pi = 3.1415926535

    x_droplet   =   zeros(droplet_bins)
    y_droplet   =   zeros(droplet_bins)

    delta_sigma     =   0.01
    sigma_range     =   200
    sdev_sigma      =   0.0

    #   Estimate standard deviation of sigma first, assuming xi is best value
        
    sigma_droplet  =   sigma_droplet_best + 0.5*float(sigma_range)*delta_sigma
    dev_flag = 0

    for idrop in range(0,sigma_range):
        sigma_droplet -= delta_sigma
        j=-1
        diff = zeros(droplet_bins)
        for i in range(0,nonzero_bins):     # Cycle through all the droplet points
            if x[i] >= x1:
                j+=1
                x_droplet[j] = x[i]
                exp_factor = 10.0**(1.5*sigma_droplet*(x[i] - xi_droplet_best))
                y_droplet[j] = cept + slope*x[i] - exp_factor
                diff[j] = y_droplet[j] - y[i]

    #
        sum_of_squares      = sum(diff[:]*diff[:])
        ratio_of_squares    = (sum_of_squares - sum_of_squares_low)/sum_of_squares_low  #   Argument of an assumed Gaussian

        sdev_sigma1  =  abs(sigma_droplet - sigma_droplet_best) + delta_sigma

    #   print 'ratio_of_squares, sigma_droplet, sigma_droplet_best, sdev_sigma1: ', ratio_of_squares, sigma_droplet, sigma_droplet_best, sdev_sigma1


        if (ratio_of_squares < 1.0 and dev_flag == 0):           #   Assume ratio_of_squares is Gaussian distributed
            sdev_sigma  =  abs(sigma_droplet - sigma_droplet_best) + delta_sigma
            dev_flag = 1

    #
    #   Now estimate standard deviation of xi, assuming sigma is best value

    x_droplet   =   zeros(droplet_bins)
    y_droplet   =   zeros(droplet_bins)

    sdev_xi     = 0.0
    delta_xi    =   0.1
    xi_range    =   30

    xi_droplet  =   xi_droplet_best + 0.5*float(xi_range)*delta_xi
    dev_flag = 0

    for idrop in range(0,xi_range):
        xi_droplet -= delta_xi
        j=-1
        diff = zeros(droplet_bins)
        for i in range(0,nonzero_bins):     # Cycle through all the droplet points
            if x[i] >= x1:
                j+=1
                x_droplet[j] = x[i]
                exp_factor = 10.0**(1.5*sigma_droplet_best*(x[i] - xi_droplet))
                y_droplet[j] = cept + slope*x[i] - exp_factor
                diff[j] = y_droplet[j] - y[i]

    #
        sum_of_squares      = sum(diff[:]*diff[:])
        ratio_of_squares    = (sum_of_squares - sum_of_squares_low)/sum_of_squares_low  #   Argument of an assumed Gaussian

        if (ratio_of_squares < 1.0 and dev_flag == 0):       #   Assume ratio_of_squares is Gaussian distributed
            sdev_xi  =  abs(xi_droplet - xi_droplet_best) + delta_xi
            dev_flag = 1

    return (sdev_xi, sdev_sigma)   

    #   Usage in calling program:  slope, cept, errs, errc, s2 = linfit(x,y)

    #.................................................................

def storeDroplet(x_droplet,y_droplet,droplet_bins, droplet_store_flag):

    x_droplet_last  =   zeros(droplet_bins)
    y_droplet_last  =   zeros(droplet_bins)

    if droplet_store_flag == 'w':
        working_file = open("Droplet_File.txt", "w")
        for i in range(0,droplet_bins):
            working_file.write("%s\t%s\n" % (x_droplet[i], y_droplet[i]))
        working_file.close()

    if droplet_store_flag == 'r':
        working_file = open("Droplet_File.txt", "r")

        j_drop = 0
        for line in working_file:
            j_drop +=  1
        working_file.close()  # Put the file pointer back at top

        x_droplet_last = zeros(j_drop)
        y_droplet_last = zeros(j_drop)

        working_file = open("Droplet_File.txt", "r")

        i = -1
        for line in working_file:
            i += 1
            items = line.strip().split()
            print i, items[0],items[1]
            x_droplet_last[i] = float(items[0])
            y_droplet_last[i] = float(items[1])

        working_file.close()
        

    return (x_droplet_last, y_droplet_last)

    #.................................................................

def weibullFit(n, bins, cum_prob, mean_eqs):

    print cum_prob[:]

#    mean_eqs = mean_eqs * 1.5

#
#     This program fits the earthquake cycle statistics to a Weibull model with parameters tau, beta
#       and returns the best fit values with the standard deviations
#
#     Definitions:
#
#         bins[i]:  Abscissas of all the bins over which the CDF is defined
#     cum_prob[i]:  Ordinates of the CDF probabilities, from 0% to 100%
#
#
    sum_of_squares_low  =   100000000.0     #   Some large number

    cum_weibull   =   np.zeros(len(bins)+1)

    delta_beta     =   0.01
    beta_range     =   300

    #   Weibull must have the same mean as the data.  So this constrains either tau or beta.  
    #       We use the mean to set tau, then Estimate standard deviation of beta.  

    #   First, find the best-fitting Weibull

    beta = 0.1

    print 'number of bins: ', len(bins)

    for j in range(0,beta_range):
        beta += delta_beta
        tau = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta)

        diff = np.zeros(len(bins)+1)

        for i in range(1,len(bins)):
            cum_weibull[i] = 1.0 - math.exp(- ((bins[i]/tau)**beta))
            diff[i] = cum_weibull[i] - 0.01*cum_prob[i]

#            print i, tau, beta, bins[i], cum_weibull[i], cum_prob[i], diff[i]

        sum_of_squares      = sum(diff[:]*diff[:])

        if sum_of_squares < sum_of_squares_low:
            sum_of_squares_low = sum_of_squares
            beta_best = beta


    print 'Low sum of squares: ', sum_of_squares_low
    print
    tau_best = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta_best)

    #   Second, find the standard deviations

    dev_flag = 0
    beta = 0.1

    for j in range(0,beta_range):
        beta += delta_beta
        tau = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta)

        diff = np.zeros(len(bins)+1)

        for i in range(1,len(bins)):
            cum_weibull[i] = 1.0 - math.exp(-((bins[i]/tau)**beta))
            diff[i] = cum_weibull[i] - 0.01*cum_prob[i]

        sum_of_squares      = sum(diff[:]*diff[:])
        ratio_of_squares    = (sum_of_squares - sum_of_squares_low)/sum_of_squares_low  #   Argument of an assumed Gaussian

        if (ratio_of_squares < 1.0 and dev_flag == 0):           #   Assume ratio_of_squares is Gaussian distributed
            sdev_beta  =  abs(beta - beta_best) + delta_beta
            dev_flag = 1

    value = (beta_best+sdev_beta)
    tau1 = mean_eqs/scipy.special.gamma(1.0 + 1.0/value)

    value = (beta_best-sdev_beta)
    tau2 = mean_eqs/scipy.special.gamma(1.0 + 1.0/value)

    sdev_tau = 0.5*abs(tau-tau1) + 0.5*abs(tau-tau2)

    return (tau_best, beta_best, sdev_tau, sdev_beta)   

    #.................................................................

def variableWeibullFit(bins, cum_prob, mean_eqs):

#
#     This program fits the earthquake cycle statistics to a Weibull model with parameters tau, beta
#       and returns the best fit values with the standard deviations
#
#     Definitions:
#
#         bins[i]:  Abscissas of all the bins over which the CDF is defined
#     cum_prob[i]:  Ordinates of the CDF probabilities, from 0% to 100%
#
#
    sum_of_squares_low  =   100000000.0     #   Some large number

    cum_weibull   =   np.zeros(len(bins)+1)

    delta_beta     =   0.01
    beta_range     =   300

    #   Weibull must have the same mean as the data.  So this constrains either tau or beta.  
    #       We use the mean to set tau, then Estimate standard deviation of beta.  

    #   -----------------------------------------

    #   First, find the best-fitting Weibull

    beta = 0.1

    for j in range(0,beta_range):
        beta += delta_beta
        tau = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta)

        diff = np.zeros(len(bins)+1)

        for i in range(1,len(bins)):
            cum_weibull[i] = 1.0 - math.exp(- ((bins[i]/tau)**beta))
            diff[i] = cum_weibull[i] - 0.01*cum_prob[i]

#            print i, tau, beta, bins[i], cum_weibull[i], cum_prob[i], diff[i]

        sum_of_squares      = sum(diff[:]*diff[:])

        if sum_of_squares < sum_of_squares_low:
            sum_of_squares_low = sum_of_squares
            beta_best = beta

        tau_best = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta_best)

    #   -----------------------------------------

    #   Second, find the standard deviations

    dev_flag = 0
    beta = 0.00

    for j in range(0,beta_range):
        beta += delta_beta
        tau = mean_eqs/scipy.special.gamma(1.0 + 1.0/beta)

        diff = np.zeros(len(bins)+1)

        for i in range(1,len(bins)):
            cum_weibull[i] = 1.0 - math.exp(-((bins[i]/tau)**beta))
            diff[i] = cum_weibull[i] - 0.01*cum_prob[i]

        sum_of_squares      = sum(diff[:]*diff[:])
        ratio_of_squares    = (sum_of_squares - sum_of_squares_low)/sum_of_squares_low  #   Argument of an assumed Gaussian

        if (ratio_of_squares < 1.0 and dev_flag == 0):           #   Assume ratio_of_squares is Gaussian distributed
            sdev_beta  =  abs(beta - beta_best) + delta_beta
            dev_flag = 1

    value = (beta_best+sdev_beta)
    tau1 = mean_eqs/scipy.special.gamma(1.0 + 1.0/value)

    value = (beta_best-sdev_beta)
    tau2 = mean_eqs/scipy.special.gamma(1.0 + 1.0/value)

    sdev_tau = 0.5*abs(tau-tau1) + 0.5*abs(tau-tau2)

    #   .......................................................................

    return (tau_best, beta_best, sdev_tau, sdev_beta)   

    #.................................................................




def change_in_latitude(km):

    #   http://www.johndcook.com/blog/2009/04/27/converting-miles-to-degrees-longitude-or-latitude/

    # Distances are measured in km.
    # Longitudes and latitudes are measured in degrees.
    # Earth is assumed to be perfectly spherical.

    earth_radius = 6371.0
    degrees_to_radians = math.pi/180.0
    radians_to_degrees = 180.0/math.pi

    #Given a distance north, return the change in latitude in degrees.

    return (km/earth_radius)*radians_to_degrees

    #.................................................................

def change_in_longitude(latitude, km):

    #   http://www.johndcook.com/blog/2009/04/27/converting-miles-to-degrees-longitude-or-latitude/

    # Distances are measured in km.
    # Longitudes and latitudes are measured in degrees.
    # Earth is assumed to be perfectly spherical.

    earth_radius = 6371.0
    degrees_to_radians = math.pi/180.0
    radians_to_degrees = 180.0/math.pi

    # Given a latitude and a distance west, return the change in longitude in degrees.
    # Find the radius of a circle around the earth at given latitude.

    r = earth_radius*math.cos(latitude*degrees_to_radians)
    return (km/r)*radians_to_degrees

def draw_big_circle(CircleCenterLat, CircleCenterLng, CircleRadius):

    # Build an array of x-y values in degrees defining a circle of the required radius

    twopi = 6.283185306

    number_points   = 100
    delta_angle     = twopi/float(number_points)

    x_circle_dg      = np.zeros(number_points)
    y_circle_dg      = np.zeros(number_points)

    angle = 0.0     #   start at due north

    for i in range(0,number_points):
        x_comp = CircleRadius * math.sin(angle) 
        y_comp = CircleRadius * math.cos(angle) 

        y_circle_dg[i] = CircleCenterLat + change_in_latitude(y_comp)
        x_circle_dg[i] = CircleCenterLng + change_in_longitude(CircleCenterLat, x_comp)

        angle += delta_angle

        x_circle_dg =   np.append(x_circle_dg, x_circle_dg[0])
        y_circle_dg =   np.append(y_circle_dg, y_circle_dg[0])


    return (x_circle_dg, y_circle_dg)
        
        
        











