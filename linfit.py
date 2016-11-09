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
#     n:     The number of data 
#
    n = int(n)
    ata = [zeros(2),zeros(2)]
    aty = [zeros(2),zeros(2)]
    atainv = [zeros(2),zeros(2)]
#    
    sumx2 = 0.
    xsum = 0.
    yxsum = 0.
    ysum = 0.

    for i in range(0,n-1):
        sumx2 = sumx2 + x[i] * x[i]
        xsum = xsum + x[i]
        yxsum = yxsum + y[i] * x[i]
        ysum = ysum + y[i]
#
    ata[0][0] = sumx2
    ata[0][1] = xsum
    ata[1][0] = xsum
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
    for i in range(0,n-1):
        s2 = s2 + (y[i] - cept - slope * x[i])**2

    s2 = s2 / (float(n) - 2.)
      
    errs = math.sqrt( float(n) * s2 / det )
    errc = math.sqrt( s2 * sumx2 / det)

    return (slope, cept, errs, errc, s2)

#   Usage in calling program (I think):  slope, cept, errs, errc, s2 = linfit(x,y,n)

