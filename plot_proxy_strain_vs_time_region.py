#!/opt/local/bin python

##############################################################################

#         Code plots a rectangular map with the earthquakes larger than 5.5 on it

#         Usage:   python plot_proxy_stress_vs_time_region.py NELat NELng SWLat SWLng MagLo
#
#         Where:    Latitude in degrees
#                   Longitude in degrees
#                   NE and SW corners of the map rectangle must be specified
#                   MagLo is the lower limit for the event magnitudes (typically 5.5)

##############################################################################

#import sys
#sys.path.reverse()

import sys
import EQMethods

##############################################################################


def main(argv=None):

    NELat       = float(sys.argv[1])
    NELng       = float(sys.argv[2])
    SWLat       = float(sys.argv[3])
    SWLng       = float(sys.argv[4])
    MagLo       = float(sys.argv[5])
    Location    = str(sys.argv[6])  

    EQMethods.proxy_strain_region(NELat, NELng, SWLat, SWLng, MagLo, Location)

#
if __name__ == "__main__": 
	sys.exit(main())



