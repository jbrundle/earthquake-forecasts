#   This file is a brief explanation of the codes for plotting earthquake data using Python
#
#   ..........................................................................
#

    Note that to use these codes, you need to have Numpy, SciPy, Matplotlib, Basemap and Python Image 
    Library (PIL) installed as part of your Python installation.

1.  EQMenu.py: is the master program and is the one that should be run from
    the command line

    Syntax:  >> python EQMenu.py

2.  EQMethods.py: is the file with the primary functions contained in it.  It
    is called by the master program.

2.  python_generate_catalog.py is the code that calls a function in EQMethods
    to generate the working data file

4.  plot_ANSS_epicenters.py: is a code that is called by EQMenu.py.  Its purpose
    is to call methods in EQMethods.py that plot earthquake epicenters on a map.

5.  plot_ANSS_seismicity.py: is a code that is called by EQMenu.py.  Its purpose
    is to call methods in EQMethods.py that plot the magnitude vs. time of large 
    earthquakes in the selected region.

6.  contour_eq_probs.py: is a code that uses the GetEarthquakeProbability API
    from www.openhazards.com/data to download and contour earthquake probabilities
    in the selected region 

7.  locations.txt:  is a csv file containing pre-defined map boundaries with a label.
    The default is California.  Each line of the file has the format:

    [ Location Name (string), S-most Lat (float), N-most Lat (float), W-most Lng (float), E-most Lng (float) ]

    To enter more pre-defined locations, edit this file with a text editor.

    If this file does not exist, you must create it with at least 1 line of data.

    *** Note that  an http error 500 will occur if values are entered that include a Lng or Lat of 0.0 ***
    *** Also, there should be no empty lines at the bottom of the location file.  If there are, you will
        see an error: "IndexError: index 1 is out of bounds for axis 0 with size 1"

PYTHON RELATED PACKAGES INSTALLED ON MY MAC VIA MACPORTS:

python py27:                python27 @2.7.9 (lang) An interpreted, object-oriented programming language
py27-matplotlib:            py27-matplotlib @1.4.2 (python, graphics, math) Matplotlib is a python plotting library
py27-matplotlib-basemap:    py27-matplotlib-basemap @1.0.7_1 (python, graphics, math) Matplotlib toolkit for plotting data on map projections
py27-pil:                   py27-pil @1.1.7_7 (python, graphics)  Python Imaging Library

Note:  A bug in python sometimes leads to the appearance of files with names such as .goutputstream*  These can be removed using the command sudo rm .goutputstream* -v

CHANGE LOG:  

New Features in Version 1.3:

- Changed the workflow so that the user generates the working data file, and can then edit it
    in the default text editor or other editor.  This allows the user to plot subsets of the complete data file.

New Features in Version 1.2:

- Added some additional error handling with input values
- Forecast time window changed from fixed 1 year interval to arbitrary future interval (user choice)
    in contour_eq_probs.py
- Added a feature to allow the user to inspect the earthquake data stream
- Added the capability of plotting on eTopo and ShadedRelief maps (py27-pil package must be 
    installed.  pil = python image library)

New Features in Version 1.1:

- Changed master code from python_menu.py to EQMenu.py
- Better method of handling pre-defined locations through the use of a location file,
    location.txt
- Renamed some of the calling codes
- More choices in parameters for forecast method, contour_eq_probs.py
- Reconfigured epicenter map so that legend is outside the map area
- Better annotations on figures
- More logical arrangement of choices in map parameters in EQMenu.py









