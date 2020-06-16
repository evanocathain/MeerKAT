
#!/usr/bin/python

import sys
from datetime import *
import ephem

if len(sys.argv)!=4:
    print ("Usage: riseset site ra dec")
    print (' ')
    print ("This program returns the rise and set LSTs for the provided RA and DEC at the requested site")
    print ("e.g. riseset mk 12:34:56 -43:21:00")
    print ("Site = MeerKAT")
    print ("Long,Lat = 21:24:40.0 -30:43:16.0")
    print ("Target Coordinates = 12:34:56.00 -43:21:00.0")
    print ("Rise Time (LST) = 6:01:30.97")
    print ("Set Time (LST)  = 19:10:28.11")
    sys.exit(0)

place=sys.argv[1]
ra   =sys.argv[2]
dec  =sys.argv[3]

site=ephem.Observer()

# NB. I got all coordinates from wikipedia. Sorry if they are incorrect!
if (place=="jbo" or place=="jodrell"):
    site.long = '-2.30715'
    site.lat = '53.23700'
    site.elevation = 0
    place = "Jodrell Bank"
elif (place=="birr" or place=="Birr" or place=="I-LOFAR" or place=="ilofar"):
    site.long = '-7.9133'
    site.lat = '53.0914'
    site.elevation = 0
    place = "Birr"
elif (place=="eff" or place=="effelsberg"):
    site.long = '6:52:58'
    site.lat = '50:31:29'
    site.elevation = 0
    place = "Effelsberg"
elif (place=="srt" or place=="sardinia"):
    site.long = '9:14:43'
    site.lat = '39:29:35'
    site.elevation = 0
    place = "Sardinia Radio Telescope"
elif (place=="wsrt" or place=="westerbork"):
    site.long = '6:36:12'
    site.lat = '52:54:53'
    site.elevation = 0
    place = "Westerbork Synthesis Radio Telescope"
elif (place=="nancay"):
    site.long = '2:12'
    site.lat = '47:23'
    site.elevation = 0
    place = "Nancay Radio Telescope"
elif (place=="gbt"):
    site.long = '-79:50:23'
    site.lat = '38:25:59'
    site.elevation = 0
    place = "Green Bank Telescope"
elif (place=="pks" or place=="parkes"):
    site.long = '148:15:46.51'
    site.lat = '-32:59:52.01'
    site.elevation = 0
    site.horizon = '30.25'
    place = "Parkes Observatory"
elif (place=="arecibo"):
    site.long = '-66:45:10'
    site.lat = '18:20:39'
    site.elevation = 0
    place = "Arecibo"
elif (place=="vla"):
    site.long = '-107.61835'
    site.lat = '34.078967'
    site.elevation = 0
    place = "Very Large Array"
elif (place=="gmrt"):
    site.long = '74:02:59.07'
    site.lat = '19:05:47.46'
    site.elevation = 0
    place = "The Giant Metrewave Radio Telescope"
elif (place=="lofar"):
    site.long = '6:52:08.18'
    site.lat = '52:54:31.55'
    site.elevation = 0
    place = "LOFAR Superterp"
elif (place=="mk" or place=="MK" or place=="meerkat" or place=="MeerKAT"):
    site.long = '21:24:40'
    site.lat = '-30:43:16'
    site.elevation = 0
    site.horizon='15.0' # Took this from SKA1 Level 0 Science Requirements. I guess it is this or better.
    place = "MeerKAT"
else:
    print ("Which telescope?! Try again!")
    sys.exit(0)

target      = ephem.FixedBody()
target._ra  = ephem.hours(ra)
target._dec = ephem.degrees(dec)
#print target._ra
#print target._dec
target.compute(site)

site.date=datetime.utcnow()
print ("Site =",place)
print ("Long,Lat =",site.long,site.lat)
#print "Current LST at Site =",site.sidereal_time()
print ("Target Coordinates =", target._ra, target._dec)
site.date = site.next_rising(target)
print ("Rise Time (LST) =", site.sidereal_time())
site.date = site.next_setting(target)
print ("Set Time (LST) =", site.sidereal_time())
#site.next_rising(target))
#print "Set Time =",site.sidereal_time(site.next_setting(target))
#print "UTC =",site.date
#print "UTC (MJD) =",site.date+15019.5
## NB pyephem records Dublin Julian Day Number not MJD
## see Table in en.wikipedia.org/wiki/Julian_day for more info

sys.exit(0)
