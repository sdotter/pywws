#!/usr/bin/env python

# pywws - Python software for USB Wireless Weather Stations
# http://github.com/jim-easterbrook/pywws
# Copyright (C) 2008-18  pywws contributors

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""conversions.py - a set of functions to convert pywws native units
(Centigrade, mm, m/s, hPa) to other popular units

"""

from __future__ import absolute_import, print_function

__docformat__ = "restructuredtext en"

import math

import pywws.localisation
import pywws.process


def scale(value, factor):
    """Multiply value by factor, allowing for None values."""
    if value is None:
        return None
    elif isinstance(value, (list, tuple)):
        return [i * factor for i in value]
    elif hasattr(value, '__iter__'): # check for iterable objects like numpy arrays or pandas series
        return type(value)(i * factor for i in value)
    else:
        return value * factor

def illuminance_wm2(lux):
    "Approximate conversion of illuminance in lux to solar radiation in W/m2"
    return scale(lux, 0.005)

def wm2_illuminance(wm2):
    """Approximate conversion of solar radiation in W/m2 to illuminance in lux"""
    return scale(wm2, 1 / 0.005)

def mph_to_kph(mph):
    ''' Convert speed from miles per hour (mph) to kilometers per hour (kph).
    '''
    return mph * 1.60934
    
def f_to_c(fahrenheit):
    ''' Convert temperature from Fahrenheit to Celsius.
    '''
    return round(((fahrenheit - 32) * 5 / 9), 2)

def inHg_to_hPa(inHg):
    ''' Convert pressure from inches of mercury (inHg) to hectopascals (hPa).
    '''
    return inHg * 33.8639

def get_dew_point_c(t_air_c, rel_humidity):
    """Compute the dew point in degrees Celsius
    :param t_air_c: current ambient temperature in degrees Celsius
    :type t_air_c: float
    :param rel_humidity: relative humidity in %
    :type rel_humidity: float
    :return: the dew point in degrees Celsius
    :rtype: float
    """
    A = 17.27
    B = 237.7
    alpha = ((A * t_air_c) / (B + t_air_c)) + math.log(rel_humidity/100.0)
    return (B * alpha) / (A - alpha)

def wind_chill(temperature, wind_speed):
    """
    Calculate the Wind Chill (feels like temperature) based on NOAA.
    Wind-chill or windchill (popularly wind chill factor) is the lowering of
    body temperature due to the passing-flow of lower-temperature air.
    Wind chill numbers are always lower than the air temperature for values
    where the formula is valid. When the apparent temperature is higher than 
    the air temperature, the heat index is used instead.
    Wind Chill Temperature is only defined for temperatures at or below
    50 F and wind speeds above 3 mph. (10C, 4.8 Km/h)

    3 Mph = 4.828 [Km/h] = 1.34 [m/s]
    50F  = (50 - 32) * 5/9 = 10C
    See:
    [1] https://en.wikipedia.org/wiki/Wind_chill
    [2] https://www.wpc.ncep.noaa.gov/html/windchill.shtml
    """

    T = temperature             # Celsius
    V = wind_speed              # Kilometer per hour

    # We should never get here...
    if T > 10 or V <= 4.8:      # if T > 50 or V <= 3:    # (F, Mph)
        e = "Wind Chill Temperature is only defined for temperatures at or below 10C and wind speeds above 4.8 Km/h."
        raise ValueError(e)

    #----------------------------------
    # WC = Wind Chill [2]
    #----------------------------------
    # (Farenheit, Mph)
    # WC = 35.74 + (0.6215 * T) - 35.75 * V**0.16 + 0.4275 * T * V**0.16
    # (Celsius, Kph)
    WC = 13.12 + (0.6215 * T) - 11.37 * V**0.16 + 0.3965 * T * V**0.16    
    
    return WC

def heat_index(temperature, humidity):
    """
    Calculate the Heat Index (feels like temperature) based on the NOAA equation.
    
    The Heat Index (HI), or humiture, or "feels like temperature", is an index that combines
    air-temperature and relative humidity in an attempt to determine the human-perceived 
    equivalent temperature.
    HI is only useful when (T > 27 C) & (RH > 40 %)
    Tf = (Tc * 9/5) + 32 
    Tc = (Tf - 32) * 5/9 
    
    See: 
    [1] https://en.wikipedia.org/wiki/Heat_index
    [2] http://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
    [3] https://github.com/geanders/weathermetrics/blob/master/R/heat_index.R
    [4] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3801457/
    """

    T = temperature     # Celsius
    H = humidity        # Relative Humidity

    # SI units (Celsius)
    c1 = -8.78469475556
    c2 = 1.61139411
    c3 = 2.33854883889
    c4 = -0.14611605
    c5 = -0.012308094
    c6 = -0.0164248277778
    c7 = 0.002211732
    c8 = 0.00072546
    c9 = -0.000003582

    #----------------------------------------------------------------
    # Maybe consider [4] in "Formula-18" by [Costanzo et al. 2006]
    # HIc   = Tc - 0.55 * (1 - 0.001 H)(Tc - 14.5)
    #       = Tc - 0.55 * (Tc - 14.5 - 0.001*H*Tc + 0.001*14.5*H 
    #----------------------------------------------------------------

    Tf = (T * 9/5) + 32 
    # Try the simplified formula from [2] first (used for HI < 80)
    HIf = 0.5 * (Tf + 61.0 + (Tf - 68.0) * 1.2 + H * 0.094)
    Tavg = (HIf + Tf)/2   # Instructions in [3] call for averaging

    if Tavg >= 80:      # [F]
        # IF (T > 27C) & (H > 40 %):
        # Use the full Rothfusz regression formula (now in Celsius)
        HI = math.fsum([
            c1,
            c2 * T,
            c3 * H,
            c4 * T * H,
            c5 * T**2,
            c6 * H**2,
            c7 * T**2 * H,
            c8 * T * H**2,
            c9 * T**2 * H**2,
        ])
    else:
        HI = (HIf - 32) * 5/9

    return HI

def feels_like(temperature, humidity, wind_speed):
    """
    Calculate the "Feels Like" temperature based on NOAA.
    Logic:
    * Wind Chill:   temperature <= 50 F and wind > 3 mph
    * Heat Index:   temperature >= 80 F
    * Temperature as is: all other cases
    
    -----------------------------------------------------
    50F  = (50F - 32) * 5/9 = 10C
    80F  = (80F - 32) * 5/9 = 26.7 C
    3 Mph = 4.828 [Km/h] = 1.34 [m/s]
    1 [m/s] =   1/[1000 m/Km] * [3600 s/hour] = 3.6 [Km/h]
    -----------------------------------------------------
    """

    T = temperature     # [C]

    #----------------------------------
    # FL = Feels Like
    #----------------------------------    
    if T <= 10 and wind_speed > 4.8:
        # Wind Chill for low temp cases (and wind)
        FL = wind_chill(T, wind_speed)
    elif T >= 26.7:
        # Heat Index for High temp cases
        FL = heat_index(T, humidity)
    else:
        FL = T

    return round(FL, 1)

def pressure_inhg(hPa):
    "Convert pressure from hectopascals/millibar to inches of mercury"
    return scale(hPa, 1 / 33.86389)

def pressure_trend_text(trend):
    """Convert pressure trend to a string, as used by the UK met
    office.

    """
    _ = pywws.localisation.translation.ugettext
    if trend > 6.0:
        return _(u'rising very rapidly')
    elif trend > 3.5:
        return _(u'rising quickly')
    elif trend > 1.5:
        return _(u'rising')
    elif trend >= 0.1:
        return _(u'rising slowly')
    elif trend < -6.0:
        return _(u'falling very rapidly')
    elif trend < -3.5:
        return _(u'falling quickly')
    elif trend < -1.5:
        return _(u'falling')
    elif trend <= -0.1:
        return _(u'falling slowly')
    return _(u'steady')

def inches_to_mm(inches):
    """
    Converts rainfall from inches to millimeters.
    
    Parameters:
    inches (float or str): Rainfall in inches. Strings will be parsed if possible.
    
    Returns:
    float: Rainfall in millimeters.
    """
    try:
        # Probeer de invoer om te zetten naar een float
        inches = float(inches)
    except ValueError:
        raise ValueError("Input must be a number or a string representing a number")
    
    return inches * 25.4

def rain_inch(mm):
    "Convert rainfall from millimetres to inches"
    return scale(mm, 1 / 25.4)

def temp_f(c):
    "Convert temperature from Celsius to Fahrenheit"
    if c is None:
        return None
    return (c * 9.0 / 5.0) + 32.0

def degrees_to_wind_direction(degrees):
    """
    Converts a wind direction in degrees to a wind direction name.
    
    Parameters:
    degrees (float): The wind direction in degrees (0-360).
    
    Returns:
    str: The corresponding wind direction name.
    """
    directions = [
        "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
        "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
    ]
    index = int(round(degrees / 22.5)) % 16  # Ensure the index is an integer
    return directions[index]

def winddir_average(data, threshold, min_count, decay=1.0):
    """Compute average wind direction (in degrees) for a slice of data.

    The wind speed and direction of each data item is converted to a
    vector before averaging, so the result reflects the dominant wind
    direction during the time period covered by the data.

    Setting the ``decay`` parameter converts the filter from a simple
    averager to one where the most recent sample carries the highest
    weight, and earlier samples have a lower weight according to how
    long ago they were.

    This process is an approximation of "exponential smoothing". See
    `Wikipedia <http://en.wikipedia.org/wiki/Exponential_smoothing>`_
    for a detailed discussion.

    The parameter ``decay`` corresponds to the value ``(1 - alpha)``
    in the Wikipedia description. Because the weather data being
    smoothed may not be at regular intervals this parameter is the
    decay over 5 minutes. Weather data at other intervals will have
    its weight scaled accordingly.

    :note: The return value is in degrees, not the 0..15 range used
        elsewhere in pywws.

    :param data: a slice of pywws raw/calib or hourly data.

    :type data: pywws.storage.CoreStore

    :param threshold: minimum average windspeed for there to be a
        valid wind direction.

    :type threshold: float

    :param min_count: minimum number of data items for there to be a
        valid wind direction.

    :type min_count: int

    :param decay: filter coefficient decay rate.

    :type decay: float

    :rtype: float
    
    """
    wind_filter = pywws.process.WindFilter()
    count = 0
    for item in data:
        wind_filter.add(item)
        if item['wind_dir'] is not None:
            count += 1
    if count < min_count:
        return None
    speed, direction = wind_filter.result()
    if speed is None or speed < threshold:
        return None
    return direction * 22.5
    
def winddir_degrees(pts):
    "Convert wind direction from 0..15 to degrees"
    return scale(pts, 22.5)

_winddir_text_array = None

def winddir_text(pts):
    "Convert wind direction from 0..15 to compass point text"
    global _winddir_text_array
    if pts is None:
        return None
    if not isinstance(pts, int):
        pts = int(pts + 0.5) % 16
    if not _winddir_text_array:
        _ = pywws.localisation.translation.ugettext
        _winddir_text_array = (
            _(u'N'), _(u'NNE'), _(u'NE'), _(u'ENE'),
            _(u'E'), _(u'ESE'), _(u'SE'), _(u'SSE'),
            _(u'S'), _(u'SSW'), _(u'SW'), _(u'WSW'),
            _(u'W'), _(u'WNW'), _(u'NW'), _(u'NNW'),
            )
    return _winddir_text_array[pts]

def wind_kmph(ms):
    "Convert wind from metres per second to kilometres per hour"
    return scale(ms, 3.6)

def wind_mph(ms):
    "Convert wind from metres per second to miles per hour"
    return scale(ms, 3.6 / 1.609344)

def wind_kn(ms):
    "Convert wind from metres per second to knots"
    return scale(ms, 3.6 / 1.852)

_bft_threshold = (
    0.3, 1.5, 3.4, 5.4, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4, 32.6)

def wind_bft(ms):
    "Convert wind from metres per second to Beaufort scale"
    if ms is None:
        return None
    for bft in range(len(_bft_threshold)):
        if ms < _bft_threshold[bft]:
            return bft
    return len(_bft_threshold)

def dew_point(temp, hum):
    """Compute dew point, using formula from
    http://en.wikipedia.org/wiki/Dew_point.

    """
    if temp is None or hum is None:
        return None
    a = 17.27
    b = 237.7
    gamma = ((a * temp) / (b + temp)) + math.log(float(hum) / 100.0)
    return (b * gamma) / (a - gamma)

def cadhumidex(temp, humidity):
    "Calculate Humidity Index as per Canadian Weather Standards"
    if temp is None or humidity is None:
        return None
    # Formulas are adapted to not use e^(...) with no appreciable
    # change in accuracy (0.0227%)
    saturation_pressure = (6.112 * (10.0**(7.5 * temp / (237.7 + temp))) *
                           float(humidity) / 100.0)
    return temp + (0.555 * (saturation_pressure - 10.0))

def usaheatindex(temp, humidity, dew=None):
    """Calculate Heat Index as per USA National Weather Service Standards

    See http://en.wikipedia.org/wiki/Heat_index, formula 1. The
    formula is not valid for T < 26.7C, Dew Point < 12C, or RH < 40%

    """
    if temp is None or humidity is None:
        return None
    if dew is None:
        dew = dew_point(temp, humidity)
    if temp < 26.7 or humidity < 40 or dew < 12.0:
        return temp
    T = (temp * 1.8) + 32.0
    R = humidity
    c_1 = -42.379
    c_2 = 2.04901523
    c_3 = 10.14333127
    c_4 = -0.22475541
    c_5 = -0.00683783
    c_6 = -0.05481717
    c_7 = 0.00122874
    c_8 = 0.00085282
    c_9 = -0.00000199
    return ((c_1 + (c_2 * T) + (c_3 * R) + (c_4 * T * R) + (c_5 * (T**2)) +
             (c_6 * (R**2)) + (c_7 * (T**2) * R) + (c_8 * T * (R**2)) +
             (c_9 * (T**2) * (R**2))) - 32.0) / 1.8

def wind_chill(temp, wind):
    """Compute wind chill, using formula from
    http://en.wikipedia.org/wiki/wind_chill

    """
    if temp is None or wind is None:
        return None
    wind_kph = wind * 3.6
    if wind_kph <= 4.8 or temp > 10.0:
        return temp
    return min(13.12 + (temp * 0.6215) +
               (((0.3965 * temp) - 11.37) * (wind_kph ** 0.16)),
               temp)

def apparent_temp(temp, rh, wind):
    """Compute apparent temperature (real feel), using formula from
    http://www.bom.gov.au/info/thermal_stress/

    """
    if temp is None or rh is None or wind is None:
        return None
    vap_press = (float(rh) / 100.0) * 6.105 * math.exp(
        17.27 * temp / (237.7 + temp))
    return temp + (0.33 * vap_press) - (0.70 * wind) - 4.00

def cloud_base(temp, hum):
    """Calculate cumulus cloud base in metres, using formula from
    https://en.wikipedia.org/wiki/Cloud_base or
    https://de.wikipedia.org/wiki/Kondensationsniveau#Konvektionskondensationsniveau
    """
    if temp is None or hum is None:
        return None
    dew_pt = dew_point(temp, hum)
    spread = float(temp) - dew_pt
    return spread * 125.0

def cloud_ft(m):
    "Convert cloud base from metres to feet."
    return scale(m, 3.28084)


def _main(argv=None):
    global _winddir_text_array
    # run some simple tests
    print('Wind speed:')
    print('%6s %8s %8s %8s %6s' % ('m/s', 'km/h', 'mph', 'knots', 'bft'))
    for ms in (0, 1, 2, 4, 6, 9, 12, 15, 18, 22, 26, 30, 34):
        print('%6g %8.3f %8.3f %8.3f %6d' % (
            ms, wind_kmph(ms), wind_mph(ms), wind_kn(ms), wind_bft(ms)))
    print('Wind direction:')
    for pts in range(16):
        print(' ' + winddir_text(pts), end='')
    print('')
    print('Wind direction, in Swedish:')
    pywws.localisation.set_translation('sv')
    _winddir_text_array = None
    for pts in range(16):
        print(' ' + winddir_text(pts), end='')
    print('')
    print('Cloud base in m and ft:')
    for hum in range(25, 75, 5):
        print("%8.3f m / %8.3f ft" % (cloud_base(15.0, hum), cloud_ft(cloud_base(15.0, hum))))
    print('')


if __name__ == "__main__":
    _main()
