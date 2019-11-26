import pytz
import tzlocal
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from astropy.time import Time
from skyfield.api import Topos, load, Loader


import dash
import dash_core_components as dcc
import dash_html_components as html

# instantiate the dash 
app = dash.Dash()

# Can use external  CSS style sheets can be used to help with the formatting. For example Bootstrap:
app.css.append_css({'external_url':'https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css'})
app.css.append_css({'external_url':'https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.0.0-beta/css/bootstrap.css'})
app.css.config.serve_locally = False


# Now the Install instructions
app.layout = html.Div([
                # OUTERMOST DIV
                html.Div([
                    html.H1(children='Track the International Space Station',
                            className='mx-auto  text-white'),
                                     ],
                    style={'background-color': '#333333'},
                    className='row rounded mx-auto  text-white'
                ),
                # END OUTERMOST DIV
            ])

stations_url = 'http://celestrak.com/NORAD/elements/stations.txt'
load = Loader('./data')
satellites = load.tle(stations_url)
satellite = satellites['ISS (ZARYA)']

dt, tm = str(datetime.utcnow()).split()
year, month, day = dt.split('-')
hour, minute, second = tm.split(':')
sec, _ = second.split('.')
    
ts = load.timescale(builtin=True)
t = ts.utc(int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(sec)
        )

# If ISS TLE is more than 7 days old, reload it
days = t - satellite.epoch

if abs(days) > 7:
    satellites = load.tle(stations_url, reload=True)
    satellite = satellites['ISS (ZARYA)']

# Time array at which to calculate the position of ISS
# Cadence set to 20 seconds, currently
t_steps = np.arange(0, (60 * 60 * 24), 20)

time_range = ts.utc(int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        t_steps
        )

# Lat Lon of Coffs Harbour
lat_lon = Topos('-30.2986 N', '153.1094 E')  

orbit = (satellite - lat_lon).at(time_range)
alt, az, _ = orbit.altaz()

# Check if sat is above the horizon, return boolean array
above_horizon = alt.degrees >= 0

# Indicies of rare times that sats are above the horizon
indicies, = above_horizon.nonzero()

# Boundary times at which the sat either rises or sets
boundaries, = np.diff(above_horizon).nonzero()

if above_horizon[0] == True:
    boundaries = [indicies[0]]+ list(boundaries)
    boundaries = np.asarray(boundaries)

if above_horizon[-1] == True:
    boundaries = list(boundaries) + [indicies[-1]]
    boundaries = np.asarray(boundaries)

# Reshape into pairs rise & set indicies
passes = boundaries.reshape(len(boundaries) // 2, 2)


def ephem_data(t_arr, pass_index, alt, az):
    '''Satellite Ephemeris Data.
    
    creates rise time, set times and alt, az arrays.
    
    Args:
        t_arr: Skyfield time array object
        pass_index: One pair of sat indicies from passed 2D array
        alt: Array of altitudes of sat at t_arr times
        az: Array of azimuths of sat at t_arr times

    Returns:
        t_rise: Rise time of sat in gps seconds (Astropy)
        t_set: Set time of sat in gps seconds (Astropy)
        sat_alt: Altitude of satellite while it is above the horizon. Array.
        sat_az: Azimuth of satellite while it is above the horizon. Array.
    '''

    i, j = pass_index
    t_rise = Time(t_arr[i].tt, format='jd').strftime("%Y-%m-%d %H:%M:%S")
    t_set = Time(t_arr[j].tt, format='jd').strftime("%Y-%m-%d %H:%M:%S")


    theta = list(az.radians)
    r = list(alt.degrees)

    sat_az = theta[i:j+1]
    sat_alt = r[i:j+1]
    
    return (t_rise, t_set, sat_alt, sat_az)



def sat_plot(alt, az, start_t, stop_t):
    '''Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    '''
    
    # Set up the polar plot.
    plt.style.use('seaborn')
    figure = plt.figure(figsize=(6,6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,10,20,30,40,50,60,70,80,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title(f"ISS rise {utc2local(start_t)}, set {utc2local(stop_t)}")
    plt.tight_layout()
    plt.plot(az, alt, '-', linewidth=1.8, alpha=0.7, color='#458323')
    plt.show()


sat_ephem = {}
sat_ephem['t_rise'] = []
sat_ephem['t_set'] = []
sat_ephem['sat_alt'] = []
sat_ephem['sat_az'] = []
    
for pass_index in passes:
    t_rise, t_set, sat_alt, sat_az = ephem_data(time_range, pass_index, alt, az)

    sat_ephem['t_rise'].append(t_rise)
    sat_ephem['t_set'].append(t_set)
    sat_ephem['sat_alt'].append(sat_alt)
    sat_ephem['sat_az'].append(sat_az)

utc_time = sat_ephem['t_rise'][0]



def utc2local(utc):
    local_timezone = tzlocal.get_localzone()
    utc_time = datetime.strptime(utc, "%Y-%m-%d %H:%M:%S")
    local_time,_ = str(utc_time.replace(tzinfo=pytz.utc).astimezone(local_timezone)).split('+')
    _, tm = local_time.split()
    return tm

#plot the next four ISS passes
#for i in range(4):
#    sat_plot(sat_ephem['sat_alt'][i], sat_ephem['sat_az'][i], sat_ephem['t_rise'][i], sat_ephem['t_set'][i])

if __name__ == '__main__':
    app.run_server(port=8080)

