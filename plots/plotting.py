import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



def read_flux(filename):
    return(pd.read_table(filename, index_col=0))

def read_surfacedensity(filename):
    pass

def plot_inwardflux(data, time, axis=None):
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
    times = data.index.values
    closest_time_index = np.argmin(abs(times-time))
    print closest_time_index
    closest_time = times[closest_time_index]
    smas = map(float, data.columns.values)
    flux = data.iloc[closest_time_index].values
    axis.semilogx(smas, flux)
    axis.set_xlabel("Distance from star [AU]")
    axis.set_ylabel("Inward pebble flux [M$_{\oplus}$/yr]")
    return(fig)


def plot_accretion(data, time, axis=None):
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
    times = data.index.values
    closest_time_index = np.argmin(abs(times-time))
    print closest_time_index
    closest_time = times[closest_time_index]
    smas = map(float, data.columns.values)
    flux = data.iloc[closest_time_index].values
    accretion = np.hstack((flux[1::]-flux[0:-1], [0]))
    axis.semilogx(smas, accretion)
    axis.set_xlabel("Distance from star [AU]")
    axis.set_ylabel("d(inward flux)/dlog(orbital distance)")
    return(fig)


