import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import pylab as py

from astropy.io import fits
from scipy import signal
from scipy.stats import f
from scipy.interpolate import interp1d

#load in data
infile = sys.argv[1] #string: name of input kepler .fits file
outfile = sys.argv[2] #string: name of master output file for this planet.
t_orbit = float(sys.argv[3]) #period, days
t_14 = (float(sys.argv[4]))/24. #transit duration, hours converted to days
t_midpt = float(sys.argv[5]) #transit ephemeris, days

#other relevant timescales
t_longcadence = 0.020434 #average elapsed time between long-cadence data points, in days
t_shortcadence = 0.0006810758 #average elapsed time between short-cadence data points, in days

if "slc" in infile:
  cadence = t_shortcadence
elif "llc" in infile:
  cadence = t_longcadence

  #open fits file and extract data from HDU 1
  fitsdata = fits.open(infile)
  data = fitsdata[1].data

  t = data.field('TIME')
  sap_flux = data.field('SAP_FLUX')
  sap_flux_err = data.field('SAP_FLUX_ERR')
  bkg = data.field('SAP_BKG')
  bkg_err = data.field('SAP_BKG_ERR')

  pdc_flux = data.field('PDCSAP_FLUX')
  pdc_flux_err = data.field('PDCSAP_FLUX_ERR')

  notnan = (~np.isnan(t) & ~np.isnan(pdc_flux) & ~np.isnan(pdc_flux_err))
  t = t[notnan]
  pdc_flux = pdc_flux[notnan]
  pdc_flux_err = pdc_flux_err[notnan]

  #begin with window = orbital period
  window = t_orbit
  tstart = np.min(t) + 0.5*window
  tend = np.max(t)

  m = 0
  j = 1

  beta_ms = []
  moving_median_fluxes = []
  while ((tstart + (j*cadence) + window) < tend):
    #print j
    times_mask = (((tstart + (j-1)*cadence) <= t) & (t < (tstart + j*cadence + window)))
    selected_times = t[times_mask]
    selected_fluxes = pdc_flux[times_mask]

    if len(selected_times) > 0:
      m += 1
      beta_ms.append((tstart + (j-0.5)*cadence))
      moving_median_fluxes.append(np.median(selected_fluxes))

    j += 1

  beta_ms = np.array(beta_ms)               #array of (tstart + (j-0.5)*cadence) for j ranging from 1 to int(obs_baseline/cadence)
                                            #(indicates position of moving median filter along time array)
  moving_median_fluxes = np.array(moving_median_fluxes)   #array of fluxes after passing a moving median filter of width (P + cadence)

  tp = np.min(beta_ms)
  tq = np.max(beta_ms)

  times_mask = ((tp <= t) & (t <= tq))

  selected_times = t[times_mask]
  selected_fluxes = pdc_flux[times_mask]
  selected_fluxerrs = pdc_flux_err[times_mask]

  beta_ms_unique, beta_ms_unique_idxs = np.unique(beta_ms, return_index=True)

  #interpolate moving median fluxes between beta_ms
  beta_ms = beta_ms[beta_ms_unique_idxs] # "xs" for interpolation
  moving_median_fluxes = moving_median_fluxes[beta_ms_unique_idxs] # "ys" for interpolation
  interpolation = interp1d(x=beta_ms, y=moving_median_fluxes)

  detrended_t = []
  detrended_flux = []
  detrended_flux_err = []
  for i in range(0, len(selected_times)):
    detrended_t.append(selected_times[i])
    detrended_flux.append(selected_fluxes[i] / interpolation(selected_times[i]))
    detrended_flux_err.append(selected_fluxerrs[i] / interpolation(selected_times[i]))

  detrended_t = np.array(detrended_t)
  detrended_flux = np.array(detrended_flux)
  detrended_flux_err = np.array(detrended_flux_err)

  detrended_t_sorted = detrended_t[np.argsort(detrended_t)]
  detrended_flux_sorted = detrended_flux[np.argsort(detrended_t)]
  detrended_flux_err_sorted = detrended_flux_err[np.argsort(detrended_t)]


  #shift time such that t0 falls at phase=0
  shifted_t = detrended_t_sorted - t_midpt
  #print shifted_t
  for i in range(0, len(detrended_t_sorted)):
    st = shifted_t[i]
    while st > 0.5*t_orbit:
      st = st - t_orbit
      shifted_t[i] = st
    while st < -0.5*t_orbit:
      st = st + t_orbit
      shifted_t[i] = st

  notnan = (~np.isnan(detrended_t_sorted) & ~np.isnan(shifted_t) & ~np.isnan(detrended_flux_sorted) & ~np.isnan(detrended_flux_err_sorted))
  detrended_t_sorted = detrended_t_sorted[notnan]
  shifted_t = shifted_t[notnan]
  detrended_flux_sorted = detrended_flux_sorted[notnan]
  detrended_flux_err_sorted = detrended_flux_err_sorted[notnan]

  if "slc" in infile:
    cadence_flag_arr = np.ones_like(detrended_t_sorted)
  elif "llc" in infile:
    cadence_flag_arr = np.zeros_like(detrended_t_sorted)

  alldata = np.vstack((detrended_t_sorted,shifted_t,detrended_flux_sorted,detrended_flux_err_sorted,cadence_flag_arr)).T

  with open(outfile,'a') as outfile_name:
      np.savetxt(outfile_name,alldata,fmt='%f %f %f %f %d')

