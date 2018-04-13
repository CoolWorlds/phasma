#!/usr/bin/env python

import numpy as np
import sys
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from astropy.stats import median_absolute_deviation
from scipy.special import erfcinv
from scipy.stats import binned_statistic


infile = sys.argv[1] #string: name of input kepler .fits file
outfile = sys.argv[2] #string: name of master output file for this planet.
binned_outfile = sys.argv[3] #string: name of binned master output file for this planet.
t_orbit = float(sys.argv[4]) #orbital period, days
t_14 = (float(sys.argv[5]))/24. #transit duration, hours converted to days
t_midpt = float(sys.argv[6]) #transit ephemeris, days
srad = float(sys.argv[7])
smass = float(sys.argv[8])
steff = float(sys.argv[9])
rp_rstar = float(sys.argv[10])
inc = float(sys.argv[11])
nplanet = int(sys.argv[12])
NEAdisp = sys.argv[13]
Kepdisp = sys.argv[14]
DVlink = "https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/" + sys.argv[15]

alldata = np.genfromtxt(infile)
#columns:
#alldata[:,0] = t
#alldata[:,1] = shifted_t
#alldata[:,2] = flux
#alldata[:,3] = flux_err
#alldata[:,4] = cadence flag

#notes on mathematica version:
# "detrend" and "lightcurve" columns: t, flux, fluxerr
# "fold" and "foldsort" columns: shifted_t, flux, fluxerr

#first pass
alldata_sort = alldata[np.argsort(alldata[:,1])]
intransit_mask = ((-0.55*t_14 < alldata_sort[:,1]) & (alldata_sort[:,1] < 0.55*t_14))
n_intransit_points = len(alldata_sort[:,1][intransit_mask])
movmed_kernelsize = np.max((np.round(n_intransit_points/100.), 10))

#make sure it's an odd number
if movmed_kernelsize %2 == 0:
	movmed_kernelsize += 1

movmed_kernelsize = int(movmed_kernelsize)
movmed_t = medfilt(alldata_sort[:,1], movmed_kernelsize)
movmed_flux = medfilt(alldata_sort[:,2], movmed_kernelsize)

movmed_interpolation = interp1d(x=movmed_t, y=movmed_flux)

ta = np.min(movmed_t)
tb = np.max(movmed_t)

times_mask = ((ta <= alldata[:,1]) & (alldata[:,1] <= tb))

selected_times = alldata[:,1][times_mask]
selected_fluxes = alldata[:,2][times_mask]
selected_fluxerrs = alldata[:,3][times_mask]

res = []
for i in range(0, len(selected_times)):
	try:
  		res.append((selected_fluxes[i] - movmed_interpolation(selected_times[i]))/selected_fluxerrs[i])
  	except ValueError:
  		pass

res = np.array(res)

fac = 1.4286 * median_absolute_deviation(res)
sig = np.sqrt(2.) * erfcinv(1./len(res))

lightcurve_idxs = []
for j in range(len(alldata)):
	try:
		if np.abs((alldata[j,2] - movmed_interpolation(alldata[j, 1]))/alldata[j,3]) < (fac*sig):
			lightcurve_idxs.append(j)
	except ValueError:
		pass

lightcurve_idxs = np.array(lightcurve_idxs)

#second pass
alldata = alldata[lightcurve_idxs]
alldata_sort = alldata[np.argsort(alldata[:,1])]
intransit_mask = ((-0.55*t_14 < alldata_sort[:,1]) & (alldata_sort[:,1] < 0.55*t_14))
n_intransit_points = len(alldata_sort[:,1][intransit_mask])
movmed_kernelsize = np.max((np.round(n_intransit_points/100.), 10))

#make sure it's an odd number
if movmed_kernelsize %2 == 0:
	movmed_kernelsize += 1

movmed_kernelsize = int(movmed_kernelsize)
movmed_t = medfilt(alldata_sort[:,1], movmed_kernelsize)
movmed_flux = medfilt(alldata_sort[:,2], movmed_kernelsize)

movmed_interpolation = interp1d(x=movmed_t, y=movmed_flux)

ta = np.min(movmed_t)
tb = np.max(movmed_t)

times_mask = ((ta <= alldata[:,1]) & (alldata[:,1] <= tb))

selected_times = alldata[:,1][times_mask]
selected_fluxes = alldata[:,2][times_mask]
selected_fluxerrs = alldata[:,3][times_mask]

res = []
for i in range(0, len(selected_times)):
	try:
  		res.append((selected_fluxes[i] - movmed_interpolation(selected_times[i]))/selected_fluxerrs[i])
  	except ValueError:
  		pass

res = np.array(res)

fac = 1.4286 * median_absolute_deviation(res)
sig = np.sqrt(2.) * erfcinv(1./len(res))

lightcurve_idxs = []
for j in range(len(alldata)):
	try:
		if np.abs((alldata[j,2] - movmed_interpolation(alldata[j, 1]))/alldata[j,3]) < (fac*sig):
			lightcurve_idxs.append(j)
	except ValueError:
		pass

lightcurve_idxs = np.array(lightcurve_idxs)


#third/final pass
alldata = alldata[lightcurve_idxs]
alldata_sort = alldata[np.argsort(alldata[:,1])]
phase = alldata_sort[:,1]/t_orbit
bins = np.linspace(-0.5,0.5,128)
#bins = np.arange(-0.5,0.502,0.002)

zeta_zs = np.empty((len(bins) - 1, 4))
zeta_lens = np.empty((len(bins) - 1))

zeta_zs[:,:] = np.nan
zeta_lens[:] = np.nan
for t in range(0, len(bins)-1):
	bin_mask = ((bins[t] < phase) & (phase < bins[t+1]))
	weights = (alldata_sort[:,3][bin_mask])**(-2)
	if len(weights) > 0:
		V1 = np.sum(weights)
		V2 = np.sum(weights**2)
		mustar = np.sum(alldata_sort[:,2][bin_mask] * weights)/np.sum(weights)
		s2 = np.sum(weights * (alldata_sort[:,2][bin_mask] - mustar)**2) / (V1 - (V2/V1))
		#print (mustar, s2)
		zeta_zs[t,0] = 0.5*(bins[t] + bins[t+1])
		zeta_zs[t,1] = mustar
		zeta_zs[t,2] = np.sqrt(s2)/np.sqrt(len(weights) - 1)
		zeta_zs[t,3] = np.sum(weights)**(-0.5)
		zeta_lens[t] = len(weights)
		#print zeta_zs[t]

zeta_zs = zeta_zs[~np.isnan(zeta_lens)]
#convert from phase to time
zeta_zs[:,0]=zeta_zs[:,0]*t_orbit

phase = alldata[:,1]/t_orbit
new_alldata = np.vstack((alldata[:,0].T,alldata[:,1].T,phase.T,alldata[:,2].T,alldata[:,3].T,alldata[:,4].T)).T

outfile_handle = open(outfile, 'a')
with outfile_handle as of:
	#of.write("# period [days]: {0},\n".format(t_orbit))
	#of.write("# stellar radius [solar radii]: {0},\n".format(srad))
	#of.write("# stellar mass [solar masses]: {0},\n".format(smass))
	#of.write("# stellar effective temp [K]: {0},\n".format(steff))

	np.savetxt(of,new_alldata,fmt="%f,%f,%f,%f,%f,%d",header="time,shifted_time,phase,pdc_flux_norm,pdc_flux_norm_err,cadence_flag")

	of.close()

#alldata[:,0] = t
#alldata[:,1] = shifted_t
#alldata[:,2] = flux
#alldata[:,3] = flux_err
#alldata[:,4] = cadence flag

binned_outfile_handle = open(binned_outfile, 'a')
with binned_outfile_handle as of:
	of.write("period[days],M[Msun],Mp[Msun],R[Rsun],Rp[Rsun],e,theta_v[rad],phi_v[rad],Teff[K],nplanet,NEA_disposition,Kepler_disposition,DV_report_link\n")
	of.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}\n".format(t_orbit,smass,0.0009543,srad,rp_rstar*srad,0,inc*(np.pi/180.),0,steff,nplanet,NEAdisp,Kepdisp,DVlink))
	np.savetxt(of,zeta_zs[:,0:3],delimiter=',',header="binned_time,binned_flux,binned_flux_err")
	#np.savetxt(of,binned_data,delimiter=',',header="binned_t,binned_flux,binned_flux_err")

	of.close()
