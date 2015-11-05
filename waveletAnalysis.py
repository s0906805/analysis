import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

__author__ = 'Evgeniya Predybaylo'


# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------

# READ THE DATA
sst = np.loadtxt('olw_InterpClim.csv', delimiter = ',')  # input SST time series

#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
variance = np.std(sst, ddof=1) ** 2
#sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
n = len(sst)
dt = 1 # units of 15 mins
time = range(0, n)  # construct time array
print sst[0], sst[3] #len(time), len(sst)

xlim = ([1, len(sst)])  # plotting range
pad = 1  # pad the time series with zeroes (recommended)
dj = 0.25  # this will do 10 sub-octaves per octave
s0 = 7* dt  # this says start at a scale of 15 mins
j1 = 16 / dj  # this says do 14 powers-of-two with dj sub-octaves each
lag1 = 0.72  # lag-1 autocorrelation for red noise background
mother = 'PAUL'

# Wavelet transform:
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
print period
# Significance levels: (variance=1 for the normalized SST)
signif = wave_signif(([1.0]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
global_ws = variance * (np.sum(power, axis=1) / n)  # time-average over all times
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother)

# Scale-average between El Nino periods of 2--8 years
avg = np.logical_and(scale >= 2, scale < 8)
Cdelta = 0.776  # this is for the MORLET wavelet
scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
scale_avg = power / scale_avg  # [Eqn(24)]
scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother)

#------------------------------------------------------ Plotting

#--- Plot time series
figure = plt.figure(figsize=(18, 9))
plt.subplot(321)
plt.plot(time, sst,alpha=0.5, rasterized=True)
plt.xlim(xlim[:])
plt.xlabel('Time (datapoint)')
plt.ylabel('W/m^2')
plt.title('a) Outgoing longwave')
plt.hold(False)


period = period/(365*96)
coi = coi/(365*96)
print np.min(power), np.max(power)
#--- Contour plot wavelet power spectrum
plt3 = plt.subplot(323)
levels = [0.125, 0.25, 0.5, 1, 2]
CS = plt.contourf(time, period, power,alpha=0.5, rasterized=True)  #*** or use 'contour'
im = plt.contourf(CS)
plt.xlabel('Time (datapoint)')
plt.ylabel('Period (years)')
plt.title('b) outgoing longwave Wavelet Power Spectrum ')
plt.xlim(xlim[:])
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
plt.hold(True)
plt.contour(time, period, sig95, [-99, 1], colors='k',alpha=0.5, rasterized=True)
# cone-of-influence, anything "below" is dubious
plt.plot(time, coi, 'k',alpha=0.5, rasterized=True)
plt.hold(False)
# format y-scale
#plt3.set_yscale('log', basey=2, subsy=None)
plt.ylim([1/366, 4])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt3.ticklabel_format(axis='y', style='plain')
plt3.invert_yaxis()
# set up the size and location of the colorbar
divider = make_axes_locatable(plt3)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
plt.colorbar(im, cax=cax, orientation='horizontal')

#--- Plot global wavelet spectrum
plt4 = plt.subplot(325)
plt.plot(global_ws, period,alpha=0.5, rasterized=True)
plt.hold(True)
plt.plot(global_signif, period, '--',alpha=0.5, rasterized=True)
plt.hold(False)
plt.xlabel('Power ')
plt.ylabel('Period (years)')
plt.title('c) Global Wavelet Spectrum')
plt.xlim([0, 1.25 * np.max(global_ws)])
# format y-scale
plt4.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt4.ticklabel_format(axis='y', style='plain')
plt4.invert_yaxis()

# READ THE DATA
sst = np.loadtxt('rsw_InterpClim.csv', delimiter = ',') 
#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
variance = np.std(sst, ddof=1) ** 2
#sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
n = len(sst)
dt = 1 # units of 15 mins
time = range(0, n)  # construct time array
print len(time), len(sst)

xlim = ([1, len(sst)])  # plotting range
pad = 1  # pad the time series with zeroes (recommended)
dj = 0.25  # this will do 10 sub-octaves per octave
s0 = 7* dt  # this says start at a scale of 15 mins
j1 = 16 / dj  # this says do 14 powers-of-two with dj sub-octaves each
lag1 = 0.72  # lag-1 autocorrelation for red noise background
mother = 'PAUL'

# Wavelet transform:
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
print period
# Significance levels: (variance=1 for the normalized SST)
signif = wave_signif(([1.0]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
global_ws = variance * (np.sum(power, axis=1) / n)  # time-average over all times
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother)

# Scale-average between El Nino periods of 2--8 years
avg = np.logical_and(scale >= 2, scale < 8)
Cdelta = 0.776  # this is for the MORLET wavelet
scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
scale_avg = power / scale_avg  # [Eqn(24)]
scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother)

#------------------------------------------------------ Plotting

#--- Plot time series
plt.subplot(322)
plt.plot(time, sst, alpha=0.5, rasterized=True)
plt.xlim(xlim[:])
plt.xlabel('Time (datapoint)')
plt.ylabel('W/m^2')
plt.title('a) reflected shortwave')
plt.hold(False)


period = period/(365*96)
coi = coi/(365*96)
print np.min(power), np.max(power)
#--- Contour plot wavelet power spectrum
plt3 = plt.subplot(324)
levels = [0.125, 0.25, 0.5, 1, 2]
CS = plt.contourf(time, period, power, alpha=0.5, rasterized=True)  #*** or use 'contour'
im = plt.contourf(CS)
plt.xlabel('Time (datapoint)')
plt.ylabel('Period (years)')
plt.title('b) reflected shortwave Wavelet Power Spectrum ')
plt.xlim(xlim[:])
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
plt.hold(True)
plt.contour(time, period, sig95, [-99, 1], colors='k',alpha=0.5, rasterized=True)
# cone-of-influence, anything "below" is dubious
plt.plot(time, coi, 'k', alpha=0.5, rasterized=True)
plt.hold(False)
# format y-scale
#plt3.set_yscale('log', basey=2, subsy=None)
plt.ylim([1/366, 4])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt3.ticklabel_format(axis='y', style='plain')
plt3.invert_yaxis()
# set up the size and location of the colorbar
divider = make_axes_locatable(plt3)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
plt.colorbar(im, cax=cax, orientation='horizontal')

#--- Plot global wavelet spectrum
plt4 = plt.subplot(326)
plt.plot(global_ws, period, alpha=0.5, rasterized=True)
plt.hold(True)
plt.plot(global_signif, period, '--', alpha=0.5, rasterized=True)
plt.hold(False)
plt.xlabel('Power')
plt.ylabel('Period (years)')
plt.title('c) Global Wavelet Spectrum')
plt.xlim([0, 1.25 * np.max(global_ws)])
# format y-scale
plt4.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt4.ticklabel_format(axis='y', style='plain')
plt4.invert_yaxis()




plt.tight_layout()


plt.savefig('plots/spec_InterpClim.pdf')

# end of code
