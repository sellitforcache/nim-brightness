#! /home/l_bergmann/anaconda/bin/python -W ignore

from pyne import mcnp, ace
import math
import pylab, numpy, sys, cPickle, progressbar, copy
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
from MCNPtools.to_energy import to_energy
from MCNPtools.to_temperature import to_temperature
from MCNPtools.to_wavelength import to_wavelength
from MCNPtools.mctal import mctal
from MCNPtools.plot import plot
import scipy.special
import numpy.linalg

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)

def make_steps(ax,bins_in,avg_in,values_in,options=['log'],label='',ylim=False,color=False,linewidth=1):
	import numpy
	assert(len(bins_in)==len(values_in)+1)

	### make copies
	bins=bins_in[:]
	values=values_in[:]
	avg=avg_in[:]
	#err=err_in[:]

	### make rectangles
	x=[]
	y=[]
	x.append(bins[0])
	y.append(0.0)
	for n in range(len(values)):
		x.append(bins[n])
		x.append(bins[n+1])
		y.append(values[n])
		y.append(values[n])
	x.append(bins[len(values)])
	y.append(0.0)

	### plot with correct scale
	if 'lin' in options:
		if 'logy' in options:
			if color:
				ax.semilogy(x,y,label=label,color=color,linewidth=linewidth)
			else:
				ax.semilogy(x,y,label=labe,linewidth=linewidthl)
		else:
			if color:
				ax.plot(x,y,label=label,color=color,linewidth=linewidth)
			else:
				ax.plot(x,y,label=label,linewidth=linewidth)
	else:   #default to log if lin not present
		if 'logy' in options:
			if color:
				ax.loglog(x,y,label=label,color=color,linewidth=linewidth)
			else:
				ax.loglog(x,y,label=label,linewidth=linewidth)
		else:
			if color:
				ax.semilogx(x,y,label=label,color=color,linewidth=linewidth)
			else:
				ax.semilogx(x,y,label=label,linewidth=linewidth)

# conversion factors
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1e3
charge_per_microamp = charge_per_amp/1e6
Na     = 6.0221409e+23  # number/mol



# tallies
std19 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarraz.mctal')
ike19 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-19IKE.mctal')
ike24 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-24IKE.mctal')

# 1d model efficiency
oned_eff = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/1deff.csv','r')
for line in f:
	nums = line.split(',')
	oned_eff[0].append(float(nums[0]))
	oned_eff[1].append(float(nums[1]))
f.close()

# mcnp calculated efficieny
mcnp_eff = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/mcnpeff.csv','r')
for line in f:
	nums = line.split(',')
	mcnp_eff[0].append(float(nums[0]))
	mcnp_eff[1].append(float(nums[1]))
f.close()


# plot efficiency
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(oned_eff[0],oned_eff[1],linewidth=2,label='1-D Model')
ax.plot(mcnp_eff[0],mcnp_eff[1],linewidth=2,label='MCNP')
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Efficiency')
ax.set_ylim([0,1])
ax.set_xlim([0,12])
ax.grid(1)
plt.legend(loc=4)
plt.show()


# final measurement
measurement = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/meas.csv','r')
for line in f:
	nums = line.split(',')
	measurement[0].append(float(nums[0]))
	measurement[1].append(float(nums[1]))
f.close()

# scale and normalize the simulation data
ene = numpy.array(ike19.tallies[5].energies[:-1])
dex = ike24.tallies[5]._hash(obj=5)
tal = ike24.tallies[5].vals[dex]['data'][:-1]
wvl = to_wavelength(ene)
widths = -1.0*numpy.diff(wvl)
sa_PD      = 2.0*numpy.pi*(1.0-102./numpy.sqrt(102.0*102.0+0.5*0.5))
sa_measure = 2.0*numpy.pi*(1.0-51./numpy.sqrt(51.0*51.0+0.5*0.5))
#sa = 6.93e-4
print sa_PD,sa_measure
tal_normed = charge_per_milliamp*numpy.divide(tal,widths*sa_measure)

measurement[1] = numpy.array(measurement[1])*6.93e-4/sa_measure

# plot spectra
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(measurement[0],measurement[1],linewidth=2,label='Measurement')
make_steps(ax,wvl,[0],tal_normed,linewidth=2,label='MCNP',options=['lin'])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
#ax.set_ylim([0,1.4e10])
ax.set_xlim([0,12])
ax.grid(1)
plt.legend(loc=1)
plt.show()
