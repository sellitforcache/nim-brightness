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
ax.grid(1)
plt.legend(loc=2)
plt.show()