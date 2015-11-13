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

#
#    function to fold a fine-binned xs with a coarse spectrum
#
def rebin_xs(xs=0,xs_e=0,spec_e=0):
	# needs to be pointwise, NOT binned
	assert(len(xs_e)==len(xs))
	#print "length of xs",len(xs)
	#print "length of spec bins",len(spec_e)

	# 
	spec_xs=[]
	for i in range(0,len(spec_e)-1):
		# get requested bin widths
		low  = spec_e[i]
		high = spec_e[i+1]
		# do logic on xs E grid
		logic_low  = xs_e < low
		logic_high = xs_e > high
		dex_low  = numpy.nonzero(numpy.diff(logic_low))[0]
		dex_high = numpy.nonzero(numpy.diff(logic_high))[0]
		# figure out edge cases
		if len(dex_low) == 0:
			if logic_low[0]:	# all ones, low is above last point
				dex_low = len(xs_e)-1
			else:				# all zeros, low is below first point
				dex_low = 0
		else:
			dex_low = dex_low[0]
		if len(dex_high) == 0:
			if logic_high[0]:   # all ones, high is below first point
				dex_high = 0
			else:				# all zeros, high is above last point
				dex_high = len(xs_e)-1
		else:
			dex_high = dex_high[0]
		#print dex_low,dex_high
		# average the pointwise data 
		if dex_low == dex_high:  # bin is within two xs points
			if dex_high == len(xs_e)-1:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high]
				a = 0.0
			else:
				e_l  = xs_e[dex_low]
				e_h  = xs_e[dex_high+1]
				xs_l = xs[  dex_low]
				xs_h = xs[  dex_high+1]
				a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg = (a/2.0)*(high*high-low*low)/(high-low)+b
		else:
			avg_vals=[]
			avg_widths=[]
			#do first bin
			e_l  = xs_e[dex_low]
			e_h  = xs_e[dex_low+1]
			xs_l = xs[  dex_low]
			xs_h = xs[  dex_low+1]
			a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg_vals.append( (a/2.0)*(e_h*e_h-low*low)/(e_h-low)+b )
			avg_widths.append(e_h-low)
			#do middle bins
			for i in range(dex_low,dex_high-1):
				e_l  = xs_e[i]
				e_h  = xs_e[i+1]
				xs_l = xs[  i]
				xs_h = xs[  i+1]
				a = (xs_h-xs_l)/(e_h-e_l)
				b = xs_l - a*e_l
				avg_vals.append( (a/2.0)*(e_h*e_h-e_l*e_l)/(e_h-e_l)+b )
				avg_widths.append(e_h-e_l)
			#do last bin
			if dex_high == len(xs_e)-1:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high]
				a=0.0
			else:
				e_l  = xs_e[dex_high]
				e_h  = xs_e[dex_high+1]
				xs_l = xs[  dex_high]
				xs_h = xs[  dex_high+1]
				a = (xs_h-xs_l)/(e_h-e_l)
			b = xs_l - a*e_l
			avg_vals.append( (a/2.0)*(high*high-e_l*e_l)/(high-e_l)+b )
			avg_widths.append(high-e_l)
			#avg by bin width and append value
			avg_widths = numpy.array(avg_widths) # normalize
			avg = numpy.average(avg_vals,weights=avg_widths)
		spec_xs.append(avg)
	# return array
	return numpy.array(spec_xs)



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


def get_view_sa(a=0.05,b=0,l=0,thk=0):
	import numpy
	### finite detector, far side of collimator
	la = l/(b/a+1.0)
	cosa = la/numpy.sqrt(la*la+a*a)
	sa_finite_d = 2.0*numpy.pi*(1.0-cosa)

	### collimator self-limit
	cosa = (thk/2.0)/numpy.sqrt(thk*thk/4.0+b*b)
	sa_collimator = 2.0*numpy.pi*(1.0-cosa)

	### return minimum
	return numpy.minimum(sa_finite_d,sa_collimator)



# conversion factors
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1e3
charge_per_microamp = charge_per_amp/1e6
Na     = 6.0221409e+23  # number/mol



# tallies
std19 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarraz.mctal')
ike19 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-19IKE.mctal')
ike24 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-24IKE.mctal')
this_tal = std19
std19.plot(tal=[5],obj=[8],options=['log','lethargy'])

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

# get new eff
mcnp_eff2 = [[],[]]
seam = mctal('/home/l_bergmann/Documents/Other Projects/ICON-brightness/detector/2mm.mctal')
ene = numpy.array(seam.tallies[1].energies[:-1])
dex =             seam.tallies[1]._hash(obj=0,cos=0)
tal_in =          seam.tallies[1].vals[dex]['data'][:-1]
dex =             seam.tallies[34]._hash(obj=0,mul=2)
tal_rr =          seam.tallies[34].vals[dex]['data'][:-1]
wvl = to_wavelength(ene)
mcnp_eff2[0] = wvl
mcnp_eff2[1] = numpy.divide(tal_rr,tal_in)

# plot efficiency
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(oned_eff[0],oned_eff[1],  linewidth=2,label='1-D Model')
ax.plot(mcnp_eff[0],mcnp_eff[1],  linewidth=2,label='MCNP - no seam')
make_steps(ax,mcnp_eff2[0],[0],mcnp_eff2[1],linewidth=2,label='MCNP - 1 mm seam',options=['lin'])
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


#### redivide by new efficicency
#old_eff = rebin_xs(xs=oned_eff[1],xs_e=oned_eff[0],       spec_e=numpy.array(measurement[0] ))
#new_eff = rebin_xs(xs=mcnp_eff2[1],xs_e=mcnp_eff2[0][1:], spec_e=numpy.array(measurement[0] ))
#measurement[1] = numpy.multiply(measurement[1][1:],numpy.divide(old_eff,new_eff))

# scale and normalize the simulation data
ene = numpy.array(this_tal.tallies[5].energies[:-1])
dex = this_tal.tallies[5]._hash(obj=8)
tal = this_tal.tallies[5].vals[dex]['data'][:-1]
wvl = to_wavelength(ene)
widths = -1.0*numpy.diff(wvl)

### solid angles
sa_collimator = get_view_sa(b=0.1,  l=10.1,   thk=10.0)
sa_cadmium    = get_view_sa(b=0.05, l= 0.1,   thk= 0.1)
sa_sapphire   = get_view_sa(b=1.0,  l=41.835, thk=25.0)
sa_steel      = get_view_sa(b=2.5,  l=784.8,  thk=30.0)
sa_aperture   = get_view_sa(b=4.0,  l=1230.8, thk=0.1)
sa_zapfen     = get_view_sa(b=4.0,  l=1577.5, thk=347.0)

sa_measure       = 2.0*numpy.pi*(1.0-101./numpy.sqrt(101.0*101.0+1.5*1.5))


print "solid angle collimator     %6.4E"%sa_collimator
print "solid angle cadmium        %6.4E"%sa_cadmium
print "solid angle sapphire       %6.4E"%sa_sapphire
print "solid angle steel          %6.4E"%sa_steel
print "solid angle aperture       %6.4E"%sa_aperture
print "solid angle zapfen         %6.4E"%sa_zapfen
print "solid angle measurement    %6.4E"%sa_measure
tal_normed = charge_per_milliamp*numpy.divide(tal,widths*sa_measure)

measurement[1] = numpy.array(measurement[1])


# CORRECTION? UNTIL DO REAL SIMULATION
target_gain = 1.0#*1.4*1.20*1.05  #  collotte, STIP, light water in loop


# plot spectra
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(measurement[0][:],measurement[1],linewidth=2,label='Measurement')
make_steps(ax,wvl,[0],tal_normed/target_gain,linewidth=2,label='MCNP',options=['lin'])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,3.2e10])
ax.set_xlim([0,12])
ax.grid(1)
plt.legend(loc=1)
plt.show()
