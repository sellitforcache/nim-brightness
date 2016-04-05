#! /home/l_bergmann/anaconda/bin/python -W ignore

from pyne import mcnp, ace
import math, os
import pylab, numpy, sys, cPickle, progressbar, copy
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec, colors
from MCNPtools.to_energy import to_energy
from MCNPtools.to_temperature import to_temperature
from MCNPtools.to_wavelength import to_wavelength
from MCNPtools.mctal import mctal
from MCNPtools.plot import plot
import scipy.special
import numpy.linalg
import matplotlib.ticker as ticker

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=16)

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

def smooth(x,window_len=11,window='flat'):
	# take from stackexchange
	import numpy

	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."

	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."


	if window_len<3:
		return x

	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


	s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
	#print(len(s))
	if window == 'flat': #moving average
		w=numpy.ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')

	y=numpy.convolve(w/w.sum(),s,mode='valid')
	return y

def coarsen(values,bins,bin_red=2):
	import numpy
	v_out=[]
	b_out=[]
	for i in range(0,len(values)/bin_red):
		v = 0.0
		for j in range(0,bin_red):
			v = v + values[i*bin_red+j]
		v_out.append(v)
		b_out.append(bins[i*bin_red])
	b_out.append(bins[-1])
	return numpy.array(v_out),numpy.array(b_out)


def make_steps(ax,bins_in,avg_in,values_in,options=['log'],linewidth=1,color=None,label='',ylim=False):
	import numpy, re
	assert(len(bins_in)==len(values_in)+1)

	### make copies
	bins=bins_in[:]
	values=values_in[:]
	avg=avg_in[:]
	#err=err_in[:]

	### smooth data?  parse format
	for opt in options:
		res = re.match('smooth',opt)
		if res:
			smooth_opts = opt.split('=')
			if len(smooth_opts)==1:
				wlen = 7
			elif len(smooth_opts)==2:
				wlen = int(smooth_opts[1])
			else:
				wlen = int(smooth_opts[1])
				print "MULTIPLE = SIGNS IN SMOOTH.  WHY?  ACCEPTING FIRST VALUE."
			if wlen%2==0:
				print "WINDOW LENGTH EVEN, ADDING 1..."
				wlen = wlen + 1
			print "smoothing %d bins..."%wlen
			#label = label + ' SMOOTHED %d BINS'%wlen
			values = smooth(numpy.array(values),window_len=wlen)
			values = values[(wlen-1)/2:-(wlen-1)/2]   # trim to original length

	### coarsen data?  parse format
	for opt in options:
		res = re.match('coarsen',opt)
		if res:
			coarsen_opts = opt.split('=')
			if len(coarsen_opts)==1:
				bin_red = 2
			elif len(coarsen_opts)==2:
				bin_red = int(coarsen_opts[1])
			else:
				bin_red = int(coarsen_opts[1])
				print "MULTIPLE = SIGNS IN SMOOTH.  WHY?  ACCEPTING FIRST VALUE."
			if len(values)%bin_red==0:
				print "Reducing bins by factor of %d ..."%bin_red
				#label = label + ' COMBINED %d BINS'%bin_red
				values,bins = coarsen(numpy.array(values),numpy.array(bins),bin_red=bin_red)
			else:
				print "DATA LENGHTH NOT EVENLY DIVISIBLE BY COARSEN FACTOR, IGNORING..."

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
			ax.semilogy(x,y,color=color,label=label,linewidth=linewidth)
		else:
			ax.plot(x,y,color=color,label=label,linewidth=linewidth)
	else:   #default to log if lin not present
		if 'logy' in options:
			ax.loglog(x,y,color=color,label=label,linewidth=linewidth)
		else:
			ax.semilogx(x,y,color=color,label=label,linewidth=linewidth)


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
	print "detector system limit",sa_finite_d," self limit ", sa_collimator
	return numpy.minimum(sa_finite_d,sa_collimator)



# conversion factors
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1e3
charge_per_microamp = charge_per_amp/1e6
Na     = 6.0221409e+23  # number/mol
total_str = 2.1355E-05  # calculated from pinhole detector


shift=-0.0

# final measurement
measurement = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/brightness_measurement.csv','r')
for line in f:
	nums = line.split(',')
	measurement[0].append(float(nums[0]))
	measurement[1].append(float(nums[1]))
f.close()
sa_measure     = 6.28E-4 #2.0*numpy.pi*(1.0-101./numpy.sqrt(101.0*101.0+1.5*1.5))
meas_edge      = numpy.array(measurement[0][:])+shift
measurement[1] = numpy.array(measurement[1])
meas_normed    = numpy.multiply(measurement[1][:],sa_measure/total_str)

# pstudy
fracs=[0.99, 0.9, 0.8, 0.762, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01]
densities=[0.162, 0.160, 0.1591, 0.1500, 0.140, 0.130, 0.120, 0.110, 0.100]


### D2 DATA COMPARE

# get sims
bin_edges={}
bin_values={}
paths={}
cases=[]
cases.append('19 K, ENDF/B-VII.1')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case028.mctal'
cases.append('19 K, Bariloche')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-bar/results/case028.mctal'
cases.append('23 K, Bariloche')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-23K-bar/results/case028.mctal'
cases.append('24 K, IKE')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-24K-ike/results/case028.mctal'


tal_num = 5
xlims = [0.75,6]

smooth_string='smooth=11'

for case in cases:
	path = paths[case]
	this_tal = mctal(path)
	wvl1 = to_wavelength(numpy.array(this_tal.tallies[tal_num].energies[:-1]))
	widths1 = -1.0*numpy.diff(wvl1)
	dex   = this_tal.tallies[tal_num]._hash(obj=0,cos=0)
	bin_edges[case]  = wvl1
	bin_values[case]  = charge_per_milliamp*numpy.divide(this_tal.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)

# plot
f=plt.figure()
ax=f.add_subplot(111)
#ax.plot(      meas_edge,          meas_normed,     linewidth=2,label='Measurement')
for case in cases:
	make_steps(ax,bin_edges[case],[0],bin_values[case],linewidth=2,label=case,options=['lin',smooth_string])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,4e11])
ax.set_xlim(xlims)
ax.grid(1)
plt.legend(loc=1)
plt.show()


### D2O DATA COMPARE

# get sims
bin_edges={}
bin_values={}
paths={}
cases=[]
cases.append(r'ENDF/B-VII.1 D$_2$O')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case028.mctal'
cases.append(r'Bariloche D$_2$O')
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std-d2o-granda/results/case028.mctal'

for case in cases:
	path = paths[case]
	this_tal = mctal(path)
	wvl1 = to_wavelength(numpy.array(this_tal.tallies[tal_num].energies[:-1]))
	widths1 = -1.0*numpy.diff(wvl1)
	dex   = this_tal.tallies[tal_num]._hash(obj=0,cos=0)
	bin_edges[case]  = wvl1
	bin_values[case]  = charge_per_milliamp*numpy.divide(this_tal.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)

# plot
f=plt.figure()
ax=f.add_subplot(111)
#ax.plot(      meas_edge,          meas_normed,     linewidth=2,label='Measurement')
for case in cases:
	make_steps(ax,bin_edges[case],[0],bin_values[case],linewidth=2,label=case,options=['lin',smooth_string])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,4e11])
ax.set_xlim(xlims)
ax.grid(1)
plt.legend(loc=1)
plt.show()



### DENSITY COMPARE

# get sims
bin_edges={}
bin_values={}
paths={}
cases=[]
cnum=4
cases.append(r'%3.0f\%% $\rho$'%(densities[0]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(1+9*(cnum-1))
cases.append(r'%3.0f\%% $\rho$'%(densities[1]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(2+9*(cnum-1))
cases.append(r'%3.0f\%% $\rho$'%(densities[2]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(3+9*(cnum-1))
cases.append(r'%3.0f\%% $\rho$'%(densities[3]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(4+9*(cnum-1))
cases.append(r'%3.0f\%% $\rho$'%(densities[4]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(5+9*(cnum-1))
cases.append(r'%3.0f\%% $\rho$'%(densities[5]/densities[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(6+9*(cnum-1))
#cases.append(r'%3.0f\%% $\rho$'%(densities[6]/densities[0]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(7+9*(cnum-1))
#cases.append(r'%3.0f\%% $\rho$'%(densities[7]/densities[0]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(8+9*(cnum-1))
#cases.append(r'%3.0f\%% $\rho$'%(densities[8]/densities[0]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(8+9*(cnum-1))

for case in cases:
	path = paths[case]
	this_tal = mctal(path)
	wvl1 = to_wavelength(numpy.array(this_tal.tallies[tal_num].energies[:-1]))
	widths1 = -1.0*numpy.diff(wvl1)
	dex   = this_tal.tallies[tal_num]._hash(obj=0,cos=0)
	bin_edges[case]  = wvl1
	bin_values[case]  = charge_per_milliamp*numpy.divide(this_tal.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)

# plot
f=plt.figure()
ax=f.add_subplot(111)
#make_steps(ax,numpy.insert(meas_edge,0.0,0) , [0], meas_normed,     linewidth=2,label='Measurement',options=['lin'])
for case in cases:
	make_steps(ax,bin_edges[case],[0],bin_values[case],linewidth=2,label=case,options=['lin',smooth_string])
ax.set_title(r'ENDF/B-VII.1, 19K, %4.3f o-D$_2$'%fracs[cnum-1])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,4e11])
ax.set_xlim(xlims)
ax.grid(1)
plt.legend(loc=1)
plt.show()



### O/P COMPARE

# get sims
bin_edges={}
bin_values={}
paths={}
cases=[]
cnum=1
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[0]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*0)
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[1]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*1)
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[2]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*2)
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[3]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*3)
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[4]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*4)
cases.append(r'%3.0f\%% o-D$_2$'%(fracs[5]*100))
paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*5)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[6]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*6)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[7]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*7)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[8]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*8)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[9]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*9)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[10]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*10)
#cases.append(r'%3.0f\%% o-D$_2$'%(fracs[11]*100))
#paths[cases[-1]]='/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/case%03d.mctal'%(cnum+9*11)

for case in cases:
	path = paths[case]
	this_tal = mctal(path)
	wvl1 = to_wavelength(numpy.array(this_tal.tallies[tal_num].energies[:-1]))
	widths1 = -1.0*numpy.diff(wvl1)
	dex   = this_tal.tallies[tal_num]._hash(obj=0,cos=0)
	bin_edges[case]  = wvl1
	bin_values[case]  = charge_per_milliamp*numpy.divide(this_tal.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)

# plot
f=plt.figure()
ax=f.add_subplot(111)
#ax.plot(      meas_edge,          meas_normed,     linewidth=2,label='Measurement')
for case in cases:
	make_steps(ax,bin_edges[case],[0],bin_values[case],linewidth=2,label=case,options=['lin',smooth_string])
ax.set_title(r'ENDF/B-VII.1, 19K, D$_2$ Density %4.3f'%(densities[cnum-1]))
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,4e11])
ax.set_xlim(xlims)
ax.grid(1)
plt.legend(loc=1)
plt.show()






#
#
# FINAL BEST GUESS
#
#
#

this_tal1 = mctal('/home/l_bergmann/repos/ICON-brightness-parametric-24K-ike/results/case030.mctal')
this_tal2 = mctal('/home/l_bergmann/repos/ICON-brightness-parametric-24K-ike/results/case033.mctal')

wvl1 = to_wavelength(numpy.array(this_tal1.tallies[tal_num].energies[:-1]))
widths1 = -1.0*numpy.diff(wvl1)
dex   = this_tal1.tallies[tal_num]._hash(obj=0,cos=0)
values1  = charge_per_milliamp*numpy.divide(this_tal1.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)

wvl2 = to_wavelength(numpy.array(this_tal2.tallies[tal_num].energies[:-1]))
widths2 = -1.0*numpy.diff(wvl2)
dex   = this_tal2.tallies[tal_num]._hash(obj=0,cos=0)
values2  = charge_per_milliamp*numpy.divide(this_tal2.tallies[tal_num].vals[dex]['data'][:-1],widths2*total_str)

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(      meas_edge    ,          meas_normed,     linewidth=2,label='Measurement')
#ax.plot(      meas_edge-0.25,         meas_normed,     linewidth=2,label='Measurement, 0.25 {\AA}  Shift')
make_steps(ax,wvl1,[0],values1,linewidth=2,label=r'MCNP 6.1, 98\% density',options=['lin',smooth_string])
make_steps(ax,wvl2,[0],values2,linewidth=2,label=r'MCNP 6.1, 80\% density',options=['lin',smooth_string])
ax.set_title(r'24 K IKE, 0.762 o-D$_2$') #0.130 g/cm$^3$ D$_2$,
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
ax.set_ylim([0,4e11])
ax.set_xlim(xlims)
ax.grid(1)
plt.legend(loc=1)
plt.show()



#
#
# total brightness from measurement
#
#
print meas_edge
print meas_normed
this_sum = 0.0
for i in range(0,len(meas_edge)-1):
	x0=meas_edge[i]
	x1=meas_edge[i+1]
	y0=meas_normed[i]
	y1=meas_normed[i+1]
	if y0>0 and y1>0:
		a=(y1-y0)/(x1-x0)
		b=y0
		this_sum = this_sum + a/2.0*(x1-x0)*(x1-x0)+b*(x1-x0)

print "TOTAL MEASURED BRIGHTNESS = %6.4E " % this_sum
print "TOTAL CALCULATED BRIGHTNESS = %6.4E " % numpy.sum(numpy.multiply(widths1[1:],values1[1:]))
