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
				ax.semilogy(x,y,label=label,linewidth=linewidthl)
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
	print "detector system limit",sa_finite_d," self limit ", sa_collimator
	return numpy.minimum(sa_finite_d,sa_collimator)



# conversion factors
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1e3
charge_per_microamp = charge_per_amp/1e6
Na     = 6.0221409e+23  # number/mol
total_str = 2.1355E-05  # calculated from pinhole detector

# final measurement
measurement = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/brightness_measurement_corrected2.csv','r')
first=True
for line in f:
	if first:
		first=False
	else:
		nums = line.split(',')
		measurement[0].append(float(nums[0]))
		measurement[1].append(float(nums[1]))
f.close()

#path = '/home/l_bergmann/repos/ICON-brightness-parametric-19K-std/results/'
path = '/home/l_bergmann/repos/ICON-brightness-parametric-24K-ike/results/'
#path = '/home/l_bergmann/repos/ICON-brightness-parametric-19K-bar/results/'
#path = '/home/l_bergmann/repos/ICON-brightness-parametric-23K-bar/results/'
cases = os.listdir(path)
cases.sort()
bin_edges={}
bin_values={}
tal_num=5
for case in cases:
	this_tal = mctal(path+case)
	wvl1 = to_wavelength(numpy.array(this_tal.tallies[tal_num].energies[:-1]))
	widths1 = -1.0*numpy.diff(wvl1)
	dex   = this_tal.tallies[tal_num]._hash(obj=0,cos=0)
	bin_edges[case]  = wvl1
	bin_values[case]  = charge_per_milliamp*numpy.divide(this_tal.tallies[tal_num].vals[dex]['data'][:-1],widths1*total_str)


shift = -0.0
sa_measure     = 6.28E-4 #2.0*numpy.pi*(1.0-101./numpy.sqrt(101.0*101.0+1.5*1.5))
meas_edge      = numpy.array(measurement[0][:])+shift
measurement[1] = numpy.array(measurement[1])
meas_normed    = numpy.multiply(measurement[1][:],sa_measure/total_str)
case_min = ''
valu_min = 99999999999999.0
norm_order = 2#numpy.inf

cut = len(meas_normed)*2.0/4.0
start_bin = 10

# plot spectra
for case in cases:
	#f=plt.figure()
	#ax=f.add_subplot(211)
	#ax2=f.add_subplot(212)
	#
	#ax.plot(      meas_edge,          meas_normed,     linewidth=2,label='Measurement')
	#make_steps(ax,bin_edges[case],[0],bin_values[case],linewidth=1,label=case[:-6],options=['lin'])
	#
	bin_edges2=bin_edges[case][1::]
	case_interp = numpy.interp(meas_edge,bin_edges2[::-1],bin_values[case][::-1])
	#ax.plot(meas_edge,case_interp)
	this_norm = numpy.linalg.norm(    numpy.divide(  (case_interp[start_bin:cut]-meas_normed[start_bin:cut]) , meas_normed[start_bin:cut])  ,ord=norm_order) / (cut-1)
	#ax2.plot(meas_edge,numpy.divide(  (case_interp-meas_normed) , meas_normed))
	print case,'CUT AT %1.2f A,  rel.err. per bin %1.1f-norm: '%(meas_edge[cut],norm_order),this_norm
	if this_norm <= valu_min:
		valu_min = this_norm
		case_min = case
	##
	#ax.set_xlabel(r'Wavelength (\AA)')
	#ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ sterad$^{-1}$)')
	#ax.set_ylim([0,5e11])
	#ax.set_xlim([0,12])
	#ax.grid(1)
	##
	#ax2.set_xlabel(r'Wavelength (\AA)')
	#ax2.set_ylabel(r'Rel.Err.')
	#ax2.set_ylim([-1,1])
	#ax2.set_xlim([0,12])
	#ax2.grid(1)
	##
	#plt.legend(loc=1)
	#plt.show()

print 'minimum rel.err. %1.1f-norm: '%norm_order,case_min,valu_min 

f=plt.figure()
ax=f.add_subplot(111)
make_steps(ax,numpy.insert(meas_edge,0.0,0) , [0], meas_normed,     linewidth=2,label='Measurement',options=['lin'])
make_steps(ax,bin_edges[case_min],[0],bin_values[case_min],linewidth=1,label=case_min[:-6],options=['lin'])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ sterad$^{-1}$)')
ax.set_ylim([0,5e11])
ax.set_xlim([0,12])
ax.grid(1)
plt.legend(loc=1)
plt.show()

#case=['MCNP - STIP,  SPHERICAL','MCNP - NOSTIP, INVERTED','MCNP - NOSTIP, INVERTED, MISALIGNED']
#mctals={}
#mctal_names={}
#areas={}
#segs={}
#tals={}
#objs={}
#dex={}
#
#mctal_names[case[0]]='/home/l_bergmann/repos/temp/ICON-center-19std-STIP.mctal'
#mctal_names[case[1]]='/home/l_bergmann/repos/temp/ICON-center-19std.mctal'
#mctal_names[case[2]]='/home/l_bergmann/repos/temp/ICON-PDarraz.mctal'
#areas[case[0]]=1.0
#areas[case[1]]=1.0
#areas[case[2]]=1.0
#segs[case[0]]=0
#segs[case[1]]=0
#segs[case[2]]=0
#tals[case[0]]=5
#tals[case[1]]=5
#tals[case[2]]=5
#objs[case[0]]=0
#objs[case[1]]=0
#objs[case[2]]=7
#dex[case[0]]=0
#dex[case[1]]=0
#dex[case[2]]=0
#
#
#f=open('PD_data.csv','w')
#header = 'Lower Bin Energy (MeV)   ,   Upper Bin Energy (MeV)    '+'Lower Bin Wavelength (A)   ,   Upper Bin Wavelength (A)    ,     '+case[0]+'      ,     '+case[1]+'     ,     '+case[2]
#f.write(header+'\n')
#
## loop over guides
#for n in case:
#
#	print "Loading "+n+" ..."
#	this_tal   = tals[   n]
#	area       = areas[  n]
#	mctal_name = mctal_names[ n]
#	this_seg   = segs[   n]
#	this_obj   = objs[   n]
#
#	tal = mctal(mctal_name)
#	tal.tallies[this_tal].energies[0] = 9e-13
#	mctals[n]=tal
#
#	dex[n] = tal.tallies[this_tal]._hash(obj=this_obj,seg=this_seg)
#
#
#ene = numpy.array(mctals[case[0]].tallies[this_tal].energies[:-1])
#ene_a = to_wavelength(ene)
#
#for i in range(0,len(ene)-1):
#
#	v0 = mctals[case[0]].tallies[tals[case[0]]].vals[dex[case[0]]]['data'][i]
#	v1 = mctals[case[1]].tallies[tals[case[1]]].vals[dex[case[1]]]['data'][i]
#	v2 = mctals[case[2]].tallies[tals[case[2]]].vals[dex[case[2]]]['data'][i]
#	e0 =   ene[i]
#	e1 =   ene[i+1]
#	a0 = ene_a[i]
#	a1 = ene_a[i+1]
#
#	f.write("% 6.4E ,  % 6.4E ,   % 6.4E ,   % 6.4E ,   % 6.4E ,   % 6.4E ,   % 6.4E\n"%(e0,e1,a0,a1,v0,v1,v2))
#
#f.close()