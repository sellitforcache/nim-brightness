#! /home/l_bergmann/anaconda/bin/python -W ignore
# /usr/local/bin/python 
#

from MCNPtools.mctal import mctal
from MCNPtools.to_wavelength import to_wavelength
import matplotlib.pyplot as plt
import re, numpy, sys, os, errno
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from subprocess import call
import xlwt

tex=True

if tex:
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rc('font', size=16)

def make_steps(ax,bins_in,avg_in,values_in,options=['log'],color=None,label='',ylim=False,linewidth=1):
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
			label = label + ' SMOOTHED %d BINS'%wlen
			values = self._smooth(numpy.array(values),window_len=wlen)
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
				label = label + ' COMBINED %d BINS'%bin_red
				values,bins = self._coarsen(numpy.array(values),numpy.array(bins),bin_red=bin_red)
			else:
				print "DATA LENGHTH NOT EVENLY DIVISIBLE BY COARSEN FACTOR, IGNORING..."

	### make rectangles
	x=[]
	y=[]
	#x.append(bins[0])
	#y.append(0.0)
	for n in range(len(values)):
		x.append(bins[n])
		x.append(bins[n+1])
		y.append(values[n])
		y.append(values[n])
	#x.append(bins[len(values)])
	#y.append(0.0)

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




tals={}
tals['specified.m']	= r'0.02+0.0007 bar'
tals['plus.m' ]		= r'0.02        bar'
tals['minus.m']		= r'0.02+0.0007 bar'

wvl={}
eff={}
eff_err={}

for talname in tals.keys():

	tal=mctal(talname)
	
	tally_number = 34
	cos=0
	dex = tal.tallies[tally_number]._hash(cos=cos,obj=0,seg=0,mul=2)
	RR_dat=tal.tallies[tally_number].vals[dex]['data'][:-1]
	RR_err=tal.tallies[tally_number].vals[dex]['err' ][:-1]
	
	
	tally_number = 1
	cos=0
	dex = tal.tallies[tally_number]._hash(cos=cos,obj=0,seg=0)
	incoming_dat=tal.tallies[tally_number].vals[dex]['data'][:-1]
	incoming_err=tal.tallies[tally_number].vals[dex]['err' ][:-1]
	
	tally_number = 1
	cos=1
	dex = tal.tallies[tally_number]._hash(cos=cos,obj=0,seg=0)
	outgoing_dat=tal.tallies[tally_number].vals[dex]['data'][:-1]
	outgoing_err=tal.tallies[tally_number].vals[dex]['err' ][:-1]
	
	energies = tal.tallies[tally_number].energies[:-1]
	this_wvl=to_wavelength(energies)
	
	this_eff=numpy.divide(RR_dat,incoming_dat)
	this_eff_err=numpy.sqrt(  numpy.add(numpy.multiply(RR_err,RR_err) , numpy.multiply(incoming_err,incoming_err)))
	
	
	wvl[talname]=this_wvl[::-1]
	eff[talname]=this_eff[::-1]
	eff_err[talname]=this_eff_err[::-1]
	
	wvl[talname]=wvl[talname][:-1]
	eff[talname]=eff[talname][:-1]
	eff_err[talname]=eff_err[talname][:-1]
	



mean_wvl	= wvl['specified.m']
mean_eff	= eff['specified.m']

wvl_avg 	= (wvl['specified.m'][1:]+wvl['specified.m'][:-1])/2
err_plus	= eff['plus.m' ] + numpy.multiply(eff['plus.m' ],eff_err['plus.m' ]) + numpy.multiply(eff['specified.m'],eff_err['specified.m'])
err_minus	= eff['minus.m'] - numpy.multiply(eff['minus.m'],eff_err['minus.m']) - numpy.multiply(eff['specified.m'],eff_err['specified.m'])

mean_rel_err 	= numpy.divide((err_plus-mean_eff),mean_eff)

# quad fits for error
x=numpy.linspace(0,12,1000)
z = numpy.polyfit(wvl_avg, mean_rel_err, 3)
p = numpy.poly1d(z)
print ""
print "ax^3 + bx^2 + cx + d = err"
print "a = % 6.4E"%z[0]
print "b = % 6.4E"%z[1]
print "c = % 6.4E"%z[2]
print "d = % 6.4E"%z[3]
print ""

# function fit for eff
def eff_func(x,a1,b1,a2,b2):
	t1=3.2
	t2=0.05
	lambda1=a1*x+b1
	lambda2=a2*x+b2
	return (1-numpy.exp(-lambda1*t1))*(numpy.exp(-lambda2*t2))
p0 = 	[ 1.6016E-01/100.,
		  1.8993E-04,
		  2.8298E-01,
		 -5.8576E-01]
popt, pcov = curve_fit(eff_func, wvl_avg, mean_eff, p0, maxfev=100000)
print ""
print "( 1-exp(-(a1*x+b1)*3.2) ) * exp(-(a2*x+b2)*0.05) = eff"
print "a1 = % 6.4E"%popt[0]
print "b1 = % 6.4E"%popt[1]
print "a2 = % 6.4E"%popt[2]
print "b2 = % 6.4E"%popt[3]
print ""


f=plt.figure()
ax=f.add_subplot(111)
ax.fill_between(wvl_avg,err_plus,err_minus, facecolor='red', linewidth=1.0, color='red', alpha=0.25,label=r'Statistical + Pressure Uncertainty')
make_steps(ax,mean_wvl,[0],mean_eff,options=['lin'],linewidth=2)
ax.plot(x,eff_func(x,*popt), linewidth=2.0, color='k', label='fit')


ax.grid(1)
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Efficiency')
ax.set_ylim([0,6e-2])
ax.set_xlim([0,12])
#ax.set_title(r'2.3 atm $^3$He, 1.2 atm Kr')
plt.legend(loc='best')
plt.tight_layout()
f.savefig('eff_0.02bar.pdf')

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(x,p(x), linewidth=2.0, color='k',label='fit')
make_steps(ax,mean_wvl,[0],mean_rel_err,options=['lin'],linewidth=2,color='g',label='rel.err. data')
ax.grid(1)
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Efficiency rel. err.')
#ax.set_ylim([0,1])
ax.set_xlim([0,12])
plt.legend(loc='best')
plt.tight_layout()
f.savefig('relerr_0.02bar.pdf')

plt.show()

f=open('0.02barHe3-1.2barKr_eff.dat','w')
f.write('wavelength (AA) lower bin boundary,  wavelength (AA) upper bin boundary,  Efficiency, rel. err.\n')
for i in range(0,len(mean_eff)-1):
	effic 	=mean_eff[i]
	err 	=mean_rel_err[i]
	lower 	=mean_wvl[i]
	upper 	=mean_wvl[i+1]
	f.write("%10.8E, %10.8E, %10.8E, %10.8E\n"%(lower,upper,effic,err))
f.close()



