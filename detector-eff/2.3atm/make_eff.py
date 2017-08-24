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
	plt.rc('font', size=12)

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


tal=mctal('mctal')

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
wvl=to_wavelength(energies)

eff=numpy.divide(RR_dat,incoming_dat)


wvl=wvl[::-1]
eff=eff[::-1]

wvl=wvl[:-1]
eff=eff[:-1]

f=plt.figure()
ax=f.add_subplot(111)
make_steps(ax,wvl,[0],eff,options=['lin'])
ax.grid(1)
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Efficiency (n,p)/incoming')
ax.set_ylim([0,1])
ax.set_title(r'2.3 atm He3, 1.2 atm Kr')

plt.show()

f=open('2.3He3-1.2Kr_eff.dat','w')
f.write('wavelength (AA) lower bin boundary,  wavelength (AA) upper bin boundary,  Efficiency\n')
for i in range(0,len(wvl)-1):
	effic=eff[i]
	lower=wvl[i]
	upper=wvl[i+1]
	f.write("%10.8E, %10.8E, %10.8E\n"%(lower,upper,effic))
f.close()



