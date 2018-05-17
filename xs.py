#! /home/l_bergmann/anaconda/bin/python -W ignore

#from pyne import mcnp
from pyne import ace

import math
import numpy, sys, cPickle, copy
import matplotlib.pyplot as plt
from MCNPtools.to_energy import to_energy
from MCNPtools.to_wavelength import to_wavelength


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=12)

#prefix = '/usr/local/LANL/MCNP6_DATA'
prefix = '/home/l_bergmann/LANL/MCNP6_DATA'


# xs
o16_lib        = ace.Library(prefix+'/xdata/endf71x/O/8016.710nc')
al27_lib       = ace.Library(prefix+'/xdata/endf71x/Al/13027.710nc')

# sab
al27_sab_lib   = ace.Library(prefix+'/xdata/ENDF71SaB/al27.22t')
sapp_sab_lib   = ace.Library(prefix+'/sapphire')


# read
o16_lib.read()
al27_lib.read()
al27_sab_lib.read()
sapp_sab_lib.read()


# get libs
o16_lib         = o16_lib.tables['8016.710nc']
al27_lib        = al27_lib.tables['13027.710nc']
al27_sab_lib    = al27_sab_lib.tables[ 'al27.22t']
osapp_sab_lib   = sapp_sab_lib.tables[ 'osapp.00t']
alsapp_sab_lib  = sapp_sab_lib.tables['alsapp.00t']




####  al27 
Md     = al27_lib.awr * 1.008664916
Na      = 6.022e23
b       = 1e-24 
density = 2.66
Nd =  density * Na * b / Md
Nd_al27 =  density * Na * b / Md
al27_xs    = al27_lib.reactions[  2].sigma * Nd
al27_a     = al27_lib.reactions[102].sigma * Nd
al27_e     = al27_lib.energy
al27_sab   = al27_sab_lib.inelastic_sigma * Nd
al27_sab_e = al27_sab_lib.inelastic_e_in
#find e dex where sab stops
al27_dex = numpy.where( al27_e >= al27_sab_e[-1] )[0][0]
# add elastic to inelastic
for i in range(0,len(al27_sab_e)):
	this_e = al27_sab_e[i]
	if this_e < al27_sab_lib.elastic_e_in[0]:
		pass
	else:
		dex = numpy.where( this_e >= al27_sab_lib.elastic_e_in )[0][-1]
		xs = al27_sab_lib.elastic_P[dex]/this_e * Nd
		al27_sab[i] = al27_sab[i] + xs
# combined
al27_total_e = numpy.concatenate((al27_sab_e, al27_e[ al27_dex:]))
al27_total   = numpy.concatenate((al27_sab,   al27_xs[al27_dex:]))+numpy.interp(al27_total_e,al27_e,al27_a)


### sapphire
al_frac = 0.4
o_frac  = 0.6
Md      = ( 26.9815385 * al_frac + 15.99491461956 * o_frac) * 1.008664916
Na      = 6.022e23
b       = 1e-24 
density = 3.98
Nd =  density * Na * b / Md
Nd_sapp = density * Na * b / Md
sapp_e     = numpy.union1d(o16_lib.energy, al27_lib.energy)
sapp_xs    = ( o_frac*numpy.interp(sapp_e,o16_lib.energy,o16_lib.reactions[2].sigma  ) + al_frac*numpy.interp(sapp_e, al27_lib.energy, al27_lib.reactions[2].sigma  )) * Nd
sapp_a     = ( o_frac*numpy.interp(sapp_e,o16_lib.energy,o16_lib.reactions[102].sigma) + al_frac*numpy.interp(sapp_e, al27_lib.energy, al27_lib.reactions[102].sigma)) * Nd
sapp_sab    = (o_frac*osapp_sab_lib.inelastic_sigma + al_frac*alsapp_sab_lib.inelastic_sigma)*Nd
sapp_sab_e  = osapp_sab_lib.inelastic_e_in
#find e dex where sab stops
sapp_dex = numpy.where( sapp_e >= sapp_sab_e[-1] )[0][0]
## add carbon to deuterium inelastic
#for i in range(0,len(sapp_sab_e)):
#	this_e =    sapp_sab_e[i]
#	if this_e < sapp_sab_lib.elastic_e_in[0]:
#		pass
#	else:
#		dex = numpy.where( this_e >= c12_lib.energy )[0][-1]
#		xs1 =  9.0*c12_lib.reactions[2].sigma[dex]* Nd
#		sapp_sab[i] = sapp_sab[i] + xs1 
# combined
sapp_total_e   = numpy.concatenate((sapp_sab_e, sapp_e[ sapp_dex:]))
sapp_total     = numpy.concatenate((sapp_sab,   sapp_xs[sapp_dex:]))+numpy.interp(sapp_total_e,sapp_e,sapp_a)



#
#   plot!
#
f=plt.figure()
ax=f.add_subplot(111)
ax.loglog(al27_total_e,      al27_total,          'b', label=r'Aluminum, Room Temp')
ax.loglog(sapp_total_e,      sapp_total,          'r', label=r'Single Crystal Sapphire, Room Temp')

ax.loglog(al27_e,      al27_a,      'b--', label='Aluminum abs')
ax.loglog(sapp_e,      sapp_a,      'r--', label='Single Crystal Sapphire abs')

plt.legend(loc=4)
ax.grid('on')
ax.set_xlabel(r'Energy (MeV)')
ax.set_ylabel(r'Macroscopic Cross Section (with S($\alpha$,$\beta$)) (cm^{-1})')
plt.show()



al27_thk = 0.7
sapp_thk = 5.0

sapp_trans = numpy.exp(-sapp_total*sapp_thk)
al27_trans = numpy.exp(-al27_total*al27_thk)
total_trans_e=numpy.union1d(al27_total_e,sapp_total_e)

total_trans = numpy.interp(total_trans_e,sapp_total_e,sapp_trans)*numpy.interp(total_trans_e,al27_total_e,al27_trans)


f=plt.figure()
ax=f.add_subplot(111)
ax.semilogx(total_trans_e,      total_trans,          'b')
ax.grid('on')
ax.set_xlabel(r'Energy (MeV)')
ax.set_ylabel(r'Transmission')
plt.show()


f=plt.figure()
ax=f.add_subplot(111)
ax.plot(to_wavelength(total_trans_e),      total_trans,          'b')
ax.grid('on')
ax.set_xlabel(r'Energy (MeV)')
ax.set_ylabel(r'Transmission')
ax.set_xlim([0,11])
plt.show()