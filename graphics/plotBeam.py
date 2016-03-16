#! /usr/bin/env python
# Script to plot proton history on the SINQ target
# Ryan M. Bergmann, Nov 17, 2014.  ryan.bergmann@psi.ch, ryanmbergmann@gmail.com

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.colors import LogNorm, PowerNorm
import sys
import numpy
import scipy as sp

### set TeX
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)

def fill_between_steps(x, y1, y2=0, h_align='post', ax=None, **kwargs):
    ''' Fills a hole in matplotlib: fill_between for step plots.

    Parameters :
    ------------

    x : array-like
        Array/vector of index values. These are assumed to be equally-spaced.
        If not, the result will probably look weird...
    y1 : array-like
        Array/vector of values to be filled under.
    y2 : array-Like
        Array/vector or bottom values for filled area. Default is 0.

    **kwargs will be passed to the matplotlib fill_between() function.

    '''
    # If no Axes opject given, grab the current one:
    if ax is None:
        ax = plt.gca()
    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = sp.repeat((x[1:] - x[:-1]), 2)
    xstep = sp.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = sp.append(xx, xx.max() + xstep[-1])

    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)#[:-1]
    if type(y2) == sp.ndarray:
        y2 = y2.repeat(2)#[:-1]

    # now to the plotting part:
    ax.fill_between(xx, y1, y2=y2, **kwargs)

    return ax


fname = sys.argv[1]
datf=open(fname)
time_str=[]
date_str=[]
current=[]
frac=[]
z=[]
time_from_zero=[]
t=0.

asp=0.02
fig=plt.figure()
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)

# recombination parameters
km=5.5e-7  #per milliamp!
b12_d2=0.9e-7
bstar=2.91e-4
c0_eq=0.955
d2_o=0.666

# get beam data
step_size = 10.0*60.0 #5.0   # seconds
index=0
measurement_index=0
for line in datf:
	items=line.split()
	time_str.append(items[5])
	date_str.append(items[4])
	if date_str[-1] == '21-07-2014' and time_str[-1] == '22:00:00':
			measurement_index = index
	current.append(float(items[2]))
	z.append(0.0)
	time_from_zero.append(t)
	t = t + step_size
	# calculate the concentration
	k = km*current[-1]/1000.0
	A = -( k + bstar*numpy.sqrt(k) + b12_d2 )
	B =  ( bstar*numpy.sqrt(k)*c0_eq + b12_d2*c0_eq + 2./3.*k )
	d2_o = (A*d2_o + B)*step_size + d2_o
	frac.append(d2_o)
	#
	index = index + 1

adjusted_time = numpy.array(time_from_zero)/(60.*60.*24.*7.)

fill_between_steps(adjusted_time,numpy.array(z),numpy.array(current),ax=ax2,color='g')
ax1.plot(          adjusted_time,numpy.array(frac),color='b',linewidth=2)

ylim=ax2.get_ylim()
xlim=ax2.get_xlim()

ax1.grid()
ax1.set_ylim([0.60,1.00])
ax1.set_xlim([0,xlim[1]])
#ax1.set_xlabel("Day of July, 2014")
ax1.set_ylabel(r"o-D$_2$ Fraction")
#ax1.set_aspect(2.0*asp/8.0)

ax2.grid()
ax2.set_ylim([ylim[0]*1.1,max(current)*1.1])
ax2.set_xlim([0,xlim[1]])
ax2.set_xlabel('Weeks Since %s'%date_str[0])
ax2.set_ylabel(r"SINQ Proton Current ($\mu$A)")
#ax2.set_aspect(2.0*asp/8.0)

ax1.annotate(r'\noindent July 21, 2014, o-D$_2$ = %4.3f'%frac[measurement_index],     
	                              xy=(adjusted_time[measurement_index],     frac[measurement_index]    ), 
	                          xytext=(adjusted_time[measurement_index]*0.95, frac[measurement_index]*1.13), 
	                          verticalalignment='bottom', horizontalalignment='right',
	                          arrowprops=dict(facecolor='black', shrink=0.05))

xlim=ax2.get_xlim()
ax2.set_xticks(numpy.arange(min(xlim), max(xlim)+1, 1.0))
xlim=ax1.get_xlim()
ax1.set_xticks(numpy.arange(min(xlim), max(xlim)+1, 1.0))
#fig.savefig('p_current.png',dpi=300)
plt.show()

