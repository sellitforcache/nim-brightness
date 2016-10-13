#! /home/l_bergmann/anaconda/bin/python -W ignore

from pyne import mcnp, ace
import math
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
	print "detector system limit",sa_finite_d," self limit ", sa_collimator
	return numpy.minimum(sa_finite_d,sa_collimator)



# conversion factors
charge_per_amp = 6.241e18
charge_per_milliamp = charge_per_amp/1e3
charge_per_microamp = charge_per_amp/1e6
Na     = 6.0221409e+23  # number/mol



# tallies
inv = mctal('/home/l_bergmann/repos/temp/ICON-PDarraz.mctal')
#ike19 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-19IKE.mctal')
#ike24 = mctal('/home/l_bergmann/Documents/nim-brightness/ICON-PDarray-24IKE.mctal')
new = mctal('/home/l_bergmann/repos/temp/ICON-center-19std.mctal')
stip= mctal('/home/l_bergmann/repos/temp/ICON-center-19std-STIP.mctal')
deg24 = mctal('/home/l_bergmann/repos/temp/ICON-sph-stip-24kden.mctal')
deg25 = mctal('/home/l_bergmann/repos/temp/ICON-sph-stip-25kden.mctal')
deg28 = mctal('/home/l_bergmann/repos/temp/ICON-sph-stip-28kden.mctal')
this_tal = stip
#std19.plot(tal=[5],obj=[8],options=['log','lethargy'])

cd_radius = 0.050

### radiography tally
rad_tal=235
rad_non_zero = 0
tol=5e-8
img_accepted=0.0
count_num=0
img = numpy.zeros((len(this_tal.tallies[rad_tal].cosines),len(this_tal.tallies[rad_tal].segments) ))
print "cosines ",this_tal.tallies[rad_tal].cosines
print "segments ",this_tal.tallies[rad_tal].segments
for y in range(0,len(this_tal.tallies[rad_tal].cosines)-1):
	for x in range(0,len(this_tal.tallies[rad_tal].segments)-1):
		dex = this_tal.tallies[rad_tal]._hash(obj=0,seg=x,cos=y)
		img[y][x]=numpy.sum(this_tal.tallies[rad_tal].vals[dex]['data'][-1])
		#
		x_val = numpy.amin(   numpy.absolute( [this_tal.tallies[rad_tal].segments[x],this_tal.tallies[rad_tal].segments[x+1]] )  )
		y_val = numpy.amin(   numpy.absolute( [this_tal.tallies[rad_tal].cosines[ y],this_tal.tallies[rad_tal].cosines[ y+1]] )  )
		#
		if img[y][x]>= tol:
			rad_non_zero = rad_non_zero + 1
		if numpy.sqrt(x_val*x_val+y_val*y_val)<=cd_radius:
			img_accepted = img_accepted + img[y][x]
			count_num = count_num+1
img_total=numpy.sum(img)
img_avg  =img_total/count_num

correction_factor_cd = img_accepted/img_total
d=19.681240460
pix_str = 2.0*numpy.pi*(1.0-d/numpy.sqrt(d*d+0.0005*0.0005))
total_str=pix_str*rad_non_zero
print "   ========  "
print "pixels >= %6.4E : %d"%(tol,rad_non_zero)
print "pixels worth approx %6.4E str"%pix_str
print "total view approx %6.4E str"%total_str
print "correction factor for Cd acceptance: %6.4E"%correction_factor_cd

#### formatter for scientific notation
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

#  plot it
f=plt.figure()
ax=f.add_subplot(111)
ax.set_title('Pinhole Image at Cd/Zapfen System Crossing Point')
ax.set_xlabel(r'x (cm)')
ax.set_ylabel(r'y (cm)')
imgax=ax.imshow(img,interpolation='none',aspect='auto',origin='lower',extent=[this_tal.tallies[rad_tal].cosines[0],this_tal.tallies[rad_tal].cosines[-1],this_tal.tallies[rad_tal].segments[0],this_tal.tallies[rad_tal].segments[-1]])
c=plt.colorbar(imgax, format=ticker.FuncFormatter(fmt))
jet = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=1)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
ax.set_axis_bgcolor(scalarMap.to_rgba(0))
c.set_label(r'Neutron Flux (n cm$^{-2}$ p$^{-1}$)')


### plot limits
# circle for Cd hole
theta  = numpy.linspace(0,2*numpy.pi,512)
circ_x = cd_radius*numpy.cos(theta)
circ_y = cd_radius*numpy.sin(theta)
ax.plot(circ_x,circ_y,color=numpy.array([155.0,155.0,155.0])/255.0,linewidth=4,linestyle='--')
# rectangle for Zapfen
la=19.681240460
b =4.0
l =1577.5
a =la*b/(l-la)
zap_x = [-a , a, a , -a, -a]
b =6.0
l =1577.5
a =la*b/(l-la)
zap_y = [-a , -a , a , a, -a]
ax.plot(zap_x,zap_y,color=numpy.array([255.0,0.0,255.0])/255.0,linewidth=4,linestyle='--')
# rectangle for steel
b =4.0
l =1510.1
a =la*b/(l-la)
steel_x = [-a , a, a , -a, -a]
b =4.0
l =1510.1
a =la*b/(l-la)
steel_y = [-a , -a , a , a, -a]
ax.plot(steel_x,steel_y,color=numpy.array([0.0,255.0,255.0])/255.0,linewidth=4,linestyle='--')
#circle for aperture
b =4.0
l =1230.8
a =la*b/(l-la)
ap_x = a*numpy.cos(theta)
ap_y = a*numpy.sin(theta)
ax.plot(ap_x,ap_y,color=numpy.array([255.0,255.0,0.0])/255.0,linewidth=4,linestyle='--')





ax.set_ylim([-0.08,0.08])
ax.set_xlim([-0.08,0.08])
plt.show()

# 1d model efficiency
oned_eff = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/1deff.csv','r')
for line in f:
	nums = line.split(',')
	oned_eff[0].append(float(nums[0]))
	oned_eff[1].append(float(nums[1]))
f.close()

# mcnp calculated efficieny
mcnp_eff = [[],[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/mcnpeff.csv','r')
for line in f:
	nums = line.split(',')
	print nums
	mcnp_eff[0].append(float(nums[0]))
	mcnp_eff[1].append(float(nums[1]))
	mcnp_eff[2].append(float(nums[2]))
f.close()

# measured efficiency
meas_eff = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/measeff.csv','r')
for line in f:
        nums = line.split(',')
        meas_eff[0].append(float(nums[0]))
        meas_eff[1].append(float(nums[1]))
f.close()


# plot efficiency
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(meas_eff[0],meas_eff[1],  linewidth=2,label='BOA Measurement')
ax.plot(oned_eff[0],oned_eff[1],  linewidth=2,label='1-D Model')
ax.plot(mcnp_eff[0],mcnp_eff[1],  linewidth=2,label='MCNP Model')
ax.fill_between(mcnp_eff[0],numpy.multiply(mcnp_eff[1],1.0+numpy.array(mcnp_eff[2])),numpy.multiply(mcnp_eff[1],1.0-numpy.array(mcnp_eff[2])), facecolor='red', linewidth=1.0, color='red', alpha=0.25,label=r'MCNP 1-$\sigma$ Error')
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Efficiency')
ax.set_ylim([0,1])
ax.set_xlim([0,12])
ax.grid(1)
plt.legend(loc=4)
plt.show()


# final measurement
measurement = [[],[]]
f=open('/home/l_bergmann/Documents/nim-brightness/brightness_measurement.csv','r')
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
tal_num=5
ene1  = numpy.array(this_tal.tallies[tal_num].energies[:-1])
dex  = this_tal.tallies[tal_num]._hash(obj=0)
tal1 = this_tal.tallies[tal_num].vals[dex]['data'][:-1]

dex  = new.tallies[tal_num]._hash(obj=0)
tal2 = new.tallies[tal_num].vals[dex]['data'][:-1]

dex  = inv.tallies[tal_num]._hash(obj=7)
tal3 = inv.tallies[tal_num].vals[dex]['data'][:-1]


ene4  = numpy.array(deg24.tallies[tal_num].energies[:-1])
dex  = deg24.tallies[tal_num]._hash(obj=0,cos=0)
tal4 = deg24.tallies[tal_num].vals[dex]['data'][:-1]

dex  = deg25.tallies[tal_num]._hash(obj=0,cos=0)
tal5 = deg25.tallies[tal_num].vals[dex]['data'][:-1]

dex  = deg28.tallies[tal_num]._hash(obj=0,cos=0)
tal6 = deg28.tallies[tal_num].vals[dex]['data'][:-1]


wvl1 = to_wavelength(ene1)
widths1 = -1.0*numpy.diff(wvl1)
wvl4 = to_wavelength(ene4)
widths4 = -1.0*numpy.diff(wvl4)
#this_tal.plot(tal=[5],obj=[7])

## compare PD and radiography
print "   ========  "
print "PD total       %6.4E"%this_tal.tallies[5].vals[dex]['data'][-1]
print "Pinhole mean   %6.4E"%(numpy.mean(img))
print "Pinhole sum    %6.4E"%(numpy.sum(img))

### solid angles
sa_collimator = get_view_sa(b=0.1,  l=10.1,   thk=10.0)
sa_cadmium    = get_view_sa(b=0.05, l= 0.1,   thk= 0.1)
sa_sapphire   = get_view_sa(b=1.0,  l=41.835, thk=25.0)
sa_steel      = get_view_sa(b=2.5,  l=784.8,  thk=30.0)
sa_aperture   = get_view_sa(b=4.0,  l=1230.8, thk=0.1)
sa_tube       = get_view_sa(b=4.0,  l=1510.1, thk=271.8)
sa_zapfen     = get_view_sa(b=4.0,  l=1577.5, thk=450.0)
sa_nozzle     = get_view_sa(b=6.9,  l=1720.0, thk=0.1)

print "  HORIZONTAL ========  "
print "solid angle collimator     %6.4E"%sa_collimator
print "solid angle cadmium        %6.4E"%sa_cadmium
print "solid angle sapphire       %6.4E"%sa_sapphire
print "solid angle steel          %6.4E"%sa_steel
print "solid angle aperture       %6.4E"%sa_aperture
print "solid angle tube           %6.4E"%sa_tube
print "solid angle zapfen         %6.4E"%sa_zapfen
print "solid angle nozzle         %6.4E"%sa_nozzle


### solid angles
sa_collimator = get_view_sa(b=0.1,  l=10.1,   thk=10.0)
sa_cadmium    = get_view_sa(b=0.05, l= 0.1,   thk= 0.1)
sa_sapphire   = get_view_sa(b=1.0,  l=41.835, thk=25.0)
sa_steel      = get_view_sa(b=2.5,  l=784.8,  thk=30.0)
sa_aperture   = get_view_sa(b=4.0,  l=1230.8, thk=150.5)
sa_tube       = get_view_sa(b=4.0,  l=1510.1, thk=271.8)
sa_zapfen     = get_view_sa(b=6.0,  l=1577.5, thk=250.0)
sa_nozzle     = get_view_sa(b=6.9,  l=1720.0, thk=0.1)

print "  VERTICAL ========  "
print "solid angle collimator     %6.4E"%sa_collimator
print "solid angle cadmium        %6.4E"%sa_cadmium
print "solid angle sapphire       %6.4E"%sa_sapphire
print "solid angle steel          %6.4E"%sa_steel
print "solid angle aperture       %6.4E"%sa_aperture
print "solid angle tube           %6.4E"%sa_tube
print "solid angle zapfen         %6.4E"%sa_zapfen
print "solid angle nozzle         %6.4E"%sa_nozzle


### solid angles
sa_collimator = get_view_sa(b=0.1,                          l=10.1,   thk=10.0)
sa_cadmium    = get_view_sa(b=0.05,                         l= 0.1,   thk= 0.1)
sa_sapphire   = get_view_sa(b=numpy.sqrt(2)*1.0,            l=41.835, thk=25.0)
sa_steel      = get_view_sa(b=numpy.sqrt(2)*2.5,            l=784.8,  thk=30.0)
sa_aperture   = get_view_sa(b=4.0,                          l=1230.8, thk=0.1)
sa_tube       = get_view_sa(b=numpy.sqrt(2)*4.0,            l=1510.1, thk=271.8)
sa_zapfen     = get_view_sa(b=numpy.sqrt(4.0*4.0+6.0*6.0),  l=1577.5, thk=450.0)
sa_nozzle     = get_view_sa(b=numpy.sqrt(2)*6.9,            l=1720.0, thk=0.1)

print "  DIAGONAL ========  "
print "solid angle collimator     %6.4E"%sa_collimator
print "solid angle cadmium        %6.4E"%sa_cadmium
print "solid angle sapphire       %6.4E"%sa_sapphire
print "solid angle steel          %6.4E"%sa_steel
print "solid angle aperture       %6.4E"%sa_aperture
print "solid angle tube           %6.4E"%sa_tube
print "solid angle zapfen         %6.4E"%sa_zapfen
print "solid angle nozzle         %6.4E"%sa_nozzle


tal1_normed = charge_per_milliamp*numpy.divide(tal1,widths1*total_str)
tal2_normed = charge_per_milliamp*numpy.divide(tal2,widths1*total_str)
tal3_normed = charge_per_milliamp*numpy.divide(tal3,widths1*total_str)

tal4_normed = charge_per_milliamp*numpy.divide(tal4,widths4*total_str)
tal5_normed = charge_per_milliamp*numpy.divide(tal5,widths4*total_str)
tal6_normed = charge_per_milliamp*numpy.divide(tal6,widths4*total_str)
sa_measure  = 6.28E-4 #2.0*numpy.pi*(1.0-101./numpy.sqrt(101.0*101.0+1.5*1.5))
measurement[1] = numpy.array(measurement[1])

# plot spectra
f=plt.figure()
ax=f.add_subplot(111)
ax.plot(measurement[0][:],numpy.multiply(measurement[1][:],sa_measure/total_str),linewidth=2,label='Measurement')
make_steps(ax,wvl1,[0],tal1_normed,linewidth=2,label='MCNP - STIP,  SPHERICAL',options=['lin'])
make_steps(ax,wvl1,[0],tal2_normed,linewidth=2,label='MCNP - NOSTIP, INVERTED',options=['lin'])
make_steps(ax,wvl1,[0],tal3_normed,linewidth=2,label='MCNP - NOSTIP, INVERTED, MISALIGNED',options=['lin'])
make_steps(ax,wvl4,[0],tal4_normed,linewidth=2,label='MCNP - STIP,  SPHERICAL, 24K density',options=['lin'])
make_steps(ax,wvl4,[0],tal5_normed,linewidth=2,label='MCNP - STIP,  SPHERICAL, 25K density',options=['lin'])
make_steps(ax,wvl4,[0],tal6_normed,linewidth=2,label='MCNP - STIP,  SPHERICAL, 28K density',options=['lin'])
ax.set_xlabel(r'Wavelength (\AA)')
ax.set_ylabel(r'Brilliance (n cm$^{-2}$ s$^{-1}$ mA$^{-1}$ \AA$^{-1}$ str$^{-1}$)')
#ax.set_ylim([0,5e11])
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
