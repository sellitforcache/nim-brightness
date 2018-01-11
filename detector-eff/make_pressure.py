#! /home/l_bergmann/anaconda/bin/python -W ignore
# /usr/local/bin/python 
#

from MCNPtools.mctal import mctal
from MCNPtools.to_wavelength import to_wavelength
from MCNPtools import calculate_materials
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

#
#  low efficiency He3 counter
#
low_eff_he3 = calculate_materials.mixture('low eff he3')
low_eff_he3.mass_density=3.485e-3
low_eff_he3.add_mixture('Kr'  , 1.2,  mode='atom')
low_eff_he3.add_mixture('He3' , 0.02, mode='atom')
low_eff_he3.finalize()
#low_eff_he3.print_material_card()

#
#  normal efficiency He3 counter
#
normal_eff_he3 = calculate_materials.mixture('normal eff he3')
normal_eff_he3.mass_density=4.468741e-3
normal_eff_he3.add_mixture('Kr'  , 1.2, mode='atom')
normal_eff_he3.add_mixture('He3' , 2.3, mode='atom')
normal_eff_he3.finalize()
#normal_eff_he3.print_material_card()

#
#  only Kr
#
#only_kr = calculate_materials.mixture('only kr det')
#only_kr.mass_density=4.468741e-3
#only_kr.add_mixture('Kr'  , 1.2, mode='atom')
#only_kr.finalize()
##only_kr.print_material_card()


print ""

pressure = 1.2+0.02-0.0007*2. #bar
temp = 20. # C
temp = temp + 273.15 # K
R = 83.144598 #cm3 bar K-1 mol-1
R_gas = R/low_eff_he3.avg_amu
low_dens= pressure/(R_gas*temp)
print "Density of 1.2 bar Kr + 0.02 bar He3 = %6.4E"%low_dens

pressure = 1.2+2.3 #bar
temp = 20. # C
temp = temp + 273.15 # K
R = 83.144598 #cm3 bar K-1 mol-1
R_gas = R/normal_eff_he3.avg_amu
normal_dens= pressure/(R_gas*temp)
print "Density of 1.2 bar Kr + 2.30 bar He3 = %6.4E"%normal_dens

#pressure = 1.2 #bar
#temp = 20. # C
#temp = temp + 273.15 # K
#R = 83.144598 #cm3 bar K-1 mol-1
#R_gas = R/only_kr.avg_amu
#only_kr_dens= pressure/(R_gas*temp)
#print "Density of 1.2 bar Kr                = %6.4E"%only_kr_dens

print ""


#
#  low efficiency He3 counter
#
low_eff_he3 = calculate_materials.mixture('low eff he3')
low_eff_he3.mass_density=low_dens
low_eff_he3.add_mixture('Kr'  , 1.2,  mode='atom')
low_eff_he3.add_mixture('He3' , 0.02, mode='atom')
low_eff_he3.finalize()
low_eff_he3.print_material_card()

#
#  normal efficiency He3 counter
#
normal_eff_he3 = calculate_materials.mixture('normal eff he3')
normal_eff_he3.mass_density=normal_dens
normal_eff_he3.add_mixture('Kr'  , 1.2, mode='atom')
normal_eff_he3.add_mixture('He3' , 2.3, mode='atom')
normal_eff_he3.finalize()
normal_eff_he3.print_material_card()

#
#  only Kr
#
#only_kr = calculate_materials.mixture('only kr det')
#only_kr.mass_density=only_kr_dens
#only_kr.add_mixture('Kr'  , 1.0, mode='atom')
#only_kr.finalize()
#only_kr.print_material_card()#