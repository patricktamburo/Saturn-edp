#The purpose of this script is to assign reconstructed SPICE .bsp files to each RSR occ
# observation. This goes through the HD finding each RSR file, then searches through the
# spk/ directory to find the matching spice file. In some cases, multiple SPICE files 
# may be chosen. It finishes by producing a text document with the correspondance.

#Import various math, science, and plotting packages.
imprt = 0
while imprt==0:
	import numpy, scipy
	from scipy.interpolate import interp1d
	import commands 
	import matplotlib
	matplotlib.use('macosx')
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import *
	from matplotlib import colors, cm
	from matplotlib.pylab import *
	from matplotlib.font_manager import FontProperties
	from mpl_toolkits.mplot3d import Axes3D
	plt.ion()
	import os, glob
	import sys, math
	from scipy import stats, signal
	from numpy import *
	from scipy import optimize
	from scipy.optimize import curve_fit
	from scipy.optimize import fsolve
	import idlsave
	import time
	import multiprocessing
	from multiprocessing import Pool
	import emcee
	import pickle
	from scipy.spatial import distance
	import itertools
	import astropy
	from astropy.io import fits
	from astropy.wcs import WCS
	from astropy.time import Time
	imprt = 1 
	
print ''
input('Wait! Do you really want to over write the RSR SPICE connection file?')


#Paths of RSRs and SPICE data
path = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Programs/'

#RSR data path
rsr_base_path = '/Volumes/PW-2TB/Cassini_Ionosphere/Data/'	
#SPICE path
spice_path = '/Volumes/PW-2TB/Cassini_Ionosphere/spk/'

#Open the output file
file = open(path+'RSR_SPICE_connection.txt','w')

#Get a list of all SPICE reconstructed files 
all_spice_bsp = sorted(glob.glob(spice_path+'*.bsp'))
all_spice_Rbsp = array([])
all_spice_Rbsp_trunc = array([])
for i in range(size(all_spice_bsp)):
	if ('R_SCPSE' in all_spice_bsp[i].split('/')[-1]) & \
		('EP1' not in all_spice_bsp[i].split('/')[-1]):
			all_spice_Rbsp = append(all_spice_Rbsp, all_spice_bsp[i].split('/')[-1])
			#Ignore the creation date for the truncated one
			all_spice_Rbsp_trunc = append(all_spice_Rbsp_trunc, \
				all_spice_bsp[i].split('/')[-1][6:])
			

#Sort the reconstructed files by starting date
all_spice_Rbsp = all_spice_Rbsp[argsort(all_spice_Rbsp_trunc)]

#Determine the start and end times for the SPICE from the headers 
spice_start, spice_end = zeros(size(all_spice_Rbsp)), zeros(size(all_spice_Rbsp))
for i in range(size(all_spice_Rbsp)):
	#Open the header
	spice_file = open(spice_path+all_spice_Rbsp[i]+'.lbl','r')
	while 1:
		#Look for the start time line
		line = spice_file.readline()
		if 'START_TIME' in line:
			#Get the corresponding time
			spice_start[i] = Time(line.split('=')[-1][1:-1],format='isot').jd
			line = spice_file.readline()
			spice_end[i] = Time(line.split('=')[-1][1:-1],format='isot').jd
			break
	spice_file.close()		

#Now begin loop over all RSR files to see which RSR file(s) they correspond to
all_cors_files = sorted(glob.glob(rsr_base_path+'cors*'))
for i in range(size(all_cors_files)):
	#All local files
	all_files = glob.glob(all_cors_files[i]+'/*')
	for j in range(size(all_files)):
		prod_id = all_files[j].split('/')[-1]
		#Skip the label files
		if ('.lbl' in prod_id) | ('.LBL' in prod_id): continue
		
		#Otherwise, get jd date of the observation and write the product id to 
		# the save file
		file.write(prod_id+',')
		prod_jd = Time(prod_id[7:11]+':'+prod_id[11:14]+':'+prod_id[15:17]+':'+\
			prod_id[17:19]+':00',format='yday').jd
		ind = where((prod_jd>spice_start)&(prod_jd<spice_end))[0][0]
		file.write(all_spice_Rbsp[ind]+',')
		
		#The rsr observation has some length. If the prod_jd is within 12 hours of the 
		# edge of a spice file, save the neighboring file as well
		if spice_start[ind+1]-prod_jd < 0.5:
			file.write(all_spice_Rbsp[ind+1]+',')
		if prod_jd - spice_end[ind-1] < 0.5:
			file.write(all_spice_Rbsp[ind-1]+',')	
		
		#Next line
		file.write(os.linesep)
file.close()