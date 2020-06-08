#The purpose of this script is to look through the raw data downloaded straight from the
# PDS and pull out as much data from the directory structure and file name is possible.
# This data will then be saved to a text file that I will import into my excel spreadsheet
# that contains the information about all the cassini rsr occs.

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

print 'stop! are you sure you want to overwrite the catalog files?'
stop = input('   ')
		
#Bodies could be 'titan', 'rioc', 'saturn', or 'enceladus'
body = 'rioc'   	

#Set the base path
raw_path = '/Volumes/PW-2TB/'+body+'_temp/atmos.nmsu.edu/PDS/data/'

#Enter the list of cors files for this body to analyze. This can be manually entered to 
# do a single one, or it will go into the raw data files and fetch all the cors numbers
cors = []
cors_files = glob.glob(raw_path+'cors*')
for i in range(size(cors_files)):
	cors.append(cors_files[i].split('_')[-1])

#No start the loop over the cors files
for l in range(size(cors)):
	#Save path
	save_path = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Misc/'+\
		'automated_data_catalog_output/v1/'+body+'/'
	save_file = open(save_path+'cors_'+cors[l]+'.txt','w')



	#Set the cors path
	cors_path = raw_path+'cors_'+cors[l]+'/'



	#First, look for the text files that contain the dataset name. Case of the text filenames
	# may be uppper or lower. Set a keyword that will be helpful later. Hopefully each cors
	# file is consistently upper or lower case
	caps = 'no'
	txt_files = glob.glob(cors_path+'*me.txt')
	if size(txt_files)==0: 
		txt_files = glob.glob(cors_path+'*ME.TXT')
		caps = 'yes'
	
	#Open the text file and search for the dataset name using readline
	file = open(txt_files[0],'r')
	while(1):
		line = file.readline()
		if ('(' in line) & (')' in line):
			dataset = line.split('(')[-1].split(')')[0]
			break
	file.close()





	#Now locate the rsr files in the occ files. There may be more than one, so a big loop must
	# begin here
	occ_files = glob.glob(cors_path+'*_*')
	for i in range(size(occ_files)):
		#Grab the directory DOY
		dir_doy = occ_files[i].split('_')[-1]
	
		#Now get down to the RSR path
		if caps=='no': rsr_path = occ_files[i]+'/rsr/'
		if caps=='yes': rsr_path = occ_files[i]+'/RSR/'
	
		#Get a list of all the files in this directory, then delete the ones that are .lbl or
		# .LBL. The capitalization varies.
		all_rsr_files = glob.glob(rsr_path+'*')
		rsr_files = []
		for j in range(size(all_rsr_files)):
			if ('.lbl' in all_rsr_files[j]) | ('.LBL' in all_rsr_files[j]): continue
			rsr_files.append(all_rsr_files[j])
	
		#Now with the rsr files of interest, begin extracting and saving info
		for k in range(size(rsr_files)):
			rsr_file = rsr_files[k].split('/')[-1]
			seq = rsr_file[1:3]
			year = rsr_file[7:11]
			doy = rsr_file[11:14]
			hr = rsr_file[15:17]
			min = rsr_file[17:19]
			uplink = rsr_file[19:22]
			band = rsr_file[22]
			pol = rsr_file[25]
			dsn = rsr_file[23:25]
			rsr = rsr_file[28:31]
			#The columns are a little different for the moons (there's an extra flyby number
			# the does not occur for saturn (at least in the spreadsheet presently).
			if (body=='titan') | (body=='enceladus'):
				save_file.write(dataset+'     '+cors[l]+'     '+dir_doy+'     '+seq+'     '+\
					'???     ???     ???     ???     ???     '+year+'     '+doy+'     '+hr+\
						'     '+min+'     '+uplink+'     '+band+'     '+pol+'     '+dsn+\
							'     ???     '+rsr+os.linesep)
			if (body=='saturn') | (body=='rioc'):
				save_file.write(dataset+'     '+cors[l]+'     '+dir_doy+'     '+seq+'     '+\
					'???     ???     ???     ???     '+year+'     '+doy+'     '+hr+\
						'     '+min+'     '+uplink+'     '+band+'     '+pol+'     '+dsn+\
							'     ???     '+rsr+os.linesep)
				
		#stop = input('stop')
	save_file.close()
	
#Lastly, create a file that has all of the occs include for easy importing into excel
save_files = glob.glob(save_path+'*.txt')
combined_file = open(save_path+'all_occs.txt','w')
for i in range(size(save_files)):
	whole_file = open(save_files[i],'r')
	combined_file.write(whole_file.read())
	whole_file.close()
combined_file.close()
	