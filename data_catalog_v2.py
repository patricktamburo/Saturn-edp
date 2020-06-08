#The purpose of this script is to look through the raw data the has been placed into the
# data directory after the final download. I will extract as much data from the header and
# filename is possible. I will import the output file into my excel spreadsheet that
# contains the information about all the cassini rsr occs.

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

#print 'stop! are you sure you want to overwrite the catalog files?'
#stop = input('   ')
		
#Set the base path
raw_path = '/Volumes/PW-2TB/Cassini_Ionosphere/Data/'
spice_path = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Programs/'

#Open the save files for each body/occ type
save_path = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Misc/'+\
	'automated_data_catalog_output/v2/'
saturn_save_file = open(save_path+'saturn.txt','w')
titan_save_file = open(save_path+'titan.txt','w')
enceladus_save_file = open(save_path+'enceladus.txt','w')
rioc_save_file = open(save_path+'rioc.txt','w')

#Enter the list of cors files for this body to analyze. This can be manually entered to 
# do a single one, or it will go into the raw data files and fetch all the cors numbers
cors = []
cors_files = sorted(glob.glob(raw_path+'cors*'))
for i in range(size(cors_files)):
	cors.append(cors_files[i].split('_')[-1])


#Now start the loop over the cors files
for l in range(size(cors)):
	#Set the cors path
	cors_path = raw_path+'cors_'+cors[l]+'/'
	
	#First, look for the label files. Case of the text filenames may be upper or lower.
	# Set a keyword that will be helpful later. Hopefully each cors file is consistently
	# upper or lower case.
	caps = 'no'
	suf = '.lbl'
	label_files = glob.glob(cors_path+'*.lbl')
	if size(label_files)==0: 
		label_files = glob.glob(cors_path+'*.LBL')
		caps = 'yes'
		suf = '.LBL'
	
	#Get a list of all the files in this directory, then delete the ones that are .lbl or
	# .LBL. Add a quick check to make sure there is an even number of files
	all_rsr_files = glob.glob(cors_path+'*')
	if mod(size(all_rsr_files),2) != 0: input('Break here, N RSR files not even!')
	rsr_files = []
	for j in range(size(all_rsr_files)):
		if suf in all_rsr_files[j]: continue
		rsr_files.append(all_rsr_files[j])
	
	
	
	#Now, loop over all the RSR only files, extracting the necessary info
	for i in range(size(rsr_files)):
		#First, get the product id
		prod_id = rsr_files[i].split('/')[-1]
	
		#Get other info out of the filename
		seq = prod_id[1:3]
		activity = prod_id[3:7]
		yyyy,doy,hh,mm = prod_id[7:11],prod_id[11:14],prod_id[15:17],prod_id[17:19]
		upband = prod_id[19]
		upstn = prod_id[20:22]
		downband = prod_id[22]
		downstn = prod_id[23:25]
		pol = prod_id[25]
		rsr = prod_id[28:31]

		#Go into the corresponding label file to get a few pieces of info
		label_file = open(rsr_files[i].split('.')[0]+suf,'r')
		while 1:
			#Read the file line by line
			line = label_file.readline()
			if 'DATA_SET_ID' in line:
				dataset = line.split('"')[-2]
			if 'START_TIME' in line:
				starttime = line.split('=')[-1][:-2]
			if 'STOP_TIME' in line:
				stoptime = line.split('=')[-1][:-2]
				break	
		label_file.close()
		
		#The last piece of information is the corresponding .bsp files. For this, search
		# the RSR SPICE connection file
		spice_file = open(spice_path+'RSR_SPICE_connection.txt','r')
		while 1:
			line = spice_file.readline()
			if prod_id in line:
				bsp_files = line.split(',')[1:-1]
				break
		spice_file.close()
		#If more than 1 spice file, join them for printing
		if size(bsp_files) > 1:
			bsp_files = ','.join(bsp_files)
		else:
			bsp_files = bsp_files[0]	
		
		#Now save everything in the right order, including the unknown info
		final_line = ','.join([dataset,cors[l],prod_id,seq,'???','???',activity,\
			'???','???',yyyy,doy,hh,mm,starttime,stoptime,upband,upstn,downband,\
				downstn,pol,rsr,bsp_files])
		
		#Determine what body/type of occ this one is for the saving
		if (prod_id[3] == 's') | (prod_id[3] == 'S'): 
			saturn_save_file.write(final_line+os.linesep)
		if (prod_id[3] == 't') | (prod_id[3] == 'T') | (prod_id[3] == 'P'): 
			titan_save_file.write(final_line+os.linesep)
		if (prod_id[3] == 'e') | (prod_id[3] == 'E'): 
			enceladus_save_file.write(final_line+os.linesep)
		if (prod_id[3] == 'r') | (prod_id[3] == 'R'): 
			rioc_save_file.write(final_line+os.linesep)
	
#Close all files
saturn_save_file.close()
titan_save_file.close()	
enceladus_save_file.close()
rioc_save_file.close()
	