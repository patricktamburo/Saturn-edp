#The purpose of this script is to look through the raw data downloaded straight from the
# PDS and relocate the RSR files into a different directory in the data storage area. This
# removes all the nested directories and puts RSR files directly in a cors directory.

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
		
#Bodies could be 'titan', 'rioc', 'saturn', or 'enceladus'
body = 'rioc'   	

#Set the base path - this will differ for the 2 TB or the SCC
raw_path = '/Volumes/PW-2TB/'+body+'_temp/atmos.nmsu.edu/PDS/data/'
save_path = '/Volumes/PW-2TB/Cassini_Ionosphere/Data/'

#Enter the list of cors files for this body to analyze. This can be manually entered to 
# do a single one, or it will go into the raw data files and fetch all the cors numbers
cors = []
cors_files = glob.glob(raw_path+'cors*')
for i in range(size(cors_files)):
	cors.append(cors_files[i].split('_')[-1])



#No start the loop over the cors files
for l in range(size(cors)):
	#Make the new cors directory for the rsr files
	new_cors_path = save_path+'cors_'+cors[l]+'/'
	os.system('mkdir '+new_cors_path)


	#Set the cors path
	cors_path = raw_path+'cors_'+cors[l]+'/'



	#First, look for the text files that contain the dataset name. Case of the text filenames
	# may be uppper or lower. Set a keyword that will be helpful later. Hopefully each cors
	# file is consistently upper or lower case
	caps = 'no'
	rsr_fmt = 'rsr/'
	txt_files = glob.glob(cors_path+'*me.txt')
	if size(txt_files)==0: 
		txt_files = glob.glob(cors_path+'*ME.TXT')
		caps = 'yes'
		rsr_fmt = 'RSR/'
	
	
	#Now locate the rsr files in the occ files. There may be more than one, so a big loop must
	# begin here
	occ_files = glob.glob(cors_path+'*_*')
	for i in range(size(occ_files)):
		
		#Move all the files in the rsr/RSR files to the new cors path
		os.system('mv '+occ_files[i]+'/'+rsr_fmt+'* '+new_cors_path)
		
		#stop = input('stop here')
		