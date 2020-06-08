#produce_PDS_products_v1.py

#Paul Dalba, Boston University

#This code is used to produce the data files that are sent to the PDS. The levels of the 
# data products are as follows:
#   Level 0: (freq)     = Frequency time series and associated variables
#   Level 1: (edp_ind)  = Individual EDP and associated variables
#   Level 2: (edp_avg)  = Averaged EDP and associated variables


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
	import ast
	imprt = 1 
close('all')

#Set paths
root_path          = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/'
path_lvl0          = root_path+'Output/'
path_lvl1          = root_path+'EDP/'
save_path_lvl0     = root_path+'PDS_data/freq/'
save_path_lvl0_pt2 = root_path+'PDS_data/other_freq/'
save_path_lvl1     = root_path+'PDS_data/edp_ind/'
save_path_lvl2     = root_path+'PDS_data/edp_avg/'
save_path_summ     = root_path+'PDS_data/'


#Level 0 save procedure: Given a body (i.e., 'Titan'), the code will search through all
# Titan RSR files and identify just those associated with the atmospheric radio occs we
# have used (i.e., not bistatic RSR files). For each appropriate RSR file, it will create
# a save file and a 'cors' directory in the appropriate save directory with the csv data.
#-----------------------------------------------------------------------------------------
if False:
	body = 'Titan'    #only option here is Titan
	cors_files = glob.glob(path_lvl0+body+'/cors_*')
	#Loop over all cors directories
	for i in range(size(cors_files)):
		cors_file = cors_files[i].split('/')[-1]
		#Get a list of all RSR directories within this cors file
		rsr_files = glob.glob(cors_files[i]+'/*')
		#Loop over all RSR directories
		for j in range(size(rsr_files)):
			rsr_file = rsr_files[j].split('/')[-1]
			#Check to see if this RSR directories was used to make an EDP. If so, there 
			# will be a flag_indices.txt file in its full_output directory.
			if size(glob.glob(rsr_files[j]+'/full_output/flag_indices.txt')) != 1: 
				continue
			
			#Otherwise, create a cors directory for this RSR file in the save path (check
			# to see if one is there already). 
			save_path = save_path_lvl0+cors_file+'/'
			if size(glob.glob(save_path)) == 0: os.system('mkdir '+save_path)
			
			#Now read in the index flag info. Since single rsr files can have both ingress
			# and egress, this file may contain either 1 or 2 lines to include ingress (0) 
			# only, egress (1) only, or both (0 and 1). 
			index_start, index_stop, geom = loadtxt(rsr_files[j]+\
				'/full_output/flag_indices.txt',skiprows=1,usecols=[0,1,2],unpack=True)
			#Convert these index arrays to integers (as they are read in as floats).
			index_start, index_end, geom = index_start.astype(int),index_stop.astype(int),\
				geom.astype(int)
			
			#Now read in the IDL data
			idl_data = idlsave.read(rsr_files[j]+'/full_output/output.sav',verbose=False)
			
			#Open the file in which to save the output data
			save_filename = '_'.join(rsr_file.split('.'))+'_freq_v01.csv'
			save_filename = save_filename.upper()
			save_file = open(save_path+save_filename,'w')
			
			#Generate a string that comprises this line in the save file
			for k in range(size(idl_data.sfduyearout)):
				s = ''
				s += '%10.i' % idl_data.sfduyearout[k] +','
				s += '%10.i' % idl_data.sfdudoyout[k] +','
				s += '%20.3f' % idl_data.sfdusecout[k] +','
				s += '%10.i' % idl_data.rftoifmhzout[k] +','
				s += '%10.i' % idl_data.ddclomhzout[k] +','
				s += '%20.12e' % idl_data.ncofreqout[k] +','
				s += '%20.12e' % idl_data.kposaaout[k] +','
				s += '%20.12e' % abs(idl_data.igcomplexaaout[k]) +','
				#Quality flag columns. Different calls based on structure of flag_indices
				# file. 
				#Only ingress in this rsr file
				if (size(geom)==1)&((geom==0).all()):
					if (k>=index_start)&(k<=index_end):
						s += '%5.1s' % '1'
					else:
						s += '%5.1s' % '0'
					s += ','
					s += '%5.1s' % '0'
				#Only egress in this rsr file
				if (size(geom)==1)&((geom==1).all()):
					s += '%5.1s' % '0'
					s += ','
					if (k>=index_start)&(k<=index_end):
						s += '%5.1s' % '1'
					else:
						s += '%5.1s' % '0'		
				#Both ingress and egress in this rsr file
				if size(geom)==2:
					ingress_ind = where(geom==0)[0][0]
					egress_ind = where(geom==1)[0][0]
					if (k>=index_start[ingress_ind])&(k<=index_end[ingress_ind]):
						s += '%5.1s' % '1'
					else:
						s += '%5.1s' % '0'
					s += ','	
					if (k>=index_start[egress_ind])&(k<=index_end[egress_ind]):
						s += '%5.1s' % '1'
					else:
						s += '%5.1s' % '0'		
				if k<(size(idl_data.sfduyearout)-1): s += '\n'		
				save_file.write(s)
			save_file.close()
#-----------------------------------------------------------------------------------------
			
			


#Level 1 save procedure: Given a body (e.g., 'Titan'), the code will comb through all
# Titan EDP directories and create the desired output file in the PDS save directory.
#-----------------------------------------------------------------------------------------
if False:
	body = 'Titan'
	occ_files = glob.glob(path_lvl1+body+'/*')
	#Loop over all occ directories
	for i in range(size(occ_files)):
		occ_file = occ_files[i].split('/')[-1]
		#Get a list of all EDP directories within this occ directory
		edp_files = glob.glob(occ_files[i]+'/'+occ_file+'*')
		#Loop over all EDP directories
		for j in range(size(edp_files)):
			edp_file = edp_files[j].split('/')[-1]
					
			#Read the IDL data for this EDP profile
			idl_data = idlsave.read(edp_files[j]+'/'+occ_file+'_'+edp_file+'_edp.sav',\
				verbose=False)
			
			#Create name of output PDS data file. This will require parsing the profile
			# file for the RSR product names that went into this particular profie, and 
			# also manipulating the edp_file variable. 
			with open(edp_files[j]+'/'+edp_file+'.profile','r') as pro:
				while 1:
					#Skip any header lines
					line = pro.readline()
					if '#' in line: continue
					split_line = line.split('=')
					if 'PRODUCT_A' in line: 
						product_a = split_line[-1]
						line = pro.readline()
						split_line = line.split('=')
						product_b = split_line[-1]
					if not line: break
			pro.close()
			#The product IDs are all the same length, make sure the prefixes are identical.
			if product_a[:19] != product_b[:19]: 
				print ''
				print 'Product A and B do not seem to match. Break here.'
				stop = input('stop')
			
			#Otherwise create and open the save file name using either product a or b
			save_filename = product_a[:19]+'_'+edp_file.split(occ_file)[-1]+'_'+body+\
				'_edp_v01.csv'
			save_filename = save_filename.upper()
			save_file = open(save_path_lvl1+save_filename,'w')
			
			#Generate a string that comprises this line in the save file
			for k in range(size(idl_data.edpsav)):
				s = ''
				s += '%20.12e' % idl_data.ettxarraysav[k] +','
				s += '%20.12e' % idl_data.etoccptarraysav[k] +','
				s += '%20.12e' % idl_data.etrxarraysav[k] +','
				s += '%30.30s' % idl_data.utctxarraysav[k] +','
				s += '%30.30s' % idl_data.utcoccptarraysav[k] +','
				s += '%30.30s' % idl_data.utcrxarraysav[k] +','
				
				s += '%20.12e' % idl_data.occptradiusarraysav[k] +','
				s += '%20.12e' % idl_data.occptlatarraysav[k] +','
				s += '%20.12e' % idl_data.occptlonarraysav[k] +','
				s += '%20.12e' % idl_data.occptszaarraysav[k] +','
				s += '%20.12e' % idl_data.occptlstarraysav[k] +','
				s += '%20.12e' % idl_data.sepanglearraysav[k] +','
				s += '%20.12e' % idl_data.epsanglearraysav[k] +','
				s += '%20.12e' % idl_data.dxdt_save[k] +','
				s += '%20.12e' % idl_data.dxdt[k] +','
				s += '%20.12e' % idl_data.newbigx[k] +','
				s += '%20.12e' % idl_data.edpsav[k] +','
				s += '%20.12e' % idl_data.electrondensity_stddevsav[k]
				if k<(size(idl_data.edpsav)-1): s += '\n'		
				save_file.write(s)
			save_file.close()
#-----------------------------------------------------------------------------------------



#Level 2 save procedure: Given a body (e.g., 'Titan'), the code will comb through all
# Titan EDP directories and create the average output EDP in the PDS save directory.
#-----------------------------------------------------------------------------------------
if False:
	body = 'Titan'				
	occ_files = glob.glob(path_lvl1+body+'/*')
	#Loop over all occ directories
	for i in range(size(occ_files)):
		occ_file = occ_files[i].split('/')[-1]
		#Get a list of all averaged EDP directories within this occ directory
		edp_files = glob.glob(occ_files[i]+'/*.pickle')
		#Loop over each EDP file in this occ
		for j in range(size(edp_files)):
			edp_file = edp_files[j].split('/')[-1].split('.')[0].split('_')[-1]
			#Revise this edp_file for saving, where the leading zero must be included for
			# flyby numbers that are below 100
			edp_file_save = edp_file[0]+edp_file[1:-1].zfill(3)+edp_file[-1]
			#Load the data with pickle
			with open(edp_files[j],'rb') as pckl:
				ettx, etoccpt, etrx, utctx, utcoccpt,utcrx, lat, lon, sza, lst, sep, \
					eps, occptradius, ed_avg, err_avg = pickle.load(pckl)
			pckl.close()
			
			#Get the directories that correspond to this edp file
			edp_directories = glob.glob(occ_files[i]+'/'+edp_file+'*')
			#Go into one of these directories (which does not matter), open the profile
			# file, and extract the sequence number of either of products. 
			profiles = glob.glob(edp_directories[0]+'/*.profile')
			with open(profiles[0],'r') as pro:
				while 1:
					#Skip any header lines
					line = pro.readline()
					if '#' in line: continue
					split_line = line.split('=')
					if 'PRODUCT_A' in line: 
						product_a = split_line[-1]
					if not line: break
			pro.close()
			
			#Create the save file
			save_filename = product_a[0:3]+'_'+edp_file_save+'_'+body+'_edp_v01.CSV'
			save_filename = save_filename.upper()
			save_file = open(save_path_lvl2+save_filename,'w')	
			
			#Generate a string that comprises this line in the save file
			for k in range(size(ed_avg)):
				s = ''
				s += '%20.12e' % ettx[k] +','
				s += '%20.12e' % etoccpt[k] +','
				s += '%20.12e' % etrx[k] +','
				s += '%30.30s' % utctx[k] +','
				s += '%30.30s' % utcoccpt[k] +','
				s += '%30.30s' % utcrx[k] +','
				
				s += '%20.12e' % occptradius[k] +','
				s += '%20.12e' % lat[k] +','
				s += '%20.12e' % lon[k] +','
				s += '%20.12e' % sza[k] +','
				s += '%20.12e' % lst[k] +','
				s += '%20.12e' % sep[k] +','
				s += '%20.12e' % eps[k] +','
				s += '%20.12e' % ed_avg[k] +','
				s += '%20.12e' % err_avg[k]
				if k<(size(ed_avg)-1): s += '\n'		
				save_file.write(s)				
			save_file.close()
#-----------------------------------------------------------------------------------------
			
			
			
			
#Summary table procedure: Do nearly the same things as for level 2, but only save the one
# line that is closest to an occ pt radius of 2575 km + 1200 km (Kliore et al. 2008) to
# a single summary table
#-----------------------------------------------------------------------------------------
if False:
	body = 'Titan'
	#Create save path
	save_filename = body+'_summary_table.csv'
	save_filename = save_filename.upper()
	save_file = open(save_path_summ+save_filename,'w')
	
	#Set the occ pt radius value to minimize distance to
	occpt_crit = 2575. + 1200.
				
	occ_files = glob.glob(path_lvl1+body+'/*')
	#Loop over all occ directories
	for i in range(size(occ_files)):
		occ_file = occ_files[i].split('/')[-1]
		#Get a list of all averaged EDP directories within this occ directory
		edp_files = glob.glob(occ_files[i]+'/*.pickle')
		#Loop over each EDP file in this occ
		for j in range(size(edp_files)):
			edp_file = edp_files[j].split('/')[-1].split('.')[0].split('_')[-1]
			#Revise this edp_file for saving, where the leading zero must be included for
			# flyby numbers that are below 100
			edp_file_save = edp_file[0]+edp_file[1:-1].zfill(3)+edp_file[-1]
			#Load the data with pickle
			with open(edp_files[j],'rb') as pckl:
				ettx, etoccpt, etrx, utctx, utcoccpt,utcrx, lat, lon, sza, lst, sep, \
					eps, occptradius, ed_avg, err_avg = pickle.load(pckl)
			pckl.close()
			
			#Get the directories that correspond to this edp file
			edp_directories = glob.glob(occ_files[i]+'/'+edp_file+'*')
			#Go into one of these directories (which does not matter), open the profile
			# file, and extract the sequence number of either of products. 
			profiles = glob.glob(edp_directories[0]+'/*.profile')
			with open(profiles[0],'r') as pro:
				while 1:
					#Skip any header lines
					line = pro.readline()
					if '#' in line: continue
					split_line = line.split('=')
					if 'PRODUCT_A' in line: 
						product_a = split_line[-1]
					if not line: break
			pro.close()
			
			#Generate a string that comprises this line in the save file
			for k in range(size(ed_avg)):
				s = ''
				s += '%10.10s' % edp_file_save +','
				s += '%30.30s' % utcoccpt[k] +','
				s += '%20.12e' % lat[k] +','
				s += '%20.12e' % lon[k] +','
				s += '%20.12e' % sza[k] +','
				s += '%20.12e' % lst[k] +','
				s += '%20.12e' % sep[k] +','
				s += '%20.12e' % eps[k] +','
				s += '%20.12e' % err_avg[k]
				s += '\n'		
				
				#Determine whether is line is the critical one to save by minimizing the 
				# distance to the critical occ pt value
				if k==0: 
					prev_occpt_dist = abs(occptradius[k]-occpt_crit)
					continue
				#If the absolute value of the distance to the critical occ pt value turns
				# around, then save the previous line
				if abs(occptradius[k]-occpt_crit) > prev_occpt_dist:
					#print edp_file_save,occptradius[k-1], occptradius[k]
					save_file.write(s_save)
					break
				#Save this string, as it may be needed after the next evaluation of s
				s_save = copy(s) 
				prev_occpt_dist = abs(occptradius[k]-occpt_crit)		
	save_file.close()
#-----------------------------------------------------------------------------------------			
			
			


#Level 0 save procedure for any body/experiment OTHER than the Titan occs that have al-
# ready been processed. It looks into the cors files and finds only those without the 
# flag_indices.txt files. For each appropriate RSR file, it will create a save file and a
# 'cors' directory in the appropriate save directory with the csv data.
#-----------------------------------------------------------------------------------------
if False:
	body = 'Saturn'   #options include Titan, Enceladus, Saturn, Rings
	cors_files = glob.glob(path_lvl0+body+'/cors_*')
	#Loop over all cors directories
	for i in range(size(cors_files)):
		cors_file = cors_files[i].split('/')[-1]
		#Get a list of all RSR directories within this cors file
		rsr_files = glob.glob(cors_files[i]+'/*')
		#Loop over all RSR directories
		for j in range(size(rsr_files)):
			rsr_file = rsr_files[j].split('/')[-1]
			#Check to see if this RSR directory successfully finished the frequency time-
			# series analysis. Some did not because their RSR files were mysteriously
			# truncated. Skip these
			if size(glob.glob(rsr_files[j]+'/full_output/')) == 0: 
				continue
			#Check to see if this RSR directories was used to make an EDP. If so, there 
			# will be a flag_indices.txt file in its full_output directory. In this case,
			# we only want the RSRs that do not have the flag_indices file.
			if size(glob.glob(rsr_files[j]+'/full_output/flag_indices.txt')) != 0: 
				continue
			#Otherwise, create a cors directory for this RSR file in the save path (check
			# to see if one is there already). 
			save_path = save_path_lvl0_pt2+cors_file+'/'
			if size(glob.glob(save_path)) == 0: os.system('mkdir '+save_path)
			
			#Now read in the IDL data
			idl_data = idlsave.read(rsr_files[j]+'/full_output/output.sav',verbose=False)
			
			#Open the file in which to save the output data
			save_filename = '_'.join(rsr_file.split('.'))+'_freq_v01.csv'
			save_filename = save_filename.upper()
			save_file = open(save_path+save_filename,'w')
			
			#Generate a string that comprises this line in the save file
			for k in range(size(idl_data.sfduyearout)):
				s = ''
				s += '%10.i' % idl_data.sfduyearout[k] +','
				s += '%10.i' % idl_data.sfdudoyout[k] +','
				s += '%20.3f' % idl_data.sfdusecout[k] +','
				s += '%10.i' % idl_data.rftoifmhzout[k] +','
				s += '%10.i' % idl_data.ddclomhzout[k] +','
				s += '%20.12e' % idl_data.ncofreqout[k] +','
				s += '%20.12e' % idl_data.kposaaout[k] +','
				s += '%20.12e' % abs(idl_data.igcomplexaaout[k]) +','
				#Quality flag columns. Since these have not been used to make an EDP, list
				# all quality flag entries as '9'. 
				s += '%5.1s' % '9' + ','
				s += '%5.1s' % '9'		
				if k<(size(idl_data.sfduyearout)-1): s += '\n'		
				save_file.write(s)
			save_file.close()
#-----------------------------------------------------------------------------------------
			
			
			
