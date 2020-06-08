"""
Name
----
batch_compiler.py

Description
-----------
This script takes several inputs from the user (cors #, band, polarization) and prepares 
to run all the matching RSR files as a batch on the SCC. It first locates all the desired
files in the Data storage directories, then it creates all the output directory structure,
the local bash file and a global bash file. It gives the user a chance to check the final
bash file before running.

Input
-----
The data directory structure must be specified, along with cors number and the band and
polarization info. Also, any qsub parameters that are not default (e.g., run time).

Output
------
For each file, a bash file that will allow it to be submitted to the SCC. Also output 
directory structure, and a bash file containing qsub commands for all selected files. 

Author
------
P. A. Dalba -- May 2017, Boston University
"""

#Import various math, science, and plotting packages.
imprt = 0
while imprt==0:
	import numpy, scipy
	from scipy.integrate import *
	from scipy.interpolate import interp1d
	from scipy.interpolate import InterpolatedUnivariateSpline
	import commands
	import os, glob
	import sys, math
	from scipy import stats, signal
	from numpy import *
	from scipy import optimize
	from scipy.optimize import curve_fit
	from scipy.optimize import fsolve
	import time
	import pickle
	from scipy.ndimage import interpolation
	from scipy.spatial import distance
	imprt = 1 
	
	
#Users input selection. user inputs, including any qsub parameters
body = 'Titan'   #Titan, Saturn, Rings, Enceladus
cors = '276'     #Must be a string with no zero in front
band = ['s','x','k']    #List of 's', 'x', or 'k'. Lower case
pol = ['r','l']      #List of 'l' or 'r'. Lower case


#Access the cors file in the data storage directory and get a list of all the data files
# there. Use this list to determine the capitalization scheme of the cors set. I assume
# that within a cors directory the capitalization does not change. 
all_cors_path = '/projectnb/cassini/Cassini_Ionosphere/Data/'
cors_path = all_cors_path+'cors_0'+cors+'/'
output_body_path = '/project/cassini/Cassini_Ionosphere/Output/'+body+'/'

#Include the path to the batch RSR proc code that will be run on all the RSR files.
program_path = '/project/cassini/Cassini_Ionosphere/Programs/'

#Determine the capitalization status using the label files suffix. Each cors file should
# have at least one label file so this should always work. If not, code will stop. This
# and the data input above initially assumes lowercase
caps = 'no'
all_LBL_files = sorted(glob.glob(cors_path+'*.LBL'))
all_lbl_files = sorted(glob.glob(cors_path+'*.lbl'))
if (size(all_LBL_files)==0) & (size(all_lbl_files)==0):
	stop = input('Break here - no label files in cors directory')
if size(all_LBL_files) > size(all_lbl_files):
	caps = 'yes'
	#Change the case if caps is 'yes'
	for i in range(size(band)):
		band[i] = band[i].upper() 
	for i in range(size(pol)):
		pol[i] = pol[i].upper() 


#Find all the files in the cors directory and isolate just the RSR files (no labels)
all_rsr_files = []
all_files = sorted(glob.glob(cors_path+'*'))
for i in range(size(all_files)):
	if ('.lbl' in all_files[i]) | ('.LBL' in all_files[i]): continue
	all_rsr_files.append(all_files[i])

#If all three bands were not listed above, then remove RSR files that should not be run.
# If the list is reduced to zero, break the code.
if size(band) != 3:
	del_list = []
	for i in range(size(all_rsr_files)):
		if all_rsr_files[i].split('/')[-1][22] not in band: del_list.append(i)
	all_rsr_files = delete(all_rsr_files,del_list)
if size(all_rsr_files)==0:
	stop = input('Cut #1 reduced the number of files to run to zero. Break here.')

#If both polarizations were not listed above, then remove RSR files that should not be run.
# If the list is reduced to zero, break the code.
if size(pol) != 2:
	del_list = []
	for i in range(size(all_rsr_files)):
		if all_rsr_files[i].split('/')[-1][25] not in pol: del_list.append(i)
	all_rsr_files = delete(all_rsr_files,del_list)
if size(all_rsr_files)==0:
	stop = input('Cut #2 reduced the number of files to run to zero. Break here.')

#Make an output directory for this cors (won't make one if one already exists)
output_cors_path = output_body_path+'cors_0'+cors+'/'
os.system('mkdir '+output_cors_path)

#Make an output directory for each RSR file (again, won't make one if it exists already)
output_rsr_paths = []
for i in range(size(all_rsr_files)):
	output_rsr_paths.append(output_cors_path+all_rsr_files[i].split('/')[-1]+'/')
	os.system('mkdir '+output_rsr_paths[i])
	
#For each RSR directory, see if a full_ouput directory exists, if it does, remove it from
# the RSR file list. It will not be run. If the list is reduced to zero, break the code.
del_list = []
for i in range(size(output_rsr_paths)):	
	full_run = size(glob.glob(output_rsr_paths[i]+'full_*'))
	if full_run == 1: del_list.append(i)
all_rsr_files = delete(all_rsr_files,del_list)
output_rsr_paths = delete(output_rsr_paths,del_list)
all_rsr_files = list(all_rsr_files)
output_rsr_paths = list(output_rsr_paths)
if size(all_rsr_files)==0:
	stop = input('Cut #3 reduced the number of files to run to zero. Break here.')


#Copy the batch RSR proc file to each of the RSR working directories
for i in range(size(all_rsr_files)):
	os.system('cp '+program_path+'rsr_proc_batch.pro '+output_rsr_paths[i])


#Create a bash script in the working directory with the correct parameters and syntax for
# each RSR file run.
for i in range(size(all_rsr_files)):
	script_file = open(output_rsr_paths[i]+'script.sh','w')
	#Start with the basic header
	script_file.write('#!/bin/bash'+os.linesep+os.linesep)
	#CD to this working directory
	script_file.write('cd '+output_rsr_paths[i]+os.linesep)
	#Compile the batch rsr code
	script_file.write('idl -IDL_CPU_TPOOL_NTHREADS 1 -e ".run rsr_proc_batch.pro"'+os.linesep)
	#Now actually run the procedure with the proper arguments and in the correct oder
	script_file.write('idl -IDL_CPU_TPOOL_NTHREADS 1 -e "rsr_proc_batch" -args '+body+\
		' '+cors+' '+all_rsr_files[i].split('/')[-1][22].upper()+' '+\
			all_rsr_files[i].split('/')[-1]+os.linesep)
	script_file.close()


#Generate a single, separate bash script in the local Programs directory that includes the
# qsub commands for all the RSR files in this batch. 
qsub_file = open(program_path+'runner.sh', 'w')
for i in range(size(all_rsr_files)):
	qsub_file.write('qsub -P cassini -m e -l h_rt=240:00:00 -e '+output_rsr_paths[i]+\
		' -o '+output_rsr_paths[i]+' '+output_rsr_paths[i]+'script.sh'+os.linesep)
qsub_file.close()


#Stop here. If the user so desires, actually submit the qsub to the cluster after the stop.
print ''
stop = raw_input('Press enter to submit the '+str(size(all_rsr_files))+' job(s) to the SCC.')  
os.system('sh runner.sh')



