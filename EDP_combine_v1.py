#The purpose of this script is to combine Titan EDP from various receivers for a given
# occultation observation. If there is a corresponding Kliore EDP on the PDS, it will
# be plotted along side the combined EDP. 

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

#Specify the Titan flyby and the observation ID
body = 'Titan'
flyby = 'T57'
obs = 'T57N'
plot_kliore = 'y'
save_data = 'n'

#Set paths
root_path = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/'
edp_path = root_path+'EDP/'

#Get all the applicable directories
obs_list = glob.glob(edp_path+body+'/'+flyby+'/'+obs+'*')	

#Create size arrays that will be needed to trim and combine EDP
tops, bottoms = zeros(size(obs_list)), zeros(size(obs_list))

#Build the individual EDP figure	
fig1 = figure(1)
ax1 = fig1.add_subplot(111)
ax1.axvline(0,c='k',ls='--')
ax1.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
ax1.set_ylabel('Altitude [km]',fontsize='x-large')
xticks(fontsize='large')
yticks(fontsize='large')
ax1.minorticks_on()
ax1.set_title('Individual EDP',fontsize='x-large')

#For each, load in all of the EDP, parse the profile files, and plot
all_IDs = []
topsides = zeros(size(obs_list))
avg_topside, std_topside, rms_topside = zeros(size(obs_list)), zeros(size(obs_list)), \
	zeros(size(obs_list))
print ''
print 'Mean          STD           RMS           Max(alt)'
for i in range(size(obs_list)):
	all_IDs += [obs_list[i].split('/')[-1]]
	alt, ed = loadtxt(obs_list[i]+'/'+flyby+'_'+all_IDs[i]+'_edp.txt',usecols=[0,1],\
		unpack=True)
	#thisrkm to Titan altitude
	alt -= 2575.  
	
	#Save start and end of alt
	tops[i], bottoms[i] = max(alt), min(alt)
	
	#Plot individual EDP
	ax1.plot(ed,alt,label=str(i)+', '+all_IDs[i])
	
	#Parse the profile file to get necessary info
	with open(obs_list[i]+'/'+all_IDs[i]+'.profile','r') as input_file:
		while 1:
			#Skip any header lines
			line = input_file.readline()
			if '#' in line: continue
			split_line = line.split('=')
			if 'TOPSIDE' in line: topsides[i] = int(split_line[-1])
			if not line: break
	input_file.close()
	
	#Print out topside information for these profiles
	topside_inds = where(alt>=topsides[i])[0]
	avg_topside[i] = mean(ed[topside_inds])
	std_topside[i] = std(ed[topside_inds])
	rms_topside[i] = sqrt(sum(ed[topside_inds]**2)/size(topside_inds))
	print avg_topside[i],std_topside[i],rms_topside[i],tops[i] 	
ax1.legend(loc='best',prop={'size':12})	


#Give the user the option to use all curves or exclude some (based on index).
#####TURN THIS OFF FOR NOW - ASSUME ALL INDIVIDUAL PROFILES WILL BE USED
print '\n'
print 'Press enter to average all EDP together'
#good_inds = raw_input('Otherwise, enter an index list of the good profiles: ')
good_inds = ''
if good_inds=='': good_inds = range(size(obs_list))
else: good_inds = ast.literal_eval(good_inds)


#I can only combine where there is overlap between all of the good EDP. 
top_cut = min(tops[good_inds])
bottom_cut = max(bottoms[good_inds])

#See if the trim produces equally sized alt arrays
alt_sizes = zeros(size(good_inds),dtype=int)
for i,j in enumerate(good_inds):
	alt, ed = loadtxt(obs_list[j]+'/'+flyby+'_'+all_IDs[j]+'_edp.txt',usecols=[0,1],\
		unpack=True)
  	#thisrkm to Titan altitude
  	alt -= 2575.
	alt_sizes[i] = size(alt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])	
if (alt_sizes != alt_sizes[0]).any():
	print ''
	stop = input('Altitude sizes are not matching up (point 1). Break here. ')


#Now, create the averaged EDP profile. Do a simple weighted average and also do the Kliore
# et al. 2011 weighted averaging method

#Simple weighted average
fig2 = figure(2)
ax2 = fig2.add_subplot(111)
ax2.axvline(0,c='k',ls='--')
ax2.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
ax2.set_ylabel('Altitude [km]',fontsize='x-large')
xticks(fontsize='large')
yticks(fontsize='large')
ax2.minorticks_on()
ax2.set_title('Simple weighted average',fontsize='x-large')

all_alt, all_ed, all_err = array([]), array([]), array([])
for i,j in enumerate(good_inds):
	alt, ed = loadtxt(obs_list[j]+'/'+flyby+'_'+all_IDs[j]+'_edp.txt',usecols=[0,1],\
		unpack=True)
  	#thisrkm to Titan altitude
  	alt -= 2575.
	all_alt = append(all_alt,alt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_ed = append(all_ed,ed[where((alt>=bottom_cut)&(alt<=top_cut))[0]])	
	all_err = append(all_err,ones(alt_sizes[i])*std_topside[j])

#Reshape arrays
all_alt = all_alt.reshape(size(good_inds),alt_sizes[0])	
all_ed = all_ed.reshape(size(good_inds),alt_sizes[0])	
all_err = all_err.reshape(size(good_inds),alt_sizes[0])	

#Run a quick check on the altitudes
if round(sum(all_alt - mean(all_alt,axis=0)),5) != 0:
	print ''
	stop = input('Altitude sizes are not matching up (point 2). Break here. ')

#Run the weighted average
alt = mean(all_alt,axis=0)
ed_avg = sum(all_ed/all_err**2,axis=0)/sum(1./all_err**2,axis=0)
err_avg = 1./sqrt(sum(1./all_err**2,axis=0))
ax2.plot(ed_avg,alt,c='k',lw=2,label='This work')
ax2.errorbar(ed_avg,alt,xerr=err_avg,fmt='none',ecolor='k',capsize=0,lw=0.5)
########################################################################################
# 
# #Kliore averaging method
# fig3 = figure(3)
# ax3 = fig3.add_subplot(111)
# ax3.axvline(0,c='k',ls='--')
# ax3.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
# ax3.set_ylabel('Altitude [km]',fontsize='x-large')
# xticks(fontsize='large')
# yticks(fontsize='large')
# ax3.minorticks_on()
# ax3.set_title('Kliore+2011 Averaging Method',fontsize='x-large')
# 
# all_alt, all_ed, all_err = array([]), array([]), array([])
# for i,j in enumerate(good_inds):
# 	alt, ed = loadtxt(obs_list[j]+'/'+flyby+'_'+all_IDs[j]+'_edp.txt',usecols=[0,1],\
# 		unpack=True)
#   	#thisrkm to Titan altitude
#   	alt -= 2575.
# 	all_alt = append(all_alt,alt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
# 	all_ed = append(all_ed,ed[where((alt>=bottom_cut)&(alt<=top_cut))[0]])	
# 	all_err = append(all_err,ones(alt_sizes[i])*std_topside[j])
# 
# #Reshape arrays
# all_alt = all_alt.reshape(size(good_inds),alt_sizes[0])	
# all_ed = all_ed.reshape(size(good_inds),alt_sizes[0])	
# all_err = all_err.reshape(size(good_inds),alt_sizes[0])	
# 
# #Run the Kliore combination methods
# alt = mean(all_alt,axis=0)
# ed_avg = mean(all_ed,axis=0)
# 
# baseline_std = sqrt(sum(std_topside[good_inds]**2))/size(good_inds)
# err_avg = sqrt(std(all_ed,axis=0)/sqrt(size(good_inds))**2 + baseline_std**2)
# 
# ax3.plot(ed_avg,alt,c='k',lw=2,label='This work')
# ax3.errorbar(ed_avg,alt,xerr=err_avg,fmt='none',ecolor='k',capsize=0,lw=0.5)
# ########################################################################################
# 

#Also, include the Kliore occs on applicable occultations

kliore_path = root_path+'Misc/Kliore2008_TitanEDP/'
csvfile = 'skip'
if obs=='T12N': csvfile = kliore_path+'S19TIIOC2006_078_0000_N_001.EDP'
if obs=='T12X': csvfile = kliore_path+'S19TIIOC2006_078_0000_X_001.EDP'
if obs=='T14N': csvfile = kliore_path+'S20TIIOC2006_140_1203_N_001.EDP'
if obs=='T14X': csvfile = kliore_path+'S20TIIOC2006_140_1203_X_001.EDP'
if obs=='T27N': csvfile = kliore_path+'S28TIIOC2007_085_0000_N_001.EDP'
if obs=='T27X': csvfile = kliore_path+'S28TIIOC2007_085_0000_X_001.EDP'
if obs=='T31N': csvfile = kliore_path+'S30TIIOC2007_148_1737_N_001.EDP'
if obs=='T31X': csvfile = kliore_path+'S30TIIOC2007_148_1737_X_001.EDP'
if (csvfile != 'skip') and (plot_kliore == 'y'):
	alt_kliore, err_kliore, ed_kliore = loadtxt(csvfile,usecols=[0,1,2],delimiter=',',\
		unpack=True)
	#Simple avg figure
	ax2.plot(ed_kliore,alt_kliore,c='r',label='Kliore+2011')
	ax2.errorbar(ed_kliore,alt_kliore,xerr=err_kliore,fmt='none',ecolor='r',capsize=0,\
		lw=0.5)	
	ax2.legend(loc='best',prop={'size':12})
	#Kliore+2011 method figure
	#ax3.plot(ed_kliore,alt_kliore,c='r',label='Kliore+2011')
	#ax3.errorbar(ed_kliore,alt_kliore,xerr=err_kliore,fmt='none',ecolor='r',capsize=0,\
	#	lw=0.5)	
	#ax3.legend(loc='best',prop={'size':12})


#Save the final EDP data and the plots in the flyby directory with the obs name
if save_data == 'y':
	save_path = edp_path+body+'/'+flyby+'/'
	fig1.savefig(save_path+'ind_edp_'+obs+'.eps',bbox_inches='tight')
	fig2.savefig(save_path+'final_edp_'+obs+'.eps',bbox_inches='tight')

	file = open(save_path+'final_edp_'+obs+'.txt','w')
	for i in range(size(ed_avg)):
		file.write('%.4f' % alt[i]+'\t\t'+'%.8f' % ed_avg[i]+'\t\t'+'%.8f' % err_avg[i]+'\n')
	file.close()





#stop=input('stop')	 	