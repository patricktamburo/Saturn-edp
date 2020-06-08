#The purpose of this script is to combine Titan EDP from various receivers for a given
# occultation observation. If there is a corresponding Kliore EDP on the PDS, it will
# be plotted along side the combined EDP.

#v2 creates the same profiles, but does so using the new IDL save files that contain all
# of the information to be saved in PDS files.

#Use keyword 'spherical' or 'ellipsoidal' to choose between different representations of Saturn.

#Import various math, science, and plotting packages.
imprt = 0
while imprt==0:
	import numpy, scipy
	from scipy.interpolate import interp1d
	#import commands
	import matplotlib
	#matplotlib.use('macosx')
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

#Specify the body, flyby and the observation ID

body = 'Saturn'
flyby = 'S68'
obs = 'S68N'
plot_kliore = 'n'
save_data = 'y'

representation = 'ellipsoidal' #Either spherical, ellipsoidal, or Schubert

#Declare one-bar semi-major and -minor axes for representation of ellipsoidal Saturn
a = 60268.
b = 54364.

#Set paths
root_path = 'C:\\Users\\tambu\\Documents\\BU\\SaturnEDP\\PD_MBP_files\\'
edp_path = root_path+'EDP\\'


#Get all the applicable directories
obs_list = glob.glob(edp_path+body+'\\'+flyby+'\\'+obs+'*')


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
print('')
print('Mean          STD           RMS           Max(alt)')


for i in range(size(obs_list)):
	all_IDs += [obs_list[i].split('\\')[-1]]

	idl_data = idlsave.read(obs_list[i]+'\\'+flyby+'_'+all_IDs[i]+'_edp_pat.sav',verbose=0)

	alt = idl_data.occptradiusarraysav
	lat = idl_data.occptlatarraysav
	ed = idl_data.edpsav

	if representation == 'ellipsoidal':
		#Calculate saturn_R, the ellipsoidal represntation of Saturn's surface, as a function of latitude
		saturn_R = a*b/(np.sqrt((b*np.cos(lat*np.pi/180.))**2.+(a*np.sin(lat*np.pi/180.))**2.))
		#Calculate occptradius - R
		r_minus_R = alt - saturn_R
		#Calculate the "representative" poinnt, which is closest to height of 2000 km (representative of
		#	ionospheric peak)
		rep_point_loc = np.where(abs(r_minus_R-2000) == np.min(abs(r_minus_R-2000)))[0]
		#Subtract off the height of representative point
		alt -= saturn_R[rep_point_loc]
		print('Representative Latitude: ',lat[rep_point_loc])
	elif representation == 'spherical':
		#thisrkm to Saturn altitude
		alt -= 60268.
	else:
		print('ERROR: representation not supported.')

	#Save start and end of alt. The round is due to T57N, which has a 0.1cm discrepancy in
	# bottom value that broke the code
	tops[i], bottoms[i] = max(alt), round(min(alt),4)

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
	print(avg_topside[i],std_topside[i],rms_topside[i],tops[i])
ax1.legend(loc='best',prop={'size':12})


#Give the user the option to use all curves or exclude some (based on index).
#####TURN THIS OFF FOR NOW - ASSUME ALL INDIVIDUAL PROFILES WILL BE USED
print('\nPress enter to average all EDP together')
print('(Averaging all EDP together)')
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
	#alt, ed = loadtxt(obs_list[j]+'/'+flyby+'_'+all_IDs[j]+'_edp.txt',usecols=[0,1],\
	#	unpack=True)
	#idl_data = idlsave.read(obs_list[i]+'/'+flyby+'_'+all_IDs[i]+'_edp.sav',verbose=0)
	idl_data = idlsave.read(obs_list[i]+'\\'+flyby+'_'+all_IDs[i]+'_edp_pat.sav',verbose=0)

	alt = idl_data.occptradiusarraysav
	lat = idl_data.occptlatarraysav
	ed = idl_data.edpsav
	#thisrkm to Titan altitude
	#alt -= 2575.
	if representation == 'ellipsoidal':
		#Calculate saturn_R, the ellipsoidal represntation of Saturn's surface, as a function of latitude
		saturn_R = a*b/(np.sqrt((b*np.cos(lat*np.pi/180.))**2.+(a*np.sin(lat*np.pi/180.))**2.))
		#Calculate occptradius - R
		r_minus_R = alt - saturn_R
		#Calculate the "representative" poinnt, which is closest to height of 2000 km (representative of
		#	ionospheric peak)
		rep_point_loc = np.where(abs(r_minus_R-2000) == np.min(abs(r_minus_R-2000)))[0]
		#Subtract off the height of representative point
		alt -= saturn_R[rep_point_loc]
	elif representation == 'spherical':
		alt -= 60268.
	else:
		print('ERROR: representation not supported.')

	alt_sizes[i] = size(alt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	if (alt_sizes != alt_sizes[0]).any():
		print('')
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

#Create the arrays of everything we want to save
all_alt, all_ed, all_err = array([]), array([]), array([])
all_ettx, all_etoccpt, all_etrx = array([]), array([]), array([])
all_utctx, all_utcoccpt, all_utcrx = array([]), array([]), array([])
all_lat, all_lon, all_sza, all_lst, all_sep, all_eps = array([]), array([]), array([]),\
	array([]), array([]), array([])

for i,j in enumerate(good_inds):
	#alt, ed = loadtxt(obs_list[j]+'/'+flyby+'_'+all_IDs[j]+'_edp.txt',usecols=[0,1],\
	#	unpack=True)
	#idl_data = idlsave.read(obs_list[i]+'/'+flyby+'_'+all_IDs[i]+'_edp.sav',verbose=0)
	idl_data = idlsave.read(obs_list[i]+'\\'+flyby+'_'+all_IDs[i]+'_edp_pat.sav',verbose=0)

	alt = idl_data.occptradiusarraysav
	ed = idl_data.edpsav
	err = idl_data.electrondensity_stddevsav
	ettx = idl_data.ettxarraysav
	etoccpt = idl_data.etoccptarraysav
	etrx = idl_data.etrxarraysav
	utctx = idl_data.utctxarraysav
	utcoccpt = idl_data.utcoccptarraysav
	utcrx = idl_data.utcrxarraysav
	lat = idl_data.occptlatarraysav
	lon = idl_data.occptlonarraysav
	sza = idl_data.occptszaarraysav
	lst = idl_data.occptlstarraysav
	sep = idl_data.sepanglearraysav
	eps = idl_data.epsanglearraysav

	if representation == 'ellipsoidal':
		saturn_R = a * b / (np.sqrt((b * np.cos(lat * np.pi / 180.)) ** 2. + (a * np.sin(lat * np.pi / 180.)) ** 2.))
		# Calculate occptradius - R
		r_minus_R = alt - saturn_R
		# Calculate the "representative" poinnt, which is closest to height of 2000 km (representative of
		#	ionospheric peak)
		rep_point_loc = np.where(abs(r_minus_R - 2000) == np.min(abs(r_minus_R - 2000)))[0]
		# Subtract off the height of representative point
		alt -= saturn_R[rep_point_loc]
	elif representation == 'spherical':
		alt -= 60268.
	else:
		print('ERROR: representation not supported.')

	all_alt = append(all_alt,alt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_ed = append(all_ed,ed[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_err = append(all_err,err[where((alt>=bottom_cut)&(alt<=top_cut))[0]])

	#All the other quantities
	all_ettx = append(all_ettx,ettx[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_etoccpt = append(all_etoccpt,etoccpt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_etrx = append(all_etrx,etrx[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_utctx = append(all_utctx,utctx[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_utcoccpt = append(all_utcoccpt,utcoccpt[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_utcrx = append(all_utcrx,utcrx[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_lat = append(all_lat,lat[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_lon = append(all_lon,lon[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_sza = append(all_sza,sza[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_lst = append(all_lst,lst[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_sep = append(all_sep,sep[where((alt>=bottom_cut)&(alt<=top_cut))[0]])
	all_eps = append(all_eps,eps[where((alt>=bottom_cut)&(alt<=top_cut))[0]])

#Reshape arrays
all_alt = all_alt.reshape(size(good_inds),alt_sizes[0])
all_ed = all_ed.reshape(size(good_inds),alt_sizes[0])
all_err = all_err.reshape(size(good_inds),alt_sizes[0])
all_ettx = all_ettx.reshape(size(good_inds),alt_sizes[0])
all_etoccpt = all_etoccpt.reshape(size(good_inds),alt_sizes[0])
all_etrx = all_etrx.reshape(size(good_inds),alt_sizes[0])
all_utctx = all_utctx.reshape(size(good_inds),alt_sizes[0])
all_utcoccpt = all_utcoccpt.reshape(size(good_inds),alt_sizes[0])
all_utcrx = all_utcrx.reshape(size(good_inds),alt_sizes[0])
all_lat = all_lat.reshape(size(good_inds),alt_sizes[0])
all_lon = all_lon.reshape(size(good_inds),alt_sizes[0])
all_sza = all_sza.reshape(size(good_inds),alt_sizes[0])
all_lst = all_lst.reshape(size(good_inds),alt_sizes[0])
all_sep = all_sep.reshape(size(good_inds),alt_sizes[0])
all_eps = all_eps.reshape(size(good_inds),alt_sizes[0])

#Run a check on the shared quantities
if round(sum(all_alt - mean(all_alt,axis=0)),5) != 0:
	print('')
	stop = input('Altitude sizes are not matching up (point 2). Break here. ')
if round(sum(all_ettx - mean(all_ettx,axis=0)),5) != 0:
	print('')
	stop = input('ET tx times are not matching up. Break here. ')
if round(sum(all_etoccpt - mean(all_etoccpt,axis=0)),5) != 0:
	print('')
	stop = input('ET occ pt times are not matching up. Break here. ')
if round(sum(all_etrx - mean(all_etrx,axis=0)),5) != 0:
	print('')
	stop = input('ET rx times are not matching up. Break here. ')
if round(sum(all_lat - mean(all_lat,axis=0)),5) != 0:
	print('')
	stop = input('Latitudes are not matching up. Break here. ')
if round(sum(all_lon - mean(all_lon,axis=0)),5) != 0:
	print('')
	stop = input('Longitudes are not matching up. Break here. ')
if round(sum(all_sza - mean(all_sza,axis=0)),5) != 0:
	print('')
	stop = input('SZAs are not matching up. Break here. ')
if round(sum(all_lst - mean(all_lst,axis=0)),5) != 0:
	print('')
	stop = input('LSTs are not matching up. Break here. ')
if round(sum(all_sep - mean(all_sep,axis=0)),5) != 0:
	print('')
	stop = input('SEP angles are not matching up. Break here. ')
if round(sum(all_eps - mean(all_eps,axis=0)),5) != 0:
	print('')
	stop = input('EPS angles are not matching up. Break here. ')


#Run the weighted average
alt = mean(all_alt,axis=0)
ed_avg = sum(all_ed/all_err**2,axis=0)/sum(1./all_err**2,axis=0)
err_avg = 1./sqrt(sum(1./all_err**2,axis=0))
ax2.plot(ed_avg,alt,c='k',lw=2,label='This work')
ax2.errorbar(ed_avg,alt,xerr=err_avg,fmt='none',ecolor='k',capsize=0,lw=0.5)

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
	fig1.savefig(save_path+'ind_edp_'+obs+'.png',bbox_inches='tight')
	fig2.savefig(save_path+'final_edp_'+obs+'.png',bbox_inches='tight')

	#Save the data as a pickle. If no errors were thrown on the shared quantities, any of
	# profiles can be used (e.g., the 0th as I do here)
	with open(save_path+'final_edp_'+obs+'_'+representation+'.pickle','wb') as pckl:
	#	pickle.dump((all_ettx[0], all_etoccpt[0], all_etrx[0], all_utctx[0], \
#			all_utcoccpt[0],all_utcrx[0], all_lat[0], all_lon[0], all_sza[0], all_lst[0],\
#				all_sep[0], all_eps[0], alt+2575., ed_avg, err_avg),pckl)
		if representation == 'spherical':
			pickle.dump((all_ettx[0], all_etoccpt[0], all_etrx[0], all_utctx[0], \
						 all_utcoccpt[0], all_utcrx[0], all_lat[0], all_lon[0], all_sza[0], all_lst[0], \
					 all_sep[0], all_eps[0], alt + 60268., ed_avg, err_avg), pckl)
		elif representation == 'ellipsoidal':
			pickle.dump((all_ettx[0], all_etoccpt[0], all_etrx[0], all_utctx[0], \
					 all_utcoccpt[0], all_utcrx[0], all_lat[0], all_lon[0], all_sza[0], all_lst[0], \
					 all_sep[0], all_eps[0], alt + saturn_R[rep_point_loc], ed_avg, err_avg, rep_point_loc,saturn_R), pckl)
	pckl.close()

	file = open(save_path+'final_edp_'+obs+'_'+representation+'.txt','w')
	for i in range(size(ed_avg)):
		file.write('%.4f' % alt[i]+'\t\t'+'%.8f' % ed_avg[i]+'\t\t'+'%.8f' % err_avg[i]+'\n')
	file.close()
	print('Saved!')





#stop=input('stop')
