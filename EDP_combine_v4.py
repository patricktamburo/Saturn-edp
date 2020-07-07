
import pdb
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
import numpy as np
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
import pdb
import pandas
import scipy.ndimage 
import datetime

close('all')

#Specify the body, flyby and the observation ID
body = 'Saturn'
flyby = 'S07'
obs_type = 'N'
obs = flyby+obs_type

sss = 's10'

#Set paths
root_path = 'C:\\Users\\tambu\\Documents\\BU\\Science Projects\\SaturnEDP\\PD_MBP_files\\'
edp_path = root_path+'EDP\\'

#Get all the applicable directories
obs_list = [edp_path+body+'\\'+flyby+'\\'+x for x in os.listdir(edp_path+body+'\\'+flyby+'\\') if (os.path.isdir(edp_path+body+'\\'+flyby+'\\'+x)) and (obs[-1] in x)]

print('Found {} individual profiles:'.format(len(obs_list)))
time.sleep(0.5)
for i in range(len(obs_list)):
    print(obs_list[i].split('\\')[-1])
    time.sleep(0.1)
print('')

#Create size arrays that will be needed to trim and combine EDP
tops, bottoms = np.zeros(size(obs_list)), np.zeros(size(obs_list))

#Build the individual EDP figure
fig1 = figure(1)
ax1 = fig1.add_subplot(111)
ax1.axvline(0,c='k',ls='--')
ax1.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
ax1.set_ylabel('Altitude [km]',fontsize='x-large')
xticks(fontsize='large')
yticks(fontsize='large')
ax1.minorticks_on()
ax1.set_title(obs+' individual EDP',fontsize='x-large')

#For each, load in all of the EDP, parse the profile files, and plot
all_IDs = []
topsides = np.zeros(np.size(obs_list))
avg_topside, std_topside, rms_topside = np.zeros(np.size(obs_list)), np.zeros(np.size(obs_list)), np.zeros(np.size(obs_list))
alt_list, ed_list, rms_list = [], [], []
occpt_list, lat_list, lon_list, sza_list, lst_list, sep_list, eps_list = [], [], [], [], [], [], []
ettx_list, etoccpt_list, etrx_list = [], [], []
utctx_list, utcoccpt_list, utcrx_list = [], [], []

for i in range(size(obs_list)):
    all_IDs += [obs_list[i].split('\\')[-1]]
    idl_data = idlsave.read(obs_list[i]+'\\'+flyby+'_'+all_IDs[i]+'_edp.sav',verbose=0)    
    edp_df = pandas.read_csv(obs_list[i]+'\\'+flyby+'_'+all_IDs[i]+'_edp.txt', names=['Radius', 'Altitude', 'ED', 'RMS']) #Read in the EDP as a pandas dataframe.
    #rad = np.array(edp_df['Radius'])
    #Get EDP info
    #alt_list.append(np.array(edp_df['Altitude']))
    #ed_list.append(np.array(edp_df['ED']))
    #ms_list.append(np.array(edp_df['RMS']))

    alt_list.append(idl_data['rminusRsav'])
    ed_list.append(idl_data['edpsav'])
    rms_list.append(idl_data['electrondensity_errorsav'])
    #Get supporting information.
    occpt_list.append(idl_data['occptradiusarraysav'])
    lat_list.append(idl_data['occptlatarraysav'])
    lon_list.append(idl_data['occptlonarraysav'])
    sza_list.append(idl_data['occptszaarraysav'])
    lst_list.append(idl_data['occptlstarraysav'])
    sep_list.append(idl_data['sepanglearraysav'])
    eps_list.append(idl_data['epsanglearraysav'])

    ettx_list.append(idl_data['ettxarraysav'])
    etoccpt_list.append(idl_data['etoccptarraysav'])
    etrx_list.append(idl_data['etrxarraysav']) #Have to reverse this?

    utctx_list.append(idl_data['utctxarraysav'])
    utcoccpt_list.append(idl_data['utcoccptarraysav'])
    utcrx_list.append(idl_data['utcrxarraysav'])

    #Find the max/min altitudes of this profile.
    tops[i] = np.max(alt_list[i])
    bottoms[i] = np.min(alt_list[i])

    #Plot it. 
    ax1.plot(ed_list[i],alt_list[i],label=str(i)+', '+all_IDs[i])

    #Parse the profile file to get necessary info.
    with open(obs_list[i]+'/'+all_IDs[i]+'.profile','r') as input_file:
        lines = input_file.readlines()
    for j in range(len(lines)):
        if 'TOPSIDE' in lines[j]:
            topsides[i] = int(lines[j].split('=')[-1])

    #Print out topside information
    topside_inds = np.where(alt_list[i] >= topsides[i])[0]
    avg_topside[i] = np.mean(ed_list[i][topside_inds])
    std_topside[i] = np.std(ed_list[i][topside_inds])
    rms_topside[i] = np.sqrt(np.sum(ed_list[i][topside_inds]**2)/np.size(topside_inds))

print('')
if len(set(topsides)) != 1:
    print('ERROR: differing topside altitudes for different bands.')
else:
    xlim_save = ax1.get_xlim()
    top_ed_extent = ed_list[i][np.where(abs(alt_list[i]-topsides[0]) == min(abs(alt_list[i]-topsides[0])))[0][0]]

plt.plot([xlim_save[0], top_ed_extent], [topsides[0], topsides[0]], color='k', linestyle='--', zorder=0)
plt.xlim(xlim_save)
plt.tight_layout()
ax1.legend(loc='best',prop={'size':12})
plt.savefig(edp_path+body+'\\'+flyby+'\\'+obs+'_individual_profiles.png')

#Combine where there is overlap between all of the good EDP.
top_cut = np.min(tops)
bottom_cut = np.max(bottoms)

#See if the trim produces equally sized alt arrays
alt_sizes = np.zeros(len(obs_list),dtype=int)

for i in range(len(obs_list)):
    alt = alt_list[i]
    alt_sizes[i] = np.size(alt[np.where((alt>bottom_cut)&(alt<top_cut))[0]])

#Break if there are unequal altitude sizes.  
if len(set(alt_sizes)) != 1:
    print('')
    stop = input('Altitude sizes are not matching up (point 1). Break here. ')


#Average EDP figure
fig2 = figure(2)
ax2 = fig2.add_subplot(111)
ax2.axvline(0,c='k',ls='--')
ax2.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
ax2.set_ylabel('Altitude [km]',fontsize='x-large')
xticks(fontsize='large')
yticks(fontsize='large')
ax2.minorticks_on()
ax2.set_title(obs+' weighted average',fontsize='x-large')

#Create the arrays of everything we want to save
all_alt = np.zeros((len(obs_list),alt_sizes[0]))
all_ed = np.zeros((len(obs_list),alt_sizes[0]))
all_err = np.zeros((len(obs_list),alt_sizes[0]))

all_occpt = np.zeros((len(obs_list),alt_sizes[0]))
all_lat = np.zeros((len(obs_list),alt_sizes[0]))
all_lon = np.zeros((len(obs_list),alt_sizes[0]))
all_sza = np.zeros((len(obs_list),alt_sizes[0]))
all_lst = np.zeros((len(obs_list),alt_sizes[0]))
all_sep = np.zeros((len(obs_list),alt_sizes[0]))
all_eps = np.zeros((len(obs_list),alt_sizes[0]))

all_ettx = np.zeros((len(obs_list),alt_sizes[0]))
all_etoccpt = np.zeros((len(obs_list),alt_sizes[0]))
all_etrx = np.zeros((len(obs_list),alt_sizes[0]))

all_utctx =  [[] for x in range(len(obs_list))]
all_utcoccpt =  [[] for x in range(len(obs_list))]
all_utcrx = [[] for x in range(len(obs_list))]

all_utctx_datetime =  [[] for x in range(len(obs_list))]
all_utcoccpt_datetime =  [[] for x in range(len(obs_list))]
all_utcrx_datetime = [[] for x in range(len(obs_list))]

for i in range(len(obs_list)):
    alt_locs = np.where((alt_list[i]>bottom_cut)&(alt_list[i]<top_cut))[0]
    all_alt[i,:] = alt_list[i][alt_locs]
    all_ed[i,:] = ed_list[i][alt_locs]
    all_err[i,:] = rms_list[i][alt_locs]
    
    all_occpt[i,:] = occpt_list[i][alt_locs]
    all_lat[i,:] = lat_list[i][alt_locs]
    all_lon[i,:] = lon_list[i][alt_locs]
    all_sza[i,:] = sza_list[i][alt_locs]
    all_lst[i,:] = lst_list[i][alt_locs]
    all_sep[i,:] = sep_list[i][alt_locs]
    all_eps[i,:] = eps_list[i][alt_locs]

    all_ettx[i,:] = ettx_list[i][alt_locs]
    all_etoccpt[i,:] = etoccpt_list[i][alt_locs]
    all_etrx[i,:] = etrx_list[i][alt_locs]

    for j in range(len(alt_locs)):
        all_utctx[i].append(utctx_list[i][alt_locs][j].decode('utf-8'))
        all_utcoccpt[i].append(utcoccpt_list[i][alt_locs][j].decode('utf-8'))
        all_utcrx[i].append(utcrx_list[i][alt_locs][j].decode('utf-8'))
        all_utctx_datetime[i].append(datetime.datetime.strptime(all_utctx[i][j], '%Y %B %d %H:%M:%S.%f'))
        all_utcoccpt_datetime[i].append(datetime.datetime.strptime(all_utcoccpt[i][j], '%Y %B %d %H:%M:%S.%f'))
        all_utcrx_datetime[i].append(datetime.datetime.strptime(all_utcrx[i][j], '%Y %B %d %H:%M:%S.%f'))

all_utctx_datetime = np.array(all_utctx_datetime)
all_utcoccpt_datetime = np.array(all_utcoccpt_datetime)
all_utcrx_datetime = np.array(all_utcrx_datetime)

pdb.set_trace()
#Run a check on most of the shared quantities (ignoring utc times for now) to make sure everything matches up. 
if np.round(sum(all_alt - np.mean(all_alt, axis=0)), 5) != 0:
    print('')
    stop = input('Altitude values are not matching up (point 2). Break here. ')

plt.figure
fig4, ax4 = plt.subplots(nrows=3, ncols=3, figsize=(10,8))

if np.round(sum(all_lat - np.mean(all_lat, axis=0)), 5) != 0:
    print('')
    stop = input('Latitude values are not matching up (point 2). Break here. ')
else:
    ax4[0,0].plot(all_lat[0])
    ax4[0,0].set_title('Latitude')

if np.round(sum(all_lon - np.mean(all_lon, axis=0)), 5) != 0:
    print('')
    stop = input('Longitude values are not matching up (point 2). Break here. ')
else:
    ax4[0,1].plot(all_lon[0])
    ax4[0,1].set_title('Longitude')

if np.round(sum(all_sza - np.mean(all_sza, axis=0)), 5) != 0:
    print('')
    stop = input('SZA values are not matching up (point 2). Break here. ')
else:
    ax4[0,2].plot(all_sza[0])
    ax4[0,2].set_title('Occpt SZA')

if np.round(sum(all_lst - np.mean(all_lst, axis=0)), 5) != 0:
    print('')
    stop = input('LST values are not matching up (point 2). Break here. ')
else:
    ax4[1,0].plot(all_lst[0])
    ax4[1,0].set_title('Occpt LST')

if np.round(sum(all_sep - np.mean(all_sep, axis=0)), 5) != 0:
    print('')
    stop = input('SEP values are not matching up (point 2). Break here. ')
else:
    ax4[1,1].plot(all_sep[0])
    ax4[1,1].set_title('SEP angle')

if np.round(sum(all_eps - np.mean(all_eps, axis=0)), 5) != 0:
    print('')
    stop = input('EPS values are not matching up (point 2). Break here. ')
else:
    ax4[1,2].plot(all_eps[0])
    ax4[1,2].set_title('EPS angle')

if np.round(sum(all_ettx - np.mean(all_ettx, axis=0)), 3) != 0:
    print('')
    stop = input('ETTX values are not matching up (point 2). Break here. ')
else:
    ax4[2,0].plot(all_ettx[0])
    ax4[2,0].set_title('ETTX')

if np.round(sum(all_etoccpt - np.mean(all_etoccpt, axis=0)), 3) != 0:
    print('')
    stop = input('ETOCCPT values are not matching up (point 2). Break here. ')
else:
    ax4[2,1].plot(all_etoccpt[0])
    ax4[2,1].set_title('ETOCCPT')

if np.round(sum(all_etrx - np.mean(all_etrx, axis=0)), 3) != 0:
    print('')
    stop = input('ETRX values are not matching up (point 2). Break here. ')
else:
    ax4[2,2].plot(all_etrx[0])
    ax4[2,2].set_title('ETRX')

#Check UTC values.
for j in range(len(all_utctx_datetime[0])):
    if (all_utctx_datetime[:,0] - all_utctx_datetime.min()).mean().total_seconds() != 0:
        print('UTCTX values are not matching up, break here.')
        pdb.set_trace()
    if (all_utcoccpt_datetime[:,0] - all_utcoccpt_datetime.min()).mean().total_seconds() != 0:
        print('UTCOCCPT values are not matching up, break here.')
        pdb.set_trace()
    if (all_utcrx_datetime[:,0] - all_utcrx_datetime.min()).mean().total_seconds() != 0:
        print('UTCRX values are not matching up, break here.')
        pdb.set_trace()

plt.tight_layout()

#Run the weighted average
alt = np.mean(all_alt,axis=0)
ed_avg = sum(all_ed/all_err**2,axis=0)/sum(1./all_err**2,axis=0)
#err_avg = 1./np.sqrt(sum(1./all_err**2,axis=0)) #Old approach: algebraic combination of individual uncertainties. 
topside_inds_avg = np.where(alt > max(topsides))[0]
rms_avg = np.zeros(len(alt)) + np.sqrt(np.mean(ed_avg[topside_inds_avg]**2)) #New approach: rms of topside. 
ax2.plot(ed_avg,alt,c='k',lw=2,label='This work')
ax2.errorbar(ed_avg,alt,xerr=rms_avg,fmt='none',ecolor='k',capsize=0,lw=0.5)
xlim_save = ax2.get_xlim()
top_ed_extent = ed_avg[np.where(abs(alt-topsides[0]) == min(abs(alt-topsides[0])))[0][0]]
ax2.plot([xlim_save[0], top_ed_extent], [topsides[0], topsides[0]], color='k', linestyle='--', zorder=0)
ax2.set_xlim(xlim_save)
plt.tight_layout()
plt.savefig(edp_path+body+'\\'+flyby+'\\'+obs+'_average_profile.png')

#Smoothed EDP figure
fig3 = figure(4)
ax3 = fig3.add_subplot(111)
ax3.axvline(0,c='k',ls='--')
ax3.set_xlabel('Electron Density [cm^-3]',fontsize='x-large')
ax3.set_ylabel('Altitude [km]',fontsize='x-large')
xticks(fontsize='large')
yticks(fontsize='large')
ax3.minorticks_on()
ax3.set_title(obs+' smoothed weighted average',fontsize='x-large')

smoothed_ed = scipy.ndimage.filters.uniform_filter(ed_avg,size=11)
#Calculate RMS of smoothed_ed, ignoring the top 6 points
topside_inds_smooth = topside_inds_avg[0:len(topside_inds_avg)-6]
rms_smooth = np.zeros(len(alt)) + np.sqrt(np.mean(smoothed_ed[topside_inds_smooth]**2))
ax3.plot(smoothed_ed,alt,c='k',lw=2,label='This work')
xlim_save = ax3.get_xlim()
top_ed_extent = smoothed_ed[np.where(abs(alt-topsides[0]) == min(abs(alt-topsides[0])))[0][0]]
ax3.plot([xlim_save[0], top_ed_extent], [topsides[0], topsides[0]], color='k', linestyle='--', zorder=0)
ax3.set_xlim(xlim_save)
plt.tight_layout()
plt.savefig(edp_path+body+'\\'+flyby+'\\'+obs+'_smoothed_average_profile.png')

#Print error information on individual profiles, averaged profile, and smoothed average profile.
print('Average topside ED: {}'.format(np.mean(ed_avg[topside_inds_avg])))
for i in range(len(obs_list)):
    print('{} topside rms: {}'.format(obs_list[i].split('\\')[-1], all_err[i][0]))
print('Averaged profile topside rms: {}'.format(rms_avg[0]))
print('Smoothed averaged profile topside rms: {}'.format(rms_smooth[0]))

#Write out aven file. If you have made it to this point, you can choose any element of the all_xxxxx arrays and they will be representative of the entire
#averaged EDP. 

#ettx, etocc, etrx, utctx, utcocc, utcrx, occptradius, alt, occptlat, occptlon, occptsza, occptlst, occptsep, occpteps, avgelecden, avgelecdenerr
if len(flyby.split('S')[1]) == 2:
    ffff = 'S0'+flyby.split('S')[1]
aven_filename = sss+'_'+ffff.lower()+obs_type.lower()+'_satur_edp_v01_r00_test.csv'

output_path = 'C:\\Users\\tambu\\Documents\\BU\\Science Projects\\SaturnEDP\\PT_files\\PDS_data\\aven\\'+obs
if not os.path.exists(output_path):
    os.mkdir(output_path)

with open(output_path+'\\'+aven_filename, 'w') as f:
    for i in range(len(all_ettx[0])):
        f.write('{:19.12e}, {:19.12e}, {:19.12e}, {}, {}, {}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}\n'.format(all_ettx[0][i], all_etoccpt[0][i], all_etrx[0][i], '  '+all_utctx[0][i], '  '+all_utcoccpt[0][i], '  '+all_utcrx[0][i].replace('500002','500000'),  all_occpt[0][i], all_alt[0][i],
            all_lat[0][i], all_lon[0][i], all_sza[0][i], all_lst[0][i], all_sep[0][i], all_eps[0][i], ed_avg[i], rms_avg[i], smoothed_ed[i], rms_smooth[i]))
pdb.set_trace()