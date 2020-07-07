import pdb
import idlsave
import os 
import pandas
import numpy as np
import datetime

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


#Output indn and aven files as csv.
for i in range(len(obs_list)):
    idl_data = idlsave.read(obs_list[i]+'\\'+flyby+'_'+obs_list[i].split('\\')[-1]+'_edp.sav',verbose=0)    
    if obs_type == 'N':
        aa = 'oi'
    elif obs_type == 'X':
        aa = 'oe'
    yyyy = str(idl_data['sfduyearoutxr'][0])
    ddd = str(idl_data['sfdudoyoutxr'][0])
    if len(ddd) == 1:
        ddd = '00'+ddd
    elif len(ddd) == 2:
        ddd = '0'+ddd
    
    hhmm =  str(idl_data['utcrxarraysav'][0]).split(':')[0][-2:] + str(idl_data['utcrxarraysav'][0]).split(':')[1]
    t = obs_list[i].split('\\')[-1].split('_')[0][-1].lower()
    bb = obs_list[i].split('\\')[-1].split('_')[1].lower()
    nn = obs_list[i].split('\\')[-1].split('_')[2]
    indn_filename = sss+'sa'+aa+yyyy+ddd+'_'+hhmm+'_'+t+'_'+bb+'_'+nn+'_satur_edp_v01_r00_test.csv'
    
    utc_tx = [idl_data['utctxarraysav'][i].decode("utf-8") for i in range(len(idl_data['utctxarraysav']))]
    utc_occ = [idl_data['utcoccptarraysav'][i].decode("utf-8") for i in range(len(idl_data['utcoccptarraysav']))]
    utc_rx = [idl_data['utcrxarraysav'][i].decode("utf-8") for i in range(len(idl_data['utcrxarraysav']))]

    output_path = 'C:\\Users\\tambu\\Documents\\BU\\Science Projects\\SaturnEDP\\PT_files\\PDS_data\\indn\\'+obs
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    with open(output_path+'\\'+indn_filename,'w') as f:
        for i in range(len(utc_tx)):
            f.write('{:19.12e}, {:19.12e}, {:19.12e}, {}, {}, {}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}, {:19.12e}\n'.format(idl_data['ettxarraysav'][i], idl_data['etoccptarraysav'][i], idl_data['etrxarraysav'][i], '  '+utc_tx[i], '  '+utc_occ[i], '  '+utc_rx[i].replace('500002','500000'),  idl_data['occptradiusarraysav'][i],idl_data['rminusRsav'][i],idl_data['occptlatarraysav'][i],
                                  idl_data['occptlonarraysav'][i],idl_data['occptszaarraysav'][i],idl_data['occptlstarraysav'][i],
                                  idl_data['sepanglearraysav'][i],idl_data['epsanglearraysav'][i],idl_data['dxdt_save'][i],
                                  idl_data['dxdt'][i],idl_data['newbigx'][i],idl_data['edpsav'][i],idl_data['electrondensity_errorsav'][i]))
