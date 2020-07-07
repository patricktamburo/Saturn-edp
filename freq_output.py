import pdb
import os
import idlsave
import numpy as np

cors = ['0106', '0107', '0108']
cors_path = 'C:\\Users\\tambu\\Documents\\BU\\Science Projects\\SaturnEDP\\PD_MBP_files\\Output\\Saturn\\'
for k in range(len(cors)):
    path = cors_path+'cors_'+cors[k]+'\\'
    sub_dirs = [x for x in os.listdir(path) if os.path.isdir(path+x)]
    for j in range(len(sub_dirs)):
        rsr_data_path = path+sub_dirs[j]+'\\full_output\\output.sav'
        idl_data = idlsave.read(rsr_data_path, verbose=0)
        sss = sub_dirs[j][0:3]
        tt = sub_dirs[j][3:5]
        aa = sub_dirs[j][5:7]
        yyyy = str(idl_data['sfduyearout'][0])
        ddd = str(idl_data['sfdudoyout'][0])
        if len(ddd) == 1:
            ddd = '00'+ddd
        elif len(ddd) == 2:
            ddd = '0'+ddd
        h = int(idl_data['sfdusecout'][0] / 3600)
        if h < 10:
            hh = '0'+str(h)
        else:
            hh = str(h)
        m = int(np.floor((idl_data['sfdusecout'][0] - (h * 3600)) / 60))
        if m < 10:
            mm = '0'+str(m)
        else:
            mm = str(m)
        xuu = sub_dirs[j][19:22]
        drr = sub_dirs[j][22:25]
        p = sub_dirs[j][25]
        q = sub_dirs[j][26]
        rcs = sub_dirs[j][28:]
        freq_filename = sss+tt+aa+yyyy+ddd+'_'+hh+mm+xuu+drr+p+q+'_'+rcs+'_freq_v01_r00_test.csv'
        output_path = 'C:\\Users\\tambu\\Documents\\BU\Science Projects\\SaturnEDP\\PT_files\\PDS_data\\freq\\'+cors[k]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        with open(output_path + '\\' + freq_filename, 'w') as f:
            for i in range(len(idl_data['sfduyearout'])):
                f.write('{:10d}, {:9d}, {:19.3f}, {:9d}, {:9d}, {:19.12e}, {:19.12e}, {:19.12e}, {:4d}, {:4d}\n'.format(idl_data['sfduyearout'][i], idl_data['sfdudoyout'][i], idl_data['sfdusecout'][i], idl_data['rftoifmhzout'][i], idl_data['ddclomhzout'][i], idl_data['ncofreqout'][i], idl_data['kposaaout'][i], abs(idl_data['igcomplexaaout'][i]), 9, 9))
