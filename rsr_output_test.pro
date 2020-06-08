PRO rsr_output_test
;restore x
RESTORE, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0233\s41saoi2008153_2120nnnx26rd.2a2_x\output.sav',/verbose
SFDUYEAROUT_X = SFDUYEAROUT
SFDUDOYOUT_X = SFDUDOYOUT
SFDUSECOUT_X = SFDUSECOUT
RFTOIFMHZOUT_X = RFTOIFMHZOUT
DDCLOMHZOUT_X = DDCLOMHZOUT
NCOFREQOUT_X = NCOFREQOUT
KPOSAAOUT_X = KPOSAAOUT
IGCOMPLEXAAOUT_X = IGCOMPLEXAAOUT

;restore k
RESTORE, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0232\s41saoi2008153_2120nnnk26rd.2b2_x\test_output_full.sav',/verbose
SFDUYEAROUT_K = SFDUYEAROUT
SFDUDOYOUT_K = SFDUDOYOUT
SFDUSECOUT_K = SFDUSECOUT
RFTOIFMHZOUT_K = RFTOIFMHZOUT
DDCLOMHZOUT_K = DDCLOMHZOUT
NCOFREQOUT_K = NCOFREQOUT
KPOSAAOUT_K = KPOSAAOUT
IGCOMPLEXAAOUT_K = IGCOMPLEXAAOUT

freqxr = (rftoifmhzout_x + ddclomhzout_x) * 1d6 - ncofreqout_x + kposaaout_x
freqkr = (rftoifmhzout_k + ddclomhzout_k) * 1d6 - ncofreqout_k + kposaaout_k
occfx = freqxr
occfk = freqkr
xdeltaf = occfx
kdeltaf = occfk
xfobs = occfx
c_mks     = 2.9979e8
me_mks    = 9.1094e-31
mp_mks    = 1.6726e-27
eps0_mks  = 8.8542e-12
e_mks     = 1.6022e-19

fdxfds = 11d/3d
fkdfds = 209d/15d
fdkfdx = 19d/5d
dxdt = (xdeltaf - kdeltaf/fdkfdx) * 8d * !pi * !pi * me_mks * eps0_mks *(xfobs) * c_mks / (e_mks^2) * 1d/(1d - fdkfdx^(-2d))
w = window(dimensions=[1500,1000])
p = plot(dxdt[0:800],/current)
STOP
END