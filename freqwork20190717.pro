pro freqwork20190717

  ; Paul Withers, Paul Dalba, 2017.09.13
  ; Astronomy Department, Boston University
  
  ;Find preliminary electron density profiles using outputs from rsr_proc.pro (run on the SCC)
  ; and simpleocc20170913.pro. Can be run for any occultation provided the correct path info
  ; below.
  
  ; Paul Dalba, 2017.11.12
  ; This version only differs from 20170913 by the input and output of data.
  ; This version streamlines the I/O.
  
  ; Paul Dalba, 2018.03.02
  ; This version only differs from 20171112 by the output of data.
  ; This version saves the indexstart(end) values to the RSR directories.
  ; This version also saves a large IDL save file that includes all pieces
  ;   of data needed to make archival files for the PDS.
  ; This version also calculates and saves the topside STDDEV of EDP.
  
  ; Pat Tamburo, 2018.10.16
  ; This version is built to handle occultations of Saturn.
  
  ; Pat Tamburo, 2019.02.05
  ; This version calculates altitudes using the Anderson-Schubert 1-bar surface for Saturn. 
  
  
itloop = 'b'
; 'a' - primary purpose is to identify time periods of good data for use by simpleocc, but also does the spherical inversion only
; 'b' - now that altitudes and pathlengths have been found by simpleocc, do all possible inversions


print, 'Careful about re-running this, as it will append the flag_indices.txt files'
print, ' unless they are commented out (they are currently commented out)'
;stop
  
;Set rootpath for the occultation data
;rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Output/'
rootpath = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\'
;;rootpath = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/'
;;;rootpath = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/'
  ;save_rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/EDP/'
save_rootpath = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\EDP\'
;;save_rootpath = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/EDP/'
;;;  save_rootpath = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/EDP/'

;Flip this switch to save figures
save_figs = 1
if save_figs eq 1 then DEVICE, DECOMPOSED=1
  
  
;You need to input the body, the flyby, the observation, the bands and the station.
; Together, these identify a unique profile. The information regarding the analysis
; of this profile are then read in from the unique profile file.
;  body =     'Titan'    ;Options are 'Titan' or 'Saturn' at present
;  flyby =    'T12'
;  obs =      'T12X'
;  bands =    'XK'       ;Options are 'SX' and 'XK' at present
;  station =  '25'
  
body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
flyby =    'S07'
obs =      'S07N'
bands =    'SX'       ;Options are 'SX' and 'XK' at present
station =  '14'

;S10N
;body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
;flyby =    'S10'
;obs =      'S10N'
;bands =    'SX'       ;Options are 'SX' and 'XK' at present
;station =  '14'

;;S56X
;body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
;flyby =    'S56'
;obs =      'S56X'
;bands =    'SX'       ;Options are 'SX' and 'XK' at present
;station =  '63'

;S51X
;body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
;flyby =    'S51'
;obs =      'S51X'
;bands =    'SX'       ;Options are 'SX' and 'XK' at present
;station =  '63'

;S54X
;body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
;flyby =    'S54'
;obs =      'S54X'
;bands =    'SX'       ;Options are 'SX' and 'XK' at present
;station =  '63'

;S75X
;body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
;flyby =    'S75'
;obs =      'S75X'
;bands =    'SX'       ;Options are 'SX' and 'XK' at present
;station =  '63'
 
baseline_dxdt = 0 ; Either 0 or 1. Specifies whether or not to 
; baseline the time series of dxdt via a linear fit. Necessary
; for some Saturn observations: S14N, S70N. 
  
;Open profile file
openr, lun, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/'+obs+'_'+bands+'_'+station+'.profile', /get_lun
;;; first set of three
;;;  openr, lun, save_rootpath+body+'/'+obs+'_'+bands+'_'+station+'/'+obs+'_'+bands+'_'+station+'.profile', /get_lun
;;; second set of three

;;  openr, lun, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\'+obs+'_'+bands+'_'+station+'.profile', /get_lun
  
;The first line is the header, throw it out
hdr = strarr(1)
readf, lun, hdr
;The next 17 lines are the data
EDP_profile = strarr(17)
readf, lun, EDP_profile
close, lun
free_lun, lun
  
;Process the input profile
cors_a = (strsplit(EDP_profile[5],'=',/extract))[1]
cors_b = (strsplit(EDP_profile[6],'=',/extract))[1]
product_a = (strsplit(EDP_profile[7],'=',/extract))[1]
product_b = (strsplit(EDP_profile[8],'=',/extract))[1]
manual_start = (strsplit(EDP_profile[9],'=',/extract))[1]
manual_end = (strsplit(EDP_profile[10],'=',/extract))[1]
bad_cut = (strsplit(EDP_profile[11],'=',/extract))[1]
trim_start = (strsplit(EDP_profile[12],'=',/extract))[1]
trim_end = (strsplit(EDP_profile[13],'=',/extract))[1]
topside = (strsplit(EDP_profile[14],'=',/extract))[1]
median_bool = (strsplit(EDP_profile[15],'=',/extract))[1]
linear_bool = (strsplit(EDP_profile[16],'=',/extract))[1]
  
;Some of these input needs to be converted
manual_start = double(manual_start)
manual_end = double(manual_end)
bad_cut = double(bad_cut)
trim_start = double(trim_start)
trim_end = double(trim_end)
topside = double(topside)
median_bool = fix(median_bool)
linear_bool = fix(linear_bool)
  
;The trim_start/end keywords are equivalent to the manual_start/end keywords and all
; .profile files have been changed to only use manual_start/end (i.e. trim_start/end is obsolete)
; If anything other than -1e5 is set for trim_start/end, break here
if trim_start ne -1e5 or trim_end ne -1e5 then begin
  print, 'Trim_start/end keyword has been called, despite being obsolete.'
  stop
endif
  
if bands eq 'SX' then begin
  scors = cors_a
  sproductid = product_a
  xcors = cors_b
  xproductid = product_b
endif
if bands eq 'XK' then begin
  xcors = cors_a
  xproductid = product_a
  kcors = cors_b
  kproductid = product_b
endif
  
;Restore the X-band .sav file from rsr_proc.pro for this data set
;;;restore, filename=rootpath+'/cors_'+xcors+'/'+xproductid+'/full_output/output.sav'
;;  restore, filename=rootpath+body+'/cors_'+xcors+'/'+xproductid+'/full_output/output.sav'
restore, filename=rootpath+body+'\cors_'+xcors+'\'+xproductid+'\full_output\output.sav'
sfduyearoutxr = sfduyearout
sfdudoyoutxr = sfdudoyout
sfdusecoutxr = sfdusecout
rftoifmhzoutxr = rftoifmhzout
ddclomhzoutxr = ddclomhzout
ncofreqoutxr = ncofreqout
kposaaoutxr = kposaaout
igcomplexaaoutxr = igcomplexaaout
  
;Calculate the sky frequency
freqxr = (rftoifmhzoutxr + ddclomhzoutxr) * 1d6 - ncofreqoutxr + kposaaoutxr
  
;Restore the S-band .sav file from rsr_proc.pro for this data set if the bands call for it
if bands eq 'SX' then begin
  ;restore, filename=rootpath+'/cors_'+scors+'/'+sproductid+'/full_output/output.sav'
  ;;    restore, filename=rootpath+body+'/cors_'+scors+'/'+sproductid+'/full_output/output.sav'
  restore, filename=rootpath+body+'\cors_'+scors+'\'+sproductid+'\full_output\output.sav'
  sfduyearoutsr = sfduyearout
  sfdudoyoutsr = sfdudoyout
  sfdusecoutsr = sfdusecout
  rftoifmhzoutsr = rftoifmhzout
  ddclomhzoutsr = ddclomhzout
  ncofreqoutsr = ncofreqout
  kposaaoutsr = kposaaout
  igcomplexaaoutsr = igcomplexaaout
  
  ;Calculate the sky frequency
  freqsr = (rftoifmhzoutsr + ddclomhzoutsr) * 1d6 - ncofreqoutsr + kposaaoutsr
endif
  
;Restore the K-band .sav file from rsr_proc.pro for this data set if the bands call for it
if bands eq 'XK' then begin
  restore, filename=rootpath+body+'/cors_'+kcors+'/'+kproductid+'/full_output/output.sav'
  ;;    restore, filename=rootpath+body+'\cors_'+kcors+'\'+kproductid+'\full_output\output.sav',/verbose
  sfduyearoutkr = sfduyearout
  sfdudoyoutkr = sfdudoyout
  sfdusecoutkr = sfdusecout
  rftoifmhzoutkr = rftoifmhzout
  ddclomhzoutkr = ddclomhzout
  ncofreqoutkr = ncofreqout
  kposaaoutkr = kposaaout
  igcomplexaaoutkr = igcomplexaaout
  
  ;Calculate the sky frequency
  freqkr = (rftoifmhzoutkr + ddclomhzoutkr) * 1d6 - ncofreqoutkr + kposaaoutkr
endif
  
;;; XXX Use X-band always to minimize places where mismatches can arise  
;Restore the SPICE data for this occultation. Choose any band, if they are the
; same receiver, it should not matter.
;;  if strmatch(bands, 'S*') then restore, rootpath+'/cors_'+scors+'/'+sproductid+'/simpleocc_output_itloop0.sav'
;;  if strmatch(bands, 'S*') then restore, rootpath+'/cors_'+scors+'/'+sproductid+'/full_output/simpleocc_output_pat.sav'
;;  if strmatch(bands, 'S*') then restore, rootpath+body+'/cors_'+scors+'/'+sproductid+'/full_output/simpleocc_output.sav'
;;  if strmatch(bands, 'S*') then restore, rootpath+body+'\cors_'+scors+'\'+sproductid+'\full_output\simpleocc_output_pat.sav',/verbose;
;;  if strmatch(bands, 'X*') then restore, rootpath+'/cors_'+xcors+'/'+xproductid+'/simpleocc_output_itloop0.sav'
;;  if strmatch(bands, 'X*') then restore, rootpath+'/cors_'+xcors+'/'+xproductid+'/full_output/simpleocc_output_pat.sav'
;;  if strmatch(bands, 'X*') then restore, rootpath+body+'/cors_'+xcors+'/'+xproductid+'/full_output/simpleocc_output.sav'
;;  if strmatch(bands, 'X*') then restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\full_output\simpleocc_output_pat.sav',/verbose

;;; XXX Use X-band always to minimize places where mismatches can arise  
;restore, rootpath+'/cors_'+xcors+'/'+xproductid+'/simpleocc_output_itloop0.sav'
restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop0.sav'

 if strmatch(obs, '*X') then begin
  ;;; DEAL WITH MY EGRESS REVERSAL
  itloop0occptradiusarray = reverse(itloop0occptradiusarray)
  itloop0occptlatarray = reverse(itloop0occptlatarray)
  itloop0occptlonarray = reverse(itloop0occptlonarray)
  itloop0occptszaarray = reverse(itloop0occptszaarray)
  itloop0occptlstarray = reverse(itloop0occptlstarray)
  itloop0sepanglearray = reverse(itloop0sepanglearray)
  itloop0epsanglearray = reverse(itloop0epsanglearray)
  itloop0ettxarray = reverse(itloop0ettxarray)
  itloop0etoccptarray = reverse(itloop0etoccptarray)
  itloop0etrxarray = reverse(itloop0etrxarray)
  itloop0utctxarray = reverse(itloop0utctxarray)
  itloop0utcoccptarray = reverse(itloop0utcoccptarray)
  itloop0utcrxarray = reverse(itloop0utcrxarray)
endif

keeptrackofutcrx = itloop0utcrxarray

;Determine UTC seconds for ERT from SPICE data
;real_ertutcsec = double(strmid(real_ertutc,12,2))*60d*60d + $
;double(strmid(real_ertutc,15,2))*60d + $
;double(strmid(real_ertutc,18,20))
  
;Use the times from SPICE and from the freq vs. time output to determine point radius for each freq value.
; The freq times and simpleocc times will match
if strmatch(bands, 'S*') then begin 
  ;if n_elements(sfdusecoutsr) ne n_elements(sfdusecoutxr) then stop 
  ;if max(abs( sfdusecoutsr - sfdusecoutxr )) ne 0. then stop
  ; if stop, then verify start times agree; that should be sufficient for index selections
  occtimesec = sfdusecoutsr
endif

if strmatch(bands, 'X*') then begin
  if n_elements(sfdusecoutxr) ne n_elements(sfdusecoutkr) then stop 
  if max(abs( sfdusecoutxr - sfdusecoutkr )) ne 0. then stop
  ;; if stop, then verify start times agree; that should be sufficient for index selections
  occtimesec = sfdusecoutxr
endif
  
; Added 3/21/19 *******
if (max(occtimesec) ge 7e4) and (min(occtimesec) lt 1e4) then begin
  aa = where(occtimesec le 4e4)
  occtimesec[aa] = occtimesec[aa] + 86400.
endif
  
;occptr = interpol(occptradiusarray,real_ertutcsec,occtimesec)
occptr = itloop0occptradiusarray
occptlat = itloop0occptlatarray

;Calculate R, the ellipsoidal represntation of Saturn's surface, as a function of latitude
a = 60268.
b = 54364.
  
saturn_R = a*b/sqrt((b*cos(itloop0occptlatarray*!pi/180.))^2.+(a*sin(itloop0occptlatarray*!pi/180.))^2.)
  
;Calculate occptradius - R
r_minus_R = itloop0occptradiusarray - saturn_R
    
;Time relative to start of observations
occtimesecadj = occtimesec - occtimesec[0]

;Define the occultation frequency values
occfx = freqxr
if bands eq 'SX' then occfs = freqsr
if bands eq 'XK' then occfk = freqkr
   
;Constants
c_mks     = 2.9979e8
me_mks    = 9.1094e-31
mp_mks    = 1.6726e-27
eps0_mks  = 8.8542e-12
e_mks     = 1.6022e-19

fdxfds = 11d/3d
fkdfds = 209d/15d
fdkfdx = 19d/5d

;The time resolution is 1 second.
dt = 1D

if strmatch(bands,'S*') then peakpower = abs(igcomplexaaoutsr)^2
if strmatch(bands,'X*') then peakpower = abs(igcomplexaaoutxr)^2
;Plot the max power at peak time series
window, 0
plot, occtimesec, peakpower, /ylog, xtitle='occtimesec', ytitle='Power of power spectrum peak'
  
;Plot the time vs. the pt radius
window, 1
plot, occtimesec, occptr-60268, xtitle='occtimesec', ytitle='occptradius - 60268'
;; XXX this occptr-60268 needs to be adjusted to oblate planet
;  hline, 0  ;the surface
;  hline, 1500, linestyle=2  ;approx peak of ionospheric ED
  
oplot, [-1,1]*1e9, [0,0]
oplot, [-1,1]*1e9, [0,0]+1500, linestyle=2

;; XXX why SX uses lower freq S, XK uses higher freq K?
if bands eq 'SX' then mixdownfreq = kposaaoutsr
if bands eq 'XK' then mixdownfreq = kposaaoutkr
  
window,7
plot,occtimesec,mixdownfreq,yr=[-1,1]*1e2

;Develop an algorithm to identify the bad portions of the data
;First, where occtimesec is below zero is not to be trusted.
if strmatch(obs, '*N') then begin
  ;This for an entry observation
  ;indexend = (where(occptr-2575. le 300.))[0] - 1
  indexend = (where(occptr-60268. le 300.))[0] - 1
  ;; XXX this occptr-60268 needs to be adjusted to oblate planet
    
  ;Now use a weighted average of the max power to estimate the start of the data
  peakpowerthresh = 10.^(0.9*max(alog10(peakpower))+0.1*min(alog10(peakpower)))
  good_locs = where(peakpower ge peakpowerthresh)
  good_locs_deriv = deriv(good_locs)
  aa = (where(ts_smooth(good_locs_deriv,3) gt 2*dt))
  indexstart = good_locs(aa[n_elements(aa)-1])
  ;;    indexstart = good_locs(   (where(ts_smooth(good_locs_deriv,3) gt 2*dt))[-1]    )
    
  ;PT added this to try and automatically identify start, fails for some occs.
  ;    if indexstart gt indexend then begin
  ;      ind = -2
  ;      while indexstart gt indexend do begin
  ;        indexstart = good_locs((where(ts_smooth(good_locs_deriv,3) gt 2*dt))[ind])
  ;        ind--
  ;      endwhile
  ;    endif
  indexstart = (where(peakpower ge peakpowerthresh))[0] ;***This is Paul D.'s original criterion!
  ;This determines the occtimesec of the ionospheric peak
  ;ionpeakocctimesec = occtimesec[(where(occptr-2575 le 1500.))[0]-1]
  ionpeakocctimesec = occtimesec[(where(occptr-60268 le 1500.))[0]-1]
  ;; XXX this occptr-60268 needs to be adjusted to oblate planet
    
  ;Manual offset
  indexend = indexend
  indexstart = indexend-650
  print, indexstart, indexend
  ;indexstart = indexstart+50 ;Paul's version
  ; indexend = indexend-300
  ; indexstart = indexstart-100
endif
   
if strmatch(obs, '*X') then begin
  ;This for an exit observation
  ;indexstart = (where(occptr-2575. le 300.))[-1] + 1
  ;;    indexstart = (where(occptr-60268. le 300.))[-1] + 1
  aa = (where(occptr-60268. le 300.))
  indexstart = aa[n_elements(aa)-1] + 1
  ;; XXX this occptr-60268 needs to be adjusted to oblate planet

  ;Now use a weighted average of the max power to estimate the end of the data
  peakpowerthresh = 10.^(0.9*max(alog10(peakpower))+0.1*min(alog10(peakpower)))
  good_locs = where(peakpower ge peakpowerthresh)
  good_locs_deriv = deriv(good_locs)
  ;If the spacing between good data locations exceeds the time resolution of the data,
  ; you've hit a bad portion of data. Make indexend the first place where this happens.
  ;indexend = good_locs((where(ts_smooth(good_locs_deriv,3) gt 2*dt))[0]) ; dodge S51X...
  if (obs ne 'S51X') and (obs ne 'S54X') and (obs ne 'S75X') then $
  indexend = good_locs((where(ts_smooth(good_locs_deriv,3) gt 2*dt))[0]) else indexend = indexstart+1000
    
  ;if indexstart gt indexend then begin
  ;  ind = 1
  ;  while indexstart gt indexend do begin
  ;    indexstart = good_locs((where(ts_smooth(good_locs_deriv,3) gt 2*dt))[ind])
  ;    ind++
  ;  endwhile
  ;endif
    
  ;;indexend = (where(peakpower ge peakpowerthresh))[-1] ***This is Paul D.'s criteria for Titan!
  ;;Paul's approach found the first good point, the last good point, and assumed everything
  ;;  in between was good data. This works for Titan, which doesn't have rings, but fails
  ;;  for Saturn. This new approach identifies the first good point, then the first bad point
  ;;  that occurs after the good point to identify the range.
    
  ;This determines the occtimesec of the ionospheric peak
  ;ionpeakocctimesec = occtimesec[(where(occptr-2575 le 1500.))[-1]+1]
  ;;    ionpeakocctimesec = occtimesec[(where(occptr-60268 le 1500.))[-1]+1]
  aa = (where(occptr-60268. le 300.))
  ionpeakocctimesec = occtimesec[aa[n_elements(aa)-1] + 1]
  ;; XXX this occptr-60268 needs to be adjusted to oblate planet
    
  ;Manual offsets
  ;indexend = indexend-50 ;This is Paul D.'s
  indexstart = indexstart - 500
  indexend = indexend - 1000

  ;; XXX HARD CODED 
  if obs eq 'S56X' then indexend = indexstart + 1000
  if obs eq 'S51X' then indexend = indexstart + 1000
    
endif
  
;Sometimes, special adjustment of the indexstart and indexend is required.
if (manual_start ne -1e5) and (manual_end ne -1e5) then begin
  indexstart = manual_start
  indexend = manual_end
endif
  
;Save the final values of indexstart and indexend in the RSR directories
; for each RSR file used in this call. Either create the file to save, or
; append it if it already exists
;product_a_flag_files = file_search(rootpath+body+'/cors_'+cors_a+'/'+product_a+'/full_output/flag_indices.txt')
;if strmatch(product_a_flag_files,'*flag_indices.txt') then begin
;  openu, lun, rootpath+body+'/cors_'+cors_a+'/'+product_a+'/full_output/flag_indices.txt', /get_lun, /append
;  if strmatch(obs,'*N') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    0'
;  if strmatch(obs,'*X') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    1'
;  close, lun
;  free_lun, lun
;endif else begin
;  openw, lun, rootpath+body+'/cors_'+cors_a+'/'+product_a+'/full_output/flag_indices.txt', /get_lun
;  printf, lun, '      indexstart            indexend            ingress (0) or egress (1)'
;  if strmatch(obs,'*N') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    0'
;  if strmatch(obs,'*X') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    1'
;  close, lun
;  free_lun, lun
;endelse
;
;product_b_flag_files = file_search(rootpath+body+'/cors_'+cors_b+'/'+product_b+'/full_output/flag_indices.txt')
;if strmatch(product_b_flag_files,'*flag_indices.txt') then begin
;  openu, lun, rootpath+body+'/cors_'+cors_b+'/'+product_b+'/full_output/flag_indices.txt', /get_lun, /append
;  if strmatch(obs,'*N') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    0'
;  if strmatch(obs,'*X') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    1'
;  close, lun
;  free_lun, lun
;endif else begin
;  openw, lun, rootpath+body+'/cors_'+cors_b+'/'+product_b+'/full_output/flag_indices.txt', /get_lun
;  printf, lun, '      indexstart            indexend            ingress (0) or egress (1)'
;  if strmatch(obs,'*N') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    0'
;  if strmatch(obs,'*X') then printf, lun, string(indexstart)+'          '+string(indexend)+'                    1'
;  close, lun
;  free_lun, lun
;endelse


print, 'indexstart = '+string(indexstart)
print, 'indexend   = '+string(indexend)

wset, 0
oplot, [1,1]*occtimesec[indexstart], [-1,1]*1e9, color=255
oplot, [1,1]*occtimesec[indexend], [-1,1]*1e9, color=255
;  vline, occtimesec[indexstart], color=255
;  vline, occtimesec[indexend], color=255
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/power_peak_time_series.jpeg', /jpeg
;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\power_peak_time_series.jpeg', /jpeg
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\power_peak_time_series_pat.jpeg', /jpeg
  
wset, 1
oplot, [1,1]*occtimesec[indexstart], [-1,1]*1e9, color=255
oplot, [1,1]*occtimesec[indexend], [-1,1]*1e9, color=255
;  vline, occtimesec[indexstart], color=255
;  vline, occtimesec[indexend], color=255
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/occptradius_vs_occtimesec.jpeg', /jpeg
;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\occptradius_vs_occtimesec.jpeg', /jpeg
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\occptradius_vs_occtimesec_pat.jpeg', /jpeg
  
wset, 7
oplot, [1,1]*occtimesec[indexstart], [-1,1]*1e9, color=255
oplot, [1,1]*occtimesec[indexend], [-1,1]*1e9, color=255
;  vline, occtimesec[indexstart], color=255
;  vline, occtimesec[indexend], color=255
    
; !!!!!!-----Manual adjustment required-----!!!!!!
;Trim portions of the freq array that contain bad data (i.e. neutral atm, planet occulted)
;if obs eq 'T12N' then bb = where((occtimesec gt 4200)and(occtimesec lt 4660))
;if obs eq 'T12X' then bb = where((occtimesec gt 5325)and(occtimesec lt 7200))
;if obs eq 'T14N' then bb = where((occtimesec gt 4.798e4)and(occtimesec lt 4.8715e4))
;if obs eq 'T14X' then bb = where((occtimesec gt 4.935e4)and(occtimesec lt 5.011e4))
;if obs eq 'T27N' then bb = where((occtimesec gt 4200)and(occtimesec lt 4862))
;if obs eq 'T27X' then bb = where((occtimesec gt 5715)and(occtimesec lt 6385))
;if obs eq 'T31N' then bb = where((occtimesec gt 6.7495e4)and(occtimesec lt 7.0405e4))
;if obs eq 'T31X' then bb = where((occtimesec gt 7.191e4)and(occtimesec lt 7.503e4))
  
;Maybe comment out...
;
;
;Plot the frequency difference versus point radius trimmed and untrimmed
;window, 1
;plot, occfs - 3d/11d * occfx, occptr, xra=[-1,1]*1e1 ;, yra=[2e3,5e3]
;oplot, occfs[bb] - 3d/11d * occfx[bb], occptr[bb]
;plot, abs(occfs - 3d/11d * occfx), occptr, /xlog ;, yra=[2e3,5e3]
;oplot, abs(occfs[bb] - 3d/11d * occfx[bb]), occptr[bb], color=255
  
;plot the point radius vs time

;window, 2
;plot, sfdusecout, occfs - 3d/11d * occfx, yra=[-1,1]*1e1;, xra=[7.49,7.51]*1e4
;stop  ;stop for inspection of bb

;Instead, make a so-called 'salt and pepper' plot (Withers 2014, Fig. 6)
;diff = occfs - 3d/11d * occfx
;window, 2
;plot, abs(diff[where(diff gt 0d)]), occptr[where(diff gt 0d)], psym=1, /xlog ;, yra=[2e3,5e3]
;oplot, abs(diff[where(diff lt 0d)]), occptr[where(diff lt 0d)], psym=1, color=255
;oplot, abs(diff[bb]), occptr[bb], psym=2, color=150

  
;Trim and re-define
keeptrackofutcrx = keeptrackofutcrx(indexstart:indexend)
xdeltaf = occfx(indexstart:indexend)
if bands eq 'SX' then begin
  sdeltaf = occfs(indexstart:indexend)
  sfobs = occfs(indexstart:indexend) ;Transmitted frequency (S)
endif
  
if bands eq 'XK' then begin
  kdeltaf = occfk(indexstart:indexend)
  xfobs = occfx(indexstart:indexend) ;Transmitted frequency (X)
endif
  
xa = occptr(indexstart:indexend) * 1D3  ;Convert point radius to meters
occtimesecsub = occtimesec(indexstart:indexend)  ;the subset of occsectime points

;Save the value of sdeltaf and xdeltaf
;save, filename='/Users/paul/Desktop/delta_freqs.sav', sdeltaf, xdeltaf, occtimesecsub
;stop
  
  
;Solve proposal equation 3 for dx/dt where x=TEC=integral(N dl)
if bands eq 'SX' then dxdt = (sdeltaf - xdeltaf/fdxfds) * 8d * !pi * !pi * me_mks * eps0_mks *(sfobs) * c_mks / (e_mks^2) * 1d/(1d - fdxfds^(-2d))
;if bands eq 'SK' then dxdt = (sdeltaf - kdeltaf/fdkfds) * 8d * !pi * !pi * me_mks * eps0_mks *(sfobs) * c_mks / (e_mks^2) * 1d/(1d - fdkfds^(-2d))
if bands eq 'XK' then dxdt = (xdeltaf - kdeltaf/fdkfdx) * 8d * !pi * !pi * me_mks * eps0_mks *(xfobs) * c_mks / (e_mks^2) * 1d/(1d - fdkfdx^(-2d))

;Save a version of dxdt before it is potentially altered below
dxdt_save = dxdt
  
;If baseline_dxdt is specified, make a linear fit to the time-series
; of dxdt, and subtract it out. This is relevant for S14N,S70N
if baseline_dxdt eq 1 then begin
  fit_end = 500
  degree = 2
  pfit = POLY_FIT(occtimesecsub[0:fit_end], dxdt[0:fit_end],degree)
  fit = dblarr(n_elements(occtimesecsub))
  for fit_ind = 0,degree do begin
    fit = fit+pfit[fit_ind]*occtimesecsub^(fit_ind)
  endfor
  ;Plot the fit
  ;;    w = window(dimensions=[2000,1500])
  ;;    p = plot(occtimesecsub,dxdt,/current,layout=[1,2,1],margin=[0.1,0.1,0.05,0.1],thick=3,title='dxdt + fit')
  ;;    p = plot(occtimesecsub[0:fit_end],fit,/overplot,color='red',layout=[1,2,1],thick=3)
  ;Subtract off baseline.
  ;;    p = plot(occtimesecsub[0:fit_end],dxdt-fit,$
  ;;      color='red',layout=[1,2,2],/current,margin=[0.1,0.1,0.05,0.1],thick=3,title='dxdt baselined')
  dxdt = dxdt - fit
endif

;Plot the dxdt vs. time and identify the ionosphere
window, 2
plot, occtimesecsub, dxdt,  xtitle='occtimesec', ytitle='dxdt', yrange=[-1,1]*1e15
oplot, occtimesecsub, smooth(dxdt,21), color=255
oplot, [1,1]*ionpeakocctimesec, [-1,1]*1e9, linestyle=2
oplot, [-1,1]*1e9, [0,0], color=255*256L
;  vline, ionpeakocctimesec,linestyle=2
;  hline, 0, color=255*256L
  
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/dxdt_vs_occtimesecsub.jpeg', /jpeg
;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\dxdt_vs_occtimesecsub.jpeg', /jpeg
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\dxdt_vs_occtimesecsub_pat.jpeg', /jpeg
  
; Test the effect of smoothing dxdt on the final EDP.
;  smooth_val = 9
;  dxdt = smooth(dxdt,smooth_val)
;  smooth_str = string(smooth_val)
  
;Plot the xa vs. dxdt and identify the ionosphere
window, 3
;plot, dxdt, xa/1e3 - 2575., xtitle='dxdt', ytitle='xa-2575', xrange=[-1,1]*5e14
plot, dxdt, xa/1e3 - 60268., xtitle='dxdt', ytitle='xa-60268', xrange=[-1,1]*2e16,yrange=[0,1.5e4]
;oplot, smooth(dxdt,21), xa/1e3 - 2575., color=255
;  xyouts, 1e16,1600,strcompress('DXDT smoothing = '+smooth_str), charsize=2
oplot, smooth(dxdt,21), xa/1e3 - 60268, color=255
print, dxdt(where(dxdt eq max(dxdt)))
print, xa(where(dxdt eq max(dxdt)))/1e3-60268
oplot, [-1,1]*1e9, [1,1]*1500, linestyle=2
oplot, [0,0], [-1,1]*1e9, color=255*256L
;;  hline, 1500, linestyle=2
;;  vline, 0, color=255*256L
;  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/xa_vs_dxdt.jpeg', /jpeg
;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\xa_vs_dxdt.jpeg', /jpeg
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\xa_vs_dxdt_pat.jpeg', /jpeg
 
;Generate a convenient monotonic array matching dxdt
;refarr = dindgen(n_elements(dxdt))

; !!!!!!-----Manual adjustment required-----!!!!!!
;Specify a baseline portion of the observations that can be used to remove any trends in the freq data set
;if obs eq 'T12N' then aa = where((refarr gt 20)and(refarr lt 400))
;if obs eq 'T12X' then aa = where((refarr gt 550)and(refarr lt 1700))
;if obs eq 'T14N' then aa = where((refarr gt 20)and(refarr lt 410))
;if obs eq 'T14X' then aa = where((refarr gt 200)and(refarr lt 600))
;if obs eq 'T27N' then aa = where((refarr gt 20)and(refarr lt 480))
;if obs eq 'T27X' then aa = where((refarr gt 100)and(refarr lt 600))
;if obs eq 'T31N' then aa = where((refarr gt 400)and(refarr lt 2600))
;if obs eq 'T31X' then aa = where((refarr gt 300)and(refarr lt 2500))
  
  
;SPECIFIC DATA REDUCTION TECHNIQUES
;Cut bad points and replace with linear interpolation between points before and after
if bad_cut ne -1e5 then begin
  stop
  bad = where(abs(dxdt) ge bad_cut)
  for i=0,n_elements(bad)-1 do begin
    dxdt[bad[i]] = interpol([dxdt[bad[i]-1],dxdt[bad[i]+1]], [occtimesec[bad[i]-1],occtimesec[bad[i]+1]], occtimesec[bad[i]])
  endfor
endif
  
if (trim_start ne -1e5) and (trim_end ne -1e5) then begin
  ;Trim dxdt more than it has already been trimmed
  keeptrackofutcrx = keeptrackofutcrx[trim_start:trim_end]
  dxdt = dxdt[trim_start:trim_end]
  xa = xa[trim_start:trim_end]
  occtimesecsub = occtimesecsub[trim_start:trim_end]
endif
  
;Specify a "topside" portion of the ionosphere that is above the plasma in the ionosphere
;topind = where((xa/1e3 - 2575.) ge topside)
topind = where((xa/1e3 - 60268.) ge topside)
;; XXX this occptr-60268 needs to be adjusted to oblate planet
; Inspect the topind range and make sure it looks okay.
;w = window(dimensions=[1800,1000])
;p = plot(occtimesecsub[topind],dxdt[topind],/current)
  
;Subtract the median dxdt value
if median_bool eq 1 then begin
  print, 'topside median=', median(dxdt[topind])
  dxdt = dxdt - median(dxdt[topind])
endif
  
;Fit a linear trend to a portion of the dxdt data
if linear_bool eq 1 then begin
  lintrend = poly_fit(occtimesecsub[topind],dxdt[topind],1)
  dxdt = dxdt - (lintrend(0) + lintrend(1)*occtimesecsub)
  print, 'b,m', lintrend[0], lintrend[1]
endif
  
;If any corrections were made, replot the dxdt curves
if (bad_cut ne -1e5) or (trim_start ne -1e5) or (trim_end ne -1e5) or (median_bool eq 1) or (linear_bool eq 1) then begin
  ;Plot the dxdt vs. time and identify the ionosphere
  window, 5
  plot, occtimesecsub, dxdt,  xtitle='occtimesec', ytitle='dxdt CORRECTED', yrange=[-1,1]*1e16, title='CORRECTED'
  oplot, occtimesecsub, smooth(dxdt,21), color=255
  vline, ionpeakocctimesec,linestyle=2
  hline, 0, color=255*256L
  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/dxdt_vs_occtimesecsub_CORRECTED.jpeg', /jpeg
  ;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\dxdt_vs_occtimesecsub_CORRECTED.jpeg', /jpeg
  ;;    if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\dxdt_vs_occtimesecsub_CORRECTED_pat.jpeg', /jpeg  
  ;Plot the xa vs. dxdt and identify the ionosphere
  window, 6
  ;plot, dxdt, xa/1e3 - 2575., xrange=[-1,1]*5e14, xtitle='dxdt CORRECTED', ytitle='xa-2575', title='CORRECTED'
  plot, dxdt, xa/1e3 - 60268., yrange=[0,1.5e4],xrange=[-1,1]*2e16, xtitle='dxdt CORRECTED', ytitle='xa-60268', title='CORRECTED'
    
  ;oplot, smooth(dxdt,21), xa/1e3 - 2575., color=255
  oplot, smooth(dxdt,21), xa/1e3 - 60268., color=255
  
  hline, 1500, linestyle=2
  vline, 0, color=255*256L
  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/xa_vs_dxdt_CORRECTED.jpeg', /jpeg
  ;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\xa_vs_dxdt_CORRECTED.jpeg', /jpeg
  ;;    if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\xa_vs_dxdt_CORRECTED_pat.jpeg', /jpeg 
endif
  
;NEW IN THIS VERSION
;  if body eq 'Saturn' then begin
;    ;Create array of Saturn radii using the 1-bar gravitational potential 
;    ; representations of Anderson-Schubert and Koskinen.
;    R_AS = dblarr(n_elements(occptlatarray))
;    R_K = dblarr(n_elements(occptlatarray))
;    for i=0,n_elements(R_AS)-1 do begin
;      R_AS[i]=sat_anderson_schubert_1bar(occptlatarray[i])
;      R_K[i] = saturn_1bar(90*!pi/180.-!pi/180*occptlatarray[i])
;    endfor
;  endif 
;stop
  
;Fit the baseline with a linear trend, then subtract it from the dxdt values
;result = poly_fit(refarr(aa), dxdt(aa), 1)
;newdxdt = dxdt  - (result(0) + result(1)*refarr) ;+ result(2)*refarr^2 + result(3)*refarr^3) + result(4)*refarr^4)
newdxdt = dxdt
  
;Show a plot of dxdt before and after the fit
;window, 3
;plot, refarr, abs(dxdt), /ylog
;plot, refarr, dxdt, yrange=[-1,1]*5e15
;vline, refarr[aa[0]], color=255
;vline, refarr[aa[-1]], color=255
;oplot, refarr, (result(0) + result(1)*refarr), linestyle=2
;oplot, refarr[aa], (result(0) + result(1)*refarr[aa])
;oplot, refarr, dxdt  - (result(0) + result(1)*refarr), color=255
;stop  ;stop for inspection of aa
  
;Calculate big X (i.e., TEC, integral(N dl))
newbigx = newdxdt*0.
i=1L
while i lt n_elements(xa) do begin
  newbigx(i) = newbigx(i-1) + 0.5*(newdxdt(i-1)+newdxdt(i))*dt
  i++
endwhile
  
if strmatch(obs, '*X') then begin
  newbigx = newbigx - newbigx[n_elements(newbigx)-1]
endif


;Find dxdr using the derivative including the points before and ahead of the current one.
; The endpoints must then be calculated separately.
; dx/dr here is domega/dx in the proposal.
newdxdr = (shift(newbigx,-1) - shift(newbigx,1)) / (shift(xa,-1)-shift(xa,1))
newdxdr(0) = 2.*newdxdr(2)-newdxdr(1)
newdxdr(n_elements(newdxdr)-1) = 2.*newdxdr(n_elements(newdxdr)-3) - newdxdr(n_elements(newdxdr)-2)
   
;Need to add a check to make sure the occ point radius is monotonic. The data trimming usually takes
; care of this, but the code should stop if it did not. If monotonic, the min of xa will either be the
; first or last element.
if where(xa eq min(xa)) ne 0 && where(xa eq min(xa)) ne n_elements(xa)-1 then begin
  print, 'Occultation radius array is not monotonic'
  stop
endif

;Sort various arrays to prepare for integration
pp = (sort(xa)) ;Returns sorted indices
reva2 = xa(pp) ; Ordered impact parameter in m
thisrkm = reva2 / 1e3 ; Ordered impact parameter in km
newrevdxdr = newdxdr(pp) ; Ordered spatial derivative of bigx in m^-3

;Set first element of newrevdxdr to zero to avoid problems upon numerical integration
newrevdxdr(n_elements(newrevdxdr)-1) = 0.

;Electrondensity is the local electron density, m-3, from corrected data
electrondensity = newrevdxdr*0.
i=0
while i lt n_elements(electrondensity) do begin
  stuff=0.d
  j=i
  while j lt n_elements(electrondensity)-1 do begin
    ;Results are very sensitive to details of integration routine - this is equation 6 of the proposal
    ;calculating integral by "discretely" summing
    stuff = stuff + 0.5 * ( alog(reva2(j)/reva2(i) + sqrt((reva2(j)/reva2(i))^(2.d)-(1.d)) ) + $
      alog(reva2(j+1)/reva2(i) + sqrt((reva2(j+1)/reva2(i))^(2.d)-(1.d)) )) * (newrevdxdr(j+1)-newrevdxdr(j))
    j++
  endwhile
  electrondensity(i) = stuff / !pi
  i++
endwhile


;; XXX GREAT BIG NEW BIT

help, keeptrackofutcrx[0]
help, keeptrackofutcrx[n_elements(keeptrackofutcrx)-1]
; use these times to select subset of data in simpleocc

window,8
plot, electrondensity, reva2, /ynoz

if itloop eq 'a' then stop

if itloop eq 'b' then begin
  restore, obs + 'pathlengthr.sav' 
  restore, obs + 'pathlengthsph.sav'
  restore, obs + 'pathlengthell.sav'
  restore, obs + 'pathlengthash.sav'
  restore, obs + 'pathlengthkos.sav'
  
  totpathlength2dm = (inpathlength2d + outpathlength2d) * 1e3 ; km to m
  altsphtotpathlength2dm = (altsphinpathlength2d + altsphoutpathlength2d) * 1e3 ; km to m
  altelltotpathlength2dm = (altellinpathlength2d + altelloutpathlength2d) * 1e3 ; km to m
  altashtotpathlength2dm = (altashinpathlength2d + altashoutpathlength2d) * 1e3 ; km to m
  altkostotpathlength2dm = (altkosinpathlength2d + altkosoutpathlength2d) * 1e3 ; km to m

  if strmatch(obs, '*N') then begin
    matrixnelec = reverse(reform(invert(transpose(totpathlength2dm)) ## newbigx))
    altsphmatrixnelec = reverse(reform(invert(transpose(altsphtotpathlength2dm)) ## newbigx))
    altellmatrixnelec = reverse(reform(invert(transpose(altelltotpathlength2dm)) ## newbigx))
    altashmatrixnelec = reverse(reform(invert(transpose(altashtotpathlength2dm)) ## newbigx))
    altkosmatrixnelec = reverse(reform(invert(transpose(altkostotpathlength2dm)) ## newbigx))
    ; best not to think too hard about the steps compressed into this line
  endif

  if strmatch(obs, '*X') then begin
    matrixnelec = reverse(reform(invert(transpose(totpathlength2dm)) ## reverse(newbigx)))
    altsphmatrixnelec = reverse(reform(invert(transpose(altsphtotpathlength2dm)) ## reverse(newbigx)))
    altellmatrixnelec = reverse(reform(invert(transpose(altelltotpathlength2dm)) ## reverse(newbigx)))
    altashmatrixnelec = reverse(reform(invert(transpose(altashtotpathlength2dm)) ## reverse(newbigx)))
    altkosmatrixnelec = reverse(reform(invert(transpose(altkostotpathlength2dm)) ## reverse(newbigx)))
    ; best not to think too hard about the steps compressed into this line
  endif

  ;;; XXX ABEL TRANSFORM USES DXDR, MATRIX METHOD USES X - SO BASELINE CORRECTION AND CONSTANT OF INTEGRATION
  ;;; NEEDS TO ENSURE X=0 AT HIGH ALTITUDE, WHEREAS PAUL'S PRIOR WORK HAD X=0 AT EARLY TIME (GOOD FOR INGRESS, BAD FOR EGRESS)

  window,0
  plot, electrondensity, reva2, /ynoz, color='FFFFFF'x ;White
  oplot, matrixnelec, reva2, color='00FF00'x ;Green
  oplot, altsphmatrixnelec, reva2, color='00FFFF'x ;Yellow
  oplot, altellmatrixnelec, reva2, color='FF00FF'x ;Purple
  ;oplot, altashmatrixnelec, reva2, color='FFFF00'x
  ;oplot, altkosmatrixnelec, reva2, color='0000FF'x

  ;window,2
  ;plot, electrondensity, (reva2-60268d3)/1e3, yra=[0,3e3], color='FFFFFF'x
  ;oplot, matrixnelec, (reva2-60268d3)/1e3, color='00FF00'x
  ;oplot, altsphmatrixnelec, (reva2-60268d3)/1e3, color='00FFFF'x
  ;oplot, altellmatrixnelec, (reva2-60268d3)/1e3, color='FF00FF'x
  ;oplot, altashmatrixnelec, (reva2-60268d3)/1e3, color='FFFF00'x
  ;oplot, altkosmatrixnelec, (reva2-60268d3)/1e3, color='0000FF'x

  ;;; XXX WARNING - THIS RESTORES A LOT OF EXISTING VARIABLES, SO OVERWRITES THEM
  ;;; JUST FIX IT IN SIMPLEOCC
  if obs eq 'S07N' then begin
    ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/' + $
    ; 'for_paul_10_06_2019/Output/cors_0107/s10sroi2005123_0230nnnx14rd.2a2/simpleocc_output_itloop2.sav'
    ;klioreprofile = read_csv('s07n_kliore.csv', n_table_header=1)
    restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'
    klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_10_06_2019\Kliore EDP\s07n_kliore.csv', n_table_header=1)
  endif
  
  if obs eq 'S10N' then begin
    ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/' + $
    ; 'for_paul_10_06_2019/Output/cors_0116/s12sroi2005177_1741nnnx14rd.2a2/simpleocc_output_itloop2.sav'
    ;klioreprofile = read_csv('s10n_kliore.csv', n_table_header=1)
    restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'
    klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_10_06_2019\Kliore EDP\s10n_kliore.csv', n_table_header=1)
  endif
  
  if obs eq 'S56X' then begin
    ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/' + $
    ; 'for_paul_10_06_2019/Output/cors_0213/s36saoc2008015_2150nnnx63rd.1a2/simpleocc_output_itloop2.sav'
    ;klioreprofile = read_csv('s56x_kliore.csv', n_table_header=1)
    restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'
    klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_10_06_2019\Kliore EDP\s56x_kliore.csv', n_table_header=1)
  endif
  
  if obs eq 'S51X' then begin
    ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0204/' + $
    ; 's34saoe2007297_0745nnnx63rd.1a2/simpleocc_output_itloop2.sav'
    ;klioreprofile = read_csv('/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/' + $
    ; 'Kliore EDP/s51x_kliore.csv', n_table_header=1)
    klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_25_06_2019\Kliore EDP\s51x_kliore.csv', n_table_header=1)
    restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'

  endif
  if obs eq 'S54X' then begin
  ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0207/' + $
  ; 's36saoe2007353_0524nnnx63rd.1a2/simpleocc_output_itloop2.sav'
  ;klioreprofile = read_csv('/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/' + $
  ; 'Kliore EDP/s54x_kliore.csv', n_table_header=1)
  restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'
  klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_25_06_2019\Kliore EDP\s54x_kliore.csv', n_table_header=1)

  endif
  if obs eq 'S75X' then begin
    ;restore, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0245/' + $
    ; 's42saoe2008189_1040nnnx63rd.1a2/simpleocc_output_itloop2.sav'
    restore, rootpath+body+'\cors_'+xcors+'\'+xproductid+'\simpleocc_output_itloop2.sav'
    ;klioreprofile = read_csv('/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/' + $
    ; 'Kliore EDP/s75x_kliore.csv', n_table_header=1)
    klioreprofile = read_csv('C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_25_06_2019\Kliore EDP\s75x_kliore.csv', n_table_header=1)
  endif

  window,1
  plot, altsphmatrixnelec, reverse(itloop2extraoccptaltspharray), color='FFFFFF'x, yra=[0,3e3] ; yra=[1.4e3,1.6e3] ;White
  oplot, altellmatrixnelec, reverse(itloop2extraoccptaltellarray), color='FF0000'x ;Blue
  ;oplot, altashmatrixnelec, reverse(itloop2extraoccptaltasharray), color='00FF00'x
  ;oplot, altkosmatrixnelec, reverse(itloop2extraoccptaltkosarray), color='0000FF'x
  oplot, klioreprofile.field2 * 1e6, klioreprofile.field1, color='00FFFF'x ;Yellow

  stop
endif


stop
 
;Specify the topside indices for thisrkm
;topind_thisrkm = where((thisrkm - 2575.) ge topside)
topind_thisrkm = where((thisrkm - 60268.) ge topside)
;; XXX this occptr-60268 needs to be adjusted to oblate planet
  
;NEW IN THIS VERSION
;topind_thisrkm = where((thisrkm - R_AS[rep_point_loc]) ge topside)

;Create an array of the standard deviation of the
electrondensity_stddev = dblarr(n_elements(electrondensity))+stddev(electrondensity[topind_thisrkm]/1e6)

;Plot the final electron density profile found here. In some cases, compare with the
; profile found on the PDS.
window, 4
;  plot, electrondensity/1e6, thisrkm-60268, yra=[0,ceil(max(thisrkm-2575.)/1e3)*1e3], xra=[-3000,3000], ytitle='thisrkm - 2575 [km]', xtitle='electrondensity [cm-3]'
  plot, electrondensity/1e6, thisrkm-60268, yra=[-5000,5000], xra=[-3000,10000], ytitle='thisrkm - 60268 [km]', xtitle='electrondensity [cm-3]'
;NEW IN THIS VERSION 
; plot, electrondensity/1e6, thisrkm-R_AS[rep_point_loc], yra=[0,5000], xra=[-3000,10000], ytitle='thisrkm - Representative Point [km]', xtitle='electrondensity [cm-3]'
oplot, [0,0], [-1,1]*1e10
  
;  xyouts, 5000,3000,strcompress('DXDT smoothing = '+smooth_str), charsize=2


;Load in EDP from PDS for comparison for the T12 data set. The others do not exist on the PDS
;kliorepath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Misc/Kliore2008_TitanEDP/'
kliorepath = 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Misc\Kliore2008_TitanEDP\'
  
csvfile = 'skip'
if obs eq 'T12N' then csvfile = kliorepath+'S19TIIOC2006_078_0000_N_001.EDP'
if obs eq 'T12X' then csvfile = kliorepath+'S19TIIOC2006_078_0000_X_001.EDP'
if obs eq 'T14N' then csvfile = kliorepath+'S20TIIOC2006_140_1203_N_001.EDP'
if obs eq 'T14X' then csvfile = kliorepath+'S20TIIOC2006_140_1203_X_001.EDP'
if obs eq 'T27N' then csvfile = kliorepath+'S28TIIOC2007_085_0000_N_001.EDP'
if obs eq 'T27X' then csvfile = kliorepath+'S28TIIOC2007_085_0000_X_001.EDP'
if obs eq 'T31N' then csvfile = kliorepath+'S30TIIOC2007_148_1737_N_001.EDP'
if obs eq 'T31X' then csvfile = kliorepath+'S30TIIOC2007_148_1737_X_001.EDP'

if csvfile ne 'skip' then begin
  ;Only plot if the selected obs has overlapping Kliore occ data
  junk =  READ_CSV(csvfile)
  ;Read in various column of edp file
  zkm = junk.(0)
  sneleccm3 = junk.(1)
  neleccm3 = junk.(2)
  latdeg = junk.(3)
  szadeg = junk.(4)
  lsthr = junk.(5)
  oplot, neleccm3, zkm, color=255
endif
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/edp.jpeg', /jpeg
;if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\edp.jpeg', /jpeg
;;  if save_figs eq 1 then saveimage, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\edp_pat.jpeg', /jpeg
  
;save, filename=path+'output_ed_profiles/'+filename, thisrkm, electrondensity, sdeltaf, xdeltaf
  
;Save the EDP as a text file
openw, lun, save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/'+flyby+'_'+obs+'_'+bands+'_'+station+'_edp.txt', /get_lun 
;openw, lun, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\'+flyby+'_'+obs+'_'+bands+'_'+station+'_edp.txt', /get_lun
;;  openw, lun, save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\'+flyby+'_'+obs+'_'+bands+'_'+station+'_edp_pat.txt', /get_lun
  
printf, lun, transpose([[thisrkm],[electrondensity/1e6],[electrondensity_stddev]])
close, lun
free_lun, lun
  
;Create a .sav file that includes all information needed to archive this EDP on the PDS.
; Note that electrondensity correspond to thisrkm, which is derived from the sorted xa values.
; For egress observations, the sort has no effect, since xa increase in time. But for ingress
; observations, electrondensity is actually reversed compared to all quantities except thisrkm.
; When saving ingress observations, electrondensity must therefore be reversed.
if strmatch(obs, '*X') then begin
  edpsav = electrondensity/1e6
  electrondensity_stddevsav = electrondensity_stddev
endif
if strmatch(obs, '*N') then begin
  edpsav = reverse(electrondensity/1e6)
  electrondensity_stddevsav = reverse(electrondensity_stddev)
endif

;;; XXX This neglects possibility of trim_start, trim_end adjusting lengths of arrays
;These are all passed straight throu
ettxarraysav = ettxarray[indexstart:indexend]
etoccptarraysav = etoccptarray[indexstart:indexend]
etrxarraysav = etrxarray[indexstart:indexend]
utctxarraysav = utctxarray[indexstart:indexend]
utcoccptarraysav = utcoccptarray[indexstart:indexend]
utcrxarraysav = utcrxarray[indexstart:indexend]
OCCPTRADIUSARRAYsav = OCCPTRADIUSARRAY[indexstart:indexend]
OCCPTLATARRAYsav = OCCPTLATARRAY[indexstart:indexend]
OCCPTLONARRAYsav = OCCPTLONARRAY[indexstart:indexend]
OCCPTSZAARRAYsav = OCCPTSZAARRAY[indexstart:indexend]
OCCPTLSTARRAYsav = OCCPTLSTARRAY[indexstart:indexend]
SEPANGLEARRAYsav = SEPANGLEARRAY[indexstart:indexend]
EPSANGLEARRAYsav = EPSANGLEARRAY[indexstart:indexend]
;NEW IN THIS VERSION
R_AS_ARRAYsav = R_AS[indexstart:indexend]
R_K_ARRAYsav = R_K[indexstart:indexend]
sav_filename = save_rootpath+body+'/'+flyby+'/'+obs+'_'+bands+'_'+station+'/'+flyby+'_'+obs+'_'+bands+'_'+station+'_edp.sav'
;;  sav_filename = save_rootpath+body+'\'+flyby+'\'+obs+'_'+bands+'_'+station+'\'+flyby+'_'+obs+'_'+bands+'_'+station+'_edp_pat.sav'
   
save, filename=sav_filename, ettxarraysav, etoccptarraysav, etrxarraysav, utctxarraysav, utcoccptarraysav, utcrxarraysav, $
  OCCPTRADIUSARRAYsav, OCCPTLATARRAYsav, OCCPTLONARRAYsav, OCCPTSZAARRAYsav, OCCPTLSTARRAYsav, SEPANGLEARRAYsav, EPSANGLEARRAYsav, $
  dxdt_save, dxdt, newbigx, edpsav, electrondensity_stddevsav, R_AS_ARRAYsav, R_K_ARRAYsav
  
;while !d.window ne -1 do wdelete, !d.window

stop
end
