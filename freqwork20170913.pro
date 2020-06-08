pro freqwork20170913

; Paul Withers, Paul Dalba, 2017.09.13
; Astronomy Department, Boston University

;Find preliminary electron density profiles using outputs from rsr_proc.pro (run on the SCC)
; and simpleocc20170913.pro. Can be run for any occultation provided the correct path info
; below.

;Set rootpath for the occultation data
rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Output/'

;Set the bands to be compared
bands = 'SX'   ;SX, XK 

;Copy the occultation profile below ---------------------------------------------------

;Specify the body. Options are 'Titan' or 'Saturn' at present
body = 'Titan'

;Specify the cors file. List as a string with a leading zero
xcors = '0145'
scors = '0145'
kcors = '0145'

;Specifiy the X and S product IDs. Make sure they are the same receiver.
xproductid = 's20tioc2006140_1317nnnx63rd.1a2'
;xproductid = 's20tioc2006140_1317nnnx55rv.1n2'

sproductid = 's20tioc2006140_1317nnns63rd.2a2'

kproductid = 's20tioc2006140_1317nnnk55rv.1n2'

;Specify observation (w/ entry (N) or exit (S))
obs = 'T14X'


;End occultation profile --------------------------------------------------------------


;Restore the X-band .sav file from rsr_proc.pro for this data set
restore, filename=rootpath+body+'/cors_'+xcors+'/'+xproductid+'/full_output/output.sav'
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


;Restore the S-band .sav file from rsr_proc.pro for this data set
restore, filename=rootpath+body+'/cors_'+scors+'/'+sproductid+'/full_output/output.sav'
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


;Restore the K-band .sav file from rsr_proc.pro for this data set
restore, filename=rootpath+body+'/cors_'+kcors+'/'+kproductid+'/full_output/output.sav'
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




;Restore the SPICE data for this occultation. Choose any band, if they are the
; same receiver, it should not matter.
if strmatch(bands, 'S*') then restore, rootpath+body+'/cors_'+scors+'/'+sproductid+'/full_output/simpleocc_output.sav'
if strmatch(bands, 'X*') then restore, rootpath+body+'/cors_'+xcors+'/'+xproductid+'/full_output/simpleocc_output.sav'

;Determine UTC seconds for ERT from SPICE data
;real_ertutcsec = double(strmid(real_ertutc,12,2))*60d*60d + $
;double(strmid(real_ertutc,15,2))*60d + $
;double(strmid(real_ertutc,18,20))



;Use the times from SPICE and from the freq vs. time output to determine point radius for each freq value.
; The freq times and simpleocc times will match 
if strmatch(bands, 'S*') then occtimesec = sfdusecoutsr 
if strmatch(bands, 'X*') then occtimesec = sfdusecoutxr 
;occptr = interpol(occptradiusarray,real_ertutcsec,occtimesec)
occptr = occptradiusarray
;Time relative to start of observations
occtimesecadj = occtimesec - occtimesec[0] 

;Define the occultaiton frequency values
occfs = freqsr
occfx = freqxr
occfk = freqkr

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
;stop

;Plot the time vs. the pt radius
window, 1
plot, occtimesec, occptr-2575., xtitle='occtimesec', ytitle='occptradius - 2575'
hline, 0  ;the surface
hline, 1500, linestyle=2  ;approx peak of ionospheric ED


;Develop an algorithm to identify the bad portions of the data
;First, where occtimesec is below zero is not to be trusted. 
if strmatch(obs, '*N') then begin
  ;This for an entry observation
  indexend = (where(occptr-2575. le 300.))[0] - 1
  
  ;Now use a weighted average of the max power to estimate the start of the data
  peakpowerthresh = 10.^(0.9*max(alog10(peakpower))+0.1*min(alog10(peakpower)))
  indexstart = (where(peakpower ge peakpowerthresh))[0]
  
  ;This determines the occtimesec of the ionospheric peak
  ionpeakocctimesec = occtimesec[(where(occptr-2575 le 1500.))[0]-1] 
  
  ;Manual offsets
  indexstart = indexstart+50
  ;indexend = indexend-300
endif



if strmatch(obs, '*X') then begin
  ;This for an exit observation
  indexstart = (where(occptr-2575. le 300.))[-1] + 1
  
  ;Now use a weighted average of the max power to estimate the end of the data
  peakpowerthresh = 10.^(0.9*max(alog10(peakpower))+0.1*min(alog10(peakpower)))
  indexend = (where(peakpower ge peakpowerthresh))[-1]
  
  ;This determines the occtimesec of the ionospheric peak
  ionpeakocctimesec = occtimesec[(where(occptr-2575 le 1500.))[-1]+1]

  ;Manual offsets
  ;indexstart = indexstart+100
  indexend = indexend-50
endif




;Sometime, special adjustment of the indexstart and indexend is required
indexstart = 1719
indexend = 2238




wset, 0
vline, occtimesec[indexstart], color=255
vline, occtimesec[indexend], color=255

wset, 1
vline, occtimesec[indexstart], color=255
vline, occtimesec[indexend], color=255

;stop



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
sdeltaf = occfs(indexstart:indexend)
xdeltaf = occfx(indexstart:indexend)
kdeltaf = occfk(indexstart:indexend)
sfobs = occfs(indexstart:indexend) ;Transmitted frequency (S)
xfobs = occfx(indexstart:indexend) ;Transmitted frequency (X)
xa = occptr(indexstart:indexend) * 1D3  ;Convert point radius to meters
occtimesecsub = occtimesec(indexstart:indexend)  ;the subset of occsectime points 

;Solve proposal equation 3 for dx/dt where x=TEC=integral(N dl)
if bands eq 'SX' then dxdt = (sdeltaf - xdeltaf/fdxfds) * 8d * !pi * !pi * me_mks * eps0_mks *(sfobs) * c_mks / (e_mks^2) * 1d/(1d - fdxfds^(-2d))
if bands eq 'SK' then dxdt = (sdeltaf - kdeltaf/fdkfds) * 8d * !pi * !pi * me_mks * eps0_mks *(sfobs) * c_mks / (e_mks^2) * 1d/(1d - fdkfds^(-2d))
if bands eq 'XK' then dxdt = (xdeltaf - kdeltaf/fdkfdx) * 8d * !pi * !pi * me_mks * eps0_mks *(xfobs) * c_mks / (e_mks^2) * 1d/(1d - fdkfdx^(-2d))



;Plot the dxdt vs. time and identify the ionosphere
window, 2
plot, occtimesecsub, dxdt,  xtitle='occtimesec', ytitle='dxdt', yrange=[-1,1]*5e14
oplot, occtimesecsub, smooth(dxdt,21), color=255
vline, ionpeakocctimesec,linestyle=2
hline, 0, color=255*256L

;Plot the xa vs. dxdt and identify the ionosphere
window, 3
plot, dxdt, xa/1e3 - 2575., xtitle='dxdt', ytitle='xa-2575', xrange=[-1,1]*5e14
oplot, smooth(dxdt,21), xa/1e3 - 2575., color=255
hline, 1500, linestyle=2
vline, 0, color=255*256L

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


specialtreatment = 1
if specialtreatment then begin
  ;Cut bad points and replace with linear interpolation between points before and after
  ;badcut = 4e14
  ;bad = where(abs(dxdt) ge badcut)
  ;for i=0,n_elements(bad)-1 do begin
  ;  dxdt[bad[i]] = interpol([dxdt[bad[i]-1],dxdt[bad[i]+1]], [occtimesec[bad[i]-1],occtimesec[bad[i]+1]], occtimesec[bad[i]])
  ;endfor
  
  ;Trim dxdt more than it has already been trimmed
  ;trimind = [56,-1]
  ;dxdt = dxdt[trimind[0]:trimind[1]]
  ;xa = xa[trimind[0]:trimind[1]]
  ;occtimesecsub = occtimesecsub[trimind[0]:trimind[1]]
  
  ;Specify a "topside" portion of the ionosphere that is above the plasma in the ionosphere
  topside = 1800.   ;km -2575. for Titan
  topind = where((xa/1e3 - 2575.) ge topside)
  
  ;Subtract the median dxdt value 
  print, 'topside median=', median(dxdt[topind])
  ;dxdt = dxdt - median(dxdt[topind])
  
  ;Fit a linear trend to a portion of the dxdt data
  lintrend = poly_fit(occtimesecsub[topind],dxdt[topind],1)
  dxdt = dxdt - (lintrend(0) + lintrend(1)*occtimesecsub)
  print, 'b,m', lintrend[0], lintrend[1]
  
  ;Plot the dxdt vs. time and identify the ionosphere
  window, 5
  plot, occtimesecsub, dxdt,  xtitle='occtimesec', ytitle='dxdt CORRECTED', yrange=[-1,1]*5e14, title='CORRECTED'
  oplot, occtimesecsub, smooth(dxdt,21), color=255
  vline, ionpeakocctimesec,linestyle=2
  hline, 0, color=255*256L
  
  ;Plot the xa vs. dxdt and identify the ionosphere
  window, 6
  plot, dxdt, xa/1e3 - 2575., xrange=[-1,1]*5e14, xtitle='dxdt CORRECTED', ytitle='xa-2575', title='CORRECTED'
  oplot, smooth(dxdt,21), xa/1e3 - 2575., color=255
  hline, 1500, linestyle=2
  vline, 0, color=255*256L
endif












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

;Find dxdr using the derivative including the points before and ahead of the current one. 
; The endpoints must then be calculated separately.
; dx/dr here is domega/dx in the proposal.
newdxdr = (shift(newbigx,-1) - shift(newbigx,1)) / (shift(xa,-1)-shift(xa,1))
newdxdr(0) = 2.*newdxdr(2)-newdxdr(1) 
newdxdr(n_elements(newdxdr)-1) = 2.*newdxdr(n_elements(newdxdr)-3) - newdxdr(n_elements(newdxdr)-2)


;Need to add a check to make sure the occ point radius is monotonic. The data trimming usually takes
; care of this, but the code shode stop if it did not. If monotonic, the min of xa will either be the 
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


;Plot the final electron density profile found here. In some cases, compare with the
; profile found on the PDS. 
window, 4
plot, electrondensity/1e6, thisrkm-2575, yra=[0,ceil(max(thisrkm-2575.)/1e3)*1e3], xra=[-3000,3000], ytitle='thisrkm - 2575 [km]', xtitle='electrondensity [cm-3]'
oplot, [0,0], [-1,1]*1e10       


;Load in EDP from PDS for comparison for the T12 data set. The others do not exist on the PDS
kliorepath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Misc/Kliore2008_TitanEDP/'
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


;save, filename=path+'output_ed_profiles/'+filename, thisrkm, electrondensity, sdeltaf, xdeltaf


stop
end