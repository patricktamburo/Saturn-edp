pro freqwork20170608

; Paul Withers, Paul Dalba, 2017.06.08
; Astronomy Department, Boston University

;Find preliminary electron density profiles using outputs from rsr_proc.pro (run on the SCC)
; and simpleocc20170608.pro. Can be run for any occultation provided the correct path info
; below.

;Set rootpath for the occultation data
rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Output/'

;Specify the body. Options are 'Titan' or 'Saturn' at present
body = 'Titan'

;Specify the cors file. List as a string with a leading zero
cors = '0140'

;Specifiy the X and S product IDs. At this point, make sure they are the same receiver.
xproductid = 's19tioc2006078_0107nnnx14rd.1a2'
sproductid = 's19tioc2006078_0107nnns14rd.2a2'



;Restore the X-band .sav file from rsr_proc.pro for this data set
restore, filename=rootpath+body+'/cors_'+cors+'/'+xproductid+'/full_output/output.sav'
sfduyearoutxr = sfduyearout
sfdudoyoutxr = sfdudoyout
sfdusecoutxr = sfdusecout
rftoifmhzoutxr = rftoifmhzout
ddclomhzoutxr = ddclomhzout
ncofreqoutxr = ncofreqout
kposaaoutxr = kposaaout

;Calculate the sky frequency
freqxr = (rftoifmhzoutxr + ddclomhzoutxr) * 1d6 - ncofreqoutxr + kposaaoutxr


;Restore the S-band .sav file from rsr_proc.pro for this data set
restore, filename=rootpath+body+'/cors_'+cors+'/'+sproductid+'/full_output/output.sav'
sfduyearoutsr = sfduyearout
sfdudoyoutsr = sfdudoyout
sfdusecoutsr = sfdusecout
rftoifmhzoutsr = rftoifmhzout
ddclomhzoutsr = ddclomhzout
ncofreqoutsr = ncofreqout
kposaaoutsr = kposaaout

;Calculate the sky frequency
freqsr = (rftoifmhzoutsr + ddclomhzoutsr) * 1d6 - ncofreqoutsr + kposaaoutsr


;Restore the SPICE data for this occultation. Choose either S or X, if they are the
; same receiver and same band, it should not matter.
restore, rootpath+body+'/cors_'+cors+'/'+sproductid+'/full_output/simpleocc_output.sav'


;Determine UTC seconds for ERT from SPICE data
real_ertutcsec = double(strmid(real_ertutc,12,2))*60d*60d + $
double(strmid(real_ertutc,15,2))*60d + $
double(strmid(real_ertutc,18,20))



;Use the times from SPICE and from the freq vs. time output to determine point radius for each freq value.
; XR and SR times from rsr_proc.pro are identical, so can use either one
occtimesec = sfdusecoutsr 
;occptr = interpol(occptradiusarray,real_ertutcsec,occtimesec)
occptr = occptradiusarray
;Time relative to start of observations
occtimesecadj = occtimesec - occtimesec[0] 

;Define the occultaiton frequency values
occfs = freqsr
occfx = freqxr


;Constants
c_mks     = 2.9979e8		  
me_mks    = 9.1094e-31		  
mp_mks    = 1.6726e-27		
eps0_mks  = 8.8542e-12            
e_mks     = 1.6022e-19		  
fdxfds = 11d/3d ; Ratio of downlink X-band frequency to downlink S-band frequency

;The time resolution is 1 second.
dt = 1D

; !!!!!!-----Manual adjustment required-----!!!!!!
;Trim portions of the freq array that contain bad data (i.e. neutral atm, planet occulted)
bb = where((occtimesec gt 4200)and(occtimesec lt 4670)) ;ingress
;bb = where((occtimesec gt 5325)and(occtimesec lt 7300))  ;egress



;Plot the frequency difference versus point radius trimmed and untrimmed
window, 0
plot, occfs - 3d/11d * occfx, occptr, xra=[-1,1]*1e1 ;, yra=[2e3,5e3]
;oplot, occfs[bb] - 3d/11d * occfx[bb], occptr[bb]
;plot, abs(occfs - 3d/11d * occfx), occptr, /xlog ;, yra=[2e3,5e3]
;oplot, abs(occfs[bb] - 3d/11d * occfx[bb]), occptr[bb], color=255

p1 = plot(occfs - 3d/11d * occfx, sfdusecout, xra=[-1,1]*1e1)
;p1.close

;Instead, make a so-called 'salt and pepper' plot (Withers 2014, Fig. 6)
diff = occfs - 3d/11d * occfx
;window, 5
;plot, abs(diff[where(diff gt 0d)]), occptr[where(diff gt 0d)], psym=1, /xlog ;, yra=[2e3,5e3]
;oplot, abs(diff[where(diff lt 0d)]), occptr[where(diff lt 0d)], psym=1, color=255
;oplot, abs(diff[bb]), occptr[bb], psym=2, color=150


;Trim and re-define
sdeltaf = occfs(bb)
xdeltaf = occfx(bb)
sfobs = occfs(bb) ;Transmitted frequency
xa = occptr(bb) * 1D3  ;Convert point radius to meters


;Solve proposal equation 3 for dx/dt where x=TEC=integral(N dl)
dxdt = (sdeltaf - xdeltaf/fdxfds) * 8d * !pi * !pi * me_mks * eps0_mks *(sfobs) * c_mks / (e_mks^2) * 1d/(1d - fdxfds^(-2d))
;dxdt = smooth(dxdt,11)


;Generate a convenient monotonic array matching dxdt
refarr = dindgen(n_elements(dxdt))

; !!!!!!-----Manual adjustment required-----!!!!!!
;Specify a baseline portion of the observations that can be used to remove any trends in the freq data set
aa = where(refarr lt 350)   ;ingress
;aa = where(refarr gt 400)   ;egress


;Fit the baseline with a linear trend, then subtract it from the dxdt values 
result = poly_fit(refarr(aa), dxdt(aa), 1)
newdxdt = dxdt  - (result(0) + result(1)*refarr) ;+ result(2)*refarr^2 + result(3)*refarr^3) + result(4)*refarr^4)

;Show a plot of dxdt before and after the fit
window, 1
plot, refarr, dxdt
oplot, refarr, (result(0) + result(1)*refarr), linestyle=2
oplot, refarr[aa], (result(0) + result(1)*refarr[aa])
oplot, refarr, dxdt  - (result(0) + result(1)*refarr), color=255

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
window, 2
plot, electrondensity/1e6, thisrkm-2575, yra=[0,5000], xra=[-500,3000]
oplot, [0,0], [-1,1]*1e10       


;Load in EDP from PDS for comparison for the T12 data set. The others do not exist on the PDS
kliorepath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/Misc/Kliore2008_TitanEDP/'
csvfile = kliorepath+'S19TIIOC2006_078_0000_N_001.EDP' ;ingress
;csvfile = kliorepath+'S19TIIOC2006_078_0000_X_001.EDP' ;egress
junk =  READ_CSV(csvfile)
;Read in various column of edp file
zkm = junk.(0)
sneleccm3 = junk.(1)
neleccm3 = junk.(2)
latdeg = junk.(3)
szadeg = junk.(4)
lsthr = junk.(5)
oplot, neleccm3, zkm, color=255



;save, filename=path+'output_ed_profiles/'+filename, thisrkm, electrondensity, sdeltaf, xdeltaf


stop
end