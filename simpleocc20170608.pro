pro simpleocc20170608

; Paul Withers, 2006.06.30
; Center for Space Physics, Boston University

; Find spacecraft occultations using SPICE

;spawn, 'date', codestarttime

; Record time program starting running
;rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/'
rootpath = 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\'
;spicerootpath = '/Volumes/PW-2TB/Cassini_Ionosphere/'
spicerootpath = 'C:\Users\tambu\Documents\BU\SaturnEDP\SPICE_kernels\'
;genspicepath = spicerootpath+'generic/'
genspicepath = spicerootpath+'generic\'
;spkspicepath = spicerootpath+'spk/'
spkspicepath = spicerootpath+'spk\'

; Set root path and SPICE paths (requires PW-2TB drive)


; Clean any lingering kernels out of memory here,
; not at end of code
cspice_ktotal, 'all', count
help, count
i=0
while i lt count do begin
cspice_kdata, 0, 'all', file, type, source, handle, found
cspice_unload, file
i=i+1
endwhile
cspice_ktotal, 'all', count
help, count

cspice_furnsh, genspicepath+'pck00010.tpc' ; planet rotational states
cspice_furnsh, genspicepath+'naif0012.tls' ; leap seconds
cspice_furnsh, genspicepath+'de430.bsp' ; SPK (trajectory kernel) for planets


;Specify the outpath path of the RSR file. Either pick a specific path... 
;loop = 'no'
;path = rootpath+'Output/Titan/cors_0140/s19tioc2006078_0107nnnx14rd.1a2/'    ;t12 xr
;path = rootpath+'Output/Titan/cors_0140/s19tioc2006078_0107nnns14rd.2a2/'    ;t12 sr
;path = rootpath+'Output/Titan/cors_0278/S51TIOI2009173_1845NNNX14RD.1A2/'    ;t57 xr
;path = rootpath+'Output/Titan/cors_0278/S51TIOI2009173_1845NNNS14RD.1B2/'    ;t57 sr
;path = rootpath+'Output/Saturn/cors_0107/s10sroi2005123_0230nnnx14rd.2a2/'   ;s07 xr
;path = rootpath+'Output/Saturn/cors_0106/s10sroi2005123_0230nnns14rd.2b2/'   ;s07 sr
;
;...or choose a cors number and run every RSR file within that cors directory.
loop = 'yes'
corsnumber = '107'
body = 'Saturn'
allrsrpaths = file_search(rootpath+'Output\'+body+'\cors_0'+corsnumber+'\*')
allrsrpaths = allrsrpaths[1]
;Loop over all rsr paths
for rsrpath = 0, n_elements(allrsrpaths)-1 do begin
path = allrsrpaths[rsrpath]+'\'



;Determine which body from the path. 's' for Saturn occultations, 't' for Titan occultations
;if strmatch(path,'*Titan*') then whichbody = 't'
;if strmatch(path,'*Saturn*') then whichbody = 's'
if body eq 'Saturn' then whichbody = 's'
if body eq 'Titan' then whichbody = 't'

;Extract the RSR filename from the path
splitpath = strsplit(path,'\',/extract)
rsrfilename = splitpath[-1]


;Check to see if this occultation already has a simpleocc_output.sav file. If so, skip it.
;if (size(file_search(path+'full_output/simpleocc_output*')))[0] eq 1 then begin
;  print, rsrfilename+' already has a simpleocc output file. Skipping.'
;  if loop eq 'no' then stop
;  if loop eq 'yes' then stop  ;replace stop with continue here to loop
;endif


;Determine the appropropriate .bsp spice files to load for this RSR file.
;Once found, that file or those files into spice
openr, lun, rootpath+'Programs/RSR_SPICE_connection.txt', /get_lun
line = ''
while not eof(lun) do begin
  readf, lun, line
  ;Check if this is the line for this RSR file
  if strmatch(line,rsrfilename+'*') then begin
    ;stop
    ;Extract the .bsp filename and furnish spice with them
    splitline = strsplit(line,',',/extract)
    spicespkfiles = splitline[1:*]
    for spk=0,n_elements(spicespkfiles)-1 do begin
      ;Ingest spacecraft trajectory kernel into memory
      cspice_furnsh, spkspicepath+spicespkfiles[spk]
    endfor
    break
  endif
endwhile
free_lun, lun




frame='j2000'
; Make SPICE work in J2000 frame as default
method='intercept'
; Used for ray-tracing, I think
; See help pages about commands that use this variable
abcorr='none'
; Don't do any aberration corrections
; See help pages about commands that use this variable

if whichbody eq 't' then rtarget = double(2575.) + 1200. ; km, Titan
if whichbody eq 's' then rtarget = double(60268.) + 2000. ; km, Saturn
; Radius of target body, km
; Values increased from reference equatorial radius since Cassini
; publications report occultation conditions at ionospheric peak

dt = double(1.)

sctstr = '-82'
sctint = -82L


;Load the output of the rsr_proc code to get the SFDU start and end times of the RSR
restore, path+'full_output/output.sav'
fakesfduyearoutstart = sfduyearout[0]
fakesfdudoyoutstart = sfdudoyout[0]
fakesfdusecoutstart = sfdusecout[0]
fakesfduyearoutend = sfduyearout[-1]
fakesfdudoyoutend = sfdudoyout[-1]
fakesfdusecoutend = sfdusecout[-1]

fakeyearstartstr = strcompress(fakesfduyearoutstart, /remove_all)
if strlen(fakeyearstartstr) ne 4 then stop

fakedoystartstr = strcompress(fakesfdudoyoutstart, /remove_all)
while strlen(fakedoystartstr) lt 3 do fakedoystartstr = '0' + fakedoystartstr
if strlen(fakedoystartstr) ne 3 then stop

fakehhstartlon = floor(fakesfdusecoutstart/ 60d^2)
fakehhstartstr = strcompress(fakehhstartlon, /remove_all)
while strlen(fakehhstartstr) lt 2 do fakehhstartstr = '0' + fakehhstartstr
if strlen(fakehhstartstr) ne 2 then stop

fakemmstartlon = floor( (fakesfdusecoutstart - fakehhstartlon*60d^2) / 60d)
fakemmstartstr = strcompress(fakemmstartlon, /remove_all)
while strlen(fakemmstartstr) lt 2 do fakemmstartstr = '0' + fakemmstartstr
if strlen(fakemmstartstr) ne 2 then stop

fakessstartdbl = fakesfdusecoutstart - fakehhstartlon*60d^2 - fakemmstartlon*60d
fakessstartstr = strcompress(fakessstartdbl, /remove_all)
if strpos(fakessstartstr, '.') eq 1 then fakessstartstr = '0' + fakessstartstr
if strpos(fakessstartstr, '.') ne 2 then stop
fakessstartstr = strmid(fakessstartstr, 0, 6)

faketimestartstr = fakeyearstartstr + '-' + fakedoystartstr + '::' + fakehhstartstr + ':' + fakemmstartstr + ':' + fakessstartstr + ' UTC'
cspice_str2et, faketimestartstr, fakeetstart


fakeyearendstr = strcompress(fakesfduyearoutend, /remove_all)
if strlen(fakeyearendstr) ne 4 then stop

fakedoyendstr = strcompress(fakesfdudoyoutend, /remove_all)
while strlen(fakedoyendstr) lt 3 do fakedoyendstr = '0' + fakedoyendstr
if strlen(fakedoyendstr) ne 3 then stop

fakehhendlon = floor(fakesfdusecoutend/ 60d^2)
fakehhendstr = strcompress(fakehhendlon, /remove_all)
while strlen(fakehhendstr) lt 2 do fakehhendstr = '0' + fakehhendstr
if strlen(fakehhendstr) ne 2 then stop

fakemmendlon = floor( (fakesfdusecoutend - fakehhendlon*60d^2) / 60d)
fakemmendstr = strcompress(fakemmendlon, /remove_all)
while strlen(fakemmendstr) lt 2 do fakemmendstr = '0' + fakemmendstr
if strlen(fakemmendstr) ne 2 then stop

fakessenddbl = fakesfdusecoutend - fakehhendlon*60d^2 - fakemmendlon*60d
fakessendstr = strcompress(fakessenddbl, /remove_all)
if strpos(fakessendstr, '.') eq 1 then fakessendstr = '0' + fakessendstr
if strpos(fakessendstr, '.') ne 2 then stop
fakessendstr = strmid(fakessendstr, 0, 6)

faketimeendstr = fakeyearendstr + '-' + fakedoyendstr + '::' + fakehhendstr + ':' + fakemmendstr + ':' + fakessendstr + ' UTC'
cspice_str2et, faketimeendstr, fakeetend

;fakeetstart = fakeetstart - 100d * 60d ; allow 100 minutes for possible light-time issues and UTC/ET differences
;fakeetend = fakeetend + 100d * 60d ; allow 100 minutes for possible light-time issues and UTC/ET differences

etstart = fakeetstart
etend = fakeetend
nsteps = ceil( (etend- etstart)/dt)
; Number of steps

if whichbody eq 't' then begin
targetint = 606L
targetstr = '606'
; SPICE label for desired target body
; 606 for Titan
endif

if whichbody eq 's' then begin
targetint = 699L
targetstr = '699'
; SPICE label for desired target body
; 699 for Saturn
endif

earthint = 399L
earthstr = '399'
; SPICE label for Earth
; Change to 10 (Sun) to get results for solar occultations

ssbint = 0L
ssbstr = '0'
; Label for solar-system barycenter, handy origin to use for state vectors

saturnint = 699L
saturnstr = '699'
; Label for Saturn, since we want to find its location for both Saturn occultations and Titan occultations

titanint = 606L
titanstr = '606'
; Label for Titan, since we want to find its location for both Saturn occultations and Titan occultations

; Variables for ingress occultations
i_et = [double(0.)] ; ET
i_sctr = i_et ; distance between spacecraft and centre of planet at occ, km
i_sctlon = i_et ; longitude of spacecraft, degE
i_sctlat = i_et ; latitude of spacecraft, degN
i_occr = i_et 
; distance between occultation point and centre of target at occ, 
; should be planetary radius to accuracy permitted by timestep, km
i_occlon = i_et ; longitude of occultation point, degE
i_occlat = i_et ; latitude of occultation point, degN
i_ssr = i_et ; distance between subsolar point and centre of planet at occ, km
i_sslon = i_et ; longitude of subsolar point, degE
i_sslat = i_et ; latitude of subsolar point, degN
i_earthsctdist = i_et ; Earth-spacecraft distance at occ, km
i_earthtargetdist = i_et ; Earth-target distance at occ, km
i_occlst = i_et ; LST of occultation point, hrs
i_occsza = i_et ; SZA of occultation point, deg
i_ls = i_et ; Ls at time of occultation, deg

; Variables for egress occultations
e_et = i_et
e_sctr = i_et
e_sctlon = i_et
e_sctlat = i_et
e_occr = i_et
e_occlon = i_et
e_occlat = i_et
e_ssr = i_et
e_sslon = i_et
e_sslat = i_et
e_earthsctdist = i_et
e_earthtargetdist = i_et
e_occlst = i_et
e_occsza = i_et
e_ls = i_et

p_et = i_et ; ET of closest approach of spacecraft to target
p_sctr = i_et ; radius of spacecraft at closest approach of spacecraft to target, km
p_lat = i_et ; latitude of spacecraft at closest approach of spacecraft to target, km
p_lon = i_et ; longitude of spacecraft at closest approach of spacecraft to target, km
p_lst = i_et ; LST of spacecraft at closest approach of spacecraft to target, km
; p_sza = i_et ; SZA of spacecraft at closest approach of spacecraft to target, km

real_et = i_et ; Time, ET, at each timestep
occptradiusarray = i_et ; distance between centre of target and
; closest approach of ray path to target, km, at each timestep
sctrarray = i_et ; distance between centre of target and spacecraft at each timestep, km
esdarray = sctrarray ; distance between Earth and spacecraft at each timestep, km
ltimetargeteartharray = i_et ; Light time between target and Earth at each timestep, UNITS
ltimescteartharray = i_et ; Light time between sct and Earth at each timestep, UNITS

et = etstart

real_et = dblarr(nsteps)
occptradiusarray = dblarr(nsteps)
sctrarray = dblarr(nsteps)
esdarray = dblarr(nsteps)
ltimetargeteartharray = dblarr(nsteps)
ltimescteartharray = dblarr(nsteps)



i=0L
while i lt nsteps do begin
;while et lt etend do begin

etrx = et ; rx means receiver, tx means transmitter
; use this time for Earth
et = 'crashifthisiseverused'
; ensure this variable doesn't play direct role in rest of code...

cspice_spkpos, earthstr, etrx, frame, 'LT', sctstr, ptarg, txrxltime
ettx = etrx - txrxltime
; use this time for events in the Saturn system, so Cassini, Saturn, Titan

cspice_spkgps, earthint, etrx, frame, ssbint, earthssb, ltimeearthssb
; Position of Earth at time etrx relative to solar system barycentre and in specified frame

stop
cspice_spkgps, targetint, ettx, frame, ssbint, targetssb, ltimetargetssb
stop
; Position of target body at time ettx relative to solar system barycentre and in specified frame
cspice_spkgps, sctint, ettx, frame, ssbint, sctssb, ltimemgsssb
; Position of spacecraft at time ettx relative to solar system barycentre and in specified frame

earthsct = earthssb - sctssb
; Position of Earth relative to spacecraft in working frame

cspice_nplnpt, earthssb, earthsct, targetssb, pnear, distoccpttarget
; pnear is nearest point on Earth-spacecraft line to target body
; Existence of occultations will depend a lot on pnear
; I loosely call position indicated by pnear the "occultation point"
; See help pages for frame in which pnear is stated
; distoccpttarget is distance between pnear and target body, km

occpttarget = pnear - targetssb
; Position of pnear relative to target body in working frame

cspice_tipbod, frame, targetint, ettx, tipm
; tipm is transformation matrix between specified inertial frame (frame, probably J2000)
; and body-fixed frame of target body
; See help pages for more on relevant frames

cspice_mxv, tipm, occpttarget, bfoccpt
; Occultation point position converted from one frame (occpttarget)
; to body-fixed frame of target body (bfoccpt)
cspice_reclat, bfoccpt, occptradius, occptlongitude, occptlatitude
stop
; Convert bfoccpt from cartesian coordinates to r/lat/lon coordinates (km/radE/radN)
occptlongitude = occptlongitude * 180./!pi
occptlatitude = occptlatitude * 180./!pi
; Convert angles from radians into degrees

scttarget = sctssb - targetssb
; Position of spacecraft relative to target body in working frame
cspice_mxv, tipm, scttarget, bfsct
; Spacecraft position transformed from one frame (scttarget)
; to body-fixed frame of target body (bfsct)
cspice_reclat, bfsct, sctradius, sctlongitude, sctlatitude
; Convert bfsct from cartesian coordinates to r/lat/lon coordinates (km/radE/radN)
sctlongitude = sctlongitude * 180./!pi
sctlatitude = sctlatitude * 180./!pi
; Convert angles from radians into degrees

;print, 'occptr/lon/lat = ', occptradius, occptlongitude, occptlatitude

earthsctdist = cspice_vdist(earthssb, sctssb)
; earthsctdist is distance between Earth and spacecraft
earthtargetdist = cspice_vdist(earthssb, targetssb)
; earthtargetdist is distance between Earth and target body
sctminustarget = earthsctdist - earthtargetdist
; Difference between these two distances
; This program examines whether an infinite line defined 
; by two points, Earth and the spacecraft, intersects the surface of the target body
; If it does, then an occultation only occurs if the spacecraft is behind the 
; target body. That is, if sctminustarget is positive.

ls = cspice_lspcn(targetstr, ettx, abcorr)
; ls (radians) is often used as a measure of seasons
; Especially common on Mars, see help pages for details

cspice_et2lst, ettx, targetint, !pi/180.*occptlongitude, 'planetocentric', $
 hr, mn, sc, time, ampm
lst = hr + (mn + sc/double(60.))/double(60.)
; lst is local solar time (0-24 hours) of occultation point at time et

cspice_et2lst, ettx, targetint, !pi/180.*sctlongitude, 'planetocentric', $
 hr, mn, sc, time, ampm
sctlst = hr + (mn + sc/double(60.))/double(60.)
; sctlst is local solar time (0-24 hours) of spacecraft at time et

;sza = cspice_vsep(double(-1.)*targetssb, scttarget)*double(180.)/!pi
; Rejected method of finding spacecraft's solar zenith angle (degrees)
; Angle between spacecraft-target vector and solar system barycentre-target vector
; Sun is not quite solar system barycentre, but it is close enough

cspice_subsol, method, targetstr, ettx, abcorr, earthstr, targetsubsolpt
; targetsubsolpt is position of sub-solar point on target body at time et
; See help pages for frame in which targetsubsolpt is stated
; Note that this may use a triaxial ellipsoid for the target body, which is not desirable for the SZA purpose

sza = cspice_vsep(targetsubsolpt, bfoccpt)*180./!pi
; Solar zenith angle (degrees) at occultation point
; Angle between subsolar point-centre of target body vector and
; occultation point-centre of target body vector

;sctsza = cspice_vsep(targetsubsolpt, subsctpoint)*180./!pi
; Solar zenith angle (degrees) at spacecraft

if i gt 0L then begin

if (occptradius lt rtarget and lastoccptradius ge rtarget and $
 sctminustarget gt double(0.)) then begin
; if occultation point was outside target radius at last timestep and
; occultation point is inside target radius at this timestep and
; spacecraft is behind target body, then
; an ingress occultation occurs

i_et = [i_et, ettx]
i_sctr = [i_sctr, sctradius]
i_sctlon = [i_sctlon, sctlongitude]
i_sctlat = [i_sctlat, sctlatitude]
i_occr = [i_occr, occptradius]
i_occlon = [i_occlon, occptlongitude]
i_occlat = [i_occlat, occptlatitude]
;i_ssr = [i_ssr, targetsubsolradius]
;i_sslon = [i_sslon, targetsubsollongitude]
;i_sslat = [i_sslat, targetsubsollatitude]
i_earthsctdist = [i_earthsctdist, earthsctdist]
i_earthtargetdist = [i_earthtargetdist, earthtargetdist]
i_occlst = [i_occlst, lst]
i_occsza = [i_occsza, sza]
i_ls = [i_ls, ls]
endif

if (occptradius ge rtarget and lastoccptradius lt rtarget and $
 sctminustarget gt double(0.)) then begin
; if occultation point was inside target radius at last timestep and
; occultation point is outside target radius at this timestep and
; spacecraft is behind target body, then
; an egress occultation occurs

e_et = [e_et, ettx]
e_sctr = [e_sctr, sctradius]
e_sctlon = [e_sctlon, sctlongitude]
e_sctlat = [e_sctlat, sctlatitude]
e_occr = [e_occr, occptradius]
e_occlon = [e_occlon, occptlongitude]
e_occlat = [e_occlat, occptlatitude]
;e_ssr = [e_ssr, targetsubsolradius]
;e_sslon = [e_sslon, targetsubsollongitude]
;e_sslat = [e_sslat, targetsubsollatitude]
e_earthsctdist = [e_earthsctdist, earthsctdist]
e_earthtargetdist = [e_earthtargetdist, earthtargetdist]
e_occlst = [e_occlst, lst]
e_occsza = [e_occsza, sza]
e_ls = [e_ls, ls]
endif

endif

if i ge 2 then begin
if (lastbutonesctr gt lastsctr) and (sctradius gt lastsctr) then begin 

; if radial distance between spacecraft and target body at last timestep
; is smaller than at preceding timestep and smaller than at current timestep,
; then last timestep was closest approach of spacecraft to target body
; Effectively a periapsis, if spacecraft is orbiting target body
p_et = [p_et, lastet]
p_sctr = [p_sctr, lastsctr]
p_lat = [p_lat, lastsctlat]
p_lon = [p_lon, lastsctlon]
p_lst = [p_lst, lastsctlst]
;p_sza = [p_sza, lastsctsza]
endif
endif
if i ge 1 then begin
lastbutonesctr = lastsctr
lastbutoneet = lastet
endif
if i ge 0 then begin
lastsctr = sctradius
lastet = etrx
lastsctlat = sctlatitude
lastsctlon = sctlongitude
lastsctlst = sctlst
;lastsctsza = sctsza
endif


lastoccptradius = occptradius
; Store occptradius for comparison at next timestep
print, i, nsteps, occptradius, (etrx-etstart)/(etend-etstart), dt

real_et[i] = etrx
occptradiusarray[i] = occptradius
sctrarray[i] = sctradius
esdarray[i] = earthsctdist
;ltimetargeteartharray[i] = ltimetargetearth
;ltimescteartharray[i] = ltimesctearth


etold = etrx
et = etrx + dt ; finally reintroduce et variable 
; Increment timestep

i++
endwhile

i_ls = i_ls*double(180.)/!pi
e_ls = e_ls*double(180.)/!pi
; Convert angles from radians to degrees

; i_xxx and e_xxx arrays were defined with a null value 
; as their first element, then augmented
; So have to remove useless first element
; But first see if any occultations occurred

if n_elements(i_et) eq 1 or n_elements(e_et) eq 1 then begin
print, 'Where are the occultations?'
print, 'Inspect variables to see if there are any ingress occs or any egress occs'
; Stopping here does not mean that no occultations occurred
; There might be an ingress occultation, but not an egress occultation
;stop
endif

if n_elements(i_et) gt 1 then begin
i_et = i_et(1:*)
i_sctr = i_sctr(1:*)
i_sctlon = i_sctlon(1:*) 
i_sctlat = i_sctlat(1:*)
i_occr = i_occr(1:*)
i_occlon = i_occlon(1:*)
i_occlat = i_occlat(1:*)
;i_ssr = i_ssr(1:*)
;i_sslon = i_sslon(1:*)
;i_sslat = i_sslat(1:*)
i_earthsctdist = i_earthsctdist(1:*)
i_earthtargetdist = i_earthtargetdist(1:*)
i_occlst = i_occlst(1:*)
i_occsza = i_occsza(1:*)
i_ls = i_ls(1:*)
endif

if n_elements(e_et) gt 1 then begin
e_et = e_et(1:*)
e_sctr = e_sctr(1:*)
e_sctlon = e_sctlon(1:*) 
e_sctlat = e_sctlat(1:*)
e_occr = e_occr(1:*)
e_occlon = e_occlon(1:*)
e_occlat = e_occlat(1:*)
;e_ssr = e_ssr(1:*)
;e_sslon = e_sslon(1:*)
;e_sslat = e_sslat(1:*)
e_earthsctdist = e_earthsctdist(1:*)
e_earthtargetdist = e_earthtargetdist(1:*)
e_occlst = e_occlst(1:*)
e_occsza = e_occsza(1:*)
e_ls = e_ls(1:*)
endif

;p_et = p_et(1:*)
;p_sctr = p_sctr(1:*)
;p_lat = p_lat(1:*)
;p_lon = p_lon(1:*)
;p_lst = p_lst(1:*)
;p_sza = p_sza(1:*)

;real_et = real_et(1:*)
;occptradiusarray = occptradiusarray(1:*)
;sctrarray = sctrarray(1:*)
;esdarray = esdarray(1:*)
;ltimetargeteartharray = ltimetargeteartharray(1:*)
;ltimescteartharray = ltimescteartharray(1:*)

help, i_et, e_et
print, i_occlat, e_occlat, i_occsza, e_occsza

cspice_et2utc, i_et, 'C', 6, i_utcstr
cspice_et2utc, e_et, 'C', 6, e_utcstr
print, i_utcstr, e_utcstr

; Kernels are not cleaned out of memory at end of code,
; instead they're cleaned out at the beginning of the code. 

;cspice_et2utc, real_et+ltimescteartharray, 'C', 6, real_ertutc
cspice_et2utc, real_et, 'C', 6, real_ertutc ; note that real_et is et at Earth, etrx


; save, filename=path+'tests07routput.sav', real_ertutc, occptradiusarray

save, filename=path+'full_output/simpleocc_output_pat_newer.sav', real_ertutc, occptradiusarray
print, 'Finishing simpleocc processing for '+rsrfilename

;This endfor ends the loop over this particular cors file.
;If running only a single file manually, comment this out
endfor 



;spawn, 'date', codeendtime
;print, codestarttime
;print, codeendtime

stop
end
