pro simpleocc20190718

; Paul Withers, 2006.09.13
; Center for Space Physics, Boston University

; Find spacecraft occultations using SPICE

itloop = 2
; 0 - use all times in RSR file, find set of CA radial distances - repeat for both bands
; 1 - use only selected subset of times for which profile will be found, find set of CA altitudes - one band only (X)
; 2 - use only selected subset of times for which profile will be found, find matrix elements - one band only (X)
; run freqwork between 0 and 1
; run freqwork after 2

; Record time program starting running
codestarttime = systime(/seconds)

;Set rootpath
rootpath = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\'

; Set root path and SPICE paths (requires PW-2TB drive)
spicerootpath = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SPICE_kernels\'
genspicepath = spicerootpath+'generic\'
spkspicepath = spicerootpath+'spk\'


; Clean any lingering kernels out of memory here, not at end of code
cspice_ktotal, 'all', count
i=0
while i lt count do begin
  cspice_kdata, 0, 'all', file, type, source, handle, found
  cspice_unload, file
  i=i+1
endwhile

cspice_ktotal, 'all', count
if count ne 0 then begin
  print, 'Error clearing out SPICE kernels, stopping.'
  stop
endif

;Some more pathing stuff.
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


;...or choose a cors number and run every RSR file within that cors directory.
loop = 'yes'
body = 'Saturn'
allcors = ['245']

print, 'Now running itloop = ',itloop,' for cors ',allcors
wait,2

;Loop over all cors files specified by allcors
for cor = 0, n_elements(allcors)-1 do begin
  corsnumber = allcors[cor]
  allrsrpaths = file_search(rootpath+'Output\'+body+'\cors_0'+corsnumber+'\*')
  ;Loop over all rsr paths
  for rsrpath = 0, n_elements(allrsrpaths)-1 do begin
    analyze_rsr = 1 ;A flag for continuing with the analysis. If certain outputs don't already exist, analysis will be skipped.

    path = allrsrpaths[rsrpath]+'\'
    
    ;Determine which body from the path. 's' for Saturn occultations, 't' for Titan occultations
    if body eq 'Saturn' then whichbody = 's'
    if body eq 'Titan' then whichbody = 't'
    
    ;Extract the RSR filename from the path
    splitpath = strsplit(path,'\',/extract)
    rsrfilename = splitpath[-1]
    
    if (itloop gt 0) and strmid(rsrfilename,22,1) ne 'x' then begin
      ;If this is itloop 1 or 2, and this is not an x-band observation, we don't care about it. 
      analyze_rsr = 0
      print, rsrfilename, ' is not x-band.'
    endif
    
    ;Determine if this rsr is associated with an ingress or egress
    norx = ''
    if strmid(rsrfilename,6,1) eq 'i' then norx = 'N'
    if strmid(rsrfilename,6,1) eq 'e' then norx = 'X
    if (body eq 'Saturn') and (corsnumber eq 213) then norx = 'X' ;Hard-coding for this particular observation. 
    if norx eq '' then begin
      print, 'ERROR: could not determine whether rsr is ingress or egress. Stopping.'
      stop
    endif
           
    ;Restore previous output if you're on loops 1 or 2. 
    if (itloop eq 1) and (analyze_rsr eq 1) then begin
      restore, filename=path+'simpleocc_output_itloop0.sav'
      extraoccptradiusarray = itloop0occptradiusarray
      extraoccptlatarray = itloop0occptlatarray
      extraoccptlonarray = itloop0occptlonarray
      extraoccptszaarray = itloop0occptszaarray
      extraoccptlstarray = itloop0occptlstarray
      extrasepanglearray = itloop0sepanglearray
      extraepsanglearray = itloop0epsanglearray
      extraettxarray = itloop0ettxarray
      extraetoccptarray = itloop0etoccptarray
      extraetrxarray = itloop0etrxarray
      extrautctxarray = itloop0utctxarray
      extrautcoccptarray = itloop0utcoccptarray
      extrautcrxarray = itloop0utcrxarray
      
      ;Check that freqwork a was performed 
      freqwork_a_bool = file_test(path+'full_output\freqwork_a_output.sav')
      if freqwork_a_bool eq 1 then begin
        restore, path+'full_output\freqwork_a_output.sav'
        carryoverutcstart = start_time_save
        carryoverutcend = stop_time_save 
      endif else begin
        analyze_rsr = 0
        ;If freqwork a output does not exist, either you need to run the offending rsr files through freqwork itloop = a, 
        ; or else this particular rsr doesn't match up with any others during the flyby (happens fairly frequently). 
        print, 'Freqwork a output does not exist for ',rsrfilename
      endelse
    endif
  
    if (itloop eq 2) and (analyze_rsr eq 1) then begin
      simpleocc_1_bool = file_test(path+'simpleocc_output_itloop1.sav')
      ;Check that simpleocc 1 output exists. 
      if simpleocc_1_bool eq 1 then begin
        restore, filename=path+'simpleocc_output_itloop1.sav'
        ;extraoccptradiusarray = occptradiusarray
        extraoccptlatarray = itloop1occptlatarray
        extraoccptlonarray = itloop1occptlonarray
        extraoccptszaarray = itloop1occptszaarray
        extraoccptlstarray = itloop1occptlstarray
        extrasepanglearray = itloop1sepanglearray
        extraepsanglearray = itloop1epsanglearray
        extraettxarray = itloop1ettxarray
        extraetoccptarray = itloop1etoccptarray
        extraetrxarray = itloop1etrxarray
        extrautctxarray = itloop1utctxarray
        extrautcoccptarray = itloop1utcoccptarray
        extrautcrxarray = itloop1utcrxarray
        
        extraoccptradiusarray = itloop1occptrarray
        extraoccptaltspharray = itloop1occptaltspharray
        extraoccptaltellarray = itloop1occptaltellarray
        extraoccptaltasharray = itloop1occptaltasharray
        extraoccptaltkosarray = itloop1occptaltkosarray
        extraoccptlatspharray = itloop1occptlatspharray
        extraoccptlatellarray = itloop1occptlatellarray
        extraoccptlatasharray = itloop1occptlatasharray
        extraoccptlatkosarray = itloop1occptlatkosarray
        
        ;Check that freqwork a was performed
        freqwork_a_bool = file_test(path+'full_output\freqwork_a_output.sav')
        if freqwork_a_bool eq 1 then begin
          restore, path+'full_output\freqwork_a_output.sav'
          carryoverutcstart = start_time_save
          carryoverutcend = stop_time_save
        endif else begin
          analyze_rsr = 0
          ;If freqwork a output does not exist, either you need to run the offending rsr files through freqwork itloop = a,
          ; or else this particular rsr doesn't match up with any others during the flyby (happens fairly frequently).
          print, 'Freqwork a output does not exist for ',rsrfilename
        endelse
      endif else begin
        ;If simpleocc 1 output doesn't exist, you either forgot to run it for this rsr, or it doesn't match with another rsr file in the 
        ; specified cors directories. 
        analyze_rsr = 0
        print, 'Simpleocc 1 output does not exist for ',rsrfilename
      endelse
    endif
  
    if analyze_rsr eq 0 then begin
      print, 'Skipping.'
      print, ''
      wait, 2
    endif else begin
      print, 'Now analyzing ',rsrfilename,', itloop = ',itloop
      wait, 2
      ;Determine the appropropriate .bsp spice files to load for this RSR file.
      ;Once found, load that file or those files into spice
      openr, lun, rootpath+'Programs\RSR_SPICE_connection.txt', /get_lun
      line = ''
      while not eof(lun) do begin
        readf, lun, line
        ;Check if this is the line for this RSR file
        if strmatch(line,rsrfilename+'*') then begin
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
      
      ; Make SPICE work in J2000 frame as default
      frame='j2000'
      
      ; Used for ray-tracing, I think
      ; See help pages about commands that use this variable
      method='intercept'
      
      ; Don't do any aberration corrections
      ; See help pages about commands that use this variable
      abcorr='none'
      
      ; Radius of target body, km
      ; Values increased from reference equatorial radius since Cassini
      ; publications report occultation conditions at ionospheric peak
      if whichbody eq 't' then rtarget = double(2575.) + 1200. ; km, Titan
      if whichbody eq 's' then rtarget = double(60268.) + 2000. ; km, Saturn
      
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
      
      fakehhstartlon = floor(fakesfdusecoutstart/ 60d^2) ;Convert seconds to hours
      fakehhstartstr = strcompress(fakehhstartlon, /remove_all)
      while strlen(fakehhstartstr) lt 2 do fakehhstartstr = '0' + fakehhstartstr
      if strlen(fakehhstartstr) ne 2 then stop
      
      fakemmstartlon = floor( (fakesfdusecoutstart - fakehhstartlon*60d^2) / 60d) ;Remainder to minutes
      fakemmstartstr = strcompress(fakemmstartlon, /remove_all)
      while strlen(fakemmstartstr) lt 2 do fakemmstartstr = '0' + fakemmstartstr
      if strlen(fakemmstartstr) ne 2 then stop
      
      fakessstartdbl = fakesfdusecoutstart - fakehhstartlon*60d^2 - fakemmstartlon*60d ;Remainder to seconds
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
      
      etstart = fakeetstart
      etend = fakeetend
      
      nsteps = ceil( (etend- etstart)/dt) ; Number of steps
      
      if itloop eq 1 or itloop eq 2 then begin ; restrict to smaller range of times
        cspice_str2et, carryoverutcstart, etstart
        cspice_str2et, carryoverutcend, etend
        
        nsteps = ceil( (etend- etstart)/dt) + 1 ; +1 seems necessary for this route
        nsteps = ceil(round((etend - etstart)/dt)) + 1 ; sometimes a tiny difference gets ceil'ed up all the way to an integer, so try this
        ; if wrong, matrix dimension incompatability when multiplying
        a = 60268d
        b = 54364d
        thisrrefsph = 60268d
        
        ndummy = 1d4
        dummylatdeg = dindgen(ndummy)/(ndummy-1d) * 180d - 90d
        dummyrrefell = a*b/sqrt((b*cos(dummylatdeg*!pi/180.))^2.+(a*sin(dummylatdeg*!pi/180.))^2.)
        dummyrrefash = sat_anderson_schubert_1bar(dummylatdeg)
        dummyrrefkos = saturn_1bar(90*!pi/180.-!pi/180*dummylatdeg)
      endif
      
      ;Flip some stuff around if an egress occ. Paul W. implementeed this to for simplicity, rather having two sets of rules for ingress/egress.
      if norx eq 'X' then begin
        blah = etstart
        etstart = etend
        etend = blah
      endif
      
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
      
      sunint = 10L
      sunstr = '10'
      ; Label for Sun, necessary for calculating angles between bodies and the probe
      
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
      i_occr = i_et ; distance between occultation point and centre of target at occ, should be planetary radius to accuracy permitted by timestep, km
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
      occptradiusarray = i_et ; distance between centre of target and closest approach of ray path to target, km, at each timestep
      sctrarray = i_et ; distance between centre of target and spacecraft at each timestep, km
      esdarray = sctrarray ; distance between Earth and spacecraft at each timestep, km
      ltimetargeteartharray = i_et ; Light time between target and Earth at each timestep, UNITS
      ltimescteartharray = i_et ; Light time between sct and Earth at each timestep, UNITS
      
      et = etstart
      
      ;Replacing arrays with definite sizes.
      real_et = dblarr(nsteps)
      occptradiusarray = dblarr(nsteps)
      sctrarray = dblarr(nsteps)
      esdarray = dblarr(nsteps)
      ltimetargeteartharray = dblarr(nsteps)
      ltimescteartharray = dblarr(nsteps)
      
      if itloop eq 1 then itloop1occptrarray = dblarr(nsteps)
      occptaltspharray = dblarr(nsteps)
      occptaltellarray = dblarr(nsteps)
      occptaltasharray = dblarr(nsteps)
      occptaltkosarray = dblarr(nsteps)
      
      occptlatspharray = dblarr(nsteps)
      occptlatellarray = dblarr(nsteps)
      occptlatasharray = dblarr(nsteps)
      occptlatkosarray = dblarr(nsteps)
      
      occptlatarray = dblarr(nsteps) ; Latitude (degN) of occultation point
      occptlonarray = dblarr(nsteps) ; Longitude (degE) of occultation point
      occptszaarray = dblarr(nsteps) ; SZA (deg) of occultation point
      occptlstarray = dblarr(nsteps) ; LST (hrs) of occultation point
      sepanglearray = dblarr(nsteps) ; Sun-Earth-Spacecraft angle (deg)
      epsanglearray = dblarr(nsteps) ; Earth-Spacecraft-Sun angle (deg)
      ettxarray = dblarr(nsteps) ; Transmission time, ET
      etoccptarray = dblarr(nsteps) ; Time as ray passes through closest approach point, ET
      etrxarray = dblarr(nsteps) ; Receiption time, ET
      ;;; still need time in various systems, ET and UTC
      
      if itloop eq 2 then begin
        inpathlength2d = dblarr(n_elements(extraoccptradiusarray),n_elements(extraoccptradiusarray))
        outpathlength2d = inpathlength2d
        
        altsphinpathlength2d = inpathlength2d
        altsphoutpathlength2d = inpathlength2d
        altellinpathlength2d = inpathlength2d
        altelloutpathlength2d = inpathlength2d
        altashinpathlength2d = inpathlength2d
        altashoutpathlength2d = inpathlength2d
        altkosinpathlength2d = inpathlength2d
        altkosoutpathlength2d = inpathlength2d
      endif
      
      i=0L
      while i lt nsteps do begin ;Begin nsteps loop
        ;while et lt etend do begin
        
        etrx = et ; rx means receiver, tx means transmitter
        ; use this time for Earth
        et = 'crashifthisiseverused'
        ; ensure this variable doesn't play direct role in rest of code...
        
        cspice_spkpos, earthstr, etrx, frame, 'LT', sctstr, ptarg, txrxltime ; choice of 'LT' means that et is time at the receiver
        ettx = etrx - txrxltime ;use this time for events in the Saturn system, so Cassini, Saturn, Titan
        
        ; insert big block of code to find etoccpt
        
        nconvergence = 1000L
        targetobject = targetstr
        sct = sctstr
        dss = earthstr
        
        ; Define the precision to solve time t-O to, units of seconds -- lower value is higher precision
        precisionto = 4d-9
        
        ; Define the precision to solve k-O to, dimensionless -- lower value is higher precision
        precisionko = 3d-16
        
        tc = ettx
        td = etrx
        
        ; Prepare initial estimates of t-OCD and k-OCD
        tocd = tc ; sct is closer to target than DSS is
        kocd = 0.
        
        ; Ensure that iterative process does not declare instant convergence
        dtocd = 10.d*precisionto
        dkocd = 10.d*precisionko
        
        itocd = 0
        while abs(dtocd) gt precisionto and abs(dkocd) gt precisionko do begin
        
          ; Values of t-OCD and k-OCD at start of loop
          ; Necessary to judge whether values at end of loop are any better or not
          oldtocd = tocd
          oldkocd = kocd
          
          ; targetssb is position of target object relative to solar system barycentre at time t-OCD, km
          cspice_spkpos, '0', tocd, frame, 'LT', targetobject, mtargetssb, ltimetargetssb
          targetssb = 0d - mtargetssb
          
          ; sctssb is position of spacecraft relative to solar system barycentre at time t-C, km
          cspice_spkpos, '0', tc, frame, 'LT', sct, msctssb, ltimesctssb
          sctssb = 0d - msctssb
          
          ; earthssb is position of DSS relative to solar system barycentre at time t-D, km
          cspice_spkpos, '0', td, frame, 'LT', dss, mearthssb, ltimeearthssb
          earthssb = 0d - mearthssb
          
          ; Position of target object relative to spacecraft in working frame, km
          targetsct = targetssb - sctssb
          ; Position of DSS relative to spacecraft in working frame, km
          earthsct = earthssb - sctssb
          
          ; Value of kocd, which is related to location of occultation point along sct-DSS line, dimensionless
          kocd = cspice_vdot(targetsct, earthsct) / cspice_vdot(earthsct, earthsct)
          
          ; occptssb is position of occultation point relative to solar system barycentre at time t-OCD
          ; Position given in previously specified frame (the variable frame, probably J2000) and in units of km
          occptssb = sctssb + kocd * (earthssb - sctssb)
          
          ; occptsct is position of occultation point at time t-OCD relative to the sct at time t-C
          ; Position given in previously specified frame (the variable frame, probably J2000) and in units of km
          occptsct = occptssb - sctssb
          
          ; Calculate new value of occultation time t-OCD using light time between position of spacecraft at time t-C and position of target object at previous estimate of time t-OCD
          tocd = tc + norm(occptsct) / cspice_clight()
          
          ; Update values of time t-O and dimensionless ko
          dtocd = tocd - oldtocd
          dkocd = kocd - oldkocd
          ;print, to, ko, dto, dko
          
          ; Increment index counter
          itocd++
          
          ; Escape from loop if convergence does not occur
          if itocd gt nconvergence then begin
            print,'STOP - Could not find accurate t-O/k-O after many iterations for iobs element: ', iobs
            stop
          endif
          
        endwhile ; while abs(dto) gt precisionto and abs(dko) gt precisionko do begin
        etoccpt = tocd
        
        ;;; end of big block of code to find etoccpt
        
        ;; XXX GREAT BIG NEW BIT
        if itloop eq 1 or itloop eq 2 then begin
        
          nxx = 1d3
          oldkocd = kocd
          xx = (dindgen(2*nxx+1) - nxx) * 1d-7 ; code is faster if nxx is smaller, but larger nxx needed to capture full extent of onion
          distarr = xx * 0d
          altsphdistarr = distarr
          altelldistarr = distarr
          altashdistarr = distarr
          altkosdistarr = distarr
          templatdegarr = distarr
          rrefspharr = distarr
          rrefellarr = distarr
          rrefasharr = distarr
          rrefkosarr = distarr
          
          cspice_tipbod, frame, targetint, etoccpt, tipm
          
          k=0L
          while k lt n_elements(xx) do begin
            kocd = oldkocd + xx[k]
            occptssb = sctssb + kocd * (earthssb - sctssb)
            distarr[k] = norm(occptssb - targetssb)
            cspice_mxv, tipm, occptssb - targetssb, bfoccpt
            cspice_reclat, bfoccpt, occptradius, occptlongitude, occptlatitude
            templatdegarr[k] = occptlatitude * 180./!pi
            k++
          endwhile
          
          rrefspharr = thisrrefsph + templatdegarr*0d
          rrefellarr = interpol(dummyrrefell,dummylatdeg,templatdegarr)
          rrefasharr = interpol(dummyrrefash,dummylatdeg,templatdegarr)
          rrefkosarr = interpol(dummyrrefkos,dummylatdeg,templatdegarr)
          ; faster than direct call to subroutines
          
          altsphdistarr = distarr - rrefspharr
          altelldistarr = distarr - rrefellarr
          altashdistarr = distarr - rrefasharr
          altkosdistarr = distarr - rrefkosarr
          
          itloop1occptrarray[i] = min(distarr)
          occptaltspharray[i] = min(altsphdistarr)
          occptaltellarray[i] = min(altelldistarr)
          occptaltasharray[i] = min(altashdistarr)
          occptaltkosarray[i] = min(altkosdistarr)
          
          aar = where(distarr eq min(distarr)) & aar = aar[0]
          aasph = where(altsphdistarr eq min(altsphdistarr)) & aasph = aasph[0]
          aaell = where(altelldistarr eq min(altelldistarr)) & aaell = aaell[0]
          aaash = where(altashdistarr eq min(altashdistarr)) & aaash = aaash[0]
          aakos = where(altkosdistarr eq min(altkosdistarr)) & aakos = aakos[0]
          
          occptlatspharray[i] = templatdegarr[aasph]
          occptlatellarray[i] = templatdegarr[aaell]
          occptlatasharray[i] = templatdegarr[aaash]
          occptlatkosarray[i] = templatdegarr[aakos]
          
          if itloop eq 2 then begin
          
            flankr = dblarr(i+1) ; assumes ingress...
            flankzsph = flankr
            flankzell = flankr
            flankzash = flankr
            flankzkos = flankr
            j=0
            while j lt i do begin
              flankr[j] = extraoccptradiusarray[i-j] + (extraoccptradiusarray[i-j-1]-extraoccptradiusarray[i-j])*0.5d
              flankzsph[j] = extraoccptaltspharray[i-j] + (extraoccptaltspharray[i-j-1]-extraoccptaltspharray[i-j])*0.5d
              flankzell[j] = extraoccptaltellarray[i-j] + (extraoccptaltellarray[i-j-1]-extraoccptaltellarray[i-j])*0.5d
              flankzash[j] = extraoccptaltasharray[i-j] + (extraoccptaltasharray[i-j-1]-extraoccptaltasharray[i-j])*0.5d
              flankzkos[j] = extraoccptaltkosarray[i-j] + (extraoccptaltkosarray[i-j-1]-extraoccptaltkosarray[i-j])*0.5d
              
              j++
            endwhile
            flankr[j] = extraoccptradiusarray[0] + (extraoccptradiusarray[0]-extraoccptradiusarray[1])*0.5d
            flankzsph[j] = extraoccptaltspharray[0] + (extraoccptaltspharray[0]-extraoccptaltspharray[1])*0.5d
            flankzell[j] = extraoccptaltellarray[0] + (extraoccptaltellarray[0]-extraoccptaltellarray[1])*0.5d
            flankzash[j] = extraoccptaltasharray[0] + (extraoccptaltasharray[0]-extraoccptaltasharray[1])*0.5d
            flankzkos[j] = extraoccptaltkosarray[0] + (extraoccptaltkosarray[0]-extraoccptaltkosarray[1])*0.5d
            
            centerr = extraoccptradiusarray[i]
            centeraltsph = extraoccptaltspharray[i]
            centeraltell = extraoccptaltellarray[i]
            centeraltash = extraoccptaltasharray[i]
            centeraltkos = extraoccptaltkosarray[i]
            
            if min(distarr) gt min(flankr) then stop
            ; the segment of the line of sight used for this ray must include all closest approach distances in the set of rays
            if max(distarr) lt max(flankr) then stop
            if min(altsphdistarr) gt min(flankzsph) then stop
            if max(altsphdistarr) lt max(flankzsph) then stop
            if min(altelldistarr) gt min(flankzell) then stop
            if max(altelldistarr) lt max(flankzell) then stop
            if min(altashdistarr) gt min(flankzash) then stop
            if max(altashdistarr) lt max(flankzash) then stop
            if min(altkosdistarr) gt min(flankzkos) then stop
            if max(altkosdistarr) lt max(flankzkos) then stop
            
            aar = where(distarr eq min(distarr)) & aar = aar[0]
            aasph = where(altsphdistarr eq min(altsphdistarr)) & aasph = aasph[0]
            aaell = where(altelldistarr eq min(altelldistarr)) & aaell = aaell[0]
            aaash = where(altashdistarr eq min(altashdistarr)) & aaash = aaash[0]
            aakos = where(altkosdistarr eq min(altkosdistarr)) & aakos = aakos[0]
            
            cakr = oldkocd + xx[aar]
            caksph = oldkocd + xx[aasph]
            cakell = oldkocd + xx[aaell]
            cakash = oldkocd + xx[aaash]
            cakkos = oldkocd + xx[aakos]
            
            inflankk = flankr * 0d
            outflankk = flankr * 0d
            altsphinflankk = flankr * 0d
            altsphoutflankk = flankr * 0d
            altellinflankk = flankr * 0d
            altelloutflankk = flankr * 0d
            altashinflankk = flankr * 0d
            altashoutflankk = flankr * 0d
            altkosinflankk = flankr * 0d
            altkosoutflankk = flankr * 0d
            
            j=0
            while j lt n_elements(flankr) do begin
              ppinr = interpol(xx[0:aar], distarr[0:aar], flankr[j])
              inflankk[j] = oldkocd + ppinr[0]
              altsphppin = interpol(xx[0:aasph], altsphdistarr[0:aasph], flankzsph[j])
              altsphinflankk[j] = oldkocd + altsphppin[0]
              altellppin = interpol(xx[0:aaell], altelldistarr[0:aaell], flankzell[j])
              altellinflankk[j] = oldkocd + altellppin[0]
              altashppin = interpol(xx[0:aaash], altashdistarr[0:aaash], flankzash[j])
              altashinflankk[j] = oldkocd + altashppin[0]
              altkosppin = interpol(xx[0:aakos], altkosdistarr[0:aakos], flankzkos[j])
              altkosinflankk[j] = oldkocd + altkosppin[0]
              
              ppoutr = interpol(xx[aar:*], distarr[aar:*], flankr[j])
              outflankk[j] = oldkocd + ppoutr[0]
              altsphppout = interpol(xx[aasph:*], altsphdistarr[aasph:*], flankzsph[j])
              altsphoutflankk[j] = oldkocd + altsphppout[0]
              altellppout = interpol(xx[aaell:*], altelldistarr[aaell:*], flankzell[j])
              altelloutflankk[j] = oldkocd + altellppout[0]
              altashppout = interpol(xx[aaash:*], altashdistarr[aaash:*], flankzash[j])
              altashoutflankk[j] = oldkocd + altashppout[0]
              altkosppout = interpol(xx[aakos:*], altkosdistarr[aakos:*], flankzkos[j])
              altkosoutflankk[j] = oldkocd + altkosppout[0]
              
              j++
            endwhile
            
            inpathlength = flankr * 0d
            outpathlength = flankr * 0d
            altsphinpathlength = flankr * 0d
            altsphoutpathlength = flankr * 0d
            altellinpathlength = flankr * 0d
            altelloutpathlength = flankr * 0d
            altashinpathlength = flankr * 0d
            altashoutpathlength = flankr * 0d
            altkosinpathlength = flankr * 0d
            altkosoutpathlength = flankr * 0d
            
            j=0
            while j lt n_elements(flankr) do begin
            
              earthsctdistance = norm(earthssb - sctssb)
              
              if j eq 0 then begin
                inpathlength[j] = abs(inflankk[j] - cakr) * earthsctdistance
                altsphinpathlength[j] = abs(altsphinflankk[j] - caksph) * earthsctdistance
                altellinpathlength[j] = abs(altellinflankk[j] - cakell) * earthsctdistance
                altashinpathlength[j] = abs(altashinflankk[j] - cakash) * earthsctdistance
                altkosinpathlength[j] = abs(altkosinflankk[j] - cakkos) * earthsctdistance
                outpathlength[j] = abs(outflankk[j] - cakr) * earthsctdistance
                altsphoutpathlength[j] = abs(altsphoutflankk[j] - caksph) * earthsctdistance
                altelloutpathlength[j] = abs(altelloutflankk[j] - cakell) * earthsctdistance
                altashoutpathlength[j] = abs(altashoutflankk[j] - cakash) * earthsctdistance
                altkosoutpathlength[j] = abs(altkosoutflankk[j] - cakkos) * earthsctdistance
              endif
              
              if j ne 0 then begin
                inpathlength[j] = abs(inflankk[j] - inflankk[j-1]) * earthsctdistance
                altsphinpathlength[j] = abs(altsphinflankk[j] - altsphinflankk[j-1]) * earthsctdistance
                altellinpathlength[j] = abs(altellinflankk[j] - altellinflankk[j-1]) * earthsctdistance
                altashinpathlength[j] = abs(altashinflankk[j] - altashinflankk[j-1]) * earthsctdistance
                altkosinpathlength[j] = abs(altkosinflankk[j] - altkosinflankk[j-1]) * earthsctdistance
                outpathlength[j] = abs(outflankk[j] - outflankk[j-1]) * earthsctdistance
                altsphoutpathlength[j] = abs(altsphoutflankk[j] - altsphoutflankk[j-1]) * earthsctdistance
                altelloutpathlength[j] = abs(altelloutflankk[j] - altelloutflankk[j-1]) * earthsctdistance
                altashoutpathlength[j] = abs(altashoutflankk[j] - altashoutflankk[j-1]) * earthsctdistance
                altkosoutpathlength[j] = abs(altkosoutflankk[j] - altkosoutflankk[j-1]) * earthsctdistance
              endif
              
              j++
            endwhile
            
            inpathlength2d[i,0:n_elements(inpathlength)-1] = reverse(inpathlength[*]) ; making sure segments match to correct N
            altsphinpathlength2d[i,0:n_elements(altsphinpathlength)-1] = reverse(altsphinpathlength[*])
            altellinpathlength2d[i,0:n_elements(altellinpathlength)-1] = reverse(altellinpathlength[*])
            altashinpathlength2d[i,0:n_elements(altashinpathlength)-1] = reverse(altashinpathlength[*])
            altkosinpathlength2d[i,0:n_elements(altkosinpathlength)-1] = reverse(altkosinpathlength[*])
            outpathlength2d[i,0:n_elements(outpathlength)-1] = reverse(outpathlength[*])
            altsphoutpathlength2d[i,0:n_elements(altsphoutpathlength)-1] = reverse(altsphoutpathlength[*])
            altelloutpathlength2d[i,0:n_elements(altelloutpathlength)-1] = reverse(altelloutpathlength[*])
            altashoutpathlength2d[i,0:n_elements(altashoutpathlength)-1] = reverse(altashoutpathlength[*])
            altkosoutpathlength2d[i,0:n_elements(altkosoutpathlength)-1] = reverse(altkosoutpathlength[*])
            
          endif
        endif  ; END GREAT BIG NEW BIT
        
        ; Position of Earth at time etrx relative to solar system barycentre and in specified frame
        cspice_spkgps, earthint, etrx, frame, ssbint, earthssb, ltimeearthssb
        
        ; Position of target body at time etoccpt relative to solar system barycentre and in specified frame
        cspice_spkgps, targetint, etoccpt, frame, ssbint, targetssb, ltimetargetssb
        
        ; Position of spacecraft at time ettx relative to solar system barycentre and in specified frame
        cspice_spkgps, sctint, ettx, frame, ssbint, sctssb, ltimemgsssb
        
        ; Position of Earth relative to spacecraft in working frame
        earthsct = earthssb - sctssb
        
        ; Position of Sun at time etrx relative to solar system barycentre and in specified frame
        cspice_spkgps, sunint, etrx, frame, ssbint, sunssb, ltimesunssb
        
        ; Position of Earth relative to Sun, not ssb, in working frame
        earthsun = earthssb - sunssb
        
        ; Sun-Earth-Spacecraft angle, redundant minus signs ensure that directions of vectors are "Sun relative to Earth" and "sct relative to Earth"
        sepangle = cspice_vsep( (-1d)*earthsun, (-1d)*earthsct ) * 180./!pi
        
        ; Earth-Spacecraft-Sun angle
        epsangle = cspice_vsep( earthsct, sunssb-sctssb) * 180./!pi
        
        ; pnear is nearest point on Earth-spacecraft line to target body
        ; Existence of occultations will depend a lot on pnear
        ; I loosely call position indicated by pnear the "occultation point"
        ; See help pages for frame in which pnear is stated
        ; distoccpttarget is distance between pnear and target body, km
        cspice_nplnpt, earthssb, earthsct, targetssb, pnear, distoccpttarget
        
        ; Position of pnear relative to target body in working frame
        occpttarget = pnear - targetssb
        
        ; tipm is transformation matrix between specified inertial frame (frame, probably J2000)
        ; and body-fixed frame of target body
        ; See help pages for more on relevant frames
        cspice_tipbod, frame, targetint, etoccpt, tipm
        
        ; Occultation point position converted from one frame (occpttarget)
        ; to body-fixed frame of target body (bfoccpt)
        cspice_mxv, tipm, occpttarget, bfoccpt
        
        ; Convert bfoccpt from cartesian coordinates to r/lat/lon coordinates (km/radE/radN)
        cspice_reclat, bfoccpt, occptradius, occptlongitude, occptlatitude
        
        ; Convert angles from radians into degrees
        occptlongitude = occptlongitude * 180./!pi
        occptlatitude = occptlatitude * 180./!pi
        
        ; Position of spacecraft relative to target body in working frame
        scttarget = sctssb - targetssb
        
        ; Spacecraft position transformed from one frame (scttarget) to body-fixed frame of target body (bfsct)
        cspice_mxv, tipm, scttarget, bfsct
        
        ; Convert bfsct from cartesian coordinates to r/lat/lon coordinates (km/radE/radN)
        cspice_reclat, bfsct, sctradius, sctlongitude, sctlatitude
        
        ; Convert angles from radians into degrees
        sctlongitude = sctlongitude * 180./!pi
        sctlatitude = sctlatitude * 180./!pi
        
        ; earthsctdist is distance between Earth and spacecraft
        earthsctdist = cspice_vdist(earthssb, sctssb)
        
        ; earthtargetdist is distance between Earth and target body
        earthtargetdist = cspice_vdist(earthssb, targetssb)
        
        ; Difference between these two distances
        ; This program examines whether an infinite line defined
        ; by two points, Earth and the spacecraft, intersects the surface of the target body
        ; If it does, then an occultation only occurs if the spacecraft is behind the
        ; target body. That is, if sctminustarget is positive.s
        sctminustarget = earthsctdist - earthtargetdist
        
        ; ls (longitude of the sun, planetocentric, radians) is often used as a measure of seasons
        ; Especially common on Mars, see help pages for details
        ls = cspice_lspcn(targetstr, etoccpt, abcorr)
        
        ; lst is local solar time (0-24 hours) of occultation point at time et
        cspice_et2lst, etoccpt, targetint, !pi/180.*occptlongitude, 'planetocentric', $
          hr, mn, sc, time, ampm
        lst = hr + (mn + sc/double(60.))/double(60.)
        
        ; sctlst is local solar time (0-24 hours) of spacecraft at time et
        cspice_et2lst, etoccpt, targetint, !pi/180.*sctlongitude, 'planetocentric', $
          hr, mn, sc, time, ampm
        sctlst = hr + (mn + sc/double(60.))/double(60.)
        
        ; targetsubsolpt is position of sub-solar point on target body at time et
        ; See help pages for frame in which targetsubsolpt is stated
        ; Note that this may use a triaxial ellipsoid for the target body, which is not desirable for the SZA purpose
        cspice_subsol, method, targetstr, etoccpt, abcorr, earthstr, targetsubsolpt
        
        ; Solar zenith angle (degrees) at occultation point
        ; Angle between subsolar point-centre of target body vector and
        ; occultation point-centre of target body vector
        sza = cspice_vsep(targetsubsolpt, bfoccpt)*180./!pi
        
        if i gt 0L then begin
          ; if occultation point was outside target radius at last timestep and
          ; occultation point is inside target radius at this timestep and
          ; spacecraft is behind target body, then
          ; an ingress occultation occurs
          if (occptradius lt rtarget and lastoccptradius ge rtarget and $
            sctminustarget gt double(0.)) then begin
            i_et = [i_et, ettx]
            i_sctr = [i_sctr, sctradius]
            i_sctlon = [i_sctlon, sctlongitude]
            i_sctlat = [i_sctlat, sctlatitude]
            i_occr = [i_occr, occptradius]
            i_occlon = [i_occlon, occptlongitude]
            i_occlat = [i_occlat, occptlatitude]
            i_earthsctdist = [i_earthsctdist, earthsctdist]
            i_earthtargetdist = [i_earthtargetdist, earthtargetdist]
            i_occlst = [i_occlst, lst]
            i_occsza = [i_occsza, sza]
            i_ls = [i_ls, ls]
          endif
          ; if occultation point was inside target radius at last timestep and
          ; occultation point is outside target radius at this timestep and
          ; spacecraft is behind target body, then
          ; an egress occultation occurs
          if (occptradius ge rtarget and lastoccptradius lt rtarget and $
            sctminustarget gt double(0.)) then begin
            e_et = [e_et, ettx]
            e_sctr = [e_sctr, sctradius]
            e_sctlon = [e_sctlon, sctlongitude]
            e_sctlat = [e_sctlat, sctlatitude]
            e_occr = [e_occr, occptradius]
            e_occlon = [e_occlon, occptlongitude]
            e_occlat = [e_occlat, occptlatitude]
            e_earthsctdist = [e_earthsctdist, earthsctdist]
            e_earthtargetdist = [e_earthtargetdist, earthtargetdist]
            e_occlst = [e_occlst, lst]
            e_occsza = [e_occsza, sza]
            e_ls = [e_ls, ls]
          endif
        endif
        
        if i ge 2 then begin
          ; if radial distance between spacecraft and target body at last timestep
          ; is smaller than at preceding timestep and smaller than at current timestep,
          ; then last timestep was closest approach of spacecraft to target body
          ; Effectively a periapsis, if spacecraft is orbiting target body
          if (lastbutonesctr gt lastsctr) and (sctradius gt lastsctr) then begin
            p_et = [p_et, lastet]
            p_sctr = [p_sctr, lastsctr]
            p_lat = [p_lat, lastsctlat]
            p_lon = [p_lon, lastsctlon]
            p_lst = [p_lst, lastsctlst]
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
        endif
        
        ; Store occptradius for comparison at next timestep
        lastoccptradius = occptradius
        ;print, i, nsteps, occptradius, (etrx-etstart)/(etend-etstart), dt
        print, i, nsteps, occptradius ;The last two really aren't necessary to print. 

        
        ;Store the output data
        real_et[i] = etrx
        occptradiusarray[i] = occptradius
        sctrarray[i] = sctradius
        esdarray[i] = earthsctdist
        occptlatarray[i] = occptlatitude
        occptlonarray[i] = occptlongitude
        occptszaarray[i] = sza
        occptlstarray[i] = lst
        sepanglearray[i] = sepangle
        epsanglearray[i] = epsangle
        ettxarray[i] = ettx
        etoccptarray[i] = etoccpt
        etrxarray[i] = etrx
        etold = etrx
        
        if norx eq 'N' then et = etrx + dt ; finally reintroduce et variable
        if norx eq 'X' then et = etrx - dt ; finally reintroduce et variable
        
        ; Increment timestep
        i++
      endwhile ;End nsteps loop
      
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
        e_earthsctdist = e_earthsctdist(1:*)
        e_earthtargetdist = e_earthtargetdist(1:*)
        e_occlst = e_occlst(1:*)
        e_occsza = e_occsza(1:*)
        e_ls = e_ls(1:*)
      endif
      
      
      help, i_et, e_et
      print, i_occlat, e_occlat, i_occsza, e_occsza
      cspice_et2utc, i_et, 'C', 6, i_utcstr
      cspice_et2utc, e_et, 'C', 6, e_utcstr
      print, i_utcstr, e_utcstr
      
      ; Kernels are not cleaned out of memory at end of code,
      ; instead they're cleaned out at the beginning of the code.
      
      cspice_et2utc, real_et, 'C', 6, real_ertutc ; note that real_et is et at Earth, etrx
      
      ;Convert all ET times to UTC times and save
      cspice_et2utc, ettxarray, 'C', 6, utctxarray
      cspice_et2utc, etoccptarray, 'C', 6, utcoccptarray
      cspice_et2utc, etrxarray, 'C', 6, utcrxarray
      
      if itloop eq 0 then begin
        itloop0occptradiusarray = occptradiusarray
        itloop0occptlatarray = occptlatarray
        itloop0occptlonarray = occptlonarray
        itloop0occptszaarray = occptszaarray
        itloop0occptlstarray = occptlstarray
        itloop0sepanglearray = sepanglearray
        itloop0epsanglearray = epsanglearray
        itloop0ettxarray = ettxarray
        itloop0etoccptarray = etoccptarray
        itloop0etrxarray = etrxarray
        itloop0utctxarray = utctxarray
        itloop0utcoccptarray = utcoccptarray
        itloop0utcrxarray = utcrxarray
        save, filename=path+'simpleocc_output_itloop0.sav', itloop0occptradiusarray, itloop0occptlatarray, itloop0occptlonarray, $
          itloop0occptszaarray, itloop0occptlstarray, itloop0sepanglearray, itloop0epsanglearray, itloop0ettxarray, itloop0etoccptarray, $
          itloop0etrxarray, itloop0utctxarray, itloop0utcoccptarray, itloop0utcrxarray
      endif
      
      if itloop eq 1 then begin
        itloop1occptradiusarray = occptradiusarray
        itloop1occptlatarray = occptlatarray
        itloop1occptlonarray = occptlonarray
        itloop1occptszaarray = occptszaarray
        itloop1occptlstarray = occptlstarray
        itloop1sepanglearray = sepanglearray
        itloop1epsanglearray = epsanglearray
        itloop1ettxarray = ettxarray
        itloop1etoccptarray = etoccptarray
        itloop1etrxarray = etrxarray
        itloop1utctxarray = utctxarray
        itloop1utcoccptarray = utcoccptarray
        itloop1utcrxarray = utcrxarray
        ; itloop1occptrarray ; already defined
        itloop1occptaltspharray = occptaltspharray
        itloop1occptaltellarray = occptaltellarray
        itloop1occptaltasharray = occptaltasharray
        itloop1occptaltkosarray = occptaltkosarray
        itloop1occptlatspharray = occptlatspharray
        itloop1occptlatellarray = occptlatellarray
        itloop1occptlatasharray = occptlatasharray
        itloop1occptlatkosarray = occptlatkosarray
        
        save, filename=path+'simpleocc_output_itloop1.sav', itloop1occptradiusarray, itloop1occptlatarray, itloop1occptlonarray, $
          itloop1occptszaarray, itloop1occptlstarray, itloop1sepanglearray, itloop1epsanglearray, itloop1ettxarray, $
          itloop1etoccptarray, itloop1etrxarray, itloop1utctxarray, itloop1utcoccptarray, itloop1utcrxarray, $
          itloop1occptrarray, itloop1occptaltspharray, itloop1occptaltellarray, itloop1occptaltasharray, itloop1occptaltkosarray, $
          itloop1occptlatspharray, itloop1occptlatellarray, itloop1occptlatasharray, itloop1occptlatkosarray
      endif
      
      if itloop eq 2 then begin
        itloop2occptradiusarray = occptradiusarray
        itloop2occptlatarray = occptlatarray
        itloop2occptlonarray = occptlonarray
        itloop2occptszaarray = occptszaarray
        itloop2occptlstarray = occptlstarray
        itloop2sepanglearray = sepanglearray
        itloop2epsanglearray = epsanglearray
        itloop2ettxarray = ettxarray
        itloop2etoccptarray = etoccptarray
        itloop2etrxarray = etrxarray
        itloop2utctxarray = utctxarray
        itloop2utcoccptarray = utcoccptarray
        itloop2utcrxarray = utcrxarray
        itloop2occptaltspharray = occptaltspharray
        itloop2occptaltellarray = occptaltellarray
        itloop2occptaltasharray = occptaltasharray
        itloop2occptaltkosarray = occptaltkosarray
        itloop2occptlatspharray = occptlatspharray
        itloop2occptlatellarray = occptlatellarray
        itloop2occptlatasharray = occptlatasharray
        itloop2occptlatkosarray = occptlatkosarray
        itloop2extraoccptradiusarray = extraoccptradiusarray
        itloop2extraoccptaltspharray = extraoccptaltspharray
        itloop2extraoccptaltellarray = extraoccptaltellarray
        itloop2extraoccptaltasharray = extraoccptaltasharray
        itloop2extraoccptaltkosarray = extraoccptaltkosarray
        itloop2extraoccptlatspharray = extraoccptlatspharray
        itloop2extraoccptlatellarray = extraoccptlatellarray
        itloop2extraoccptlatasharray = extraoccptlatasharray
        itloop2extraoccptlatkosarray = extraoccptlatkosarray
        
        save, filename=path+'simpleocc_output_itloop2.sav', itloop2occptradiusarray, itloop2occptlatarray, itloop2occptlonarray, $
          itloop2occptszaarray, itloop2occptlstarray, itloop2sepanglearray, itloop2epsanglearray, itloop2ettxarray, itloop2etoccptarray, $
          itloop2etrxarray, itloop2utctxarray, itloop2utcoccptarray, itloop2utcrxarray, itloop2occptaltspharray, itloop2occptaltellarray, $
          itloop2occptaltasharray, itloop2occptaltkosarray, itloop2occptlatspharray, itloop2occptlatellarray, $
          itloop2occptlatasharray, itloop2occptlatkosarray, itloop2extraoccptradiusarray, itloop2extraoccptaltspharray, $
          itloop2extraoccptaltellarray, itloop2extraoccptaltasharray, itloop2extraoccptaltkosarray, itloop2extraoccptlatspharray, $
          itloop2extraoccptlatellarray, itloop2extraoccptlatasharray, itloop2extraoccptlatkosarray
      endif
      
      print, 'Finished simpleocc processing for '+rsrfilename
      print, ''
      if itloop eq 2 then begin
        ;; RAW MATERIALS FOR THE MATRIX INVERSION
        ;Saves to rsr path. 
        save, filename=path + 'pathlengthr.sav',         inpathlength2d,       outpathlength2d
        save, filename=path + 'pathlengthsph.sav', altsphinpathlength2d, altsphoutpathlength2d
        save, filename=path + 'pathlengthell.sav', altellinpathlength2d, altelloutpathlength2d
        save, filename=path + 'pathlengthash.sav', altashinpathlength2d, altashoutpathlength2d
        save, filename=path + 'pathlengthkos.sav', altkosinpathlength2d, altkosoutpathlength2d
      endif
    endelse ;Ends analyze_rsr conditional 
  endfor ;Ends loop over particular cors file
endfor ;End cors loop 

print, 'Finished processing cors ', allcors, ', itloop = ',itloop
codeendtime = systime(/seconds)
runtime = codeendtime-codestarttime
print, 'Run time (minutes): '
print, runtime/60.,  format='(d20.2)'

stop
end
