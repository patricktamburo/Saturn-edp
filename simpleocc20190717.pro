pro simpleocc20190717

; Paul Withers, 2006.09.13
; Center for Space Physics, Boston University

; Find spacecraft occultations using SPICE

itloop = 2
; 0 - use all times in RSR file, find set of CA radial distances - repeat for both bands
; 1 - use only selected subset of times for which profile will be found, find set of CA altitudes - one band only (X)
; 2 - use only selected subset of times for which profile will be found, find matrix elements - one band only (X)
; run freqwork between 0 and 1
; run freqwork after 2

;spawn, 'date', codestarttime
; Record time program starting running

;rootpath = '/Users/paul/Documents/Boston_University/Research/Cassini_Ionosphere/'
;;rootpath = 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\'
;spicerootpath = '/Volumes/PW-2TB/Cassini_Ionosphere/'
spicerootpath = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SPICE_kernels\'
;genspicepath = spicerootpath+'generic/'
genspicepath = spicerootpath+'generic\'
;spkspicepath = spicerootpath+'spk/'
;;spkspicepath = spicerootpath+'spk\'

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

;cspice_furnsh, '../pck00010.tpc' ; planet rotational states, Paul's machine.
;cspice_furnsh, '../naif0012.tls' ; leap seconds
;cspice_furnsh, '../de430.bsp' ; SPK (trajectory kernel) for planets

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
;loop = 'yes'
;body = 'Saturn'
;;allcors = ['139','140','145','175','184','185','186','259','260','276','277','278','493','494','501','505','536','538','544','545']
;;allcors = ['544','545']
;allcors = ['122']

allcors =0 ; null value to ensure one loop only

for cor = 0, n_elements(allcors)-1 do begin
  ;;  corsnumber = allcors[cor]
    
    ;allrsrpaths = file_search(rootpath+'Output/'+body+'/cors_0'+corsnumber+'/*')
  ;;  allrsrpaths = file_search(rootpath+'Output\'+body+'\cors_0'+corsnumber+'\*')
    ;Loop over all rsr paths
  ;;  for rsrpath = 0, n_elements(allrsrpaths)-1 do begin
      ;path = allrsrpaths[rsrpath]+'/'
  ;;    path = allrsrpaths[rsrpath]+'\'
      
  ;;; XXX HARD-CODED
  ;S07N
  ;cspice_furnsh, '../050606R_SCPSE_05114_05132.bsp' ; S07N
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0106/s10sroi2005123_0230nnns14rd.2b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0107/s10sroi2005123_0230nnnx14rd.2a2/'
  ;thisocc = 'S07N'
  ;carryoverutcstart = '2005 MAY 03 05:46:30.500002'
  ;carryoverutcend = '2005 MAY 03 05:57:20.500002'
  ;norx = 'N' ; 'N' or 'X', need to flip some things
  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SPICE_kernels\spk\050606R_SCPSE_05114_05132.bsp' ;S07N, Pat
  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0106\s10sroi2005123_0230nnns14rd.2b2\' ;S07N S-band
  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0107\s10sroi2005123_0230nnnx14rd.2a2\' ;S07N X-band
  thisocc = 'S07N'
  carryoverutcstart = '2005 MAY 03 05:46:30.500002'
  carryoverutcend = '2005 MAY 03 05:57:20.500002'
  norx = 'N' ; 'N' or 'X', need to flip some things
  
  
  ;S10N
  ;cspice_furnsh, '../050802R_SCPSE_05169_05186.bsp' ; S10N
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0116/s12sroi2005177_1741nnns14rd.2b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0116/s12sroi2005177_1741nnnx14rd.2a2/'
  ;thisocc = 'S10N'
  ;carryoverutcstart = '2005 JUN 26 20:18:56.500003'
  ;carryoverutcend = '2005 JUN 26 20:29:46.500003'
  ;norx = 'N' ; 'N' or 'X', need to flip some things
;  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SPICE_kernels\spk\050802R_SCPSE_05169_05186.bsp' ;S10N, Pat
;  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0116\s12sroi2005177_1741nnns14rd.2b2\' ;S10N S-band
;  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0116\s12sroi2005177_1741nnnx14rd.2a2\' ;S10N X-band
;  thisocc = 'S10N'
;  carryoverutcstart = '2005 JUN 26 20:18:56.500003'
;  carryoverutcend = '2005 JUN 26 20:29:46.500003'
;  norx = 'N' ; 'N' or 'X', need to flip some things

  
  
  ;S56X
  ;cspice_furnsh, '../080327R_SCPSE_07365_08045.bsp' ; S56X
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0213/s36saoc2008015_2150nnns63rd.1b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190610/for_paul_10_06_2019/Output/cors_0213/s36saoc2008015_2150nnnx63rd.1a2/'
  ;thisocc = 'S56X'
  ;carryoverutcstart = '2008 JAN 15 22:52:40.500001'
  ;carryoverutcend = '2008 JAN 15 23:09:20.500001'
  ;norx = 'X' ; 'N' or 'X', need to flip some things
;  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SPICE_kernels\spk\080327R_SCPSE_07365_08045.bsp' ;S56X, Pat
;  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0213\s36saoc2008015_2150nnns63rd.1b2\' ;S56X S-band
;  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0213\s36saoc2008015_2150nnnx63rd.1a2\' ;S56X X-band
;  thisocc = 'S56X'
;  carryoverutcstart = '2008 JAN 15 22:52:40.500001'
;  carryoverutcend = '2008 JAN 15 23:09:20.500001'
;  norx = 'X' ; 'N' or 'X', need to flip some things
;  
  
  ;S51X
  ;cspice_furnsh, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Spice Files/' + $
  ; '080117R_SCPSE_07262_07309.bsp' ; S51X
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0203/' + $
  ; 's34saoe2007297_0745nnns63rd.1b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0204/' + $
  ; 's34saoe2007297_0745nnnx63rd.1a2/'
  ;thisocc = 'S51X'
  ;carryoverutcstart = '2007 OCT 24 08:04:52.500001'
  ;carryoverutcend = '2007 OCT 24 08:21:32.500001'
  ;norx = 'X' ; 'N' or 'X', need to flip some things
;  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\SaturnEDP_backup_6_17_19\PD_MBP_files\Programs\checks\ed\spice_data\080117R_SCPSE_07262_07309.bsp' ;S51X, Pat
;  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0203\s34saoe2007297_0745nnns63rd.1b2\' ;S51X S-band
;  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0204\s34saoe2007297_0745nnnx63rd.1a2\' ;S51X X-band
;  thisocc = 'S51X'
;  carryoverutcstart = '2007 OCT 24 08:04:52.500001'
;  carryoverutcend = '2007 OCT 24 08:21:32.500001'
;  norx = 'X' ; 'N' or 'X', need to flip some things
  
  ;S54X
  ;cspice_furnsh, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Spice Files/' + $
  ; '080307R_SCPSE_07345_07365.bsp' ; S54X
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0207/' + $
  ; 's36saoe2007353_0524nnns63rd.1b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0207/' + $
  ; 's36saoe2007353_0524nnnx63rd.1a2/'
  ;thisocc = 'S54X'
  ;carryoverutcstart = '2007 DEC 19 05:51:34.500002'
  ;carryoverutcend = '2007 DEC 19 05:59:54.500001'
  ;norx = 'X' ; 'N' or 'X', need to flip some things
;  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_25_06_2019\Spice Files\080307R_SCPSE_07345_07365.bsp' ;S54X, Pat
;  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0207\s36saoe2007353_0524nnns63rd.1b2\' ;S54X S-band
;  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0207\s36saoe2007353_0524nnnx63rd.1a2\' ;S54X X-band
;  thisocc = 'S54X'
;  carryoverutcstart = '2007 DEC 19 05:51:34.500002'
;  carryoverutcend = '2007 DEC 19 05:59:54.500001'
;  norx = 'X' ; 'N' or 'X', need to flip some things
  
  ;S75X
  ;cspice_furnsh, '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Spice Files/' + $
  ; '080819R_SCPSE_08141_08206.bsp' ; S75X
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0245/' + $
  ; 's42saoe2008189_1040nnns63rd.1b2/'
  ;path = '/deimos/deimos_backup/olddeimos00/spice/icy/cassiniocc/tamburo/20190625/for_paul_25_06_2019/Output/cors_0245/' + $
  ; 's42saoe2008189_1040nnnx63rd.1a2/'
;  cspice_furnsh, 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\for_paul_25_06_2019\Spice Files\080819R_SCPSE_08141_08206.bsp' ;S75X, Pat
;  ;path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0245\s42saoe2008189_1040nnns63rd.1b2\' ;S75X S-band
;  path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PD_MBP_files\Output\Saturn\cors_0245\s42saoe2008189_1040nnnx63rd.1a2\' ;S75X X-band
;  thisocc = 'S75X'
;  carryoverutcstart = '2008 JUL 07 11:14:48.499999'
;  carryoverutcend = '2008 JUL 07 11:23:08.499999'
;  norx = 'X' ; 'N' or 'X', need to flip some things
  
  
  body = 'Saturn'
  
  if itloop eq 1 then begin
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
  endif
  
  if itloop eq 2 then begin
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
  endif
  
  ;Determine which body from the path. 's' for Saturn occultations, 't' for Titan occultations
  ;if strmatch(path,'*Titan*') then whichbody = 't'
  ;if strmatch(path,'*Saturn*') then whichbody = 's'
  if body eq 'Saturn' then whichbody = 's'
  if body eq 'Titan' then whichbody = 't'
  ;if body eq
      
  ;Extract the RSR filename from the path
  ;splitpath = strsplit(path,f'/',/extract)
  ;;    splitpath = strsplit(path,'\',/extract)
  
  ;;    rsrfilename = splitpath[-1]
      
  ;Check to see if this occultation already has a simpleocc_output.sav file. If so, skip it.
  ;if (size(file_search(path+'full_output/simpleocc_output*')))[0] eq 1 then begin
  ;  print, rsrfilename+' already has a simpleocc output file. Skipping.'
  ;  if loop eq 'no' then stop
  ;  if loop eq 'yes' then stop  ;replace stop with continue here to loop
  ;endif
      
      
  ;Determine the appropropriate .bsp spice files to load for this RSR file.
  ;Once found, that file or those files into spice
  ;openr, lun, rootpath+'Programs//deimos/withers/homedir/latex/felicidust2019RSR_SPICE_connection.txt', /get_lun
  ;;    openr, lun, rootpath+'Programs\RSR_SPICE_connection.txt', /get_lun
  
  ;;    line = ''
  ;;    while not eof(lun) do begin
  ;;      readf, lun, line
  ;;      ;Check if this is the line for this RSR file
  ;;      if strmatch(line,rsrfilename+'*') then begin
  ;;        ;Extract the .bsp filename and furnish spice with them
  ;;        splitline = strsplit(line,',',/extract)
  ;;        spicespkfiles = splitline[1:*]
  ;;        for spk=0,n_elements(spicespkfiles)-1 do begin
  ;;          ;Ingest spacecraft trajectory kernel into memory
  ;;          cspice_furnsh, spkspicepath+spicespkfiles[spk]
  ;;        endfor
  ;;        break
  ;;      endif
  ;;    endwhile
  ;;    free_lun, lun
  
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
  ;;    restore, path+'full_output\output.sav'
  fakesfduyearoutstart = sfduyearout[0]
  fakesfdudoyoutstart = sfdudoyout[0]
  fakesfdusecoutstart = sfdusecout[0]
  ;;;fakesfduyearoutend = sfduyearout[-1]
  ;;;fakesfdudoyoutend = sfdudoyout[-1]
  ;;;fakesfdusecoutend = sfdusecout[-1]
  ;;;
  fakesfduyearoutend = sfduyearout[n_elements(sfduyearout)-1]
  fakesfdudoyoutend = sfdudoyout[n_elements(sfdudoyout)-1]
  fakesfdusecoutend = sfdusecout[n_elements(sfdusecout)-1]
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
      
  ;fakeetstart = fakeetstart - 100d * 60d ; allow 100 minutes for possible light-time issues and UTC/ET differences
  ;fakeetend = fakeetend + 100d * 60d ; allow 100 minutes for possible light-time issues and UTC/ET differences
  
  etstart = fakeetstart
  etend = fakeetend
  
  nsteps = ceil( (etend- etstart)/dt) ; Number of steps
  
  if itloop eq 1 or itloop eq 2 then begin ; restrict to smaller range of times
    ;;; XXX HARD-CODED
    ;; VALUES OBTAINED FROM FREQWORK AND HARD-CODED HERE (AND ABOVE)
    cspice_str2et, carryoverutcstart, etstart
    cspice_str2et, carryoverutcend, etend
    
    ;cspice_str2et, '2005 MAY 03 05:46:30.500002', etstart ; S07N
    ;cspice_str2et, '2005 MAY 03 05:57:20.500002', etend
    ;cspice_str2et, '2005 JUN 26 20:18:56.500003', etstart ; S10N
    ;cspice_str2et, '2005 JUN 26 20:29:46.500003', etend
    ;cspice_str2et, '2008 JAN 15 22:52:40.500001', etstart ; S56X
    ;cspice_str2et, '2008 JAN 15 23:09:20.500001', etend
    
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
  
  if norx eq 'N' then begin
  endif
  
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
  
  ;;;
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
    ettx = etrx - txrxltime
    ; use this time for events in the Saturn system, so Cassini, Saturn, Titan
    
    ;;;
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
  
      ;  a = 60268d
      ;  b = 54364d
      ;thisrrefsph = 60268d
      cspice_tipbod, frame, targetint, etoccpt, tipm
  
      k=0L 
      while k lt n_elements(xx) do begin
        kocd = oldkocd + xx[k]
        occptssb = sctssb + kocd * (earthssb - sctssb)
        distarr[k] = norm(occptssb - targetssb)
        cspice_mxv, tipm, occptssb - targetssb, bfoccpt
        cspice_reclat, bfoccpt, occptradius, occptlongitude, occptlatitude
        templatdegarr[k] = occptlatitude * 180./!pi
        
        ;rrefspharr[k] = thisrrefsph
        ;rrefellarr[k] = a*b/sqrt((b*cos(templatdegarr[k]*!pi/180.))^2.+(a*sin(templatdegarr[k]*!pi/180.))^2.)
        ;rrefasharr[k] = sat_anderson_schubert_1bar(templatdegarr[k])
        ;rrefkosarr[k] = saturn_1bar(90*!pi/180.-!pi/180*templatdegarr[k])
        
        k++
      endwhile
  
      rrefspharr = thisrrefsph + templatdegarr*0d
      rrefellarr = interpol(dummyrrefell,dummylatdeg,templatdegarr)
      rrefasharr = interpol(dummyrrefash,dummylatdeg,templatdegarr)
      rrefkosarr = interpol(dummyrrefkos,dummylatdeg,templatdegarr)
      ; hopefully faster than direct call to subroutines
      
      ;rrefspharr = thisrrefsph + templatdegarr*0d
      ;rrefellarr = a*b/sqrt((b*cos(templatdegarr*!pi/180.))^2.+(a*sin(templatdegarr*!pi/180.))^2.)
      ;rrefasharr = sat_anderson_schubert_1bar(templatdegarr)
      ;rrefkosarr = saturn_1bar(90*!pi/180.-!pi/180*templatdegarr)
  
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
  
            ;inpathlength[j] = norm( (inflankk[j] - cakr) * (earthssb - sctssb))
            ;altsphinpathlength[j] = norm( (altsphinflankk[j] - caksph) * (earthssb - sctssb))
            ;altellinpathlength[j] = norm( (altellinflankk[j] - cakell) * (earthssb - sctssb))
            ;altashinpathlength[j] = norm( (altashinflankk[j] - cakash) * (earthssb - sctssb))
            ;altkosinpathlength[j] = norm( (altkosinflankk[j] - cakkos) * (earthssb - sctssb))
            ;outpathlength[j] = norm( (outflankk[j] - cakr) * (earthssb - sctssb))
            ;altsphoutpathlength[j] = norm( (altsphoutflankk[j] - caksph) * (earthssb - sctssb))
            ;altelloutpathlength[j] = norm( (altelloutflankk[j] - cakell) * (earthssb - sctssb))
            ;altashoutpathlength[j] = norm( (altashoutflankk[j] - cakash) * (earthssb - sctssb))
            ;altkosoutpathlength[j] = norm( (altkosoutflankk[j] - cakkos) * (earthssb - sctssb))
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
  
            ;inpathlength[j] = norm( (inflankk[j] - inflankk[j-1]) * (earthssb - sctssb))
            ;altsphinpathlength[j] = norm( (altsphinflankk[j] - altsphinflankk[j-1]) * (earthssb - sctssb))
            ;altellinpathlength[j] = norm( (altellinflankk[j] - altellinflankk[j-1]) * (earthssb - sctssb))
            ;altashinpathlength[j] = norm( (altashinflankk[j] - altashinflankk[j-1]) * (earthssb - sctssb))
            ;altkosinpathlength[j] = norm( (altkosinflankk[j] - altkosinflankk[j-1]) * (earthssb - sctssb))
            ;outpathlength[j] = norm( (outflankk[j] - outflankk[j-1]) * (earthssb - sctssb))
            ;altsphoutpathlength[j] = norm( (altsphoutflankk[j] - altsphoutflankk[j-1]) * (earthssb - sctssb))
            ;altelloutpathlength[j] = norm( (altelloutflankk[j] - altelloutflankk[j-1]) * (earthssb - sctssb))
            ;altashoutpathlength[j] = norm( (altashoutflankk[j] - altashoutflankk[j-1]) * (earthssb - sctssb))
            ;altkosoutpathlength[j] = norm( (altkosoutflankk[j] - altkosoutflankk[j-1]) * (earthssb - sctssb))
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
  
    cspice_spkgps, earthint, etrx, frame, ssbint, earthssb, ltimeearthssb
    ; Position of Earth at time etrx relative to solar system barycentre and in specified frame
        
    ;;;cspice_spkgps, targetint, ettx, frame, ssbint, targetssb, ltimetargetssb
    ; Position of target body at time ettx relative to solar system barycentre and in specified frame
    cspice_spkgps, targetint, etoccpt, frame, ssbint, targetssb, ltimetargetssb
        
    ;Testing against older version of the code. 
    ; restore, 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Programs\time_test.sav'
    ; etoccpt = ettx
    ; cspice_spkgps, targetint, etoccpt, frame, ssbint, targetssb, ltimetargetssb

    ;stop
    ; Position of target body at time etoccpt relative to solar system barycentre and in specified frame
    ;;stop
    cspice_spkgps, sctint, ettx, frame, ssbint, sctssb, ltimemgsssb
    ; Position of spacecraft at time ettx relative to solar system barycentre and in specified frame
    
    earthsct = earthssb - sctssb
    ; Position of Earth relative to spacecraft in working frame
        
    cspice_spkgps, sunint, etrx, frame, ssbint, sunssb, ltimesunssb
    ; Position of Sun at time etrx relative to solar system barycentre and in specified frame
    earthsun = earthssb - sunssb
    ; Position of Earth relative to Sun, not ssb, in working frame
    sepangle = cspice_vsep( (-1d)*earthsun, (-1d)*earthsct ) * 180./!pi
    ; Sun-Earth-Spacecraft angle, redundant minus signs ensure that directions of vectors are "Sun relative to Earth" and "sct relative to Earth"
    epsangle = cspice_vsep( earthsct, sunssb-sctssb) * 180./!pi
    ; Earth-Spacecraft-Sun angle
        
    cspice_nplnpt, earthssb, earthsct, targetssb, pnear, distoccpttarget
    ; pnear is nearest point on Earth-spacecraft line to target body
    ; Existence of occultations will depend a lot on pnear
    ; I loosely call position indicated by pnear the "occultation point"
    ; See help pages for frame in which pnear is stated
    ; distoccpttarget is distance between pnear and target body, km

    occpttarget = pnear - targetssb
    ; Position of pnear relative to target body in working frame
        
    ;cspice_tipbod, frame, targetint, ettx, tipm
    cspice_tipbod, frame, targetint, etoccpt, tipm
    ; tipm is transformation matrix between specified inertial frame (frame, probably J2000)
    ; and body-fixed frame of target body
    ; See help pages for more on relevant frames
        
    cspice_mxv, tipm, occpttarget, bfoccpt
    ; Occultation point position converted from one frame (occpttarget)
    ; to body-fixed frame of target body (bfoccpt)
    cspice_reclat, bfoccpt, occptradius, occptlongitude, occptlatitude
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
        
    ;;;ls = cspice_lspcn(targetstr, ettx, abcorr)
    ls = cspice_lspcn(targetstr, etoccpt, abcorr)
    ; ls (longitude of the sun, planetocentric, radians) is often used as a measure of seasons
    ; Especially common on Mars, see help pages for details
    
    ;;;cspice_et2lst, ettx, targetint, !pi/180.*occptlongitude, 'planetocentric', $
    cspice_et2lst, etoccpt, targetint, !pi/180.*occptlongitude, 'planetocentric', $
     hr, mn, sc, time, ampm
    lst = hr + (mn + sc/double(60.))/double(60.)
    ; lst is local solar time (0-24 hours) of occultation point at time et
    
    ;;;cspice_et2lst, ettx, targetint, !pi/180.*sctlongitude, 'planetocentric', $
    cspice_et2lst, etoccpt, targetint, !pi/180.*sctlongitude, 'planetocentric', $
     hr, mn, sc, time, ampm
    sctlst = hr + (mn + sc/double(60.))/double(60.)
    ; sctlst is local solar time (0-24 hours) of spacecraft at time et
        
    ;sza = cspice_vsep(double(-1.)*targetssb, scttarget)*double(180.)/!pi
    ; Rejected method of finding spacecraft's solar zenith angle (degrees)
    ; Angle between spacecraft-target vector and solar system barycentre-target vector
    ; Sun is not quite solar system barycentre, but it is close enough
    
    ;;;cspice_subsol, method, targetstr, ettx, abcorr, earthstr, targetsubsolpt
    cspice_subsol, method, targetstr, etoccpt, abcorr, earthstr, targetsubsolpt
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
        
    occptlatarray[i] = occptlatitude     
    occptlonarray[i] = occptlongitude
    occptszaarray[i] = sza
    occptlstarray[i] = lst
    sepanglearray[i] = sepangle
    epsanglearray[i] = epsangle
    ettxarray[i] = ettx
    etoccptarray[i] = etoccpt
    etrxarray[i] = etrx
    ;Store the output data
        
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
  
  ;Convert all ET times to UTC times and save
  cspice_et2utc, ettxarray, 'C', 6, utctxarray
  cspice_et2utc, etoccptarray, 'C', 6, utcoccptarray
  cspice_et2utc, etrxarray, 'C', 6, utcrxarray
      
  ; save, filename=path+'tests07routput.sav', real_ertutc, occptradiusarray
      
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
    ;save, filename=path+'simpleocc_output_itloop0.sav', occptradiusarray, occptlatarray, occptlonarray, occptszaarray, occptlstarray, $
    ; sepanglearray, epsanglearray, ettxarray, etoccptarray, etrxarray, utctxarray, utcoccptarray, utcrxarray
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
    ;save, filename=path+'simpleocc_output_itloop1.sav', occptradiusarray, occptlatarray, occptlonarray, occptszaarray, occptlstarray, $
    ; sepanglearray, epsanglearray, ettxarray, etoccptarray, etrxarray, utctxarray, utcoccptarray, utcrxarray, $
    ; itloop1occptrarray, occptaltspharray, occptaltellarray, occptaltasharray, occptaltkosarray
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
    ;save, filename=path+'simpleocc_output_itloop2.sav', occptradiusarray, occptlatarray, $
    ; occptlonarray, occptszaarray, occptlstarray, sepanglearray, epsanglearray, $
    ; ettxarray, etoccptarray, etrxarray, utctxarray, utcoccptarray, utcrxarray, $
    ; occptaltspharray, occptaltellarray, occptaltasharray, occptaltkosarray, $
    ; extraoccptradiusarray, extraoccptaltspharray, extraoccptaltellarray, extraoccptaltasharray, extraoccptaltkosarray 
  endif
  
  ;       save, filename=path+'full_output/simpleocc_output.sav', occptradiusarray, occptlatarray, occptlonarray, occptszaarray, occptlstarray, sepanglearray, epsanglearray, $
      ; ettxarray, etoccptarray, etrxarray, utctxarray, utcoccptarray, utcrxarray
  
  ;    save, filename=path+'simpleocc_output_test.sav', occptradiusarray, occptlatarray, $
  ;;    save, filename=path+'full_output\simpleocc_output_pat.sav', occptradiusarray, occptlatarray, $
  ;     occptlonarray, occptszaarray, occptlstarray, sepanglearray, epsanglearray, $
  ;     ettxarray, etoccptarray, etrxarray, utctxarray, utcoccptarray, utcrxarray, $
  ; occptaltspharray, occptaltellarray, occptaltasharray, occptaltkosarray  
  ;;    print, 'Finishing simpleocc processing for '+rsrfilename
  print, 'Finishing simpleocc processing for '
       
  if itloop eq 2 then begin
    ;; RAW MATERIALS FOR THE MATRIX INVERSION
    save, filename=thisocc + 'pathlengthr.sav', inpathlength2d, outpathlength2d
    save, filename=thisocc + 'pathlengthsph.sav', altsphinpathlength2d, altsphoutpathlength2d
    save, filename=thisocc + 'pathlengthell.sav', altellinpathlength2d, altelloutpathlength2d
    save, filename=thisocc + 'pathlengthash.sav', altashinpathlength2d, altashoutpathlength2d
    save, filename=thisocc + 'pathlengthkos.sav', altkosinpathlength2d, altkosoutpathlength2d
  endif
  
  ;save, filename = 'coolstuff.sav', occptradiusarray, occptlatarray, occptlonarray, occptszaarray, occptlstarray, sepanglearray, epsanglearray, $
  ; ettxarray, etoccptarray, etrxarray
  
  ;This endfor ends the loop over the list of cors file.
  ;If running only a single file manually, comment this out
endfor ;End cors loop 

;spawn, 'date', codeendtime
;print, codestarttime
;print, codeendtime

stop
end
