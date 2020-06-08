PRO pat_paul_comp_saturn
  rootpath = 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\'
  body = 'Saturn'
  corsnumber='107'
  allrsrpaths = file_search(rootpath+'Output\'+body+'\cors_0'+corsnumber+'\*')
  ;allrsrpaths = allrsrpaths[1]
  ;stop
  allrsrpaths = allrsrpaths[1]

  for i = 0, n_elements(allrsrpaths)-1 do begin
    ans = strsplit(allrsrpaths[i],'\',/extract)
    file = ans[10]
    print,file
    RESTORE,'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Output\'+body+'\cors_0'+corsnumber+'\'+file+'\full_output\simpleocc_output.sav'
    paul_occptradiusarray = occptradiusarray
    paul_utcrxarray    = real_ertutc
    ;RESTORE,'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Output\'+body+'\cors_0'+corsnumber+'\'+file+'\full_output\simpleocc_output_pat_old.sav'
    RESTORE,'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Output\'+body+'\cors_0'+corsnumber+'\'+file+'\full_output\simpleocc_output_pat_newer.sav'

    pat_occptradiusarray = occptradiusarray
    pat_utcrxarray    = real_ertutc
    
    print,'Max. occpt    diff: ', max(abs(pat_occptradiusarray-paul_occptradiusarray[6000:21150]))
    print,''
    print, 'Pat start/end times:  ',pat_utcrxarray[0],', ',pat_utcrxarray[-1]
    print, 'Paul start/end times: ',paul_utcrxarray[0],', ',paul_utcrxarray[-1]

    wait,1.5
  endfor
  print,'Done!'
  stop
END