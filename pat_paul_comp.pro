PRO pat_paul_comp
;cors = 'cors_0139\'
;file = 's19tibi2006077_2337nnnx26rv.1n2'
rootpath = 'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\'
body = 'Titan'
corsnumber='538'
allrsrpaths = file_search(rootpath+'Output\'+body+'\cors_0'+corsnumber+'\*')
for i = 0, n_elements(allrsrpaths)-1 do begin
  ans = strsplit(allrsrpaths[i],'\',/extract)
  file = ans[10]
  print,file
  RESTORE,'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Output\Titan\cors_0'+corsnumber+'\'+file+'\full_output\simpleocc_output.sav'
  paul_occptradiusarray = occptradiusarray
  paul_occptlatarray = occptlatarray
  paul_occptlonarray = occptlonarray
  paul_occptszaarray = occptszaarray
  paul_occptlstarray = occptlstarray
  paul_sepanglearray = sepanglearray
  paul_epsanglearray = epsanglearray
  paul_ettxarray     = ettxarray
  paul_etoccptarray  = etoccptarray
  paul_etrxarray     = etrxarray
  paul_utctxarray    = utctxarray
  paul_utcoccptarray = utcoccptarray
  paul_utcrxarray    = utcrxarray
  RESTORE,'C:\Users\tambu\Documents\BU\SaturnEDP\PD_MBP_files\Output\Titan\cors_0'+corsnumber+'\'+file+'\full_output\simpleocc_output_pat.sav'
  pat_occptradiusarray = occptradiusarray
  pat_occptlatarray = occptlatarray
  pat_occptlonarray = occptlonarray
  pat_occptszaarray = occptszaarray
  pat_occptlstarray = occptlstarray
  pat_sepanglearray = sepanglearray
  pat_epsanglearray = epsanglearray
  pat_ettxarray     = ettxarray
  pat_etoccptarray  = etoccptarray
  pat_etrxarray     = etrxarray
  pat_utctxarray    = utctxarray
  pat_utcoccptarray = utcoccptarray
  pat_utcrxarray    = utcrxarray
  
  print,'Max. occpt    diff: ', max(abs(pat_occptradiusarray-paul_occptradiusarray))
  print,'Max. occptlat diff: ', max(abs(pat_occptlatarray-paul_occptlatarray))
  print,'Max. occptlon diff: ', max(abs(pat_occptlonarray-paul_occptlonarray))
  print,'Max. occptsza diff: ', max(abs(pat_occptszaarray-paul_occptszaarray))
  print,'Max. occptlst diff: ', max(abs(pat_occptlstarray-paul_occptlstarray))
  print,'Max. sepangle diff: ', max(abs(pat_sepanglearray-paul_sepanglearray))
  print,'Max. epsangle diff: ', max(abs(pat_epsanglearray-paul_epsanglearray))
  print,'Max. ettx     diff: ', max(abs(pat_ettxarray-paul_ettxarray))
  print,'Max. etoccpt  diff: ', max(abs(pat_etoccptarray-paul_etoccptarray))
  print,'Max. etrx     diff: ', max(abs(pat_etrxarray-paul_etrxarray))
  print,''
  wait,1.5
endfor
print,'Done!'
END