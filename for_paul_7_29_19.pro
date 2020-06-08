pro for_paul_7_29_19
restore, 'for_paul_7_29_19.sav'
for i = 0,n_elements(full_thing)-1 do begin
  ;print, occtimesecsub[i] - 57835.500
;print, occtimesecsub[i] - 57835.500;, occfs_save[i] -    2.2983678000D+009, occfx_save[i] - 8.4273478000D+009, full_thing[i], format='(10E30)'
;print,format='(10E30)', occfs_save[i] -    2.2983678000D+009
;print,format='(10E30)', occfx_save[i] - 8.4273478000D+009
;print,format='(10E30)', full_thing[i]

print, occtimesecsub[i], occfs_save[i], occfx_save[i], occfs_save[i]-(3d/11d)*occfx_save[i], format='(10E30)'
endfor

write_csv, 's14n_sx14_new.csv', occtimesecsub, occfs_save, occfx_save, occfs_save-(3d/11d)*occfx_save,header=['occtimesec','fs','fx','fs-3/11*fx']
stop
end