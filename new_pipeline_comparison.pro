PRO new_pipeline_comparison
;Compares output from Paul W.'s code with Pat T.'s.
occ = 'S75X'
;fname = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PT_files\pat_paul_comparisons\'+'output' + occ + '_pat_new.sav'
;save,ALTSPHMATRIXNELEC, ALTELLMATRIXNELEC, ALTASHMATRIXNELEC, ALTKOSMATRIXNELEC, ITLOOP2EXTRAOCCPTALTSPHARRAY, ITLOOP2EXTRAOCCPTALTELLARRAY, ITLOOP2EXTRAOCCPTALTASHARRAY, ITLOOP2EXTRAOCCPTALTKOSARRAY, filename=fname 
root_path = 'C:\Users\tambu\Documents\BU\Science Projects\SaturnEDP\PT_files\pat_paul_comparisons\'
occ_path = root_path + occ + '\'
paul_file = occ_path + 'output' + occ + '_pat_new.sav'
pat_file  = occ_path + 'output' + occ + '_pat.sav'

restore,paul_file,/verbose
paul_elec_sph = ALTSPHMATRIXNELEC
paul_elec_ell = ALTELLMATRIXNELEC
paul_elec_ash = ALTASHMATRIXNELEC
paul_elec_kos = ALTKOSMATRIXNELEC
paul_alt_sph  = ITLOOP2EXTRAOCCPTALTSPHARRAY
paul_alt_ell  = ITLOOP2EXTRAOCCPTALTELLARRAY
paul_alt_ash  = ITLOOP2EXTRAOCCPTALTASHARRAY
paul_alt_kos  = ITLOOP2EXTRAOCCPTALTKOSARRAY

restore,pat_file,/verbose
pat_elec_sph = ALTSPHMATRIXNELEC
pat_elec_ell = ALTELLMATRIXNELEC
pat_elec_ash = ALTASHMATRIXNELEC
pat_elec_kos = ALTKOSMATRIXNELEC
pat_alt_sph  = ITLOOP2EXTRAOCCPTALTSPHARRAY
pat_alt_ell  = ITLOOP2EXTRAOCCPTALTELLARRAY
pat_alt_ash  = ITLOOP2EXTRAOCCPTALTASHARRAY
pat_alt_kos  = ITLOOP2EXTRAOCCPTALTKOSARRAY

print,'Max. alt. diff., spherical: ',max(abs(pat_alt_sph-paul_alt_sph))
print,'Max. alt. diff., ellipsoidal: ',max(abs(pat_alt_ell-paul_alt_ell))
print,'Max. alt. diff., A + S: ',max(abs(pat_alt_ash-paul_alt_ash))
print,'Max. alt. diff., Koskinen: ',max(abs(pat_alt_kos-paul_alt_kos))
print,'Max. e- diff., spherical: ',max(abs(pat_elec_sph-paul_elec_sph))
print,'Max. e- diff., ellipsoidal: ',max(abs(pat_elec_ell-paul_elec_ell))
print,'Max. e- diff., A + S: ',max(abs(pat_elec_ash-paul_elec_ash))
print,'Max. e- diff., Koskinen: ',max(abs(pat_elec_kos-paul_elec_kos))

w = window(dimensions=[1800,1800])
p = plot(paul_elec_sph,reverse(paul_alt_sph),/current,thick=2,layout=[2,2,1],title='Spherical',name='Withers',margin=[0.23,0.2,0.1,0.2],$
  xtitle='e$^-$ m$^{-3}$',ytitle='Altitude (km)')
p2 = plot(pat_elec_sph,reverse(pat_alt_sph),/overplot,thick=2,linestyle='--',color='r',name='Tamburo')
l = legend(target=[p,p2],position=[0.6,0.52],shadow=0)

p = plot(paul_elec_ell,reverse(paul_alt_ell),/current,thick=2,layout=[2,2,2],title='Ellipsoidal',name='Withers',margin=[0.23,0.2,0.1,0.2],$
  xtitle='e$^-$ m$^{-3}$',ytitle='Altitude (km)')
p2 = plot(pat_elec_ell,reverse(pat_alt_ell),/overplot,thick=2,linestyle='--',color='r',name='Tamburo')

p = plot(paul_elec_ash,reverse(paul_alt_ash),/current,thick=2,layout=[2,2,3],title='Anderson Schubert',name='Withers',margin=[0.23,0.2,0.1,0.2],$
  xtitle='e$^-$ m$^{-3}$',ytitle='Altitude (km)')
p2 = plot(pat_elec_ash,reverse(pat_alt_ash),/overplot,thick=2,linestyle='--',color='r',name='Tamburo')

p = plot(paul_elec_kos,reverse(paul_alt_kos),/current,thick=2,layout=[2,2,4],title='Koskinen',name='Withers',margin=[0.23,0.2,0.1,0.2],$
  xtitle='e$^-$ m$^{-3}$',ytitle='Altitude (km)')
p2 = plot(pat_elec_kos,reverse(pat_alt_kos),/overplot,thick=2,linestyle='--',color='r',name='Tamburo')

t = text(0.5,0.96,occ,font_size=16)

STOP
END