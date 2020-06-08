function sat_anderson_schubert_1bar,lat,r0=r0,rot=rot

  ;Here we reproduce the best fit reference geoid of 
  ;Anderson and Schubert (2007)
  ;==================================================
  ;Argument: latitude in deg
  ;Optional arguments: 
  ;r0  - alternative equatorial radius (m)
  ;rot - alternative rotation rate (s-1)
  
  G   = 6.674215d-11
  GMS = 37931207.7d9
  MS  = GMS/G
  Ref = 60330.d3
  Req = 60268.d3
  J2  = 16290.73d-6
  J4  = -935.5d-6
  J6  = 85.3d-6
  J8  = -10.d-6
  J10 = 2.d-6
  om0 = 1.655430196d-4
  
  rg = lat*0.d0
  rg(*) = Req
  
  if(keyword_set(r0)) then Req  = r0
  if(keyword_set(rot)) then om0 = rot 
  
  shape = geoid(Req,rg,lat,Ref,om0,J2,J4,J6,J8,MS)*1.d-3
  
  return,shape

end