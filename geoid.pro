function geoid,rref,rg,lats,req,om0,J2,J4,J6,J8,Mp,lref=lref

  ;This function generates a reference geoid for uniform rotation
  ;===================================================================
  ;rtest - initial altitude guesses for each latitude point
  ;lats  - array of latitudes (deg)
  ;req   - equatorial reference radius for J coeffs (lower atmosphere)
  ;rref  - fitted equatorial radius
  ;om0   - angular rotation rate
  ;J2    - zonal gravity harmonic
  ;Mp    - planet mass
  ;lref  - possible latitude reference (deg)
  
  if(keyword_set(lref)) then $
  Uref = grav_pot(Mp,rref,lref*!DtoR,J2,req,om0,J4=J4,J6=J6,J8=J8) $
  else $
  Uref = grav_pot(Mp,rref,0.d0,J2,req,om0,J4=J4,J6=J6,J8=J8)
  outs = lats*0.d0
  
  for n = 0,n_elements(lats)-1 do begin
    rtest = rg(n)
    dr = rtest
    while ( (abs(dr)/rtest) ge 1.d-14) do begin 
      U1 = grav_pot(Mp,rtest,lats(n)*!DtoR,J2,req,om0,J4=J4,J6=J6,J8=J8)
      g1 = grav_pot(Mp,rtest,lats(n)*!DtoR,J2,req,om0,J4=J4,J6=J6,J8=J8,/gradU)
      dU = Uref-U1
      dr = -dU/g1
      rtest = rtest+dr
    endwhile
    outs(n) = rtest
  endfor
  return,outs

end