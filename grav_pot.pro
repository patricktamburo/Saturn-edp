function grav_pot,Mp,r,lat,J2,req,omega,J4=J4,J6=J6,J8=J8,gradU=gradU
 
  ;Note - lat must be in radians!!!
  
  G = 6.674215d-11
  
  ;Calculate Legendre polynomials
  const = G*Mp/r 
  mu    = sin(lat)
  P2    = (3.d0*(mu^2.d0)-1.d0)/2.d0
  if(keyword_set(J4)) then P4 = (35.d0*(mu^4.d0)-30.d0*(mu^2.d0)+3.d0)/8.d0 $
  else P4 = 0.d0
  if(keyword_set(J6)) then P6 = (231.d0*(mu^6.d0)-315.d0*(mu^4.d0) $ 
                              + 105.d0*(mu^2.d0)-5.d0)/16.d0 $
  else P6 = 0.d0
  if(keyword_set(J8)) then P8 = (6435.d0*(mu^8.d0)-12012.d0*(mu^6.d0) $
                              + 6930.d0*(mu^4.d0)-1260.d0*(mu^2.d0)+35.d0)/128.d0 $
  else P8 = 0.d0
  
  if(keyword_set(gradU)) then begin
    ;Calculate perpendicular gravity
    first  = 1.d0
    second = -3.d0*J2*((req/r)^2.d0)*P2
    if(P4 ne 0.d0) then third = -5.d0*J4*((req/r)^4.d0)*P4 $
    else third = 0.d0
    if(P6 ne 0.d0) then fourth = -7.d0*J6*((req/r)^6.d0)*P6 $
    else fourth = 0.d0
    if(P8 ne 0.d0) then fifth = -9.d0*J8*((req/r)^8.d0)*P8 $
    else fifth = 0.d0
    rot    = -(cos(lat)^2.d0)*(omega^2.d0)*(r^3.d0)/(G*Mp)
    pot    = const*(first+second+third+fourth+fifth+rot)/r
  endif else begin   
    first  = 1.d0
    second = -J2*((req/r)^2.d0)*P2
    if(P4 ne 0.d0) then third  = -J4*((req/r)^4.d0)*P4 $
    else third = 0.d0
    if(P6 ne 0.d0) then fourth = -J6*((req/r)^6.d0)*P6 $
    else fourth = 0.d0
    if(P8 ne 0.d0) then fifth  = -J8*((req/r)^8.d0)*P8 $
    else fifth = 0.d0 
    rot    = 0.5d0*(cos(lat)^2.d0)*(omega^2.d0)*(r^3.d0)/(G*Mp)
    pot    = const*(first+second+third+fourth+fifth+rot)
  endelse
  
  ;if(keyword_set(gradU)) then begin
  ;  ;Calculate perpendicular gravity
  ;  const  = G*Mp/(r^2.d0)
  ;  first  = 1.d0
  ;  second = -1.5d0*((req/r)^2.d0)*J2*(3.d0*(sin(lat)^2.d0)-1.d0)
  ;  rot    = -(cos(lat)^2.d0)*(omega^2.d0)*(r^3.d0)/(G*Mp)
  ;  pot    = const*(first+second+rot)
  ;endif else begin   
  ;  const  = G*Mp/r
  ;  first  = 1.d0
  ;  second = -0.5d0*((req/r)^2.d0)*J2*(3.d0*(sin(lat)^2.d0)-1.d0) 
  ;  ;mu    = sin(lat)
  ;  ;P4    = (35.d0*(mu^4.d0)-30.d0*(mu^2.d0)+3.d0)/8.d0
  ;  ;third = J4*((req/r)^4.d0)*P4    
  ;  rot   = 0.5d0*(cos(lat)^2.d0)*(omega^2.d0)*(r^3.d0)/(G*Mp)
  ;  pot   = const*(first+second+rot)
  ;endelse
  
  return,pot

end