function saturn_1bar,colat,ref=ref
  
  ;This is the recommended 1 bar level shape from the Saturn Atmospheric Modeling 
  ;Working Group (SAMWG): EMBARGOED, NOT FOR OFFICIAL USE.  
  ;==============================================================================
  ;Input:
  ;colat - colatitude in radians
  ;ref   - reference radius (km) for alterations
  
  if(keyword_set(ref)) then s = ref else s = 58113.501d0
  s1  = 1.8444684d-5
  s2  = -0.069578323d0
  s3  = -3.0122526d-5
  s4  = 0.0056924209d0
  s5  = -5.7401771d-5
  s6  = -0.00075962418d0
  s7  = 1.2933112d-5
  s8  = 0.00034638233d0
  s9  = 1.2565650d-5
  s10 = -0.00021321543d0 
  x   = cos(colat)
  P1  = x
  P2  = 0.5d0*(3.d0*(x^2.d0)-1.d0)
  P3  = 0.5d0*(5.d0*(x^3.d0)-3.d0*x)
  P4  = (35.d0*(x^4.d0)-30.d0*(x^2.d0)+3.d0)/8.d0
  P5  = (63.d0*(x^5.d0)-70.d0*(x^3.d0)+15.d0*x)/8.d0
  P6  = (231.d0*(x^6.d0)-315.d0*(x^4.d0)+105.d0*(x^2.d0)-5.d0)/16.d0 
  P7  = (429.d0*(x^7.d0)-693.d0*(x^5.d0)+315.d0*(x^3.d0)-35.d0*x)/16.d0
  P8  = (6435.d0*(x^8.d0)-12012.d0*(x^6.d0)+6930.d0*(x^4.d0)-1260.d0*(x^2.d0) $
      + 35.d0)/128.d0
  P9  = (12155.d0*(x^9.d0)-25740.d0*(x^7.d0)+18018.d0*(x^5.d0)-4620.d0*(x^3.d0) $
      + 315.d0*x)/128.d0
  P10 = (46189.d0*(x^10.d0)-109395.d0*(x^8.d0)+90090.d0*(x^6.d0)-30030.d0*(x^4.d0) $
      + 3465.d0*(x^2.d0)-63.d0)/256.d0     
  fac   = 1.d0+s1*P1+s2*P2+s3*P3+s4*P4+s5*P5+s6*P6+s7*P7+s8*P8+s9*P9+s10*P10
  shape = s*fac
  return,shape

end