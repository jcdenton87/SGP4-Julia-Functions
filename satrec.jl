
"""
#This type file defines the satrec type
#Fields are used by various parts of the SGP4 module
#Initialize to default values as needed, pay attention to integers
#and float (meaning 0 versus 0.0)
"""

type satrec
  error
  satnum
  epochyr
  epochdays
  ndot
  nddot
  bstar
  inclo
  nodeo
  ecco
  argpo
  mo
  no
  a
  alta
  altp
  jdsatepoch
  isimp
  method
  aycof
  con41
  cc1
  cc4
  cc5
  d2
  d3
  d4
  delmo
  eta
  argpdot
  omgcof
  sinmao
  t
  t2cof
  t3cof
  t4cof
  t5cof
  x1mth2
  x7thm1
  mdot
  nodedot
  xlcof
  xmcof
  nodecf
  irez
  d2201
  d2211
  d3210
  d3222
  d4410
  d4422
  d5220
  d5232
  d5421
  d5433
  dedt
  del1
  del2
  del3
  didt
  dmdt
  dnodt
  domdt
  e3
  ee2
  peo
  pgho
  pho
  pinco
  plo
  se2
  se3
  sgh2
  sgh3
  sgh4
  sh2
  sh3
  si2
  si3
  sl2
  sl3
  sl4
  gsto
  xfact
  xgh2
  xgh3
  xgh4
  xh2
  xh3
  xi2
  xi3
  xl2
  xl3
  xl4
  xlamo
  zmol
  zmos
  atime
  xli
  xni
  init
end