"""
# -----------------------------------------------------------------------------
#
#                              procedure sgp4init
#
#   this procedure initializes variables for sgp4.
#
# Author:
#   Jeff Beck <--Converted to Julia by Jonathan Denton (jonathan.denton.1@us.af.mil)
#   beckja@alumni.lehigh.edu
#   1.0 (aug 7, 2006) - update for paper dav
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600   28 jun 2005
#
#   inputs        :
#     satn        - satellite number
#     bstar       - sgp4 type drag coefficient              kg/m2er
#     ecco        - eccentricity
#     epoch       - epoch time in days from jan 0, 1950. 0 hr
#     argpo       - argument of perigee (output if ds)
#     inclo       - inclination
#     mo          - mean anomaly (output if ds)
#     no          - mean motion
#     nodeo      - right ascension of ascending node
#
#   outputs       :
#     satrc      - common values for subsequent calls
#     return code - non-zero on error.
#                    1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
#                    2 - mean motion less than 0.0
#                    3 - pert elements, ecc < 0.0  or  ecc > 1.0
#                    4 - semi-latus rectum < 0.0
#                    5 - epoch elements are sub-orbital
#                    6 - satellite has decayed
#
#   locals        :
#     CNODM  , SNODM  , COSIM  , SINIM  , COSOMM , SINOMM
#     Cc1sq  , Cc2    , Cc3
#     Coef   , Coef1
#     cosio4      -
#     day         -
#     dndt        -
#     em          - eccentricity
#     emsq        - eccentricity squared
#     eeta        -
#     etasq       -
#     gam         -
#     argpm       - argument of perigee
#     ndem        -
#     inclm       - inclination
#     mm          - mean anomaly
#     nm          - mean motion
#     perige      - perigee
#     pinvsq      -
#     psisq       -
#     qzms24      -
#     rtemsq      -
#     s1, s2, s3, s4, s5, s6, s7          -
#     sfour       -
#     ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
#     sz1, sz2, sz3
#     sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
#     tc          -
#     temp        -
#     temp1, temp2, temp3       -
#     tsi         -
#     xpidot      -
#     xhdot1      -
#     z1, z2, z3          -
#     z11, z12, z13, z21, z22, z23, z31, z32, z33         -
#
#   coupling      :
#     getgravconst
#     initl       -
#     dscom       -
#     dpper       -
#     dsinit      -
#     sgp4        -
#
#   references    :
#     hoots, roehrich, norad spacetrack report #3 1980
#     hoots, norad spacetrack report #6 1986
#     hoots, schumacher and glover 2004
#     vallado, crawford, hujsak, kelso  2006
#  ----------------------------------------------------------------------------
"""

function sgp4init(whichconst, satrc, xbstar, xecco, epoch,
         xargpo, xinclo, xmo, xno, xnodeo);

   # /* ------------------------ initialization --------------------- */
   # /* ----------- set all near earth variables to zero ------------ */
   satrc.isimp   = 0;   satrc.method = 'n'; satrc.aycof    = 0.0;
   satrc.con41   = 0.0; satrc.cc1    = 0.0; satrc.cc4      = 0.0;
   satrc.cc5     = 0.0; satrc.d2     = 0.0; satrc.d3       = 0.0;
   satrc.d4      = 0.0; satrc.delmo  = 0.0; satrc.eta      = 0.0;
   satrc.argpdot = 0.0; satrc.omgcof = 0.0; satrc.sinmao   = 0.0;
   satrc.t       = 0.0; satrc.t2cof  = 0.0; satrc.t3cof    = 0.0;
   satrc.t4cof   = 0.0; satrc.t5cof  = 0.0; satrc.x1mth2   = 0.0;
   satrc.x7thm1  = 0.0; satrc.mdot   = 0.0; satrc.nodedot = 0.0;
   satrc.xlcof   = 0.0; satrc.xmcof  = 0.0; satrc.nodecf  = 0.0;

   # /* ----------- set all deep space variables to zero ------------ */
   satrc.irez  = 0;   satrc.d2201 = 0.0; satrc.d2211 = 0.0;
   satrc.d3210 = 0.0; satrc.d3222 = 0.0; satrc.d4410 = 0.0;
   satrc.d4422 = 0.0; satrc.d5220 = 0.0; satrc.d5232 = 0.0;
   satrc.d5421 = 0.0; satrc.d5433 = 0.0; satrc.dedt  = 0.0;
   satrc.del1  = 0.0; satrc.del2  = 0.0; satrc.del3  = 0.0;
   satrc.didt  = 0.0; satrc.dmdt  = 0.0; satrc.dnodt = 0.0;
   satrc.domdt = 0.0; satrc.e3    = 0.0; satrc.ee2   = 0.0;
   satrc.peo   = 0.0; satrc.pgho  = 0.0; satrc.pho   = 0.0;
   satrc.pinco = 0.0; satrc.plo   = 0.0; satrc.se2   = 0.0;
   satrc.se3   = 0.0; satrc.sgh2  = 0.0; satrc.sgh3  = 0.0;
   satrc.sgh4  = 0.0; satrc.sh2   = 0.0; satrc.sh3   = 0.0;
   satrc.si2   = 0.0; satrc.si3   = 0.0; satrc.sl2   = 0.0;
   satrc.sl3   = 0.0; satrc.sl4   = 0.0; satrc.gsto  = 0.0;
   satrc.xfact = 0.0; satrc.xgh2  = 0.0; satrc.xgh3  = 0.0;
   satrc.xgh4  = 0.0; satrc.xh2   = 0.0; satrc.xh3   = 0.0;
   satrc.xi2   = 0.0; satrc.xi3   = 0.0; satrc.xl2   = 0.0;
   satrc.xl3   = 0.0; satrc.xl4   = 0.0; satrc.xlamo = 0.0;
   satrc.zmol  = 0.0; satrc.zmos  = 0.0; satrc.atime = 0.0;
   satrc.xli   = 0.0; satrc.xni   = 0.0;

   # sgp4fix - note the following variables are also passed directly via satrc.
   # it is possible to streamline the sgp4init call by deleting the "x"
   # variables, but the user would need to set the satrc.* values first. we
   # include the additional assignment in case twoline2rv is not used.
   satrc.bstar      = xbstar;
   satrc.ecco       = xecco;
   satrc.argpo      = xargpo;
   satrc.inclo      = xinclo;
   satrc.mo         = xmo;
   satrc.no         = xno;
   satrc.nodeo      = xnodeo;

   #     /* -------------------- wgs-72 earth constants ----------------- */
   #     // sgp4fix identify constants and allow alternate values
   global tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2
   (tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2) = getgravc( whichconst );

   ss     = 78.0 / radiusearthkm + 1.0;
   qzms2t = ((120.0 - 78.0) / radiusearthkm)^4;
   x2o3   =  2.0 / 3.0;
   # sgp4fix divisor for divide by zero check on inclination
   # the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
   # 1.5 e-12, so the threshold was changed to 1.5e-12 for consistancy
   temp4    =   1.5e-12;

   satrc.init = 'y';
   satrc.t    = 0.0;

   (ainv,  ao,     satrc.con41,   con42,  cosio,  cosio2, einv,   eccsq,
           satrc.method,  omeosq, posq,   rp,     rteosq, sinio,
           satrc.gsto,    satrc.no) = initl(satrc.ecco,    epoch,  satrc.inclo,   satrc.no,
           satrc.satnum);

   satrc.error = 0;

   # sgp4fix remove this check as it is unnecessary
   # the mrt check in sgp4 handles decaying satellite cases even if the starting
   #if (rp < 1.0)
   #   printf("# *** satn#d epoch elts sub-orbital ***\n", satn);
   #    satrc.error = 5;
   #end

   if ((omeosq >= 0.0 ) | ( satrc.no >= 0.0))
       satrc.isimp = 0;
       if (rp < (220.0 / radiusearthkm + 1.0))
           satrc.isimp = 1;
       end
       sfour  = ss;
       qzms24 = qzms2t;
       perige = (rp - 1.0) * radiusearthkm;

       # /* - for perigees below 156 km, s and qoms2t are altered - */
       if (perige < 156.0)
           sfour = perige - 78.0;
           if (perige < 98.0)
               sfour = 20.0;
           end
           qzms24 = ((120.0 - sfour) / radiusearthkm)^4.0;
           sfour  = sfour / radiusearthkm + 1.0;
       end
       pinvsq = 1.0 / posq;

       tsi  = 1.0 / (ao - sfour);
       satrc.eta  = ao * satrc.ecco * tsi;
       etasq = satrc.eta * satrc.eta;
       eeta  = satrc.ecco * satrc.eta;
       psisq = abs(1.0 - etasq);
       coef  = qzms24 * tsi^4.0;
       coef1 = coef / psisq^3.5;
       cc2   = coef1 * satrc.no * (ao * (1.0 + 1.5 * etasq + eeta *
           (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satrc.con41 *
           (8.0 + 3.0 * etasq * (8.0 + etasq)));
       satrc.cc1   = satrc.bstar * cc2;
       cc3   = 0.0;
       if (satrc.ecco > 1.0e-4)
           cc3 = -2.0 * coef * tsi * j3oj2 * satrc.no * sinio / satrc.ecco;
       end
       satrc.x1mth2 = 1.0 - cosio2;
       satrc.cc4    = 2.0* satrc.no * coef1 * ao * omeosq *
           (satrc.eta * (2.0 + 0.5 * etasq) + satrc.ecco *
           (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
           (-3.0 * satrc.con41 * (1.0 - 2.0 * eeta + etasq *
           (1.5 - 0.5 * eeta)) + 0.75 * satrc.x1mth2 *
           (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrc.argpo)));
       satrc.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
           (etasq + eeta) + eeta * etasq);
       cosio4 = cosio2 * cosio2;
       temp1  = 1.5 * j2 * pinvsq * satrc.no;
       temp2  = 0.5 * temp1 * j2 * pinvsq;
       temp3  = -0.46875 * j4 * pinvsq * pinvsq * satrc.no;
       satrc.mdot     = satrc.no + 0.5 * temp1 * rteosq * satrc.con41 +
           0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
       satrc.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *
           (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
           temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
       xhdot1            = -temp1 * cosio;
       satrc.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
           2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
       xpidot            =  satrc.argpdot+ satrc.nodedot;
       satrc.omgcof   = satrc.bstar * cc3 * cos(satrc.argpo);
       satrc.xmcof    = 0.0;
       if (satrc.ecco > 1.0e-4)
           satrc.xmcof = -x2o3 * coef * satrc.bstar / eeta;
       end
       satrc.nodecf = 3.5 * omeosq * xhdot1 * satrc.cc1;
       satrc.t2cof   = 1.5 * satrc.cc1;

       # // sgp4fix for divide by zero with xinco = 180 deg
       if (abs(cosio+1.0) > 1.5e-12)
          satrc.xlcof   = -0.25 * j3oj2 * sinio *
              (3.0 + 5.0 * cosio) / (1.0 + cosio);
       else
         satrc.xlcof   = -0.25 * j3oj2 * sinio *
              (3.0 + 5.0 * cosio) / temp4;
       end
       satrc.aycof   = -0.5 * j3oj2 * sinio;
       satrc.delmo   = (1.0 + satrc.eta * cos(satrc.mo))^3;
       satrc.sinmao  = sin(satrc.mo);
       satrc.x7thm1  = 7.0 * cosio2 - 1.0;

       # /* --------------- deep space initialization ------------- */
       if ((2*pi / satrc.no) >= 225.0)
           satrc.method = 'd';
           satrc.isimp  = 1;
           tc    =  0.0;
           inclm = satrc.inclo;

           (sinim,cosim,sinomm,cosomm,snodm,cnodm,day,satrc.e3,satrc.ee2,
               em,emsq,gam,satrc.peo,satrc.pgho,satrc.pho,satrc.pinco,
               satrc.plo,rtemsq,satrc.se2,satrc.se3,satrc.sgh2,
               satrc.sgh3,satrc.sgh4,satrc.sh2,satrc.sh3,satrc.si2,
               satrc.si3,satrc.sl2,satrc.sl3,satrc.sl4,s1,s2,s3,s4,s5,
               s6,s7,ss1,ss2,ss3,ss4,ss5,ss6,ss7,sz1,sz2,sz3,sz11,sz12,
               sz13,sz21,sz22,sz23,sz31,sz32,sz33,satrc.xgh2,satrc.xgh3,
               satrc.xgh4,satrc.xh2,satrc.xh3,satrc.xi2,satrc.xi3,
               satrc.xl2,satrc.xl3,satrc.xl4,nm,z1,z2,z3,z11,z12,z13,
               z21,z22,z23,z31,z32,z33,satrc.zmol,satrc.zmos) = dscom(epoch,
               satrc.ecco,satrc.argpo,tc,satrc.inclo,
               satrc.nodeo,satrc.no);

           (satrc.ecco,satrc.inclo,satrc.nodeo,satrc.argpo,satrc.mo)= dpper(satrc.e3,
               satrc.ee2,satrc.peo,satrc.pgho,
               satrc.pho,satrc.pinco,satrc.plo,satrc.se2,satrc.se3,
               satrc.sgh2,satrc.sgh3,satrc.sgh4,satrc.sh2,satrc.sh3,
               satrc.si2,satrc.si3,satrc.sl2,satrc.sl3,satrc.sl4,
               satrc.t,satrc.xgh2,satrc.xgh3,satrc.xgh4,satrc.xh2,
               satrc.xh3,satrc.xi2,satrc.xi3,satrc.xl2,satrc.xl3,
               satrc.xl4,satrc.zmol,satrc.zmos,inclm,satrc.init,
               satrc.ecco,satrc.inclo,satrc.nodeo,satrc.argpo,satrc.mo);

           argpm  = 0.0;
           nodem  = 0.0;
           mm     = 0.0;

           (em,argpm,inclm,mm,nm,nodem,satrc.irez,satrc.atime,
               satrc.d2201,satrc.d2211,satrc.d3210,satrc.d3222,
               satrc.d4410,satrc.d4422,satrc.d5220,satrc.d5232,
               satrc.d5421,satrc.d5433,satrc.dedt,satrc.didt,
               satrc.dmdt,dndt,satrc.dnodt,satrc.domdt,satrc.del1,
               satrc.del2,satrc.del3,
                #ses,sghl,sghs,sgs,shl,shs,sis,sls,theta,
               satrc.xfact,satrc.xlamo,satrc.xli,satrc.xni)=dsinit(cosim,emsq,satrc.argpo,s1,s2,s3,s4,s5,sinim,ss1,ss2,ss3,
               ss4,ss5,sz1,sz3,sz11,sz13,sz21,sz23,sz31,sz33,satrc.t,tc,
               satrc.gsto,satrc.mo,satrc.mdot,satrc.no,satrc.nodeo,
               satrc.nodedot,xpidot,z1,z3,z11,z13,z21,z23,z31,z33,em,
               argpm,inclm,mm,nm,nodem,satrc.ecco,eccsq);
       end

       # /* ----------- set variables if not deep space ----------- */
       if (satrc.isimp != 1)
           cc1sq          = satrc.cc1 * satrc.cc1;
           satrc.d2    = 4.0 * ao * tsi * cc1sq;
           temp           = satrc.d2 * tsi * satrc.cc1 / 3.0;
           satrc.d3    = (17.0 * ao + sfour) * temp;
           satrc.d4    = 0.5 * temp * ao * tsi *
               (221.0 * ao + 31.0 * sfour) * satrc.cc1;
           satrc.t3cof = satrc.d2 + 2.0 * cc1sq;
           satrc.t4cof = 0.25 * (3.0 * satrc.d3 + satrc.cc1 *
               (12.0 * satrc.d2 + 10.0 * cc1sq));
           satrc.t5cof = 0.2 * (3.0 * satrc.d4 +
               12.0 * satrc.cc1 * satrc.d3 +
               6.0 * satrc.d2 * satrc.d2 +
               15.0 * cc1sq * (2.0 * satrc.d2 + cc1sq));
       end
   end # // if omeosq = 0

   # /* finally propogate to zero epoch to initialise all others. */
   # sgp4fix take out check to let satellites process until they are actually below earth surface
   #if(satrc.error == 0)

       (satrc, r, v) = sgp4(satrc, 0.0);
   #end

   satrc.init = 'n';

   return satrc

end
