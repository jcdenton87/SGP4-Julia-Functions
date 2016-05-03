"""
# This is a vectorized version of the MATLAB sgp4.m function provided with
# AIAA 2006-6753, "Revisiting Spacetrack Report #3," available at
# http://celestrak.com/publications/AIAA/2006-6753/.  Instead of producing
# only a single position and velocity vector at time tsince, it now
# produces a matrix of position and velocity vectors for an input time
# vector.

# For the vectorized version to fully function, the following original
# files must be replaced:
# dpper.m -> dpper_vectorized.m
# dspace.m -> dspace_vectorized.m
# sgp4.m -> sgp4_vectorized.m
# sgp4init.m -> sgp4init_vectorized.m (minor change, simply points to
#               vectorized versions so that originals may be deleted from
#               working directory)
# twoline2rv.m -> twoline2rv_simple.m (initializes given input TLEs)

# Any questions or comments, contact Matthew Schmunk at
# matthew.schmunk@afit.edu.

# -----------------------------------------------------------------------------
#
#                              procedure sgp4
#
#  this procedure is the sgp4 prediction model from space command. this is an
#    updated and combined version of sgp4 and sdp4, which were originally
#    published separately in spacetrack report #3. this version follows the
#    methodology from the aiaa paper (2006) describing the history and
#    development of the code.
#
# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu
#   1.0 (aug 7, 2006) - update for paper dav
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600   28 jun 2005
#
#   inputs        :
#     satrec    - initialised structure from sgp4init() call.
#     tsince    - time eince epoch (minutes)
#
#   outputs       :
#     r           - position vector                     km
#     v           - velocity                            km/sec
#     return code - non-zero on error.
#                    1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
#                    2 - mean motion less than 0.0
#                    3 - pert elements, ecc < 0.0  or  ecc > 1.0
#                    4 - semi-latus rectum < 0.0
#                    5 - epoch elements are sub-orbital
#                    6 - satellite has decayed
#
#   locals        :
#     am          -
#     axnl, aynl        -
#     betal       -
#     COSIM   , SINIM   , COSOMM  , SINOMM  , Cnod    , Snod    , Cos2u   ,
#     Sin2u   , Coseo1  , Sineo1  , Cosi    , Sini    , Cosip   , Sinip   ,
#     Cosisq  , Cossu   , Sinsu   , Cosu    , Sinu
#     Delm        -
#     Delomg      -
#     Dndt        -
#     Eccm        -
#     EMSQ        -
#     Ecose       -
#     El2         -
#     Eo1         -
#     Eccp        -
#     Esine       -
#     Argpm       -
#     Argpp       -
#     Omgadf      -
#     Pl          -
#     R           -
#     RTEMSQ      -
#     Rdotl       -
#     Rl          -
#     Rvdot       -
#     Rvdotl      -
#     Su          -
#     T2  , T3   , T4    , Tc
#     Tem5, Temp , Temp1 , Temp2  , Tempa  , Tempe  , Templ
#     U   , Ux   , Uy    , Uz     , Vx     , Vy     , Vz
#     inclm       - inclination
#     mm          - mean anomaly
#     nm          - mean motion
#     nodem      - longi of ascending node
#     xinc        -
#     xincp       -
#     xl          -
#     xlm         -
#     mp          -
#     xmdf        -
#     xmx         -
#     xmy         -
#     nodedf     -
#     xnode       -
#     nodep      -
#     np          -
#
#   coupling      :
#     getgravconst
#     dpper
#     dspace
#
#   references    :
#     hoots, roehrich, norad spacetrack report #3 1980
#     hoots, norad spacetrack report #6 1986
#     hoots, schumacher and glover 2004
#     vallado, crawford, hujsak, kelso  2006
#  ----------------------------------------------------------------------------
"""

function sgp4_vectorized(satrec,tsince)

   # /* ------------------ set mathematical constants --------------- */
   #twopi = 2.0 * pi; #all instances replaced with 2*pi()
   #x2o3  = 2.0 / 3.0; #replaced with 2/3
   #temp4    =   1.0 + cos(pi-1.0e-9); #subbed directly into equation

   #  // sgp4fix identify constants and allow alternate values
   #global tumin mu radiusearthkm xke j2 j3 j4 j3oj2

   #ALTERED by M. Schmunk - fixing these values saves enormous time on long
   #runs.  The examples distributed with the 2006 "Revisited" document uses
   #WGS-72 constants, so those are used here.
   #ALTERED again by J. Denton to use globals...

  global tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2
   vkmpersec     = radiusearthkm * xke/60.0; #only for velocity

   # /* --------------------- clear sgp4 error flag ----------------- */
   satrec.t     = tsince;
   tsize = length(tsince);
   satrec.error = 0;

   # /* ------- update for secular gravity and atmospheric drag ----- */
   xmdf    = ones(tsize)*satrec.mo + satrec.mdot * satrec.t;
   argpdf  = ones(tsize)*satrec.argpo + satrec.argpdot * satrec.t;
   nodedf  = ones(tsize)*satrec.nodeo + satrec.nodedot * satrec.t;
   argpm   = argpdf;
   mm      = xmdf;
   t2      = satrec.t.^2;
   nodem   = nodedf + satrec.nodecf*t2;
   tempa   = ones(tsize) - satrec.cc1 * satrec.t;
   tempe   = satrec.bstar * satrec.cc4 * satrec.t;
   templ   = satrec.t2cof * t2;

   if (satrec.isimp != 1)
       delomg = satrec.omgcof * satrec.t;
       delm   = satrec.xmcof *
          ((ones(tsize) + satrec.eta * cos(mm)).^3 -
          ones(tsize)*satrec.delmo);
       temp   = delomg + delm;
       mm     = xmdf + temp;
       argpm  = argpdf - temp;
       t3     = t2 .* satrec.t;
       t4     = t3 .* satrec.t;
       tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 -
           satrec.d4 * t4;
       tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -
           ones(tsize)*satrec.sinmao);
       templ  = templ + satrec.t3cof * t3 + t4 .* (ones(tsize)*satrec.t4cof +
          satrec.t * satrec.t5cof);
   end

   nm    = satrec.no;
   em    = satrec.ecco;
   inclm = satrec.inclo;
   if (satrec.method == "d")
       tc = satrec.t;
       satrec.atime,em,argpm,inclm,satrec.xli,mm,satrec.xni,nodem,dndt,nm = dspace_vectorized(
           satrec.d2201,satrec.d2211,satrec.d3210,
           satrec.d3222,satrec.d4410,satrec.d4422,
           satrec.d5220,satrec.d5232,satrec.d5421,
           satrec.d5433,satrec.dedt,satrec.del1,
           satrec.del2,satrec.del3,satrec.didt,
           satrec.dmdt,satrec.dnodt,satrec.domdt,
           satrec.irez,satrec.argpo,satrec.argpdot,satrec.t,
           tc,satrec.gsto,satrec.xfact,satrec.xlamo,satrec.no,
           satrec.atime,em,argpm,inclm,satrec.xli,mm,
           satrec.xni,nodem,nm);
   end # // if method = d

   if (minimum(nm) <= 0.0) == 1
#       @printf(1,'# error nm #f\n', nm);
       @printf("\nError, object %d: object has nonpositive mean motion (Code 2).\n",
                satrec.satnum);
       satrec.error = 2;
   end
   am = (xke ./ nm).^(2/3) .* tempa.^2;

    if (maximum(am) < 0.95) == 1
    #       @printf(1,'# error em #f\n', em);
        @printf("\nError, object %d: object has probably reentered (Code 1).\n",
            satrec.satnum);
        satrec.error = 1;
       r = zeros(3,tsize)';
       v = zeros(3,tsize)'; #only for velocity
       return #leave this record
    end

   nm = xke ./ am.^1.5; #only for velocity
   em = ones(tsize).*em - tempe;

   if (maximum(em) >= 1.0) == 1
#       @printf(1,'# error em #f\n', em);
       @printf("\nError, object %d: orbit is parabolic or hyperbolic (Code 1).\n",
           satrec.satnum);
       satrec.error = 1;
       r = zeros(3,tsize)';
       v = zeros(3,tsize)'; #only for velocity
       return #leave this record
   end

   if (maximum(em) < -0.001) == 1
       #       @printf(1,'# error em #f\n', em);
       @printf("\nError, object %d: negative (<-0.001) eccentricity occurs (Code 1).\n",
           satrec.satnum);
       satrec.error = 1;
   end

   if (minimum(em) < 0.0) == 1;
       em[em .< 0.0]  = 1.0e-6;
   end

   mm     = mm + satrec.no * templ;
   xlm    = mm + argpm + nodem;
   #emsq   = em * em; #apparently never used, probably left in by mistake
   #temp   = 1.0 - emsq; #likewise
   nodem  = rem(nodem, 2*pi);
   argpm  = rem(argpm, 2*pi);
   xlm    = rem(xlm, 2*pi);
   mm     = rem(xlm - argpm - nodem, 2*pi);

   # /* ----------------- compute extra mean quantities ------------- */
   if length(inclm) != tsize
       sinim = ones(tsize)*sin(inclm);
       cosim = ones(tsize)*cos(inclm);
   else
       sinim = sin(inclm);
       cosim = cos(inclm);
   end

   # /* -------------------- add lunar-solar periodics -------------- */
   ep     = em;
   if length(inclm) != tsize
       xincp  = ones(tsize)*inclm;
   else
       xincp = inclm;
   end
   argpp  = argpm;
   nodep  = nodem;
   mp     = mm;
   sinip  = sinim;
   cosip  = cosim;

   if (satrec.method == "d")
       ep,xincp,nodep,argpp,mp = dpper_vectorized(
           satrec.e3,satrec.ee2,satrec.peo,
           satrec.pgho,satrec.pho,satrec.pinco,
           satrec.plo,satrec.se2,satrec.se3,
           satrec.sgh2,satrec.sgh3,satrec.sgh4,
           satrec.sh2,satrec.sh3,satrec.si2,
           satrec.si3,satrec.sl2,satrec.sl3,
           satrec.sl4,satrec.t,satrec.xgh2,
           satrec.xgh3,satrec.xgh4,satrec.xh2,
           satrec.xh3,satrec.xi2,satrec.xi3,
           satrec.xl2,satrec.xl3,satrec.xl4,
           satrec.zmol,satrec.zmos,satrec.inclo,
           satrec.init,ep,xincp,nodep,argpp,mp);

       if (minimum(xincp) < 0.0) == 1
           nodeptemp = nodep + ones(tsize)*pi;
           nodep[xincp .< 0.0] = nodeptemp[xincp .< 0.0];
           argpptemp = argpp - ones(tsize)*pi;
           argpp[xincp .< 0.0]  = argpptemp[xincp .< 0.0];
           xincp[xincp .< 0.0]  = -xincp[xincp .< 0.0];
       end

       if (maximum(ep) < 0.0) ==1
#           @printf(1,'# error ep #f\n', ep);
           satrec.error = 3;
           @printf("\nError, object %d: (Code 3).\n",
                satrec.satnum);
       end
       if (maximum(ep) > 1.0) == 1
#           @printf(1,"# error ep #f\n", ep);
           @printf("\nError, object %d: (Code 3).\n",
                satrec.satnum);
           satrec.error = 3;
       end
   end # // if method = d

   # /* -------------------- long period periodics ------------------ */
   if (satrec.method == "d")
       sinip =  sin(xincp);
       cosip =  cos(xincp);
       satrec.aycof = -0.5*j3oj2.*sinip;
       # // sgp4fix for divide by zero with xinco = 180 deg
       if max(abs(cosip+ones(tsize)) > 1.5e-12) == 1
           satrec.xlcof = -0.25 * j3oj2 * sinip .*
               ((3.0*ones(tsize) + 5.0 * cosip) ./ (ones(tsize)+cosip));
       else
           satrec.xlcof = -0.25 * j3oj2 * sinip .*
               ((3.0*ones(tsize) + 5.0 * cosip) ./ (ones(tsize) + ones(tsize)*cos(pi-1.0e-9)));
       end;
   end

   axnl = ep .* cos(argpp);
   temp = 1.0 ./ (am .* (ones(tsize) - ep.^2));
   aynl = ep .* sin(argpp) + temp .* satrec.aycof;
   xl   = mp + argpp + nodep + satrec.xlcof .* temp .* axnl;

#    disp('matching?')
#    [tsince;axnl;aynl;xl]

   # /* --------------------- solve kepler's equation --------------- */
   u    = rem(xl - nodep, 2*pi);
   eo1  = copy(u);
   tem5 = 9999.9*ones(tsize);
   ktr = ones(tsize);

   sineo1 = sin(eo1); #initialize
   coseo1 = cos(eo1);

   # //   sgp4fix for kepler iteration
   # //   the following iteration needs better limits on corrections
   while (maximum(abs(tem5)) >= 1.0e-12 )== 1 && maximum(ktr) <= 10
       tem5orig = abs(tem5) .>= 1.0e-12; #we care about updating records *originally* meeting conditions
       sineo1[tem5orig] = sin(eo1[tem5orig]);
       coseo1[tem5orig] = cos(eo1[tem5orig]);
       tem5temp1 = ones(tsize) - coseo1 .* axnl - sineo1 .* aynl;
       tem5[tem5orig]  = tem5temp1[tem5orig];
       tem5[tem5orig]   = (u[tem5orig] - aynl[tem5orig] .*coseo1[tem5orig] + axnl[tem5orig] .*sineo1[tem5orig] - eo1[tem5orig])./tem5[tem5orig];
       tem5temp2 = tem5;
       if (maximum(abs(tem5temp2)) > 0.95) == 1
           if maximum((abs(tem5temp2) .> 0.95) + (tem5temp2 .> 0) .> 1) == 1
               tem5temp2[((abs(tem5temp2) .> 0.95) + (tem5temp2 .> 0)) .> 1] = 0.95;
           end
           if (maximum((abs(tem5temp2) .> 0.95) + (tem5temp2 .<= 0)) .> 1) == 1
               tem5temp2[((abs(tem5temp2) .> 0.95) + (tem5temp2 .<= 0)) .> 1] = -0.95;
           end
       end
       tem5[tem5orig] = tem5temp2[tem5orig];
       eo1[tem5orig] = eo1[tem5orig] + tem5[tem5orig];
       ktrtemp = ktr + 1;
       ktr[tem5orig] = ktrtemp[tem5orig];
   end

   # /* ------------- short period preliminary quantities ----------- */
   ecose = axnl.*coseo1 + aynl.*sineo1;
   esine = axnl.*sineo1 - aynl.*coseo1;
   el2   = axnl.^2 + aynl.^2;
   pl    = am.*(1.0-el2);
   if maximum(pl .< 0.0) == 1
#       @printf(1,"# error pl #f\n", pl);
       @printf("\nError, object %d: (Code 4).\n",
                satrec.satnum);
       satrec.error = 4;
       r = zeros(tsize)';
       v = zeros(tsize)'; #only for velocity
   else
       rl     = am .* (ones(tsize) - ecose);
       rdotl  = sqrt(am) .* esine./rl; #only for velocity
       rvdotl = sqrt(pl) ./ rl;
       betal  = sqrt(ones(tsize) - el2);
       temp   = esine ./ (ones(tsize) + betal);
       sinu   = am ./ rl .* (sineo1 - aynl - axnl .* temp);
       cosu   = am ./ rl .* (coseo1 - axnl + aynl .* temp);
       su     = atan2(sinu, cosu);
       sin2u  = (cosu + cosu) .* sinu;
       cos2u  = ones(tsize) - 2.0*ones(tsize) .* sinu.^2;
       temp   = 1.0 ./ pl;
       temp1  = 0.5 * j2 * temp;
       temp2  = temp1 .* temp;

       # /* -------------- update for short period periodics ------------ */
       if (satrec.method == "d")
           cosisq = cosip.^2;
           satrec.con41  = 3.0*cosisq - ones(tsize);
           satrec.x1mth2 = ones(tsize) - cosisq;
           satrec.x7thm1 = 7.0*cosisq - ones(tsize);
       end

       mrt   = rl .* (ones(tsize) - 1.5 * temp2 .* betal .* satrec.con41) +
           0.5 * temp1 .* satrec.x1mth2 .* cos2u;
       su    = su - 0.25 * temp2 .* satrec.x7thm1 .* sin2u;
       xnode = nodep + 1.5 * temp2 .* cosip .* sin2u;
       xinc  = xincp + 1.5 * temp2 .* cosip .* sinip .* cos2u;
       mvt   = rdotl - nm .* temp1 .* satrec.x1mth2 .* sin2u ./ xke;
       rvdot = rvdotl + nm .* temp1 .* (satrec.x1mth2 .* cos2u +
           1.5 * satrec.con41) ./ xke; #only for velocity

       # /* --------------------- orientation vectors ------------------- */
       sinsu =  sin(su);
       cossu =  cos(su);
       snod  =  sin(xnode);
       cnod  =  cos(xnode);
       sini  =  sin(xinc);
       cosi  =  cos(xinc);
       xmx   = -snod .* cosi;
       xmy   =  cnod .* cosi;
       ux    =  xmx .* sinsu + cnod .* cossu;
       uy    =  xmy .* sinsu + snod .* cossu;
       uz    =  sini .* sinsu;
       vx    =  xmx .* cossu - cnod .* sinsu; #only for velocity
       vy    =  xmy .* cossu - snod .* sinsu;
       vz    =  sini .* cossu;

       # /* --------- position and velocity (in km and km/sec) ---------- */
       r = radiusearthkm*[(mrt .* ux)' ;(mrt .* uy)'; (mrt .* uz)'];
       v = vkmpersec*[(mvt .* ux + rvdot .* vx)';
           (mvt .* uy + rvdot .* vy)'; (mvt .* uz + rvdot .* vz)'];
   end # // if pl > 0

       # // sgp4fix for decaying satellites
        if maximum(mrt .< 1.0) == 1
   #         printf("# decay condition #11.6f \n",mrt);
            @printf("\nError, object %d: reenters during run (Code 6) ",
                satrec.satnum);
            satrec.error = 6;
            zeroout = (mrt < 1.0)';
            r(zeroout,:) = zeros(size(r(zeroout,:)));
            v(zeroout,:) = zeros(size(v(zeroout,:))); #only for velocity
        end

#    global idebug dbgfile  #removed for timemms
#    if idebug
#        debug7;
#    end

   return satrec, r, v
 end;
