  whichconst=72;
  typerun="v";
  typeinput="e";
  
  global tumin, radiusearthkm, xke, j2, j3, j4, j3oj2

    xpdotp   =  1440.0 / (2.0*pi);   # 229.1831180523293;  # [rev/day]/[rad/min]

    revnum = 0;
    elnum  = 0;
    year   = 0;
    sat=SGP4.satrec(0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0); #default values initialization for satrec...

#     // set the implied decimal points since doing a formated read
#     // fixes for bad input data values (missing, )

    # parse first line
    carnumb = parse(Int,longstr1[1]);
    sat.satnum = parse(Int,longstr1[3:7]);
    classification = longstr1[8];
    intldesg = longstr1[10:17];
    sat.epochyr = parse(Int,longstr1[19:20]);
    sat.epochdays = parse(Float64,longstr1[21:32]);
    sat.ndot = parse(Float64,longstr1[34:43]);
    sat.nddot = parse(Float64,longstr1[45:50]);
    nexp = parse(Int,longstr1[51:52]);
    sat.bstar = parse(Float64,longstr1[54:59]);
    ibexp = parse(Int,longstr1[60:61]);
    numb = parse(Int,longstr1[63]);
    elnum = parse(Float64,longstr1[65:68]);

    # parse second line
    if (typerun == "v")
        cardnumb = parse(Int,longstr2[1]);
        sat.satnum = parse(Float64,longstr2[3:7]);
        sat.inclo = parse(Float64,longstr2[8:16]);
        sat.nodeo = parse(Float64,longstr2[17:25]);
        sat.ecco = parse(Float64,"."*longstr2[27:33]);
        sat.argpo = parse(Float64,longstr2[34:42]);
        sat.mo = parse(Float64,longstr2[43:51]);
        sat.no = parse(Float64,longstr2[52:63]);
        revnum = parse(Float64,longstr2[64:68]);
        startmfe = parse(Float64,longstr2[70:81]);
        stopmfe  = parse(Float64,longstr2[83:96]);
        deltamin = parse(Float64,longstr2[97:105]);
    else
        cardnumb = parse(Int,longstr2[1]);
        sat.satnum = parse(Float64,longstr2[3:7]);
        sat.inclo = parse(Float64,longstr2[8:16]);
        sat.nodeo = parse(Float64,longstr2[17:25]);
        sat.ecco = parse(Float64,"."*longstr2[27:33]);
        sat.argpo = parse(Float64,longstr2[34:42]);
        sat.mo = parse(Float64,longstr2[43:51]);
        sat.no = parse(Float64,longstr2[52:63]);
        revnum = parse(Float64,longstr2[64:68]);
    end

#     // ---- find no, ndot, nddot ----
    sat.no   = sat.no / xpdotp; #//* rad/min
    sat.nddot= sat.nddot * 10.0^nexp;
    sat.bstar= sat.bstar * 10.0^ibexp;

#     // ---- convert to sgp4 units ----
    sat.a    = (sat.no*tumin)^(-2/3);                # [er]
    sat.ndot = sat.ndot  / (xpdotp*1440.0);          # [rad/min^2]
    sat.nddot= sat.nddot / (xpdotp*1440.0*1440);     # [rad/min^3]

#     // ---- find standard orbital elements ----
    sat.inclo = deg2rad(sat.inclo);
    sat.nodeo = deg2rad(sat.nodeo);
    sat.argpo = deg2rad(sat.argpo);
    sat.mo    = deg2rad(sat.mo);

    sat.alta = sat.a*(1.0 + sat.ecco) - 1.0;
    sat.altp = sat.a*(1.0 - sat.ecco) - 1.0;

#     // ----------------------------------------------------------------
#     // find sgp4epoch time of element set
#     // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
#     // and minutes from the epoch (time)
#     // --------------------------------------------------------------

#     // ------------- temp fix for years from 1957-2056 ----------------
#     // ------ correct fix will occur when year is 4-digit in 2le ------
     if (sat.epochyr < 57)
         year= sat.epochyr + 2000;
       else
         year= sat.epochyr + 1900;
     end;

     (mon,day,hr,minute,sec) = SGP4.days2mdh( year,sat.epochdays );
     sat.jdsatepoch = SGP4.jday( year,mon,day,hr,minute,sec );

#     // input start stop times manually
     if ((typerun != "v") && (typerun != "c"))
         # ------------- enter start/stop ymd hms values --------------------
           if (typeinput == 'e')
               print("input start year "); startyear = parse(Int,readline(STDIN));
               print("input start mon "); startmon  = parse(Int,readline(STDIN));
               print("input start day "); startday  = parse(Int,readline(STDIN));
               print("input start hr "); starthr   = parse(Int,readline(STDIN));
               print("input start min "); startmin  = parse(Int,readline(STDIN));
               print("input start sec "); startsec  = parse(Int,readline(STDIN));
               jdstart = jday( startyear,startmon,startday,starthr,startmin,startsec );

               print("input stop year "); stopyear = parse(Int,readline(STDIN));
               print("input stop mon "); stopmon  = parse(Int,readline(STDIN));
               print("input stop day "); stopday  = parse(Int,readline(STDIN));
               print("input stop hr "); stophr   = parse(Int,readline(STDIN));
               print("input stop min "); stopmin  = parse(Int,readline(STDIN));
               print("input stop sec "); stopsec  = parse(Int,readline(STDIN));
               jdstop = jday( stopyear,stopmon,stopday,stophr,stopmin,stopsec );

               startmfe = (jdstart - sat.jdsatepoch) * 1440.0;
               stopmfe  = (jdstop - sat.jdsatepoch) * 1440.0;
               deltamin = print("input time step in minutes"); parse(Int,readline(STDIN));
           end;
           # -------- enter start/stop year and days of year values -----------
           if (typeinput == "d")
               print("input start year "); startyear    = parse(Int,readline(STDIN));
               print("input start dayofyr "); startdayofyr = parse(Int,readline(STDIN));
               print("input stop year "); stopyear     = parse(Int,readline(STDIN));
               print("input stop dayofyr "); stopdayofyr  = parse(Int,readline(STDIN));

               (mon,day,hr,minute,sec) = days2mdh( startyear,startdayofyr);
               jdstart = jday( startyear,mon,day,hr,minute,sec);
               (mon,day,hr,minute,sec) = days2mdh( stopyear,stopdayofyr);
               jdstop = jday( stopyear,mon,day,hr,minute,sec);

               startmfe = (jdstart - sat.jdsatepoch) * 1440.0;
               stopmfe  = (jdstop - sat.jdsatepoch) * 1440.0;
               deltamin = print("input time step in minutes"); parse(Int,readline(STDIN));
           end;
           # ------------------ enter start/stop mfe values -------------------
           if (typeinput == "m")
               print("input start mfe"); startmfe = parse(Int,readline(STDIN));
               print("input stop mfe"); stopmfe  = parse(Int,readline(STDIN));
               print("input time step in minutes"); deltamin = parse(Int,readline(STDIN));
           end;
       end;
#     // perform complete catalog evaluation
     if (typerun == "c")
         startmfe =  -1440.0;
         stopmfe  =  1440.0;
         deltamin = 20.0;
     end;

#     // ------------- initialize the orbit at sgp4epoch --------------
     sgp4epoch = sat.jdsatepoch - 2433281.5; # days since 0 Jan 1950
	 sat = SGP4.sgp4init(whichconst, sat, sat.bstar, sat.ecco, sgp4epoch,
         sat.argpo, sat.inclo, sat.mo, sat.no, sat.nodeo);

return sat, startmfe, stopmfe, deltamin,longstr1,longstr2