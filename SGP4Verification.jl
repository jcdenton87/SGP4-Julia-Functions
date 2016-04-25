include("SGP4Module.jl")

tic()
#Verification Test Case #1
longstr1="1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
longstr2="2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.00      4320.0        360.00 ";

sat, startmfe, stopmfe, deltamin,longstr1,longstr2=SGP4.twoline2rv(72,longstr1,longstr2,"v","e");
sat,r,v=SGP4.sgp4(sat,60); #in fractions of a minute
toc()

output1=zeros(13,20);
for ii=0:360:4320
sat, r, v = SGP4.sgp4(sat,ii);
year, month,day,hour,minute,sec= SGP4.invjday( sat.jdsatepoch + SGP4.jday(2000,1,1,0,ii,0) -2451544.5 );
mu = 3.986004418000000e+05;
p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = SGP4.rv2coe(r,v, mu);
output1[div(ii,360)+1,1:20]=[ii,r[1],r[2],r[3],v[1],v[2],v[3],a, ecc, rad2deg(incl),rad2deg(omega),rad2deg(argp),rad2deg(nu),rad2deg(m), year,month,day,hour,minute,sec]'
end


tic()
#Verification Test Case #2
longstr1="1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955";
longstr2="2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145  -5184.0     -4896.0        120.00 ";

sat, startmfe, stopmfe, deltamin,longstr1,longstr2=SGP4.twoline2rv(72,longstr1,longstr2,"v","e");
sat,r,v=SGP4.sgp4(sat,60); #in fractions of a minute
toc()

output2=zeros(5,20);
tempcount=1;
for ii=[0 -5184 -5064 -4944 -4896]
sat, r, v = SGP4.sgp4(sat,ii);
year, month,day,hour,minute,sec= SGP4.invjday( sat.jdsatepoch + SGP4.jday(2000,1,1,0,ii,0) -2451544.5 );
mu = 3.986004418000000e+05;
p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = SGP4.rv2coe(r,v, mu);
output2[tempcount,1:20]=[ii,r[1],r[2],r[3],v[1],v[2],v[3],a, ecc, rad2deg(incl),rad2deg(omega),rad2deg(argp),rad2deg(nu),rad2deg(m), year,month,day,hour,minute,sec]'
tempcount=tempcount+1;
end

tic()
#Verification Test Case #3
longstr1="1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985";
longstr2="2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774      0.0      2880.0        120.00 ";

sat, startmfe, stopmfe, deltamin,longstr1,longstr2=SGP4.twoline2rv(72,longstr1,longstr2,"v","e");
sat,r,v=SGP4.sgp4(sat,60); #in fractions of a minute
toc()

output3=zeros(25,20);
for ii=0:120:2880
sat, r, v = SGP4.sgp4(sat,ii);
year, month,day,hour,minute,sec= SGP4.invjday( sat.jdsatepoch + SGP4.jday(2000,1,1,0,ii,0) -2451544.5 );
mu = 3.986004418000000e+05;
p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = SGP4.rv2coe(r,v, mu);
output3[div(ii,120)+1,1:20]=[ii,r[1],r[2],r[3],v[1],v[2],v[3],a, ecc, rad2deg(incl),rad2deg(omega),rad2deg(argp),rad2deg(nu),rad2deg(m), year,month,day,hour,minute,sec]'
end

println("Output1")
println(showall(output1))
println("Output2")
println(showall(output2))
println("Output3")
println(showall(output3))
