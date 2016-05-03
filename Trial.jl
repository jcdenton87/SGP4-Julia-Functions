cd("C:\\Users\\jonathan.denton\\Documents\\Julia\\FirstCode\\")

include("SGP4Module.jl")

tic()
line01="1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
line02="2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";#     0.00      4320.0        360.00";

satrec=SGP4.twoline2rv(72,line01,line02,"m","e");
satrec,r,v=SGP4.sgp4(satrec,60); #in fractions of a minute
toc()

starttime=0;
stoptime=4320;
increment=360;
output1=zeros(13,20);
for ii=0:360:4320
satrec, r, v = SGP4.sgp4(satrec,ii);
year, month,day,hour,minute,sec= SGP4.invjday( satrec.jdsatepoch + SGP4.jday(2000,1,1,0,ii,0) -2451544.5 );
mu = 3.986004418000000e+05;
p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = SGP4.rv2coe(r,v, mu);
output1[div(ii,360)+1,1:20]=[ii,r[1],r[2],r[3],v[1],v[2],v[3],a, ecc, rad2deg(incl),rad2deg(omega),rad2deg(argp),rad2deg(nu),rad2deg(m), year,month,day,hour,minute,sec]'
end
