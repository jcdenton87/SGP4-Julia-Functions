cd("$(homedir())/Documents/Julia/FirstCode")

include("SGP4Module.jl")

tic()
#Verification Test Case #1


line01="1 41435U 16021B   16090.96201042  .00007039  00000-0  10449-2 0  9996";
line02="2 41435  55.0041  77.5955 7274051 174.8379 202.5445  2.30404667    34     0.00      4320.0        360.00 ";

sat,startmfe,stopmfe,deltamin,line01,line02=SGP4.twoline2rv(72,line01,line02,"v","c");
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
