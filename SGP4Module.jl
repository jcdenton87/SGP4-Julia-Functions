"""
#Author: Jonathan Denton (jonathan.denton.1@us.af.mil)
#Date: 31 Mar 16
#This file includes all of the constants and functions in the SGP4 module
#They can be accessed via SGP4."modulename"
#Function names include:
#dspace,satrec,getgravc,dscom,gstime,initl,dpper,dsinit,jday,invjday,sgp4init,
#days2mdh,twoline2rv,newtonnu,rv2coe,sgp4,and angl
#Use ?SGP4."name" to read markdown for each function
#Note: I'm not a programmer so coding is bare minimum to get MatLab version running in Julia...
"""


module SGP4

include("dspace.jl")
include("satrec.jl")
include("constants.jl")
include("getgravc.jl")
include("dscom.jl")
include("gstime.jl")
include("initl.jl")
include("dpper.jl")
include("dsinit.jl")
include("jday.jl")
include("invjday.jl")
include("sgp4init.jl")
include("days2mdh.jl")
include("newtonnu.jl")
include("rv2coe.jl")
include("sgp4.jl")
include("angl.jl")
include("twoline2rv.jl")
include("tleformat.jl")
include("sgp4init_vectorized.jl")
include("sgp4_vectorized.jl")
include("twoline2rv_vectorized.jl")
end
