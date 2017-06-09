using Thermo
using Unitful
using Unitful: K, °C, hPa, Pa
using Base.Test

T=300 # K
p=800 # hPa

# water vapor saturation methods
Lv()
Lv(T-273.15)
Lv(T*u"K")
Lv((300:303)u"K")
Lv([300 301]u"K")

@show es(T-273.15,p)
@show es(T*u"K",p*u"hPa")
# weird mixed units/unitless
@show es(273.15u"K",800)
@show es(0.0,800u"hPa")
@show es(T*u"K",p)
@show es((T-273.15)u"°C",p*u"hPa")
@show es((T-273.15),p*u"hPa")

@show qs(T-273.15,p)
@show qs(T*u"K",p*u"hPa")
@show qs(°C.(T*K),p*u"hPa")

# potential temperature methods
# nondimensional methods
@show theta(T,p)
@show theta(T,p,1000)
# weird cases with mixed unitless and units should be errors
@show theta(T*u"K",p)
@show theta((T-273.15)*u"°C",p)
#
@show theta(T,p*u"hPa")
@show theta(T,p*100*u"Pa")
@show theta(T*u"K",p*u"hPa")
@show theta(T*u"K",p*u"hPa")
@show theta(T*u"K",p*100*u"Pa")
@show theta((T-273.15)*u"°C",p*u"hPa")
@show theta((T-273.15)*u"°C",p*u"hPa")
@show theta((T-273.15)*u"°C",p*100*u"Pa")
# unit conversions are easy with Unitful; dispatch is hard

# test broadcasting arrays
T=collect(273:3:303)
@show theta(T,p)
p=collect(100:100:1000)
nd=min(length(p),length(T))
T,p = T[1:nd],p[1:nd]
@show theta(T,p)
@show theta(T,p.')
