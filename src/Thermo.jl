module Thermo

using Unitful
using Unitful: K, °C, Pa, hPa
using Unitful: NoUnits, Temperature, Pressure, Amount
#using Unitful: Length, Mass
#using Unitful: Units
#=
# extend Unitful.unit function to be able to find the units of arrays
Unitful.unit(a::AbstractArray)=Unitful.unit(a[1])
Unitful.dimension(a::AbstractArray)=Unitful.dimension(a[1])
=#

# Practical extensions to Unitful:
# Define a nondimensional number to include nondimensional unitless builtin Number types but not dimensional Quantities.
# Define a union of NondimensionalNumber and DimensionlessQuantity types.
# Includes dimensionless ratios of Pressure, Mass, Temperature quantities.
# It's easier to dispatch on Unions than typeintersections.
typealias NondimensionalNumber Union{Real, Complex} # not including Quantities
typealias NondimensionalQuantity Union{NondimensionalNumber, DimensionlessQuantity}
NdN = NondimensionalNumber # for short
NdQ = NondimensionalQuantity

# Extend uconvert to facilitate consistent interfaces for functions of Unitful and/or NondimensionalNumber inputs
# by converting to expected units, dimensional or nondimensional, for all inputs as they come in. E.g:
# most specific method calculates:
#   theta(T,::Temperature, p::Pressure, p0::Pressure, kappa::NondimensionalNumber) = K(T)/Exner(p, p0, kappa)
# one catch-all method handles all unit conversions:
#   theta(T, p, p0, kappa) = theta(uconvert(K,T), uconvert(hPa,p), uconvert(hPa,p0), kappa)
"""
ucon(args...)=Unitful.uconvert(args...)
ucon(a::Units, x::NondimensionalNumber) = a*x
extends and shadows Unitful.uconvert to add units to NondimensionalNumbers
x::NondimensionalNumber assumed to be in a::Units.
examples:
ucon(K, x*°C)      # Unitful.uconvert standard method
ucon(K, x+273.15)  # extended to NondimensionalNumbers
"""
ucon(args...)=Unitful.uconvert(args...)
ucon(a::Unitful.Units, x::NondimensionalNumber) = a*x
ucon{T<:NondimensionalNumber}(a::Unitful.Units, x::AbstractArray{T}) = ucon.(a,x)
#ucon(a::Unitful.Units, x::NondimensionalNumber, sourceunits::Unitful.Units) = ucon(a, x*sourceunits) # superfluous

"""
Nondimensional unit convert
undconvert(a::Unitful.Units, x::Unitful.Quantity) = NoUnits(ucon(a,x)/a)
undconvert(a::Unitful.Units, x::NondimensionalNumber) = x
outputs nondimensional number assumed in the requested units. Assumes nondimensional x already is in units of a.
"""
undconvert(a::Unitful.Units, x::Unitful.Quantity) = NoUnits(ucon(a,x)/a)
undconvert(a::Unitful.Units, x::NondimensionalNumber) = x # == NoUnits(ucon(a,a*x)/a)

export theta, ev,qv,es,qs, Tlcl, Lv, dqsdT, moistad

include("ThermoConstants.jl")
using .constants_unitful
using .constants_unitful: pref, Rd, Rv, Cp, Cpv, Cw
using .constants_unitful: EarthGravity

"T/Θ = Exner(p, p0=1000hPa, kappa=Rd/Cp)  Exner function"
Exner(p::NdN  , p0::NdN, kappa=NoUnits(Rd/Cp)) = (p/p0)^kappa # both unitless
Exner(p::Pressure, p0::Pressure, kappa=NoUnits(Rd/Cp)) = NoUnits(p/p0)^kappa # both pressure, but may be of different units
#Exner(p::Pressure, p0::Pressure=hPa(pref), kappa=NoUnits(Rd/Cp)) = NoUnits(p/p0)^kappa # define intersecting specific case to avoid method ambiguity
# weird cases where one of p,p0 has units but the other doesn't probably should be errors, but here we convert them to have same units as the unitful one:
Exner(p::NdN, p0::Pressure=hPa(pref), kappa=NoUnits(Rd/Cp)) = (p/NoUnits(p0/unit(p0)))^kappa
Exner(p::Pressure, p0::NdN=NoUnits(pref/hPa), kappa=NoUnits(Rd/Cp)) = (NoUnits(p/unit(p))/p0)^kappa
# extend to broadcasting arrays, iterables with .
Exner(p::AbstractArray, p0, kappa)=Exner.(p, p0, kappa)

"Exner exponent kappa(qv[kg/kg])=(Rd/Cp)*(1-0.28*qv) varies with the specific humidity."
kappa(qv::NdQ)=NoUnits(Rd/Cp*(1-0.28*qv))

"""
theta(T,p, p0=1000hPa, kappa=Rd/Cp)
calculates potential temperature from
temperature T (default K), pressure p (default hPa).
Quantities with other compatible units may be specified
by the Unitful module.
"""
# unitless input assumes T in K and returns unitless T values corresponding to K
theta(T::NdN,      p, p0=hPa(pref), kappa=NoUnits(Rd/Cp)) =   T /Exner(p, p0,kappa)
theta(T::Temperature, p, p0=hPa(pref), kappa=NoUnits(Rd/Cp)) = K(T)/Exner(p, p0,kappa)
# output unitful quantities with units K, automatically converts Celsius to K
theta(T::AbstractArray,p,p0=hPa(pref),kappa=NoUnits(Rd/Cp))=theta.(T,p,p0,kappa) # extend to arrays
# theta(T,p::AbstractArray,p0=hPa(pref),kappa=NoUnits(Rd/Cp))=theta.(T,p,p0,kappa)
# multiple dispatch of p is handed by Exner

"""
qv(p/ev)
qv(p::Pressure,ev::Pressure)
qv(p::Amount,ev::Amount)
specific humidity [kg/kg]
"""
#qv(p::Pressure,ev::Pressure) = NoUnits(Rd/(Rv*(p/ev + (Cp/Rv-1.))))
#qv(p::Amount,ev::Amount) = NoUnits(Rd/(Rv*(p/ev + (Cp/Rv-1.))))
qv(poev::Number) = ucon(u"kg/kg",Rd/(Rv*(poev + (Cp/Rv-1.))))
qv(p::Number,ev::Number) = qv(p/ev)
#qv(p::Union{AbstractArray,Number},ev::Union{AbstractArray,Number})=qv.(p,ev) # handle abstractarrays

"""
es(T,p) [hPa] computes saturation vapor pressure based on Wexler's formula,
with enhancement factor for moist air rather than water vapor.
The enhancement factor requires a pressure.
T [degrees C], p [hPa] (note the reversed input order), es [hPa]
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
# procedure: convert to expected units, nondimensionalize, calculate
# methods "diagonally dispatch" nondimensionalizing one argument at a time
es(T::NdN,p::NdN=NoUnits(pref/hPa)) = 6.1121hPa*(1.0007 + 3.46e-8*p).*exp((17.502*T)./(240.97 + T))
es(T::Temperature,p::Pressure=pref) = es(NoUnits(°C(T)/°C),NoUnits(p/hPa)) # resolves ambiguity below
# weird corner cases probably should be errors, but are allowed for now:
es(T::Temperature,p::NdN) = es(NoUnits(°C(T)/°C),p)
es(T::NdN,   p::Pressure) = es(T,NoUnits(p/hPa))
# NoUnits form expects T in °C, p in hPa
# uses Buck 1981 numerical constants
#es(T::AbstractArray,p               )=es.(T,p) # broadcasts
#es(T::AbstractArray,p::AbstractArray)=es.(T,p)
#es(T               ,p::AbstractArray)=es.(T,p)

"""
qs(T[°C],p[hPa]) [kg/kg]
Saturation specific humidity based on Wexler's formula for es
with enhancement factor. See es(T,p).
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
qs(T::Temperature,p::Pressure=pref)=qv(p,es(T,p))
qs(T::NdN,p::NdN=NoUnits(pref/hPa))=qv(p,NoUnits(es(T,p)/hPa))
#qs(T::AbstractArray,p::AbstractArray)=qs.(T,p)  # broadcast
#qs(T::AbstractArray,p)=qs.(T,p)
#qs(T,p::AbstractArray)=qs.(T,p)

"ev(p,qv) vapor pressure in units of p"
ev(p::Pressure,qv::NdQ) = p/NoUnits(Rd/(Rv*qv) +1)
ev(p::NdN,qv::NdQ) = p/NoUnits(Rd/(Rv*qv) +1)
#ev(p,qv)=ev.(p,qv)
#  ev = p.*qv./(Rd/Rv + qv)
#  Tl = 2840./(3.5*log(T) - log(.01*ev) - 4.805) + 55

"""
moistad(T[°C], p[hPa])
Moist adaibatic lapse rate of _potential_ temperature [K/m]
from Wood and Bretherton 2006
(also Rogers and Yau, but with mixing ratio instead of specific humidity).
(c) Simon de Szoeke, 2011
"""
moistad(T::Temperature,p=pref)=moistad(NoUnits(°C(T)/°C), p)
moistad(T::NdN,p::Pressure=pref)=moistad(T, NoUnits(p/hPa))
function moistad(T::NdN,p::NdN=NoUnits(pref/hPa))
    L=Lv(T)
    qsat=qs(T,p)
    (EarthGravity/Cp) * (1 - (1+L*qsat/(Rd*(T+273.15))) / (1+((L/(T+273.15))^2 *qsat/(Cp*Rv))) )
end

"Lv(T) [J/kg] Latent heat of vaporization of water as a function of temperature [°C]."
Lv(T::NdN=0.0) = (2.501e6 + (Cpv-Cw)*T/unit(Cpv))*u"J/kg"
Lv(T::Temperature=0.0°C) = (2.501e6u"J/kg" + (Cpv-Cw)*°C(T)) |> u"J/kg" # more specific method called first
Lv{Y<:Union{Range,AbstractArray}}(T::Y)=Lv.(T)
# Lv(T::AbstractArray)=Lv.(T)

"dqsdT(T[K],p[hPa]) = dqs/dT = qs(p,T)*Lv(T)/(Rv*T^2) [1/K] using Wexler and Clausius-Clapeyron formulas"
dqsdT(T::Temperature,p::Pressure=pref) = qs(T,p)*Lv(T)/(Rv*(K(T)^2))
dqsdT(T::Temperature,p::NdN) = dqsdT(T,p*hPa)
dqsdT(T::NdN,p::Pressure=pref) = dqsdT(T*K,p)
dqsdT(T::NdN,p::NdN) = dqsdT(T*K,p*hPa)
dqsdT(T::AbstractArray,p=pref) = dqsdT.(T,p) # broadcasts
dqsdT(T,p::AbstractArray) = dqsdT.(T,p)


"""
Tlcl(T[K], p[hPa], qv[kg/kg])
Tlcl(T[K], p[hPa], ev[Pa])
Temperature at the LCL [K]. From Bolton, 1980, MWR, 108, 1046-1053.
"""
Tlcl(T::Temperature,p,qv) = Tlcl(NoUnits(K(T)/K),p,vapor)
Tlcl(T::NdN,p::Pressure,qv) = Tlcl(T,NoUnits(p/Pa),vapor)
Tlcl(T::NdN,p::NdN,ev::Pressure) = (2840./(3.5*log(T) - log(NoUnits(ev/hPa)           ) - 4.805) + 55.)*K
Tlcl(T::NdN,p::NdN,qv::NdQ) =   (2840./(3.5*log(T) - log(p/(Rd/(Rv*NoUnits(qv)) +1)) - 4.805) + 55.)*K
Tlcl(T,p,q)=Tlcl.(T,p,q)

"theta_e(T[K],p[hPa],qv[kg/kg]) from eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053."
theta_e(T::Temperature,p::Pressure,qv::NdQ) = T/Exner(p,p0,kappa(qv)) * exp((3376./Tlcl(T,p,qv) - 2.54)*qv*(1 + 0.81*qv))
theta_e(T,p,qv) = theta_e(ucon(K,T),ucon(hPa,p),NoUnits(ucon(u"kg/kg",qv)))

include("theta_w.jl")

"Blackbody radiation bb(T)"
bb(t::Temperature)=ucon(u"W/m^2", Unitful.σ*K(t)^4) # specific method calculates
bb(t)=bb(ucon(K, t)) # wrapper converts arguments' units

end
