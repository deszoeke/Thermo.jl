"Useful constants for atmospheric science on Earth. Uses Unitful constants."
module constants_unitful

using Unitful

#export UniversalGasConst, AvogadroNumber, BoltzmannConstant, SpeedOfLight
#export PermittivityOfVacuum, PlanckConstant
#export AirDryMolecularWeight, AirDrySpecificHeatConstPressure, AirDrySpecificHeatConstVolume
#export AirDryGasConst, AirDryThermalConductivity
#export Cp, Cpv, Cw
#export H2OMolecularWeight, H2OGasConst, LatentHeatVapor, LatentHeatFusion
#export EarthGravity, EarthAngularVelocity, EarthSunSurfaceDistance
#export EarthRadius, SolarFluxTOA
#export p0, KelvinCelsius

# physical and chemical constants
UniversalGasConst=Unitful.R # 8314.3 # J/K/kmol
AvogadroNumber=Unitful.Na # 6.022e26 # Avogadro's number
BoltzmannConstant=Unitful.k # 1.381e-23 # J/K/molecule Boltzmann's constant
SpeedOfLight=Unitful.c0 #2.998e8 # m/s
PermittivityOfVacuum=Unitful.ϵ0 # 8.85e-12  # C^-2 N^-1 m^-2 permittivity of vacuum
PlanckConstant=Unitful.h # 6.6262e-34 # J*s

# properties of dry air
AirDryMolecularWeight=28.97u"kg/kmol"
AirDrySpecificHeatConstPressure=1004.0u"J/K/kg" # Wallace and Hobbs
AirDrySpecificHeatConstVolume=717.0u"J/K/kg"
#AirDryGasConst=287.0u"J/K/kg"
AirDryThermalConductivity=2.40e-2u"J/m/s/K" # at 0C
Rd=287.04u"J/K/kg" # J/K/kg Bolton
Cp=1005.7u"J/K/kg" # J/K/kg Bolton

# properties of water
H2OMolecularWeight=18.016u"g/mol"
H2OGasConst=461u"J/K/kg"
LatentHeatVapor=2.5e6u"J/kg" # at 0 C
LatentHeatFusion=3.34e5u"K/kg"
Rv=461.5u"J/kg/K" # Bolton
# specific heats of vapor and liquid
Cpv=1870u"J/kg/K" # Bolton
Cw=4190u"J/kg/K"  # Bolton

# properties of Earth
EarthGravity=9.8u"m/s/s"            # m/s/s on gravity on Earth surface
EarthAngularVelocity=7.292e-5u"s^-1"
EarthSunSurfaceDistance=1.50e11u"m"
EarthRadius=6.37e6u"m"
SolarFluxTOA=1.38e3u"W/m^2"

# reference values
#hPa=100u"Pa"
pref=1.e5u"Pa"
#KelvinCelsius=273.15 # K
end # module constants_unitful

module constants_unitless

#export UniversalGasConst, AvogadroNumber, BoltzmannConstant
#export PermittivityOfVacuum, PlanckConstant, SpeedOfLight
#export AirDryMolecularWeight, AirDrySpecificHeatConstPressure, AirDrySpecificHeatConstVolume
#export AirDryThermalConductivity
#export Rd, Cp, Cpv, Cw
#export H2OMolecularWeight, H2OGasConst, LatentHeatVapor, LatentHeatFusion
#export EarthGravity, EarthAngularVelocity, EarthRadius
#export EarthSunSurfaceDistance, SolarFluxTOA

# physical and chemical constants
UniversalGasConst=8314.31 #J/K/kmol
AvogadroNumber=6.022e26 # Avogadro's number
BoltzmannConstant=1.381e-23 # J/K/molecule Boltzmann's constant
SpeedOfLight=2.998e8 # m/s
PermittivityOfVacuum=8.85e-12  # C^-2 N^-1 m^-2 permittivity of vacuum
PlanckConstant=6.6262e-34 # J/s

# properties of dry air
AirDryMolecularWeight=28.97 # kg/kmol
AirDrySpecificHeatConstPressure=1004.  # J/K/kg, Wallace and Hobbs
AirDrySpecificHeatConstVolume=717.   # J/K/kg
AirDryGasConst=287.   # J/K/kg
AirDryThermalConductivity=2.40e-2 # J/m/s/K at 0C
Rd=287.04 # J/K/kg Bolton
Cp=1005.7 # J/K/kg Bolton

# properties of water
H2OMolecularWeight=18.016 # kg/kmol
H2OGasConst=461        # J/K/kg
LatentHeatVapor=2.5e6  # J kg-1 at 0 C
LatentHeatFusion=3.34e5 # K kg-1
Rv=461.5 # Bolton
# specific heats of vapor and liquid
Cpv=1870 # Bolton
Cw=4190  # Bolton

# properties of Earth
EarthGravity=9.8                    # m/s/s on gravity on Earth surface
EarthAngularVelocity=7.292e-5       # s-1
EarthSunSurfaceDistance=1.50e11     # m
EarthRadius=6.37e6                  # m
SolarFluxTOA=1.38e3                 # W/m2

# reference values
pref=1.e5   # hPa
KelvinCelsius=273.15 # K
end # constants unitless
