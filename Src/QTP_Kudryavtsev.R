#====================================================================#
# Modified Kudryavtsev model: compute the depth of seasonal freezing #
# or thawing                                                         #
#====================================================================#
##### Input Data #####
# Aa: Annual amplitude of air temperate (°C) 
# Ta (MAAT): Mean annual air temperature (°C)
# Ags: Annual amplitude of the ground surface temperature (°C)
# Tgs: Mean annual surface temperature (°C)
# Tps (TTOP): Mean annual temperature at the bottom of active layer (°C)
##### Inuput Parameter #####
# L: Volumetric latent heat (J kg-1)
# λsn: Thermal conductivity of snow  (W m-1 °C-1)
# λf: Thermal conductivity of frozen ground (W m-1 °C-1) 
# λt: Thermal conductivity of thawed ground (W m-1 °C-1)
# Wvol: Volumetric water content (fraction of 1) 
# Csn: Volumetric heat capacity of snow cover (J kg-1 °C-1)
# C: Volumetric heat capacity of soil skeleton (J kg-1 °C-1)
# ρ: Density of soil skeleton (kg m-3)
# Cf: Volumetric heat capacity of frozen ground (J kg-1 °C-1)
# Ct: Volumetric heat capacity of thawed ground (J kg-1° C-1)

# Z: the depth of seasonal freezing or thawing (m)
# 温度振幅计算；潜热？
# Air temperature amplitude
AirTemp_amp_daily <- (Air_Ground$Tmax - Air_Ground$Tmin)
Aa <- 0.1*by(AirTemp_amp_daily, Air_Ground$Year, mean)

Ta <- 0.1*by(Air_Ground$Temperaute, Air_Ground$Year, mean)

# Ground temperature amplitude
GroundTemp_amp_daily <- (Air_Ground$GT_0_MAX - Air_Ground$GT_0_MIN)
Tgs <- 0.1*by(GroundTemp_amp_daily, Air_Ground$Year, mean)

Tgs <- 0.1*by(Air_Ground$GT, Air_Ground$Year, mean)

# TTOP
# GT_ann: Mean annual surface temperature
GT_ann <- 0.1*abs(tapply(Air_Ground$GT, Air_Ground$Year, mean, na.rm=T))

numerator <- (0.5*GT_ann*(λt + λf) + 
   (Tgs*(λf - λt)/pi)*(GT_ann*asin(GT_ann/Tgs)/Tgs + 1 - sqrt(pi^2/Tgs^2)))

if (is.na(numerator)||numerator == NA||anyNA(numerator)){
  numerator < 0
  TTOP_K <- NA
} 

if(numerator < 0) {
  TTOP_K <- numerator/λf
}
if(numerator > 0){
  TTOP_K <- numerator/λt 
}

