#====================================================================#
# Spatial Stefan, please reference the QTP_Stefan.R                  #
# Stephen solution:  compute freezing and thawing depth using        #
# accumulated ground surface degree-day total I (either the freezing #
# index DDF or thawing index DDT)                                    #
#====================================================================#
# Date: 06/10/2015
# Author: Lihui Luo, luolh@lzb.ac.cn

#---------------------------------------------------------------------------------------------------
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

#---------------------------------------------------------------------------------------------------
library(RNetCDF)

# Df: freezing depth; Dt: thawing depth; Da: Active layer depth
# W为融化时土的总含水量(%);Wu为冻土中未冻水含量(%)
# L为冰的融化热（kJ·m-3）;γ为土的干容重（kg·m-3）
# λt为融化时的导热系数(W·m-1·K-1); λf为冻结时的导热系数;
# 冻结时的导热系数大于融化时的导热系数
# 不同的植被和土壤类型参数差异化
W <- 0.19
Wu <- 0.05
L <- 3.3e5
γ <- 1240
λt <- 1.2*86400
λf <- 1.5*86400

stefan_nc <- create.nc("Data/QTEC_Stefan.nc")
dim.def.nc(stefan_nc, "longitude", 51)
dim.def.nc(stefan_nc, "latitude", 74)
dim.def.nc(stefan_nc, "time", unlim=TRUE)

var.def.nc(stefan_nc, "longitude", "NC_DOUBLE", "longitude")
var.def.nc(stefan_nc, "latitude", "NC_DOUBLE", "latitude")
var.def.nc(stefan_nc, "time", "NC_DOUBLE", "time")

var.def.nc(stefan_nc, "Df", "NC_DOUBLE", c("longitude", "latitude", "time"))
var.def.nc(stefan_nc, "Dt", "NC_DOUBLE", c("longitude", "latitude", "time"))

att.put.nc(stefan_nc, "longitude", "long_name", "NC_CHAR", "longitude")
att.put.nc(stefan_nc, "longitude", "units", "NC_CHAR", "degrees_east")
att.put.nc(stefan_nc, "longitude", "axis", "NC_CHAR", "X")

att.put.nc(stefan_nc, "latitude", "long_name", "NC_CHAR", "longitude")
att.put.nc(stefan_nc, "latitude", "units", "NC_CHAR", "degrees_north")
att.put.nc(stefan_nc, "latitude", "axis", "NC_CHAR", "Y")

att.put.nc(stefan_nc, "time", "long_name", "NC_CHAR", "time")
att.put.nc(stefan_nc, "time", "units", "NC_CHAR", "hours since 1900-01-01 00:00:00")
att.put.nc(stefan_nc, "time", "calendar", "NC_CHAR", "standard")

#att.put.nc(stefan_nc, "Df", "_FillValue", "NC_FLOAT", "-3.4e+38")
#att.put.nc(stefan_nc, "Dt", "_FillValue", "NC_FLOAT", "-3.4e+38")

att.put.nc(stefan_nc, "Df", "missing_value", "NC_FLOAT", "-9999")
att.put.nc(stefan_nc, "Dt", "missing_value", "NC_FLOAT", "-9999")

att.put.nc(stefan_nc, "Df", "projection", "NC_CHAR", "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
att.put.nc(stefan_nc, "Dt", "projection", "NC_CHAR", "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

att.put.nc(stefan_nc, "Df", "projection_format", "NC_CHAR", "PROJ.4")
att.put.nc(stefan_nc, "Dt", "projection_format", "NC_CHAR", "PROJ.4")

QTEC_Forcing_Att <- open.nc("Data/QTEC_buffer30km_Forcing_Daymean/QTEC_buffer30_Forcing_1980_Daymean.nc", 
                                write=FALSE)
year <- c( 710037, 718797, 727557, 736317, 745101, 753861, 762621, 
           771381, 780165, 788925, 797685, 806445, 815229, 823989, 832749, 841509, 
           850293, 859053, 867813, 876573, 885357, 894117, 902877, 911637, 920421, 
           929181, 937941, 946701, 955485, 964245, 973005)

latitude <- var.get.nc(QTEC_Forcing_Att, "latitude")
longitude <- var.get.nc(QTEC_Forcing_Att, "longitude")

var.put.nc(stefan_nc, "latitude", latitude)
var.put.nc(stefan_nc, "longitude", longitude)
var.put.nc(stefan_nc, "time", year)

close.nc(QTEC_Forcing_Att)

Years <- c(1980:2010)
for (Y in Years) {
  QTEC_Forcing_Daymean <- open.nc(paste("Data/QTEC_buffer30km_Forcing_Daymean/QTEC_buffer30_Forcing_", 
                                        Y ,"_Daymean.nc", sep=""), write=FALSE)
  ALL_Days <- dim.inq.nc(QTEC_Forcing_Daymean,"time")$length
  print(ALL_Days)
  AIR_TEMP <- var.get.nc(QTEC_Forcing_Daymean, "temp")
  Air_Temp_less_0 <- AIR_TEMP[,,1:ALL_Days]-273.15
  Air_Temp_upper_0 <- AIR_TEMP[,,1:ALL_Days]-273.15
  
  Air_Temp_less_0[Air_Temp_less_0[,,1:ALL_Days] > 0] =0
  Air_Temp_upper_0[Air_Temp_upper_0[,,1:ALL_Days] < 0] =0
  
  DDTa <- rowSums(Air_Temp_upper_0, dims=2)
  DDFa <- rowSums(Air_Temp_less_0, dims=2)
  
  # DDFa <- 0.1*abs(tapply(Air_Ground$Temperaute[Air_Ground$Temperaute<=0], 
  #         Air_Ground$Year[Air_Ground$Temperaute<=0], sum, na.rm=T))
  # DDTa <- 0.1*tapply(Air_Ground$Temperaute[Air_Ground$Temperaute>0], 
  #         Air_Ground$Year[Air_Ground$Temperaute>0], sum, na.rm=T)
  
  Dt <- sqrt((10*λt*DDTa)/(L*γ*(W-Wu)))
  Df <- sqrt((-0.5*λf*DDFa)/(L*γ*(W-Wu)))
  # Df[Df(,,)==NA] = -3.4e+38
  # Dt[Dt(,,)==NA] = -3.4e+38
  
  YY <- which(Years==Y)
  print(YY)
  var.put.nc(stefan_nc, "Df", Df, c(1, 1, YY), c(51, 74, 1), na.mode=2)
  var.put.nc(stefan_nc, "Dt", Dt, c(1, 1, YY), c(51, 74, 1), na.mode=2)

  sync.nc(stefan_nc)
  close.nc(QTEC_Forcing_Daymean)
}

# Global attribution
att.put.nc(stefan_nc, "NC_GLOBAL", "Title", "NC_CHAR", "Compute Spatial Stefan from QTEC_Forcing")
att.put.nc(stefan_nc, "NC_GLOBAL", "Author", "NC_CHAR", "Lihui Luo email:luolh@lzb.ac.cn")

close.nc(stefan_nc)
