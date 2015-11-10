
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows

# The script requires the following R packages
#require(sp)
#require(rgdal)
#require(raster)

# output all figure in postscript format
# postscript("images/EQTEC_OUTPUT_all_sites.ps", width=6, height=5)

#!!!! 用高程和经纬度计算的地温与实际的地温进行比较后，优选一种地温计算算法，
#!!!! 然后再扩展到面上。

#====================================================================#
# TTOP model:  refers to the temperature at the top of the           #
# permafrost or at the bottom of the seasonally frozen layer.        #
#====================================================================#
# Smith and Riseborough, 1996, 2002, 2012;Juliussen and Humlum, 2007; 
# Riseborough,Shiklomanov, et al., 2008;
#!!!! Snow depth and Potential Incoming Solar Radiation
QZ_sites <- read.csv("Data/QZ_sites.csv", sep = ",", header=T)

TTOP_1980_2010 <- array(NA, dim=c(length(QZ_sites$SID), 31))
DA_1980_2010 <- array(NA, dim=c(length(QZ_sites$SID), 31))
Dt_Df_1980_2010 <- array(NA, dim=c(length(QZ_sites$SID), 31))
Thermal_Offset_1980_2010 <- array(NA, dim=c(length(QZ_sites$SID), 31))
GT_ann_1980_2010 <- array(NA, dim=c(length(QZ_sites$SID), 31))
#colors <- c("black", "red", "blue", "green", "pink", "grey")
#plot(1980:2010, seq(-5, -1,length.out = 31), type="n", xlab="Year", ylab="Offset")
#legend("topleft", legend=c("Wudaoliang", "Tuotuohe", "Bange", "Anduo", "Naqu", "Dangxiong"), col=colors, 
#       ncol=3, lty=1, bty="n")

# Please define the site ID !!!
#station <- c(52707, 52737, 52818, 52908, 56004, 55279, 55294, 55299, 55493, 55578, 55591)

#station <- c(52908, 56004, 55279, 55294, 55299, 55493)
#station <- 55299
#write.table(t(c("SID", "TTOP", "Lon", "Lat", "Ele")), file="TTOP.txt", col.names= F, 
#            row.names = F, sep=',')

station <- QZ_sites$SID

for (SID in station) {

position <- which(station==SID)

#fromYear <- QZ_sites$FromYear[QZ_sites$SID == SID]
#toYear <- QZ_sites$ToYear[QZ_sites$SID == SID]
SNAME <- QZ_sites$SNameEN[QZ_sites$SID == SID]
fromYear <- 1980
toYear <- 2010
Lon <- QZ_sites$Lon[QZ_sites$SID == SID]
Lat <- QZ_sites$Lat[QZ_sites$SID == SID]
Ele <- QZ_sites$Ele[QZ_sites$SID == SID]

# 32766 is the NULL Value for Air_Data
Ground_Data <- read.table(paste("Data/QZ-ATM/E_GT",SID,".txt", sep=""), header=F, 
                          na.strings=c("32766", "NA"), 
                          col.names = c("Year","Mon","Day","evaporation", "GT_0_AVG", 
                                        "GT_0_MAX", "GT_0_MIN", "GT_5_AVG", "GT_10_AVG", 
                                        "GT_15_AVG", "GT_20_AVG", "GT_40_AVG", "GT_50_AVG", 
                                        "GT_80_AVG", "GT_160_AVG", "GT_320_AVG"))

# 888888 & 32766 are the NULL Value for Ground_Data
Air_Data <- read.table(paste("Data/QZ-ATM/TRPFUPS_daily_",SID,".dat", sep=""), 
                       header=T, na.strings=c(c("888888","32766"), "NA"))

# joint Air_Data and Ground_Data in a data frame
Air_Ground <- merge(Air_Data, Ground_Data, by = c("Year","Mon","Day"), all=T)

# P: Period (365 days)
P <- 365
# atmospheric pressure = -0.011*elevation + 100.8
# + 1 Kpa, - 0.0075℃  Temperature of the freezing point
Atm_Pressure <- -0.011*Ele + 100.8

# Ts: Temperature of the ground surface (℃)
# Ta: Temperature of the air (℃)
# Tf: Temperature of the freezing point (0℃)
Tf <- (101.3-Atm_Pressure)*0.0075/101.3
Tf <- rep(Tf, length(Air_Ground$GT))

# Kt: Thermal conductivity of ground in thawed state
# Kf: Thermal conductivity of ground in frozen state
# rk: Thermal conductivity ratio; An rk-value equal to 1 implies that MAGST = TTOP, and
# rk< 1 indicates a thermal offset. Smith and Riseborough (1996, 2002) expect values in 
# the range 0.6–0.9 for mineral soil and as low as 0.3 for organic soils, 
# depending on the amount of soil moisture. 
# Where soil moisture content is low,rk > 1 can be expected

# rk <- Kt/Kf
rk <- 0.30

# uts: Period of the thawing season for the ground surface (days)
# uta: Period of the thawing season for the air (days)
# ufs: Period of the freezing season for the ground surface (days)
# ufa: Period of the freezing season for the air (days)
# DDTa: Thawing degree-days of air (Celsius degree-days)
# DDTs: Thawing degree-days of ground surface (Celsius degree-days)
# DDFs: Freezing degree-days of ground surface (Celsius degree-days)
# DDFa: Freezing degree-days of air (Celsius degree-days) 
# 数值需要*0.1
DDFa <- 0.1*abs(tapply(Air_Ground$Temperaute[Air_Ground$Temperaute<=0], Air_Ground$Year[Air_Ground$Temperaute<=0], sum, na.rm=T))
DDTa <- 0.1*tapply(Air_Ground$Temperaute[Air_Ground$Temperaute>0], Air_Ground$Year[Air_Ground$Temperaute>0], sum, na.rm=T)

DDFs <- 0.1*abs(tapply(Air_Ground$GT[Air_Ground$GT<=0], Air_Ground$Year[Air_Ground$GT<=0], sum, na.rm=T))
DDTs <- 0.1*tapply(Air_Ground$GT[Air_Ground$GT>0], Air_Ground$Year[Air_Ground$GT>0], sum, na.rm=T)

#plot(names(DDFa), DDFa, type="o", ylim=c(min(DDFa, DDFs, DDTa, DDTs)*0.8, max(DDFa, DDFs, DDTa, DDTs)*1.2), 
#     xlab="Year", ylab="DDF&T", cex=0.6, main=SNAME)
#lines(names(DDFs), DDFs, type="o", col="red", cex=0.6)
#lines(names(DDTa), DDTa, type="o", col="green", cex=0.6)
#lines(names(DDTs), DDTs, type="o", col="blue", cex=0.6)
#legend("topleft", legend=c("DDFa", "DDFs", "DDTa", "DDTs"), col=c("black","red", "green", "blue"), ncol=2, lty=1, bty="n")

# Because DDFa, DDFs, DDTa and DDTs have different length, operation should have the same length
DDFs4 <- which(names(DDFs)==fromYear)
DDFs2 <- which(names(DDFs)==toYear)
DDFa4 <- which(names(DDFa)==fromYear)
DDFa2 <- which(names(DDFa)==toYear)

DDTs4 <- which(names(DDTs)==fromYear)
DDTs2 <- which(names(DDTs)==toYear)
DDTa4 <- which(names(DDTa)==fromYear)
DDTa2 <- which(names(DDTa)==toYear)

# nt: Thawing n-factor. Scaling factor between DDTa and DDTs
# 计算这两个值的长度，然后比较，奇异值的处理？
nt <- DDTs[DDTs4:DDTs2]/DDTa[DDTa4:DDTa2]
# nf: Freezing n-factor. Scaling factor between DDFa and DDFs
nf <- DDFs[DDFs4:DDFs2]/DDFa[DDFa4:DDFa2]
#lines(names(nt), nt, type="l", col=colors[position], pch=position, cex=0.6)

#plot(names(nt), nt, type="o", ylim=c(min(nt,nf)*0.8, max(nt,nf)*1.2), 
#     xlab="Year", ylab="nt & nf", cex=0.6, main=SNAME)
#lines(names(nf), nf, type="o", col="red", cex=0.6)
#legend("topleft", legend=c("nt: Thawing n-factor", "nf: Freezing n-factor"), 
#       col=c("black","red"), ncol=2, lty=1, bty="n")

# TTOP = MAAT + Surface_Offset + Thermal_Offset
# MAAT: Mean annual air temperature (℃)
# MAGT: mean annual ground temperature (depth unspecified) (℃)
# MAGST: Mean annual ground surface temperature (℃)
MAAT <- (DDTa[DDTa4:DDTa2]-DDFa[DDFa4:DDFa2])/P
# MAGST <- ((nt*DDTa[DDTa4:DDTa2])-(nf*DDFa[DDFa4:DDFa2]))/P
MAGST <- (DDTs[DDTs4:DDTs2]-DDFs[DDFs4:DDFs2])/P
TTOP <- ((rk*nt*DDTa[DDTa4:DDTa2])-(nf*DDFa[DDFa4:DDFa2]))/P

TTOP_1980_2010[position, 1:31] <- TTOP

#write.table(t(c(SID, TTOP[length(TTOP)], Lon, Lat, Ele)), file="TTOP.txt",
#            col.names= F, row.names = F, sep=',', append=T)

# 是否需要比较这些值的大小，然后打印显示冻土处于何种状况之中？
if(MAAT<0||MAGST<0||TTOP<0) {
  MMT_min <- min(MAAT, MAGST, TTOP)*1.2
}else 
  MMT_min <- min(MAAT, MAGST, TTOP)*0.8

# 不同的土壤有不同的导热系数
# rk=1, MASGT = TTOP; 0<rk<1, MAGST > TTOP; rk>1, MAGST < TTOP
#plot(names(MAAT), MAAT, type="o", ylim=c(MMT_min, max(MAAT, MAGST, TTOP)*1.2), 
#     xlab="Year", ylab=paste("MAAT & MAGST & TTOP"), cex=0.6, main=SNAME)
#lines(names(MAGST), MAGST, type="o", col="red", cex=0.6)
#lines(names(TTOP), TTOP, type="o", col="green", pch=position, cex=0.6)
#legend("topleft", legend=c("MAAT", "MAGST", "TTOP"), col=c("black","red", "green"), ncol=3, lty=1, bty="n")

# Surface Offset is equal to the difference between the MAGST and MAAT
Surface_Offset <- (DDFa[DDTa4:DDTa2]*(1-nf)/P)-(DDTa[DDTa4:DDTa2]*(1-nt)/P)

# Thermal Offset: the thermal offset is equal to (TTOP - MAGST).Between the ground surface 
# and the top of permafrost, heat transfer by conduction varies seasonally between frozen 
# and thawed states, since the thermal conductivity of ice is four times that of water.
Thermal_Offset <- (nt*DDTa[DDTa4:DDTa2]*(rk-1))/P
Thermal_Offset_1980_2010[position, 1:31] <- Thermal_Offset

# if Thermal offset > Annual mean ground temperature, exists permafrost
GT_ann <- 0.1*abs(tapply(Air_Ground$GT, Air_Ground$Year, mean, na.rm=T))
GT_ann4 <- which(names(GT_ann)==fromYear)
GT_ann2 <- which(names(GT_ann)==toYear)
GT_ann <- GT_ann[GT_ann4:GT_ann2]
GT_ann_1980_2010[position, 1:31] <- GT_ann

# Nival Offset: The first term (Surface_Offset) on the right-hand-side is positive, 
# and represents the elevation of MAGST over MAAT due to the insulating effect of 
# winter snow cover (nival offset) 
Nival_Offset <- DDFa[DDFa4:DDFa2]*(1-nf)/P

# vegetation offset: The second term (Surface_Offset) is negative, and represents the 
# reduction in MAGST due to vegetation effects in summer (vegetation offset)
Vegetation_Offset <- DDTa[DDTa4:DDTa2]*(1-nt)/P

# Because Vegetation_Offset exist negative value
if(Surface_Offset<0||Thermal_Offset<0||Nival_Offset<0||Vegetation_Offset<0){
  STNV_min <- min(Surface_Offset, Thermal_Offset, Nival_Offset, Vegetation_Offset)*1.2
}else 
  STNV_min <- min(Surface_Offset, Thermal_Offset, Nival_Offset, Vegetation_Offset)*0.8

#plot(names(Surface_Offset), Surface_Offset, type="o", 
#     ylim=c(STNV_min, 
#            max(Surface_Offset, Thermal_Offset, Nival_Offset, Vegetation_Offset)*1.3), 
#     xlab="Year", ylab="Different Offset", cex=0.6, main=SNAME)
#lines(names(Surface_Offset), Surface_Offset, type="o", col=colors[position], cex=0.6)
#lines(names(Thermal_Offset), Thermal_Offset, type="o", col="red", cex=0.6)
#lines(names(Nival_Offset), Nival_Offset, type="o", col="green", cex=0.6)
#lines(names(Vegetation_Offset), Vegetation_Offset, type="o", col=colors[position], cex=0.6)
#legend("topleft", legend=c("Surface Offset", "Thermal Offset", "Nival Offset", "Vegetation Offset"), 
#       col=c("black","red", "green", "blue"), ncol=2, lty=1, bty="n")

#高程和MAAT、MAGST、TTOP的关系

