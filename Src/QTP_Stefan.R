#====================================================================#
# Stephen solution:  compute freezing and thawing depth using        #
# accumulated ground surface degree-day total I (either the freezing #
# index DDF or thawing index DDT);                                   #
# is widely used for spatial active-layer characterisation by        #
# estimating soil properties (‘edaphic parameters’) empirically,     #
# using summer air temperature records and active-layer data         #
# obtained from representative locations                             #
#====================================================================#
# Lunardini, 1981; Nelson et al., 1997; Shiklomanov and Nelson, 2003; Zhang T. et al., 2005.
# Df: freezing depth; Dt: thawing depth; Da: Active layer depth
# W为融化时土的总含水量(%);Wu为冻土中未冻水含量(%)
# L为冰的融化热（kJ·m-3）;γ为土的干容重（kg·m-3）
# λt为融化时的导热系数(W·m-1·K-1); λf为冻结时的导热系数;
# 冻结时的导热系数大于融化时的导热系数
W <- 0.19
Wu <- 0.05
L <- 3.3e5
γ <- 1240
λt <- 1.3*86400
λf <- 1.8*86400
# 导热系数单位换算 1 W·m-1·K-1 = 1 J·cm-1·s-1·℃-1 = 1 * 60 * 60 * 24 J·m-1d-1℃
# The quantity tT/86,400 defines the ‘thawing index’ DDT (°C days), a time–temperature 
# integral usually calculated by summing mean daily temperatures above 0°C.
# The n-factors are identified based on vegetation type:
# forest land, 2.30; grassland and shrub, 1.89; grassland, 1.39; mixed area
# of meadow and steppe meadow, 1.60; and no vegetation, 2.55.
# Dt: Depth of Thaw
Dt <- sqrt((2*λt*DDTa[DDTa4:DDTa2])/(L*γ*(W-Wu)))
Df <- sqrt((2*λf*DDFa[DDFa4:DDFa2])/(L*γ*(W-Wu)))
Da <- sqrt((λt*(nt+nf)*DDTa[DDTa4:DDTa2])/(L*γ*(W-Wu)))

DA_1980_2010[position, 1:31] <- Da
Dt_Df_1980_2010[position, 1:31] <- Dt-Df

#plot(names(Da), Da, type="o", pch=0, cex=0.6, xlab="Year", ylim=c(min(Dt,Da,Df)*0.9, 
#                                                                  max(Dt,Da,Df)*1.1), ylab="Depth (m)")
#lines(names(Df), Df, type="l", col=colors[position], pch=position, cex=0.6)
#lines(names(Dt), Dt, type="o", pch=1, cex=0.6, col="red")
#lines(names(Df), Df, type="o", pch=2, cex=0.6, col="green")
#legend("topleft", legend=c("Active layer depth", "Thawing depth", "Freezing depth"), 
#       col=c("black","red", "green"), ncol=1, lty=1, pch=c(0:2), bty="n", cex=0.8)

#Da.lm <- lm(Da ~ names(Da))
#abline(Da.lm)

}
