#====================================================================#
#   Statistical Analysis and Mapping for all sites                   #                              
#====================================================================#
require(animation)
require(ggmap)
require(ggplot2)

ani.options(
  convert = shQuote('C:/Program Files/ImageMagick-6.9.0-Q16/convert.exe'), outdir = getwd()
)

year <- rep(fromYear:toYear,1)
# gif  animation function
# TTOP
draw = function(i,title, col){  
  ggmap(get_googlemap(center=c(lon = 90, lat = 35),zoom=5, maptype='terrain')) +
    geom_point(data=as.data.frame(TTOP_1980_2010[,i]),aes(x=QZ_sites$Lat,y=QZ_sites$Lon),
               size=TTOP_1980_2010[,i]+7, colour=col,alpha=0.6) +
    ggtitle(paste('TTOP ', year[i])) +
    theme(legend.position=c(80,40))
} 

saveGIF( for(i in 1:31) print(draw(i, "TTOP", "blue")), 
         movie.name = "TTOP_1980_2010+7.gif")

# Active-layer Thickness
draw1 = function(i,var, title, col){  
  ggmap(get_googlemap(center=c(lon = 90, lat = 35),zoom=5, maptype='terrain')) +
    geom_point(data=as.data.frame(var[,i]),aes(x=QZ_sites$Lat,y=QZ_sites$Lon),colour=col,
               size=var[,i]*3,alpha=0.6) +
    ggtitle(paste(title, " ", year[i]))
} 

saveGIF( for(i in 1:31) print(draw1(i,DA_1980_2010, "Active-layer Thickness (m)", "red")), 
         movie.name = "DA_1980_2010.gif")

# Thawing Depth - Freezing Depth
quartile_min <- quantile(Dt_Df_1980_2010)[1]
quartile_max <- quantile(Dt_Df_1980_2010)[5]

draw1 = function(i, title){  
  ggmap(get_googlemap(center=c(lon = 90, lat = 35),zoom=5, maptype='terrain')) +
    geom_point(data=as.data.frame(Dt_Df_1980_2010[,i]),
               aes(x=QZ_sites$Lat,y=QZ_sites$Lon, colour=Dt_Df_1980_2010[,i]),
               size=(Dt_Df_1980_2010[,i]+3)*1.5) +
    scale_colour_gradientn(
      name = "Dt - Df",
      limits = c(quartile_min, quartile_max),
      colours = c("#0000FF","#00FF00","#FF0000"),
      guide = "colourbar") +
    ggtitle(paste(title, " ", year[i]))
} 

saveGIF( for(i in 1:31) print(draw1(i, "Thawing Depth - Freezing Depth (m)")), 
         movie.name = "Dt_Df_1980_2010_color.gif")

# How to determine whether a permafrost exist
Permafrost_0_1 <- GT_ann_1980_2010 - abs(Thermal_Offset_1980_2010)
quartile_min <- quantile(Permafrost_0_1)[1]
quartile_max <- quantile(Permafrost_0_1)[5]

draw1 = function(i, title){  
  ggmap(get_googlemap(center=c(lon = 90, lat = 35),zoom=5, maptype='terrain')) +
    geom_point(data=as.data.frame(Permafrost_0_1[,i]),
               aes(x=QZ_sites$Lat,y=QZ_sites$Lon, colour=Permafrost_0_1[,i]),
               size=Permafrost_0_1[,i]+5) +
    scale_colour_gradientn(
      name = "Permafrost Exist",
      limits = c(quartile_min, quartile_max),
      colours = c("#0000FF","#00FF00","#FF0000"),
      guide = "colourbar") +
    ggtitle(paste(title, " ", year[i]))
} 

saveGIF( for(i in 1:31) print(draw1(i, "Annual Mean Ground Temperature - Thermal Offset")), 
         movie.name = "Permafrost_1980_2010_color.gif", convert = "convert")

#bm.files = sprintf("Rplot%d.png", 1:31)
#im.convert(files = bm.files, output = "Permafrost_1980_2010_color.gif")
#convert -delay 100 Rplot1.png Rplot2.png Rplot3.png Rplot4.png Rplot5.png Rplot6.png Rplot7.png 
#Rplot8.png Rplot9.png Rplot10.png Rplot11.png Rplot12.png Rplot13.png Rplot14.png Rplot15.png 
#Rplot16.png Rplot17.png Rplot18.png Rplot19.png Rplot20.png Rplot21.png Rplot22.png Rplot23.png 
#Rplot24.png Rplot25.png Rplot26.png Rplot27.png Rplot28.png Rplot29.png Rplot30.png Rplot31.png 
#Permafrost_1980_2010_color.gif
