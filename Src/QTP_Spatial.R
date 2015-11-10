rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows

# The script requires the following R packages
require(sp)
require(rgdal)
require(raster)
require(rgl)
require(rasterVis)
require(maptools)
require(plotrix)

# output all figure in postscript format
postscript("EQTEC_OUTPUT.ps")

#!!!! 用高程和经纬度计算的地温与实际的地温进行比较后，优选一种地温计算算法，
#!!!! 然后再扩展到面上。

#====================================================================#
# Functions: Compute the gorund surface temperature                  #
#====================================================================#

# Function: Compute average annual ground surface temperature 
# ZHANG Zhong-qiong，WU Qing-bai. 2012. 
comp_ground_temp_LLE <- function(lon, lat, elev){
  ground_temp_LLE <- 62.05 - 0.11*lon - 0.74*lat - 0.006*elev
  return (ground_temp_LLE)
}

# Function: Compute average annual ground surface temperature with slope aspect
# LU Jiahao, NIU Fujun, Cheng Hua,et al. 2013.
comp_ground_temp_LEA <- function(lat, elev, asp){
  ground_temp_LEA <- 65.461 - 1.222*lat - 0.005*elev - 0.299*cos(asp)
  return (ground_temp_LEA)
}

#====================================================================#
# The project (highway, railway, power line, bridge)                 #
#====================================================================#

# use ogrinfo to get the layer after install gdal-bin
QZ_highway <- readOGR("data/Highway.kml", "Highway")
highway <- list("sp.lines", as(QZ_highway, "SpatialLines"), col="white")

QZ_railway <- readOGR("data/Railway.kml", "Railway")
railway <- list("sp.lines", as(QZ_railway, "SpatialLines"), col="blue")

QZ_powerline <- readOGR("data/PowerGrid.kml", "PowerGrid")
powerline <- list("sp.lines", as(QZ_powerline, "SpatialLines"), col="red")

#====================================================================#
# Compute the gorund surface temperature                             #
#====================================================================#

# "SpatialGridDataFrame"
study_area_dem <- readGDAL("data/srtm_studyarea_dem.tif")
lat_lon <- coordinates(study_area_dem)
elevation <- study_area_dem@data[[1]]

# compute the gorund surface temperature with lon, lat and elevlation
ground_temp_LLE <- comp_ground_temp_LLE(lat_lon[,1], lat_lon[,2], elevation)
ground_temp_data = data.frame(x=lat_lon[,1], y=lat_lon[,2], col=ground_temp_LLE)
ground_temp_LLE_grid = SpatialPixelsDataFrame(points=ground_temp_data[,1:2], 
                                              data=data.frame(ground_temp_LLE))
proj4string(ground_temp_LLE_grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs 
                                         +ellps=WGS84 +towgs84=0,0,0")

# Output the tif map with the gorund surface temperature
writeGDAL(ground_temp_LLE_grid, fname = "data/annualGroundTemp_LLE.tif", 
          drivername = "GTiff", type = "Float32")

# summary the gorund surface temperature
study_area_groundtemp_LLE <- readGDAL("data/annualGroundTemp_LLE.tif")
summary(study_area_groundtemp_LLE)

# compute the slope aspect
dem <- raster("data/srtm_studyarea_dem.tif")
aspectdem = terrain(dem, opt='aspect', unit='degrees')
slopedem = terrain(dem, opt='slope', unit='degrees')

# Output the tif map with the aspect
aspectdem_grid <- as(aspectdem, "SpatialPixelsDataFrame")
proj4string(aspectdem_grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs 
                                   +ellps=WGS84 +towgs84=0,0,0")
writeGDAL(aspectdem_grid, fname = "data/study_area_aspect.tif", 
          drivername = "GTiff", type = "Float32")

aspectdem_data <- readGDAL("data/study_area_aspect.tif")

# compute the gorund surface temperature with lat, elevlation and aspect
ground_temp_LEA <- comp_ground_temp_LEA(lat_lon[,2], elevation, aspectdem)
ground_temp_LEA_grid <- as(ground_temp_LEA, "SpatialPixelsDataFrame")
proj4string(ground_temp_LEA_grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs 
                                         +ellps=WGS84 +towgs84=0,0,0")

# Output the tif map with the gorund surface temperature
writeGDAL(ground_temp_LEA_grid, fname = "data/annualGroundTemp_LEA.tif", 
          drivername = "GTiff", type = "Float32")

# summary the gorund surface temperature
study_area_groundtemp_LEA <- readGDAL("data/annualGroundTemp_LEA.tif")
summary(study_area_groundtemp_LEA)

# figure attribution
trellis.par.set(sp.theme())
north <- list("SpatialPolygonsRescale", layout.north.arrow(),
              offset = c(92.85000, 35.25000), scale = 0.04)
scale <- list("SpatialPolygonsRescale", layout.scale.bar(),
              offset = c(93.15000, 34.68000), scale = 0.1, 
              fill = c("transparent", "black"))
txt1 <- list("sp.text", c(93.15000, 34.69500), "0")
txt2 <- list("sp.text", c(93.25000, 34.69500), "1 km")
dem.layout <- list(north, scale, txt1, txt2)
#dem.layout <- list(highway, railway, powerline, north, scale, txt1, txt2)

# plot dem
ele <- spplot(study_area_dem, scales=list(draw=T), sp.layout = dem.layout, 
              main="Elevation of Qinghai-Tibet Highway (m)")

plot(ele)

#legend("bottomright", legend=c("High Way", "Railway", "Power Line"),
#       title="Legend", bty="n", inset=0.05,
#       lty=c( 1,-1,-1), pch=c(-1,15, 1),
#       col=c("white", "blue", "red"))

# plot 3D dem
plot3D(dem)

# plot aspect
print(spplot(aspectdem_data, scales=list(draw=T), sp.layout = dem.layout, 
             main=expression("Aspect of Qinghai-Tibet Highway ( " * degree * " )")))

# plot ground surface temperature
print(spplot(study_area_groundtemp_LLE, scales=list(draw=T), sp.layout = dem.layout, 
             ylab=expression('Ground temperature (' * degree * 'C)'), main="LLE"))
print(spplot(study_area_groundtemp_LEA, scales=list(draw=T), sp.layout = dem.layout, 
             ylab=expression('Ground temperature (' * degree * 'C)'), main="LEA"))

# Boxplots
opar<-par(cex=1,mfrow = c(1, 2))
boxplot(ground_temp_LLE, ylim=c(-6,1), 
        ylab=expression('Ground temperature (' * degree * 'C)'), main="LLE")
boxplot(ground_temp_LEA, ylim=c(-6,1), 
        ylab=expression('Ground temperature (' * degree * 'C)'), main="LEA")
par(opar)

# compute highway Roadbed Trend (azimuth angle)
highway <- raster("data/QZ_highway_raster.tif")
highway_pix <- as(highway, "SpatialPixelsDataFrame")
highway_azimuth <- trackAzimuth(highway_pix@coords, type="snyder_sphere")

# plot highway trend
V <- length(highway_azimuth)
for (l in 1:V){
  if (highway_azimuth[l] < 0){
    highway_azimuth[l] <- -highway_azimuth[l] + 180
  }
}

highway_azimuth[highway_azimuth[1:length(highway_azimuth)] < 0] <- -highway_azimuth + 180 
opar<-par(cex=1,mfrow = c(1, 2))
polar.plot(highway_azimuth/36, names(highway_azimuth), start=90,clockwise=TRUE,
           main=expression('The Azimuth of Qinghai-Tibet Highway (' * degree * ')'),
           lwd=3,line.col=4) 
boxplot(highway_azimuth)
par(opar)

# write KML format for WebGIS
writeOGR(ground_temp_LLE_grid, "data/annualGroundTemp_LLE.kml", 
         layer = "LLE", driver = "KML", overwrite_layer = TRUE)
writeOGR(ground_temp_LEA_grid, "data/annualGroundTemp_LEA.kml", 
         layer = "LEA", driver = "KML", overwrite_layer = TRUE)
writeOGR(aspectdem_grid, "data/study_area_aspect.kml", 
         layer = "aspect", driver = "KML", overwrite_layer = TRUE)

#====================================================================#
# Define the land cover type (Plant Function Type)                   #
# Please refer to Globcover2009_Legend.csv                           #
#====================================================================#

study_area_LUCC <- readGDAL("data/LUCC.tif")
lat_lon_lucc <- coordinates(study_area_LUCC)
landcover <- study_area_LUCC@data[[1]]

# 植被状况对灾害的影响？
lucc <- raster("data/LUCC.tif")
# NA, 11, 14, 20, 30, 140, 150, 200, 210 were filled.

length(values(lucc))

#====================================================================#
# Comopute the critical embamkment height design value               #
#====================================================================#


# ice content
# highway Roadbed Height
# Highway Road runtime
# highway Roadbed Construction (engineering measure)

dev.off()

#====================================================================#
# 1. 考虑LUCC、DEM、土壤质地、地温、气温、公路路基、运行时间的灾害分析
# 2. 假设情景，以突发事件作为驱动，如热融滑塌、泥石流、滑坡等；
# 3. 以IPCC未来100年气候预测结果作为气候驱动来模拟
#====================================================================#


QZ_highway <- readShapeSpatial("data/青藏公路/QT_road.shp")
QZ_railway <- readShapeSpatial("data/青藏公路/QZ_railway.shp")
ggmap(get_googlemap(center=c(lon = 90, lat = 35),zoom=5, maptype='terrain')) +
  geom_line(QZ_highway)+
  geom_line(QZ_railway)