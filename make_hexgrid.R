#make hex_grid
library(rgdal)
library(pracma)

#OUR CHICKENS ARE CALLED KAREN, ETHEL and LINDA

#point geometry shapefile(s)
GPS <- readOGR(".", "cmax_lightgeo_wgs84")

#project to equal area projection (https://spatialreference.org/ref/epsg/3035/proj4/) / https://epsg.io/3035
GPS  <- spTransform(GPS, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

#write projected shapefile (to aid validation in ArcGIS if required)
writeOGR(GPS, ".", "cmax_lightgeo_EPSG3035", driver="ESRI Shapefile", overwrite=TRUE)

#specify area of polygon (in sqkm; given we are now working in projected coordinate systems)
PArea <- 250

#FUNCTION: basic hexagon (create single hexagon given specified centroid [cx], [cy] and hexagon height [in metres])
hexagon    <- function(cx, cy, h){
  y <- 0.5 * h
  p <- y/sin(deg2rad(60))
  x <- p * cos(deg2rad(60))
  x1 <- cx
  y1 <- cy
  x2 <- x1 - x
  y2 <- y1 + y
  x3 <- x1
  y3 <- y1 + h
  x4 <- x1 + (x*2)
  y4 <- y1 + h
  x5 <- x1 + (x*3)
  y5 <- y1 + y
  x6 <- x1 + (x*2)
  y6 <- y1
  coords <- data.frame(x=c(x1, x2, x3, x4, x5, x6, x1), y=c(y1, y2, y3, y4, y5, y6, y1))
  return(coords)
}

#FUNCTION: tessellation of hexagons (using SpatialPointsDataframe and required hexagon area)
tessellate <- function(spdf, ha){

  #calculate spatial extents of final pattern
  xmin <- min(spdf@coords[,1])
  xmax <- max(spdf@coords[,1])
  ymin <- min(spdf@coords[,2])
  ymax <- max(spdf@coords[,2])
  xrange <- xmax-xmin
  yrange <- ymax-ymin
  offset <- 0.025
  xmin <- floor(xmin-(xrange*offset))
  xmax <- ceiling(xmax+(xrange*offset))
  ymin <- floor(ymin-(yrange*offset))
  ymax <- ceiling(ymax+(yrange*offset))

  #estimate dimensions of hexagon given the requested area (in km2)
  num <- 4 * ha
  den <- 6 * sqrt(3)
  h  <- (2 * sqrt(num/den)) * 1000 #convert to metres

  #calculate tessallation offsets
  yc <- h * 0.5
  p  <- yc/sin(deg2rad(60))
  xc <- p * cos(deg2rad(60))

  xv1 <- seq(xmin,xmax,xc+(xc*2))
  yv1 <- seq(ymin,ymax,h)
  yv2 <- yv1 + h/2

  #create hexagons
  k <- 1
  HexPoly <- list()
  for (i in 1:numel(xv1)){
    for (j in 1:numel(yv1)){

      if (rem(i,2)){
        ph <- hexagon(xv1[i], yv1[j], h);
      }
      else {
        ph <- hexagon(xv1[i], yv2[j], h);
      }
      Hex <- Polygon(ph)
      HexPoly[[k]] <- Polygons(list(Hex), k)
      k <- k + 1
    }
  }
  sps     <- SpatialPolygons(HexPoly)
  data    <- data.frame(ID=seq(length(sps)),row.names=seq(length(sps)))
  spolydf <- SpatialPolygonsDataFrame(sps,data)
  proj4string(spolydf) <- proj4string(spdf)

  return(spolydf)
}

#request tessellated hexagon grid by providing a point pattern in a spatial points dataframe and specified hexagon area (sq.km)
HTes <- tessellate(GPS, PArea)

#count the point data occuring within a polygonal hexagon tessallation
res    <- over(GPS, HTes)
B      <- as.data.frame(table(res$ID))
ix     <- is.element(HTes$ID, B$Var1)
HTes$NLoc[ix] <- B$Freq

#write
writeOGR(HTes, ".", "cmax_hexgrid_250sqkm_EPSG3035", driver="ESRI Shapefile", overwrite=TRUE)

#plot
plot(HTes, col=HTes@data$NLoc)