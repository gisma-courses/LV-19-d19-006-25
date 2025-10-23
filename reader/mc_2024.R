#------------------------------------------------------------------------------
# Author: creuden@gmail.com
# Description:  cleans climate data from ecowitt stations and interpolates the air temp
# Copyright:GPL (>= 3)  Date: 2024-09-28 
#------------------------------------------------------------------------------

# 0 ---- project setup ----

# load packages (if not installed please install via install.packages())
require("pacman")
# packages installing if necessary and loading
pacman::p_load(mapview, mapedit, tmap, tmaptools, raster, terra, stars, gdalcubes,
               sf,webshot, dplyr,CDSE,webshot, downloader, tidyverse,RStoolbox,
               rprojroot, exactextractr, randomForest, ranger, e1071, caret, 
               link2GI, rstac, OpenStreetMap,colorspace,ows4R,httr,
               lwgeom,readxl,highfrequency,tibble,xts,data.table,gstat)

# create a string containing the current working directory
wd=paste0(find_rstudio_root_file(),"/reader/data/")

# define time period to aggregate temp dat
time_period = 3

# multiplication factor for blowing up the Copernicus DEM
blow_fac = 15

# reference system as proj4 string for old SP package related stuff
crs = raster::crs("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
sfcrs <- st_crs("EPSG:32633")

# Copernicus DEM (https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1)
fnDTM = paste0(wd,"copernicus_DEM.tif")  

# Weather Data adapt if you download a new file 
# https://www.ecowitt.net/home/index?id=20166)
# https://www.ecowitt.net/home/index?id=149300
fn_dataFC29 = paste0(wd,"all_GW1000A-WIFIFC29.xlsx")
fn_dataDB2F =paste0(wd,"all_GW1000A-WIFIDB2F.xlsx")

# station data as derived by the field group
fn_pos_data= paste0(wd,"stations_prelim.shp")

# arbitrary plot borders just digitized for getting a limiting border of the plot area
fn_area =paste0(wd,"plot.shp")

# rds file for saving the cleaned up weather data
cleandata = paste0(wd,"climdata.RDS")

# 1 ---- read data ----
# read_sf("data/de_nuts1.gpkg") |> st_transform(crs) -> de
# read DEM data
DTM = terra::rast(fnDTM) # DTM.
# increase resolution by 15
DTM=disagg(DTM, fact=c(blow_fac, blow_fac)) 
#rename layer to altitude
names(DTM)="altitude"
r=DTM*0

# read station position data
pos=st_read(fn_pos_data)
# read station position data
area=st_read(fn_area)
# reproject the dataset to the project crs
area=st_transform(area,crs)
# read temperature data we need to skip row 1 due to excel format
clim_dataFC29 = as_tibble(readxl::read_excel(fn_dataFC29, skip = 1)) 
clim_dataDB2F = as_tibble(readxl::read_excel(fn_dataDB2F, skip = 1))

# select the required cols
tempFC29 = clim_dataFC29 %>% dplyr::select(c(1,2,32,36,40,44,48))
tempDB2F = clim_dataDB2F %>% dplyr::select(c(1,25,29,33,37,41,45,49,53))
# rename header according to the pos file names and create a merge field time
names(tempDB2F) = c("time","ch1_r","ch2_r","ch3_r","ch4_r","ch5_r","ch6_r","ch7_r","ch8_r")
names(tempFC29) = c("time","base","ch1","ch2","ch3","ch4","ch5")
#merge files
temp=merge(tempFC29,tempDB2F)
# convert datum which is a string to date format
temp$time=as.POSIXct(temp$time)
# aggregate timeslots according to the value in time_period
temp3h = highfrequency::aggregateTS(xts::as.xts(temp), alignBy = "hours",dropna = T,alignPeriod = time_period)
# add the datum colum (which is now a pointer of the timeseries) as first col in the dataset
temp_fin=as_tibble(temp3h) %>% add_column(time = zoo::index(temp3h), .before = 1)
# transpose and combine the table
temp_fin=as_tibble(cbind(nms = names(temp_fin), t(temp_fin)))
# delete first row 
names(temp_fin) = temp_fin[1,]
temp_fin=temp_fin[-1,]
# replace names specially time by stationid
names(temp_fin)[names(temp_fin) == 'time'] = 'stationid'
# extract altitudes for positions
pos$altitude= exactextractr::exact_extract(DTM,st_buffer(pos,1),"mean")
# merge positions and values via id
m=merge(pos,temp_fin)
# make the var name working for gstat by replacing all patterns
n= gsub(x = names(m),pattern = "-",replacement = "")
n= gsub(x = n,pattern = " ",replacement = "")
n= gsub(x = n,pattern = ":",replacement = "")
n= gsub(x = n,pattern = "2023",replacement = "A2023")
# and rename couse this as new names
names(m)=n
m= st_transform(m,sfcrs)

saveRDS(m,cleandata)

# grep the varnames for an interpolation loop
vars=grep(glob2rx("A2023*"), n, value = TRUE)
vars
# convert final sf vector to terra vector
temperature_vect = vect(m)
temperature_vect 
# create table containing x, y, value (A20230829220000) to interpolate this values in space
xyz=cbind(geom(temperature_vect)[,3],geom(temperature_vect)[,4],as.numeric(temperature_vect$A20230829220000))
# convert to data frame and name header
xyz=data.frame(xyz)
names(xyz) =c("x","y","temp")
xyz

# the same just for x,y
xy=cbind(geom(temperature_vect)[,3],geom(temperature_vect)[,4])
#the same just for z
z=as.numeric(temperature_vect$A20230829220000)

# -terra package
# Voronoi Segmentation
p = vect(xyz, geom=c("x", "y")) 
voronoi = voronoi(p)
v = sf::st_as_sf(p)
sf::st_crs(v)= sfcrs


# Nearest neighbor interpolation
interpNN = interpNear(r, as.matrix(xyz),radius=100)

# Inverse Distance interpolation
interpIDW = interpIDW(r, as.matrix(xyz), radius=300, power=2, smooth=1, maxPoints=3)
# ploting
plot(interpIDW)

# -gstat package
# Inverse Distance interpolation
idw_gstat <- gstat::idw(A20230829220000~1, m, st_as_stars(r),nmin = 3, maxdist = 100, idp = 2.0)
# ploting
plot(idw_gstat)



# kriging

# first convert terra spatraster to stars  and set name to altitude
dtm = setNames(st_as_stars(DTM),"altitude")

# create an autovariogram model
vm.auto = automap::autofitVariogram(formula = as.formula(paste("altitude", "~ 1")),
                                    input_data = m)
k <- krige(A20230829220000 ~ altitude, m, dtm,
           vm.auto$var_model)
plot(k)

# map it
mapview(raster(interpNN) ,col=rainbow(25)) + 
  mapview( raster(interpIDW),col=rainbow(25)) + 
  mapview(k,col=rainbow(25))+ v

## universal kriging with gstat 
# we use the terrain model for prediction


# for all time slots
for (var in vars[7:8]){
  # autofit variogramm for kriging 
  vm.auto = automap::autofitVariogram(formula = as.formula(paste("altitude", "~ 1")),
                                      input_data = m)
  plot(vm.auto)
  
  #   # kriging   
  print(paste0("kriging ", var))
  k <- gstat::krige( as.formula(paste(var, " ~ altitude")), m, dtm,
              vm.auto$var_model)
  plot(k)
  # save to geotiff
  stars::write_stars(k,paste0(wd,var,"v_interpol.tif"),overwrite=TRUE)
  
}

