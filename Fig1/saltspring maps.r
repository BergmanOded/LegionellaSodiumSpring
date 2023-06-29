library("ggplot2")
theme_set(theme_bw())
library("sf")
library(rnaturalearth)
library(rnaturalearthdata)
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library(ggspatial)
library(mapdata)
library(mapproj)
library(raster)
library(maps)
library(extrafont)

#copy the zip file with the data (this is a shape file, to be added as a layer to the map):
LakeZipFile <- "http://biogeo.ucdavis.edu/data/diva/wat/ISR_wat.zip"
#download with the path and file name:
download.file(LakeZipFile, destfile = "/home/oded/Documents/post/map/Lake/Lakes.zip")
unzip("/home/oded/Documents/post/map/Lake/Lakes.zip")
#check all the files areb there
list.files("/home/oded/Documents/post/map/Lake")
Lake <-sf::read_sf("/home/oded/Documents/post/map/Lake/ISR_water_areas_dcw.shp")
Lake1 <-sf::read_sf("/home/oded/Documents/post/map/Lake/ISR_water_lines_dcw.shp")
Lake
Lake1

ElevationZipFile <- "http://biogeo.ucdavis.edu/data/diva/alt/ISR_alt.zip"
#download with the path and file name:
download.file(ElevationZipFile, destfile = "/home/oded/Documents/post/map/Elevation/Elevation.zip")
unzip("/home/oded/Documents/post/map/Elevation/Elevation.zip")
#check all the files areb there
list.files("/home/oded/Documents/post/map/Elevation")
Elevationgrd <- raster("/home/oded/Documents/post/map/Elevation/ISR_alt.grd")
Elevationgrd
str(Elevationgrd)
dim(Elevationgrd)
plot(Elevationgrd)



#first, we creat a map of the world with all contries----
world <- ne_countries(scale = "medium", returnclass = "sf")
head(world)
class(world)
site_A <- data.frame(longitude = c(35.615), latitude = c(32.82)) # coordinates are converted and thus differ from old articles
#longitude + latitude = WGS84 system, can be taken from google earth
#for more then one site: data.frame(longitude = c(-80.144005, -80.109), latitude = c(26.479005, 
str(world)

# map of israel (for station A article)----
mapIl <- ggplot(data = world) +
  geom_sf() + geom_sf(data=Lake, color = 'deepskyblue3')+ #geom_sf(data = Lake1)+ #Lake1 adds river layer
  layer_spatial(Elevationgrd, aes(alpha = after_stat(band1)), fill = "darkgray") +
  scale_alpha_continuous(na.value = 0) +
  geom_rect(xmin = 35.5, xmax = 35.69, ymin = 32.69, ymax = 32.92, 
            fill = NA, colour = "black", size = 0.7) +
  coord_sf(xlim = c(34.2, 36), ylim = c(29.4, 33.5), expand = FALSE)+ 
  annotate(geom = "text", x = 34.65, y = 33.3, label = "Israel", 
           fontface = "plain", color = "black", size = 5, family='Arial') +
  annotation_scale(pad_x = unit(3.2,"cm"), pad_y = unit(0.4, "cm"))+
  theme(legend.position = "none")+
  theme(axis.ticks = element_blank())+
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

mapIl


# create map of Mediterranean sea (with Israel - for springs Legionella article)---- 
world_points<- st_point_on_surface(world)
world_points <- cbind(world, st_coordinates(st_point_on_surface(world$geometry)))
mapME <- ggplot(data = world) +
  geom_sf() + geom_sf(data=Lake, color = "black") + #geom_sf(data = Lake1)+ #Lake1 adds river layer
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "darkblue", check_overlap = TRUE, size = 3, family='Arial') + 
  #zoom
  coord_sf(xlim = c(5, 57), ylim = c(20, 50), expand = FALSE)+ 
  annotate(geom = "text", x = 31.1, y = 33, label = "Israel", 
           fontface = "bold.italic", color = "red", size = 5, family='Arial') +
  annotate(geom = "text", x = 35, y = 27, label = "Red", 
           fontface = "bold.italic", color = "red", size = 3, family='Arial') +
  annotate(geom = "text", x = 35, y = 26.5, label = "Sea", 
           fontface = "bold.italic", color = "red", size = 3, family='Arial') +
  annotate(geom = "text", x = 24.9, y = 34.25, label = "Mediterranean", 
           fontface = "bold.italic", color = "red", size = 4, family='Arial') +
  annotate(geom = "text", x = 24.9, y = 33.4, label = "Sea", 
           fontface = "bold.italic", color = "red", size = 4, family='Arial') +
  geom_rect(xmin = 34, xmax = 36, ymin = 29.5, ymax = 33.5, 
            fill = NA, colour = "black", size = 0.7) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(8, "in"), pad_y = unit(5.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  annotation_scale(pad_x = unit(0.5,"cm"), pad_y = unit(0.4, "cm"))+
  theme(axis.title = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.text = element_blank())+
  theme(text=element_text(family="Arial"))
mapME

#zoom Kinneret - springs+ station A+ Haon_Drill+ Barbutim----
(site_Tabgha <- data.frame(longitude = c(35.55), latitude = c(32.868)))
(site_Fuliya <- data.frame(longitude = c(35.538), latitude = c(32.808)))
(site_Tiberias <- data.frame(longitude = c(35.567), latitude = c(32.768)))
(Haon_Drill <- data.frame(longitude = c(35.621), latitude = c(32.724)))
(Barbutim <- data.frame(longitude = c(35.55), latitude = c(32.85)))

mapzoom_kinneret <- ggplot(data = world, fill="white") +
  geom_sf() + geom_sf(data=Lake)+ #geom_sf(data = Lake1)+ #Lake1 adds river layer
  layer_spatial(Elevationgrd, aes(alpha = stat(band1)), fill = "darkgray") +
  scale_alpha_continuous(na.value = 0) +
  geom_point(data = site_Tabgha, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Fuliya, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Tiberias, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = Haon_Drill, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = Barbutim, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  annotate(geom = "text", x = 35.522, y = 32.8685, label = "Tabgha", 
           color = "black", size = 4, family='Arial') +
  annotate(geom = "text", x = 35.513, y = 32.808, label = "Fuliya", 
           color = "black", size = 4, family='Arial') +
  annotate(geom = "text", x = 35.538, y = 32.768, label = "THS", 
           color = "black", size = 4, family='Arial') +
  annotate(geom = "text", x = 35.57, y = 32.72, label = "Haon Borehole", 
           color = "black", size = 4, family='Arial') +
  annotate(geom = "text", x = 35.517, y = 32.85, label = "Barbutim", 
           color = "black", size = 4, family='Arial') +
  annotate(geom = "text", x = 35.6, y = 32.86, label = "Sea of", 
           color = "black", size = 5.5, family='Arial') +
  annotate(geom = "text", x = 35.6, y = 32.85, label = "Galilee", 
           color = "black", size = 5.5, family='Arial') +
  coord_sf(xlim = c(35.45, 35.7), ylim = c(32.7, 32.92), expand = FALSE)+ #total map coordinates 
  scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL))+
  annotation_scale(pad_x = unit(0.3,"cm"), pad_y = unit(0.3, "cm"))+
  theme(axis.title = element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text(angle=45, hjust=1, size=8))+
  theme(legend.position = "none")+
  theme(text=element_text(family="Arial"))
mapzoom_kinneret


# zoom tabgha (springs Legionella article)----
(site_Ein7 <- data.frame(longitude = c(35.55), latitude = c(32.866)))
(site_Maayan1 <- data.frame(longitude = c(35.5497), latitude = c(32.8665)))
(site_Nur <- data.frame(longitude = c(35.5505), latitude = c(32.867)))
(site_GesherHana <- data.frame(longitude = c(35.548), latitude = c(32.8655)))
(site_Barbutim <- data.frame(longitude = c(35.5505), latitude = c(32.862)))

mapzoom_tabgha <- ggplot(data = world, fill="white") +
  geom_sf() + geom_sf(data=Lake)+ #geom_sf(data = Lake1)+ #Lake1 adds river layer
  layer_spatial(Elevationgrd, aes(alpha = stat(band1)), fill = "darkgray") +
  scale_alpha_continuous(na.value = 0) +
  geom_point(data = site_Ein7, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Maayan1, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Nur, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_GesherHana, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Barbutim, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  annotate(geom = "text", x = 35.55, y = 32.8695, label = "Tabgha stations", 
           color = "darkblue", size = 3.5, family='Arial') +
  annotate(geom = "text", x = 35.5487, y = 32.8666, label = "Nur", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5493, y = 32.8673, label = "Maayan1", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5489, y = 32.866, label = "Ein7", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5488, y = 32.8652, label = "Gesher-Hana", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5498, y = 32.862, label = "Barbutim", 
           color = "darkblue", size = 3, family='Arial') +
  coord_sf(xlim = c(35.547, 35.553), ylim = c(32.86, 32.87), expand = FALSE)+ #total map coordinates 
  scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL))+
  annotation_scale(pad_x = unit(0.3,"cm"), pad_y = unit(0.3, "cm"))+
  theme(axis.title = element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text(angle=45, hjust=1, size=8))+
  theme(legend.position = "none")+
  theme(text=element_text(family="Arial"))
mapzoom_tabgha 


# zoom Tiberias (spring Legionella article)----
(site_Tiberias_Pool <- data.frame(longitude = c(35.568), latitude = c(32.765)))
(site_Tiberias_Head <- data.frame(longitude = c(35.567), latitude = c(32.767)))
(site_Tiberias_Pump <- data.frame(longitude = c(35.5685), latitude = c(32.7645)))
mapzoom_Tiberias <- ggplot(data = world, fill="white") +
  geom_sf() + geom_sf(data=Lake)+ #geom_sf(data = Lake1)+ #Lake1 adds river layer
  layer_spatial(Elevationgrd, aes(alpha = stat(band1)), fill = "darkgray") +
  scale_alpha_continuous(na.value = 0) +
  geom_point(data = site_Tiberias_Head, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Tiberias_Pump, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Tiberias_Pool, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  annotate(geom = "text", x = 35.565, y = 32.7674, label = "Tiberias Head", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5665, y = 32.764, label = "Tiberias Pump", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.56595, y = 32.7654, label = "Tiberias Pool", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.565, y = 32.7687, label = "THS stations", 
           color = "darkblue", size = 3.5, family='Arial') +
  coord_sf(xlim = c(35.563, 35.577), ylim = c(32.762, 32.769), expand = FALSE)+ #total map coordinates 
  scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL))+
  annotation_scale(pad_x = unit(0.3,"cm"), pad_y = unit(0.3, "cm"))+
  theme(axis.title = element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text(angle=45, hjust=1, size=8))+
  theme(legend.position = "none")+
  theme(text=element_text(family="Arial"))
 # theme(plot.margin = unit(c(6,0,0,0), "lines"))

mapzoom_Tiberias

# zoom Fuliya (Legionella artile)----
(site_Fuliya1 <- data.frame(longitude = c(35.5427), latitude = c(32.8056)))
(site_Fuliya2 <- data.frame(longitude = c(35.5425), latitude = c(32.806)))
(site_Fuliya3 <- data.frame(longitude = c(35.542), latitude = c(32.8063)))
(site_Fuliya4 <- data.frame(longitude = c(35.54), latitude = c(32.807)))
(site_Fuliya5 <- data.frame(longitude = c(35.5385), latitude = c(32.8085)))
mapzoom_Fulia <- ggplot(data = world, fill="white") +
  geom_sf() + geom_sf(data=Lake)+ #geom_sf(data = Lake1)+ #Lake1 adds river layer
  layer_spatial(Elevationgrd, aes(alpha = stat(band1)), fill = "darkgray") +
  scale_alpha_continuous(na.value = 0) +
  geom_point(data = site_Fuliya1, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Fuliya2, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Fuliya3, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Fuliya4, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  geom_point(data = site_Fuliya5, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = "lightblue", stroke = 1) +
  annotate(geom = "text", x = 35.5415, y = 32.80485, label = "Fuliya 1", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5407, y = 32.8056, label = "Fuliya 2", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5405, y = 32.8063, label = "Fuliya 3", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.5385, y = 32.8071, label = "Fuliya 4", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.537, y = 32.8092, label = "Fuliya 5", 
           color = "darkblue", size = 3, family='Arial') +
  annotate(geom = "text", x = 35.541, y = 32.8138, label = "Fuliya stations", 
           color = "darkblue", size = 3.5, family='Arial') +
  coord_sf(xlim = c(35.535, 35.545), ylim = c(32.802, 32.815), expand = FALSE)+ #total map coordinates 
  scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL))+
  annotation_scale(pad_x = unit(0.3,"cm"), pad_y = unit(0.3, "cm"))+
  theme(axis.title = element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), 
        axis.text.y= element_text(angle=45, hjust=1, size=8))+
  theme(legend.position = "none")+
  theme(text=element_text(family="Arial"))

mapzoom_Fulia


library(gridExtra)
library(gridtext)
library(grid)

#map_ME_IL_Kenneret <- grid.arrange(mapME, mapIl, mapzoom_kinneret, nrow=1) 
map_ME_IL_Kenneret <- cowplot::plot_grid(mapME, mapIl, mapzoom_kinneret, nrow = 1, rel_widths = c(2.7, 1, 2))
map_stations <- cowplot::plot_grid(NULL, mapzoom_Tiberias, mapzoom_tabgha, mapzoom_Fulia, nrow = 1, rel_widths = c(0.4, 1.2, 1.5, 1.1))

png("/home/oded/Documents/post/map/spring Legionella/new/01.13.22_map_ME_IL_Kenneret.png",  width = 20000, height = 12000, res = 1200)
print(map_ME_IL_Kenneret)
dev.off()
    

png("/home/oded/Documents/post/map/spring Legionella/new/01.13.22_map_stations_final.png",  width = 20000, height = 12000, res = 1200)
print(map_stations)
dev.off()


