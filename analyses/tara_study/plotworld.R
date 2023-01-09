library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
require(sf)

# options(stringsAsFactors = F)         # no automatic data transformation
# options("scipen" = 100, "digits" = 4) # suppress math annotation
# op <- options(gvis.plot.tag='chart')  # set gViz options
# install libraries
# install.packages(c("OpenStreetMap", "DT", "RColorBrewer", "mapproj", "sf", "RgoogleMaps", 
#                    "scales", "rworldmap", "maps", "tidyverse", "rnaturalearth", 
#                    "rnaturalearthdata", "rgeos", "ggspatial", "maptools", "leaflet", "sf", 
#                    "tmap", "here", "rgdal", "scales"))
# # install package from github
# devtools::install_github("dkahle/ggmap", ref = "tidyup")
# 




## MYD AD
study="Tara"
input_dir =  paste0('~/Documents/FEASTX/',study,'Files/')
code_dir = '~/Documents/FEASTX/code/pipeline/'


# load data
require(rgeos)
world <- ne_countries(scale = "medium", returnclass = "sf")

meta_by_accession = read.csv(paste0(input_dir,"metadata_merge.csv"),stringsAsFactors = FALSE)
unique_meta = meta_by_accession[!duplicated(meta_by_accession$station), ]
unique_meta$longitude = unique_meta$Longitude_East
unique_meta$latitude = unique_meta$Latitude_North
unique_meta$coord = paste0(unique_meta$longitude,",", unique_meta$latitude)
replace_number = as.integer(gsub('TARA_','',unique_meta$station))
replace_number[is.na(replace_number)] = "148b"
unique_meta$station_number = replace_number
unique(unique_meta$sea_area)
# gene world map
library(rworldmap)
# get map
# worldmap <- getMap(resolution = "coarse")
# plot(worldmap, col = "lightgrey", 
#      fill = T, border = NA,
#      xlim = c(-180, 180), ylim = c(-90, 90),
#      bg = "aliceblue",
#      asp = 1, wrap=c(-180,180)) 
# 
# points(unique_meta$rounded_longitude, unique_meta$rounded_latitude, 
#        col = "red", cex = .1)


require(wesanderson)

mypalette = c("#2c6192","#fd954f","#6916A9","#94c95e","#80abcb","#c74e50","#008043","#000000","#b0aeae")
names(mypalette) = c("NPO","SPO","NAO","SAO","MS","RS","IO","SO","Unknown")
barplot(c(1:length(mypalette)), axes=FALSE, col=mypalette)

mylightpalette = c("#95B7D6","#9ACAB3","#EBEBEB","#B3DBFF","#E3B7F5","#FA8D8F","#CEE8B4","#FCD8C0","#969696")
names(mylightpalette) = c("IO","NAO","Unknown","MS","NPO","RS","SAO","SPO","SO")

unique_meta$sea_area = factor(unique_meta$sea_area,levels = c("NPO","SPO","NAO","SAO","MS","RS","IO","SO","Unknown"))

require(ggrepel)
# ggplot(unique_meta, aes(x=rounded_longitude, y= rounded_latitude,color=sea_area)) +   
#   borders("world", colour=NA, fill="lightgrey")  +
#   geom_point(alpha = 1, size = 1) +
#   scale_fill_manual(values = mypalette) + scale_color_manual(values = mypalette) + 
#   theme(panel.background = element_rect(fill = "azure1", colour = "azure1")) +
#   geom_label(aes(x=rounded_longitude, y= rounded_latitude, label=station_number),label.size = 0.25,size=3)
# 
single_ocean_label = data.frame(ocean = c("Indian Ocean", "Red Sea", "Mediterreanean Sea", "South Atlantic Ocean",
                                          "North Atlantic Ocean", "South Pacific Ocean", "North Pacific Ocean",
                                          "South Ocean"),
                                sea_area = c("IO", "RS", "MS", "SAO",
                                             "NAO", "SPO", "NPO",
                                             "SO"),
                                longitude = c(100, 30,60, -10,
                                              -20,-120,-140,
                                              -20),
                                latitude = c(-25,10,50,-50,
                                             20,-50, 50,
                                             -70))

                            
p <- ggplot(unique_meta, aes(x=longitude, y= latitude,
                        color=sea_area, fill=sea_area)) +   
  borders("world", colour=NA, fill="lightgrey")  +
  geom_point(alpha = 1, size = 1) +
  scale_color_manual(values = mypalette) + 
  scale_fill_manual(values = mypalette) + 
  theme(panel.background = element_rect(fill = "azure1", colour = "azure1")) +
  geom_label_repel(color = "white",
                  min.segment.length = 0,
                  aes(x=longitude, y= latitude, label=station_number),
                  box.padding = 0.1) +
  geom_text(data = single_ocean_label, aes(x=longitude, y= latitude, label=ocean),
            color = "black", check_overlap = T, size = 5) + guides(fill = guide_legend(override.aes = aes(label = "")))
  

p
  
# +
#   geom_text(aes(x=rounded_longitude, y= rounded_latitude, label=station_number),
#             color = "gray20", check_overlap = T, size = 3,
#             position = position_dodge(0.5),vjust = 0.5,hjust = -2.3)

ggsave(plot=p,paste0(input_dir,"worldmap.pdf"),width = 13, height = 8) 
# 
