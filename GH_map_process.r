# GH_map_process
library(gtools)
library(ggplot2)
library(rgdal)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(BMS)

# library(ggrepel)
# library(ggmap)
# library(reshape2)
# require(rgdal) 
# require(maptools)
# require(plyr)
# require(broom)
# require(rgeos)
# require(raster)
# require(leaflet)

setwd("../GIS")
HGlnd_sp = readOGR("HGland_sm.shp", layer="HGland_sm")
ggplot(data = HGlnd_sp, aes(x=long,y=lat,group=group)) + 
  geom_polygon() + coord_equal(ratio=1)
HGlnd = HGlnd_sp %>% fortify()
HGblk_sp = readOGR("HG_Blocks_sm.shp", layer="HG_Blocks_sm")
HGblk_sf = HGblk_sp %>% fortify(region = 'BlockID')
HGblk_sf  = merge(HGblk_sf, HGblk_sp@data, by.x = 'id', by.y = 'BlockID')
HGblkpts = readOGR("GH_Blocks_pts.shp", layer="GH_Blocks_pts")
HGblkpts = data.frame(HGblkpts)

# Merge with values dataframe (use same approach for combining with density data)
val = as.numeric(unique(unique(HGblk_sf$id)))
valF = as.factor(val)
values = data.frame(BlockID = valF,id=val,
                    value = val)
HGblk = merge(HGblk_sf, values, by.x='id', by.y = 'id')
rm(HGblk_sf)

c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
# cols = c24[sample(length(c24),nrow(values),replace = T)]
cols = rep(c24,3)
# Plot of blocks with labeled numbers (save this for plotting)
mapplot = ggplot() + 
  geom_polygon(data=HGblk, aes(x = long, y = lat, fill = BlockID, color=BlockID, group = group),
               alpha = 1,size = 1) +
  scale_fill_manual(values = cols, guide = FALSE) + 
  scale_color_manual(values = cols, guide = FALSE) + 
  geom_polygon(data = HGlnd, aes(x=long,y=lat,fill=piece,group=group),
               color="wheat4", fill="cornsilk1",size = 0.1) +      
  geom_text_repel(data=HGblkpts, aes(x=Xcoord, y=Ycoord, label=BlockID), 
                  force=4, point.padding = 0.1, box.padding = 0.5,
                  min.segment.length = .05,
                  nudge_x = -2, nudge_y = 1) +
  scale_x_continuous(name = "East-west (m)") + # , breaks = NULL, labels = NULL
  scale_y_continuous(name = "North-south (m)") + # , breaks = NULL, labels = NULL
  north(HGlnd,location = "topright") +
  scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
           transform = FALSE, location = "bottomleft") +  
  ggtitle("Haida Gwaii, coastal sections for sea otter model") +
  coord_equal(ratio=1) + theme_minimal()
print(mapplot)
# Example of plotting continuous values (e.g. replace value with densities)
ggplot() + 
  geom_polygon(data=HGblk, aes(x = long, y = lat, fill = value, color=value, group = group),
               alpha = 1,size = 1) +
  scale_fill_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_x_continuous(name = "East-west (m)") + # , breaks = NULL, labels = NULL
  scale_y_continuous(name = "North-south (m)") + # , breaks = NULL, labels = NULL
  geom_polygon(data = HGlnd, aes(x=long,y=lat,fill=piece,group=group),
               color="wheat4", fill="cornsilk1",size = 0.1) +   
  north(HGlnd,location = "topright") +
  scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
           transform = FALSE, location = "bottomleft") +
  coord_equal(ratio=1) + theme_minimal()

# Save mapplot and HGlnd and HGblk and HGblkpts
setwd("../GHmetapop")
setwd("./data")
save(mapplot,HGlnd,HGblk,HGblkpts,file="GISdata.rdata")
setwd("..")


