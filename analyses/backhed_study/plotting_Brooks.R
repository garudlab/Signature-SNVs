require(dplyr)
library(tidyr)
require(ggbeeswarm)
library(ggpubr)
study = "Brooks"
input_dir = paste0('~/Documents/FEASTX/',study,'Files/snv_feast_results/')
out_dir = paste0('~/Documents/FEASTX/',study,'Files/pretty_plots/')
#out_dir = paste0('~/Documents/FEASTX/',study,'Files/snv_feast_results/')

dir.create(out_dir)
# seed = 9
# sink ="sink_swabs" #Infant" # "sink_swabs (9)" # #"other_surfaces_wipes"

# seed = 7
# sink = "other_surfaces_wipes"
# # # file_name = paste0("FEAST_box_input_snps_only_start_1_seed_7_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_dropout_0_sink_", sink)
# 
# seed = 8
# sink = "touched_surfaces_swabs"
# # # file_name = paste0("FEAST_box_input_snps_only_start_1_seed_8_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_dropout_0_sink_", sink)


seed = 10
sink = "Infant"
#file_name = paste0("FEAST_box_input_snps_only_start_1_seed_10_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_dropout_0_sink_", sink)

#"Intended_Seed_7"

inclusion_criteria = "snps_only"

file_name = paste0("FEAST_box_input_",inclusion_criteria,"_start_1_seed_",seed,"_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_dropout_0_sink_",sink)
#file_name = "Intended_Seed_10"
data = read.csv(paste0(input_dir,file_name,".csv"),stringsAsFactors = FALSE,row.names = 1)



new_rownames = sapply(row.names(data),function(x){
  strung_out = unlist(strsplit(x,"_")[[1]])
  paste(strung_out[1:(length(strung_out)/2)],collapse = "_")})

# remap rowbnames
remap= list()
remap[["touched_surfaces_swabs"]] = "Touched Surface"
remap[["sink_swabs" ]] = "Sink Basin"

remap[["other_surfaces_wipes"]] = "Floor and Isolette Top"
remap[["rand_sink_swabs" ]] = "Sink Basin (diff. room)"

remap[["rand_touched_surfaces_swabs"]] = "Touched Surface (diff. room)"
remap[["rand_other_surfaces_wipes"]] ="Floor and Isolette Top (diff. room)"
remap[["Infant"]] = "Infant"
remap[["Unknown"]] = "Unknown"
new_rownames = sapply(new_rownames, function(x){remap[[x]]})
row.names(data) = new_rownames

to_plot = data %>% tibble::rownames_to_column(var = "Source")%>% pivot_longer(-Source )

to_plot$Source <- factor(to_plot$Source, 
                         levels = c("Infant",
                                    "Touched Surface",
                                    "Sink Basin",
                                    "Floor and Isolette Top",
                                    "Touched Surface (diff. room)",
                                    "Sink Basin (diff. room)",
                                    "Floor and Isolette Top (diff. room)",
                                    "Unknown"))

mypalette = c("#0073C2","#6B8E23","#EEAD0E","#FF1493", "#CCE698","#D6BD63","#FFA0D3","#CDC9C9")
names(mypalette) =  c("Infant","Touched Surface","Sink Basin","Floor and Isolette Top",
                      "Touched Surface (diff. room)","Sink Basin (diff. room)",
                      "Floor and Isolette Top (diff. room)","Unknown")


# require(wesanderson)
# require(RColorBrewer)
# colors <- wes_palette("GrandBudapest2" ,n = 3)
# colors <- colorRampPalette(colors)(8)
# 
# color_list = 
#length(unique(to_plot$name))
# c(settings[1],settings[3]),
# c(settings[1],settings[4]) ,
if( seed == 7){
  # sink is isolette
  title_color = mypalette[["Floor and Isolette Top"]]
  settings =   c("Infant", "Touched Surface" , "Sink Basin",
                 "Floor and Isolette Top (diff. room)",  "Touched Surface (diff. room)" ,
                 "Sink Basin (diff. room)","Unknown"                   )
  my_comparisons <- list(c("Floor and Isolette Top (diff. room)","Infant"),
                         c("Floor and Isolette Top (diff. room)","Touched Surface"),
                         c("Floor and Isolette Top (diff. room)","Sink Basin")) #list( c(settings[1],settings[2]))#
  
}else if (seed == 10){
  title_color = mypalette[["Infant"]]
  #s ink is infant
  settings =   c("Floor and Isolette Top", "Touched Surface" , "Sink Basin",
                 "Floor and Isolette Top (diff. room)",  "Touched Surface (diff. room)" ,
                 "Sink Basin (diff. room)","Unknown"                   )
  my_comparisons <- list(c("Floor and Isolette Top","Floor and Isolette Top (diff. room)"),
                         c("Touched Surface","Touched Surface (diff. room)"),
                         c("Sink Basin","Sink Basin (diff. room)")) #list( c(settings[1],settings[2]))#
  
}else if( seed == 8){
  #sink is touched surface
  title_color = mypalette[["Touched Surface"]]
  settings =   c("Infant", "Floor and Isolette Top" , "Sink Basin",
                 "Floor and Isolette Top (diff. room)",  "Touched Surface (diff. room)" ,
                 "Sink Basin (diff. room)","Unknown"                   )
  my_comparisons <- list(c("Touched Surface (diff. room)" ,"Infant"),
                         c("Touched Surface (diff. room)","Floor and Isolette Top"),
                         c("Touched Surface (diff. room)","Sink Basin")) #list( c(settings[1],settings[2]))#
  
}else if( seed == 9){
  #sink swab
  title_color = mypalette[["Sink Basin"]]
  settings =   c("Infant", "Floor and Isolette Top" , "Touched Surface" ,
                 "Floor and Isolette Top (diff. room)",  "Touched Surface (diff room)" ,
                 "Sink Basin (diff. room)","Unknown"                   )
  my_comparisons <- list(c("Sink Basin (diff. room)","Infant"),
                         c("Sink Basin (diff. room)","Floor and Isolette Top"),
                         c("Sink Basin (diff. room)","Touched Surface")) #list( c(settings[1],settings[2]))#
  
}
range(to_plot$value)

to_plot$value = 100* as.numeric(to_plot$value)
# g <- ggplot_build(p)
# unique(g$data[[1]]["colour"])
# color_list = list()
# color_list[["sink_swabs"]] = "#CD534C"
# color_list[["other_surfaces_wipes"]] = "#AF8D00"
# color_list[["touched_surfaces_swabs"]] = "#CD534C"
# color_list[["Infant"]] = "#0073C2"
#EFC000FF


 

p <- ggboxplot(to_plot, x = "Source", y = "value",
               color = "Source", palette = "jco",add =  "jitter", size = 1) + ylab("% Contribution from Source") +
  
  theme(text = element_text(size=14),legend.position="none",
        aspect.ratio=4/3,plot.margin = margin(0, 0, 0, 2, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(color = title_color))+ylim(0,130) + 
  ggtitle(paste0("Sink: ",remap[[sink]])) + 
  scale_color_manual(values =mypalette)+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",
                     hide.ns=FALSE)
p
ggsave(plot=p,paste0(out_dir, file_name,".png"),
       width = 6, height = 7) 
 
