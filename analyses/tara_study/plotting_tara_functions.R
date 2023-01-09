metadata_maker <- function(input_dir,seed){
  
  # travels2[c("TARA_018","TARA_023","TARA_025"),c("TARA_018","TARA_023","TARA_025")]
  meta_by_accession = read.csv(paste0(input_dir,"Metadata/metadata_merge.csv"),stringsAsFactors = FALSE)
  
  #"ERR599144", "ERR599052",ERR598963,ERR315858
  unique_meta = meta_by_accession[!duplicated(meta_by_accession$station), ]
  unique_meta$rounded_longitude = round(unique_meta$Longitude_East,1)
  unique_meta$rounded_latitude = round(unique_meta$Latitude_North,1)
  unique_meta$coord = paste0(unique_meta$rounded_longitude,",", unique_meta$rounded_latitude)

  row.names(meta_by_accession) = meta_by_accession$run_accession
  meta_by_pangeaid = meta_by_accession
  row.names(meta_by_pangeaid) = meta_by_accession$PANGAEA.sample.identifier
  
  
  
  fam_reference = read.csv(paste0(input_dir,"sink_source_config/Pipe_ChosenFamily_seed",seed,".csv"),stringsAsFactors = FALSE,row.names = 1)
  mapping = data.frame(shorthand = c(colnames(fam_reference),"Unknown"), acc = c(unlist(fam_reference[1,]),"Unknown"), 
                     station = c(meta_by_accession[unlist(fam_reference[1,]),"station"],"Unknown"),
                     id = c(meta_by_accession[unlist(fam_reference[1,]),"PANGAEA.sample.identifier"],"Unknown"))
  newrn = row.names(mapping)
  newrn[length(newrn)] = "Unknown"
  row.names(mapping) = newrn
  
  return(list(meta_by_accession,meta_by_pangeaid,unique_meta,mapping, fam_reference))
}


load_data <- function(args,onesource_TF){
  analysis_type =args[1] # # 
  seed = args[2]
  sink= args[3]
  one_source = args[4]
  distance_metric =args[5]
  inclusion_criteria =args[6]
  
  
  #"combo" ###"combo"##
  uniqueness = 1
  threshold = "1.0"
  familythreshold = 0.1
  start =1 # NA
  min_reads = 5
  
  dropout = 0
  
  file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                       familythreshold,"_minreads_",min_reads)
  if(onesource_TF){
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_one_source_",one_source,"_dropout_", dropout,"_sink_",sink,".csv")
    
  }else{
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".csv")
    
  }
  print(ALLtimes_path)
  
  all_results = read.csv(ALLtimes_path,stringsAsFactors = FALSE,row.names = 1)
  row.names(all_results) = sapply(row.names(all_results),function(x){unlist(strsplit(x,"_"))[1]})
  
  return(all_results)
}
plotting_dots <- function(args,all_results,input_dir,labels_on_or_off,plot_always){
  
  analysis_type =args[1] # # 
  seed = args[2]
  sink= args[3]
  one_source = args[4]
  distance_metric =args[5]
  inclusion_criteria =args[6]
  
  all_results = all_results[1:(nrow(all_results)-1),]
  if(distance_metric ==  "CurrentsTravelTime"){
    travels = read.csv(paste0(input_dir,"travel_times_all_stations.csv"),stringsAsFactors = FALSE,row.names=1,header=TRUE)
    
  }
  if(distance_metric ==  "GeosphereDistance"){
    travels = read.csv(paste0(input_dir,"geo_distance.csv"),stringsAsFactors = FALSE,row.names=1,header=TRUE)
    
  }
  
  

  metadata_frames = metadata_maker(input_dir = input_dir,seed = seed)
  meta_by_accession = metadata_frames[[1]]
  meta_by_pangeaid = metadata_frames[[2]]
  unique_meta = metadata_frames[[3]]
  mapping = metadata_frames[[4]]
  fam_reference = metadata_frames[[5]]
  
  ## GEtting the SINK DATA
  source_accessions = unlist(fam_reference[1,row.names(all_results)])
  source_stations = meta_by_accession[source_accessions,"station"]
  
  ## GETTING THE SOURCE DATA
  sink_accessions = sapply(colnames(all_results),function(x){
    paste0("TARA",unlist(strsplit(x,"TARA"))[2])
  })
  
  sink_stations = meta_by_pangeaid[sink_accessions,"station"]
  colnames(all_results) = sink_stations
  
  if(analysis_type == "correlation_of_sites"){
    correlation_matrix = cor(all_results,all_results)
    names_matrix = correlation_matrix
    for(col in 1:ncol(names_matrix)){
      names_matrix[,col] = paste0(row.names(names_matrix),"_vs_",colnames(names_matrix)[col])
    }
    travel_matrix = travels[row.names(correlation_matrix),colnames(correlation_matrix)]
    
    all_results = data.frame(correlation = correlation_matrix[upper.tri(correlation_matrix)],
                             sites_compared = names_matrix[upper.tri(correlation_matrix)],
                             distance = travel_matrix[upper.tri(correlation_matrix)])
    to_plot = all_results
    
    
    # STATS FOR PLOTS
    
    if(plot_always){
      y_max = max(to_plot$correlation,na.rm=TRUE)
      lm = lm( correlation ~ distance,data=to_plot)
      stats_vec = summary(lm)
      if(labels_on_or_off == "on"){
        label_to_plot = all_results$sites_compared
      }else{
        label_to_plot = NA
      }
      
      
      ggp <- ggplot(to_plot, aes(x=distance, y= correlation,label =label_to_plot)) + 
        geom_point() + theme_bw() + geom_text(size=2.7,vjust= -1,hjust=0.6) +
        xlab(distance_metric)+
        geom_smooth(method = "lm", linetype = "dashed", size = 1,se=F) + 
        annotate(geom="text", 
                 x=min(to_plot$distance) + 0.1*(max(to_plot$distance)- min(to_plot$distance)), y=1.15*y_max, 
                 label=paste("Slope = ",signif(stats_vec$coefficients["distance","Estimate"],digits=2),
                             ", Adj R2 = ", signif(stats_vec$adj.r.squared,digits=2),
                             ", p = ",signif(stats_vec$coefficients["distance",4],digits=2)),
                 color = "#00A08A",size = 3,hjust=0) + ggtitle(paste0(one_source, " -> ", sink))
      
    }
    
    
  }else if(analysis_type == "contribution_of_ocean"){
    #
    colnames(all_results) = paste0(colnames(all_results),"o",1:ncol(all_results))
    all_results$source = row.names(all_results)
    all_results$source_station = source_stations
    
    all_results_long = gather(all_results, key = "sink_station",value = "contribution",1:(ncol(all_results)-2))
    nrow(all_results_long)
    travel_matrix = as.vector(as.matrix(travels[source_stations,sink_stations ]))
    all_results_long$distance = travel_matrix
    to_plot = all_results_long
    to_plot$contribution = 100*to_plot$contribution
    
    
    # STATS FOR PLOTS
    if(plot_always){
      y_max = max(to_plot$contribution,na.rm=TRUE)
      lm = lm( contribution ~ distance,data=to_plot)
      stats_vec = summary(lm)
      
      if(labels_on_or_off == "on"){
        label_to_plot = paste0(source_station,"->", sink_station)
      }else{
        label_to_plot = NA
      }
      
      ggp <- ggplot(to_plot, aes(x=distance, y= contribution,label = label_to_plot)) + 
        geom_point() + theme_bw() + geom_text(size=2.7,vjust= -1,hjust=0.6) +
        xlab(distance_metric)+
        geom_smooth(method = "lm", linetype = "dashed", size = 1,se=F) +
        annotate(geom="text", 
                 x=min(to_plot$distance) + 0.1*(max(to_plot$distance)- min(to_plot$distance)), y=1.15*y_max, 
                 label=paste("Slope = ",signif(stats_vec$coefficients["distance","Estimate"],digits=2),
                             ", Adj R2 = ", signif(stats_vec$adj.r.squared,digits=2),
                             ", p = ",signif(stats_vec$coefficients["distance",4],digits=2)),
                 color = "#00A08A",size = 3,hjust=0) + 
        ylim(0,1.15*y_max) + ggtitle(paste0(one_source, " -> ", sink))+
        theme(text = element_text(size=15),aspect.ratio=1,plot.title = element_text(size=15)) 
      
      
      
      
      
      
    }
    
  }
  
  
  
  
  
  if(plot_always){
    ggsave(plot=ggp  ,paste0(out_dir,"PLOT_lm_",analysis_type,"_",distance_metric, "_",study,"_sink_",sink, "_one_source_",
                             one_source,"_",inclusion_criteria,  ".pdf"),
           height = 5.5, width = 5.5) 
  }
 
  return(to_plot)
}

plotting_box <- function(sink1, sink2,source1, source2, onesource_TF,all_results1, all_results2,input_dir,labels_on_or_off,plot_always){
  
  if(!onesource_TF){
    all_results1 = all_results1[grepl(source1,row.names(all_results1)) |grepl("Unknown",row.names(all_results1)) ,]
    
    all_results2 = all_results2[grepl(source2,row.names(all_results2)) |grepl("Unknown",row.names(all_results2)) ,]
    
  }
  
  
  to_plot1 = all_results1 %>% tibble::rownames_to_column(var = "Source")%>% pivot_longer(-Source )
  to_plot1$Sink = paste0(source1, "->", sink1)
  to_plot2 = all_results2 %>% tibble::rownames_to_column(var = "Source")%>% pivot_longer(-Source )
  to_plot2$Sink = paste0(source2, "->", sink2)
  
  to_plot1$Direction = "Forward"
  to_plot2$Direction = "Backward"
  to_plot = rbind(to_plot1,to_plot2)
 
  my_comparisons <- list(c(paste0(source1, "->", sink1),paste0(source2, "->", sink2)))
  if(plot_always){
    p <- ggboxplot(to_plot%>% filter(Source != "Unknown"), x = "Sink", y = "value",
                   palette = "jco",add =  "jitter") + ylab("Proportion") +
      theme(text = element_text(size=10),legend.position="right",
            aspect.ratio=4/3,plot.margin = margin(0, 0, 0, 2, "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle(paste0("Sinks: ",sink1, " ", sink2)) + 
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",hide.ns = FALSE)
    return(p)
  }else{
    return(to_plot%>% filter(Source != "Unknown"))
  }
  

}


  
  