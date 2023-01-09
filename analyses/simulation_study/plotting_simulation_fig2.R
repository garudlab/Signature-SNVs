require(tidyr)
require(dplyr)
seed = 27
sink= "B"
one_source = ""
study = "SimulationD"
code_dir = '~/Documents/FEASTX/VAST/pipeline/scripts_source_tracking/'
input_dir = paste0('~/Documents/FEASTX/',study,'Files/')
out_dir = paste0('~/Documents/FEASTX/',study,'Files/pretty_plots/')
dir.create(out_dir)


inclusion_criterias =c("snps_only","otus_only")##
uniqueness = 1
threshold = "1.0"
familythreshold = 0.1
start =1 # NA
min_reads = 10
dropout = -1 # 0 #
robustness_testing =FALSE

all_trials = list()
for(trial in 1:length(inclusion_criterias )){
  inclusion_criteria = inclusion_criterias[trial]

  
  
  file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                       familythreshold,"_minreads_",min_reads)
  
  if(one_source != ""){
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_one_source_",one_source,"_dropout_", dropout,"_sink_",sink,".csv")
    
  }else{
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".csv")
    
  }
  
  
  all_results = read.csv(ALLtimes_path,stringsAsFactors = FALSE,row.names = 1)
  row.names(all_results) = sapply(row.names(all_results),function(x){unlist(strsplit(x,"_"))[1]})
  colnames(all_results) = gsub("fam_","",colnames(all_results))
  
  meta = read.csv(paste0(input_dir,"metadata/metadata_merge_seed", seed,".csv"),stringsAsFactors = FALSE)
  meta$cohort = sapply(meta$cohort, function(x){
    if(x == "UNK1"){
      return("Unknown")
    }else{
      return(x)
    }
  })
  
  if((seed != 24) & (seed != 21)){
    backup = meta %>% select(study_id,cohort,True_prop,complexity ) 
    colnames(backup)[1:2] = c("sink","sources")
  }

 
  true_props= meta %>% select(study_id,cohort,True_prop ) %>% spread(cohort,True_prop)
  row.names(true_props) = true_props$study_id
  true_props = t(true_props[,3:ncol(true_props)])
  
  #row.names(true_props)[nrow(true_props)] = "Unknown"
  
  same_trials = intersect(colnames(true_props),colnames(all_results))
  all_results =all_results[,same_trials]
  true_props = data.frame(true_props[,same_trials])
  true_props = true_props[row.names(all_results),]
  
  all_results$sources = row.names(all_results)
  true_props$sources = row.names(true_props)
  #head(true_props)
  all_results_long = gather(all_results, key = "sink",value = "est_contribution",1:(ncol(all_results)-1))
  true_props_long = gather(true_props, key = "sink",value = "true_contribution",1:(ncol(all_results)-1))
  if((seed != 24) & (seed != 21)){
    true_props_long_with_complexity = merge(true_props_long, backup, by=c("sink","sources"))
    to_plot = merge(all_results_long , true_props_long_with_complexity, by=c("sources","sink"))
    
  }else{
    to_plot = merge(all_results_long , true_props_long, by=c("sources","sink"))
  }
  to_plot$datatype = inclusion_criteria
  
  
  # complexity = sapply(to_plot$sink,function(x){
  #   meta %>% filter(study_id == x) %>% select(Complexity)
  # })
  
  all_trials[[trial]] = to_plot
}
all_trial_frame = do.call(rbind,all_trials)
if((seed == 24) | (seed == 21)){
  instrain_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/inStrain_compare_tables/"
  
  write.table(all_trial_frame,paste0(instrain_dir ,"seed", seed,"_data.csv"),sep = ",",
              col.names = TRUE,
              row.names= TRUE,quote = FALSE)
}



for(comp in c("Simple","Complex")){
  #comp = "Simple"
  #for(comp in c("Simple")){
  
  to_plot  = all_trial_frame  %>% filter(complexity == comp)
  
  to_plot$est_contribution = 100 * to_plot$est_contribution
  to_plot$true_contribution = 100 * to_plot$true_contribution
  
  cor_stats = list()
  RMSE = list()
  stats_vec = list()
  for(u in inclusion_criterias){
    #head(all_trial_frame)
    #u = "snps_only"
    sub_data  = to_plot %>% filter(datatype ==u)
    cor_stats[[u]] = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = "pearson")
    RMSE[[u]] = sqrt(sum(( sub_data$est_contribution - sub_data$true_contribution)^2)/nrow(sub_data))
    stats_vec[[u]] = c( round(cor_stats[[u]]$estimate,3),
                         signif(cor_stats[[u]]$p.value, digits=3),round(RMSE[[u]],3))
  }
  
  #head(all_trial_frame)
  require("wesanderson")
  colors = wes_palette("Darjeeling1",n = length(inclusion_criterias))
  #colors = c("#868686FF","#0073C2FF","#EFC000FF","#A73030FF","#003C67FF","#7AA6DCFF","#8F7700FF")
  #"#00A08A" "#FF0000"
  head(to_plot)
  to_plot$datatype = sapply(to_plot$datatype,function(x){
    if(x == "otus_only"){
      return("Species-FEAST")
    }else{
      return("SNV-FEAST")
    }
  })
  #factor(to_plot$datatypes)
  to_plot$datatype = factor(to_plot$datatype, levels = c("Species-FEAST","SNV-FEAST"))
  to_plot$unknown_or_not = sapply(to_plot$sources, function(x){
    if(x == "Unknown"){
      return(x)
    }else{
      return("Known")
    }
  })
  to_plot$unknown_or_not = factor(to_plot$unknown_or_not, levels = c("Unknown","Known"))
  p <- ggplot(to_plot, aes(x = true_contribution, y =est_contribution,
                           color=datatype)) +
    annotate(geom="text", 
             x=0.0*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=100, 
             label=paste("R = ", stats_vec[["otus_only"]][1],", p= ",stats_vec[["otus_only"]][2], 
                         ", RMSE = ",stats_vec[["otus_only"]][3]),
             color = colors[1],size = 3,hjust=0)+
    
    annotate(geom="text", 
             x=0.0*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=96, 
             label=paste("R = ", stats_vec[["snps_only"]][1],", p= ",stats_vec[["snps_only"]][2], 
                         ", RMSE = ",stats_vec[["snps_only"]][3]),
              color = colors[2], size=3,hjust=0)+
    xlim(c(0,100)) +
    geom_smooth(method=lm, aes(fill=datatype),show.legend = FALSE,size=0.5) +
    scale_color_manual("Method",values = colors)+
    geom_point(alpha = 0.9, aes(shape = as.factor(unknown_or_not)))+ geom_abline() +
    scale_shape_manual("Source Type",values=c(1,16))+
    ggtitle( paste0("Community Complexity: ", comp)) +
    theme_minimal() + xlab(paste0("True Percentage")) +
    ylab(paste0("Estimated Percentage")) + 
    theme(text = element_text(size=13),plot.title = element_text(size=13)) 
  
  
  p
  
  


  ggsave(plot=p,paste0(out_dir,"PLOT_lm_", "_",study,"_seed_",
                     seed, "_complexity_",comp, ".pdf"),
         width = 6.5, height = 5) 
  
  # CORRELATION STATISTIC COMPARISON
  
  cor_stats_per_sample = list()
  head(to_plot)
  unique_sinks = unique(to_plot$sink)
  for(u in c("Species-FEAST","SNV-FEAST")){
    for(ind in unique_sinks){
      to_plot %>% filter(sink == ind)
      sub_data  = to_plot %>% filter(datatype ==u,sink == ind)
      if(nrow(sub_data)> 1){
        cor_stats = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = "pearson")
        cor_stats_per_sample[[u]] = c(cor_stats_per_sample[[u]],cor_stats$estimate)
        
      }else{
        cor_stats_per_sample[[u]] = c(cor_stats_per_sample[[u]],NA)
      }
    }
    
   
    #u = "snps_only"
    
  }
  print(comp)
  print("Wilcox Paired Test")
  result_wilcox = wilcox.test(cor_stats_per_sample[["SNV-FEAST"]],cor_stats_per_sample[["Species-FEAST"]],paired=TRUE)
  print(result_wilcox)
  print(result_wilcox$statistic)
  print(result_wilcox$p.value)
  #head(all_trial_frame)
}
  
