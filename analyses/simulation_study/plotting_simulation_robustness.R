require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)
all_trial_frame_list = list()
seeds = c(21,24)
style_plot= "groupedbar"
code_dir = '~/Documents/FEASTX/VAST/pipeline/scripts_source_tracking/'
source(paste0(code_dir,"robust_calculations.R"))

for( s in 1:length(seeds)){
  #s = 1
  
  seed = seeds[s]
  sink= "B"
  one_source = ""
  study = "SimulationD"
  
  input_dir = paste0('~/Documents/FEASTX/',study,'Files/')
  out_dir = paste0('~/Documents/FEASTX/',study,'Files/pretty_plots/')
  dir.create(out_dir)
  
  
  inclusion_criteria ="snps_only" #"combo" ###"combo"##
  uniqueness = 1
  threshold = "1.0"
  familythreshold = 0.1
  start =1 # NA
  min_reads = 10
  
  robustness_testing = TRUE
  trial_parameters = list()
  first_param = c("1","1","1","0.5","0.5","0.5","0.1","0.1","0.1")
  second_param = c("1","0.5","0.1","1","0.5","0.1","1","0.5","0.1")
  third_param = c(1,2,3)
  count = 0
  
  for(i in 1:length(first_param)){
    for(j in 1:length(third_param)){
      count = count + 1
      trial_parameters[[count]] = c(first_param[i],second_param[i],third_param[j])
      
    }
  }
  all_trials = list()
  for(trial in 1:length(trial_parameters)){
    #trial = 1
    proportion_species = trial_parameters[[trial]][1]
    proportion_snvs =  trial_parameters[[trial]][2]
    robustness_testing_seed = trial_parameters[[trial]][3]
    
    dropout = 0
    
    file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                         familythreshold,"_minreads_",min_reads)
    
    if(robustness_testing){
      file_string = paste0(file_string,"_robustnesstest", robustness_testing_seed, "_species_",proportion_species,"_snv_",proportion_snvs,"_")
    }
    
    if(one_source != ""){
      ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_one_source_",one_source,"_dropout_", dropout,"_sink_",sink,".csv")
      
    }else{
      ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".csv")
      
    }
    
    
    all_results = read.csv(ALLtimes_path,stringsAsFactors = FALSE,row.names = 1)
    row.names(all_results) = sapply(row.names(all_results),function(x){unlist(strsplit(x,"_"))[1]})
    colnames(all_results) = gsub("fam_","",colnames(all_results))
    
    meta = read.csv(paste0(input_dir,"metadata/metadata_merge_seed", seed,".csv"),stringsAsFactors = FALSE)
  
    true_props= meta %>% select(study_id,cohort,True_prop) %>% spread(cohort,True_prop)
    row.names(true_props) = true_props$study_id
    true_props = t(true_props[,3:ncol(true_props)])
    row.names(true_props)[nrow(true_props)] = "Unknown"
    
    same_trials = intersect(colnames(true_props),colnames(all_results))
    all_results =all_results[,same_trials]
    true_props = data.frame(true_props[,same_trials])
    true_props = true_props[row.names(all_results),]
    
    all_results$sources = row.names(all_results)
    true_props$sources = row.names(true_props)
    all_results_long = gather(all_results, key = "sink",value = "est_contribution",1:(ncol(all_results)-1))
    true_props_long = gather(true_props, key = "sink",value = "true_contribution",1:(ncol(all_results)-1))
    to_plot = merge(all_results_long , true_props_long, by=c("sources","sink"))
    to_plot$robust_trial = robustness_testing_seed 
    to_plot$prop_species = proportion_species
    to_plot$prop_snvs = proportion_snvs
    
    all_trials[[trial]] = to_plot
  }
  
  
  all_trial_frame = do.call(rbind,all_trials)
  head(all_trial_frame)
  all_trial_frame$params = paste0("species_",all_trial_frame$prop_species,"_snvs_",all_trial_frame$prop_snvs)
  all_trial_frame_list[[seed]] = all_trial_frame
  
  
}


full_data = do.call(rbind,all_trial_frame_list)

source_count_dict = list()
unique_sinks = unique(full_data$sink)
for(us in unique_sinks){
  #us = unique_sinks[1]
  sub_data = full_data %>% filter(sink == us)
  sub_data_prop = sub_data[!duplicated(sub_data$sources), "true_contribution"]
  source_count_dict[[us]] = sum(sub_data_prop > 0)
}
full_data$source_count = sapply(full_data$sink, function(x){
  source_count_dict[[x]]
})

unique_params = c("species_1_snvs_1" , "species_1_snvs_0.5"  , "species_1_snvs_0.1",
                  "species_0.5_snvs_1" ,"species_0.5_snvs_0.5" ,"species_0.5_snvs_0.1" ,
                  "species_0.1_snvs_1" ,"species_0.1_snvs_0.5","species_0.1_snvs_0.1")


plot_all = FALSE
full_data$params = factor(full_data$params ,levels = unique_params)
head(full_data)

if(plot_all){


  subset = c(0,0.25)
  all_trial_frame = full_data %>% filter(true_contribution > subset[1],true_contribution < subset[2])
  
  stats_vec  = calculate_metrics(all_trial_frame,unique_params)
  
  all_trial_frame$params = factor(all_trial_frame$params ,levels = unique_params)
  
  head(all_trial_frame)
  #require("wesanderson")
  require(RColorBrewer)
  colors = brewer.pal(9, "Paired")
  #colors = wes_palette("Zissou1" ,n = 9)
  
  
  head(all_trial_frame)
  unique(all_trial_frame$sources)
  all_trial_frame$unknown_or_not = sapply(all_trial_frame$sources, function(x){
    if(x == "Unknown"){
      return(x)
    }else{
      return("Known")
    }
  })
  all_trial_frame$unknown_or_not = factor(all_trial_frame$unknown_or_not, levels = c("Unknown","Known"))
  
  to_plot= all_trial_frame
  to_plot$params = factor(  to_plot$params ,levels = unique_params)
  
  
  p <- ggplot(all_trial_frame, aes(x = true_contribution, y =est_contribution,color=params)) +
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=1.0, 
             label=paste("R = ", stats_vec[[unique_params[1]]][1],", p= ",stats_vec[[unique_params[1]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[1]]][3]),
             color = colors[1],size = 3,hjust=0)+
    
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.96, 
             label=paste("R = ", stats_vec[[unique_params[2]]][1],", p= ",stats_vec[[unique_params[2]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[2]]][3]),
             color = colors[2], size=3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.92, 
             label=paste("R = ", stats_vec[[unique_params[3]]][1],", p= ",stats_vec[[unique_params[3]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[3]]][3]),
             color = colors[3] , size=3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.88, 
             label=paste("R = ", stats_vec[[unique_params[4]]][1],", p= ",stats_vec[[unique_params[4]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[4]]][3]),
             color = colors[4],size = 3,hjust=0) + 
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution) - min(to_plot$true_contribution)), y=0.84, 
             label=paste("R = ", stats_vec[[unique_params[5]]][1],", p= ",stats_vec[[unique_params[5]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[5]]][3]),
             color = colors[5],size = 3,hjust=0) +
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.8, 
             label=paste("R = ", stats_vec[[unique_params[6]]][1],", p= ",stats_vec[[unique_params[6]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[6]]][3]),
             color = colors[6],size = 3,hjust=0) + 
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.76, 
             label=paste("R = ", stats_vec[[unique_params[7]]][1],", p= ",stats_vec[[unique_params[7]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[7]]][3]),
             color = colors[7],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.76, 
             label=paste("R = ", stats_vec[[unique_params[8]]][1],", p= ",stats_vec[[unique_params[8]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[8]]][3]),
             color = colors[8],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=0.76, 
             label=paste("R = ", stats_vec[[unique_params[9]]][1],", p= ",stats_vec[[unique_params[9]]][2], 
                         ", RMSE = ",stats_vec[[unique_params[9]]][3]),
             color = colors[9],size = 3,hjust=0)+
    
    
    geom_smooth(method=lm, show.legend = FALSE,size=0.5,se=FALSE) + 
    scale_color_manual("Subsampling type",values = colors)+
    geom_point(alpha = 0.5, aes(shape = as.factor(unknown_or_not)))+ geom_abline() +
    scale_shape_manual("Source Type",values=c(1,16))+
    ggtitle("Trials") +
    theme_minimal() + xlab(paste0("True Percentage")) +
    ylab(paste0("Estimated Percentage")) + 
    theme(text = element_text(size=13),plot.title = element_text(size=13)) 
  
  
  p
  

  
  ggsave(plot=p,paste0(out_dir,"PLOT_lm_", "_",study,"_robustness_",
                       inclusion_criteria,  ".pdf"),
         width = 7, height = 5) 
}else if(style_plot == "lines"){
  
  subset = c(0,1.0)
  full_data = full_data %>% filter(true_contribution > subset[1],true_contribution < subset[2])
  
  mean_coefficient_list = list()
  mean_RMSE_list = list()
  mean_intercept_list = list()
  mean_corr_list = list()
  
  
  coefficient_list = list()
  RMSE_list = list()
  intercept_list = list()
  corr_list = list()
  for(u in unique_params){
    
    coefficient_vec = c()
    RMSE_vec = c()
    intercept_vec = c()
    corr_vec = c()
    
    for(trial_t in c(1:3)){
      
      
      sub_data = full_data %>% filter( 
                                             params == u, robust_trial == trial_t)
      
    
      model_t = lm(  est_contribution ~ true_contribution ,data = sub_data )
      intercept_vec = c(intercept_vec, model_t$coefficients[1])
      
      coefficient_vec = c(coefficient_vec, model_t$coefficients[2])
      cor_stats = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = "spearman")
     
      corr_vec = c( corr_vec, round(cor_stats$estimate,3) )
      RMSE_vec = c(RMSE_vec, sqrt(sum(( sub_data$est_contribution - sub_data$true_contribution)^2)/nrow(sub_data)))
     
    }
    mean_coefficient_list[[u]] = round(mean(coefficient_vec),3)
    mean_RMSE_list[[u]] =  round(mean(RMSE_vec),3)
    mean_intercept_list[[u]] =  round(mean(intercept_vec),3)
    mean_corr_list[[u]] =  round(mean(corr_vec),3)
    
    coefficient_list[[u]] = coefficient_vec
    RMSE_list[[u]] =  RMSE_vec
    intercept_list[[u]] =  intercept_vec
    corr_list[[u]] =  corr_vec
  }
    
  all_trial_frame = full_data
  head(all_trial_frame)
  #require("wesanderson")
  require(RColorBrewer)
  
  #colors = c(brewer.pal(6,"Blues")[4:6], brewer.pal(6,"Reds")[4:6], brewer.pal(6,"Greens")[4:6])
  
  blues   = brewer.pal(6,"Blues")[4:6]
  reds = brewer.pal(6,"Reds")[4:6]
  greens = brewer.pal(6,"Greens")[4:6]
  colors = c(blues[1], reds[1], greens[1],blues[2], reds[2], greens[2] , blues[3], reds[3], greens[3])
  
    #brewer.pal(9, "Paired")
  #colors = wes_palette("Zissou1" ,n = 9)
  
  
  head(all_trial_frame)
  unique(all_trial_frame$sources)
  all_trial_frame$unknown_or_not = sapply(all_trial_frame$sources, function(x){
    if(x == "Unknown"){
      return(x)
    }else{
      return("Known")
    }
  })
  all_trial_frame$unknown_or_not = factor(all_trial_frame$unknown_or_not, levels = c("Unknown","Known"))
  
  to_plot = all_trial_frame
  
  
  to_plot$params = factor(  to_plot$params ,levels = unique_params)
  
  p <- ggplot( to_plot , aes(x = true_contribution, y =est_contribution,color=params)) +
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=1.0, 
             label=paste("Corr = ", mean_corr_list[[unique_params[1]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[1]]]),
             color = colors[1],size = 3,hjust=0)+
    
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.96, 
             label=paste("Corr = ", mean_corr_list[[unique_params[2]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[2]]]),
             color = colors[2],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.92, 
             label=paste("Corr = ", mean_corr_list[[unique_params[3]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[3]]]),
             color = colors[3],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.88, 
             label=paste("Corr = ", mean_corr_list[[unique_params[4]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[4]]]),
             color = colors[4],size = 3,hjust=0)+
    
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.84, 
             label=paste("Corr = ", mean_corr_list[[unique_params[5]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[5]]]),
             color = colors[5],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.80, 
             label=paste("Corr = ", mean_corr_list[[unique_params[6]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[6]]]),
             color = colors[6],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.76, 
             label=paste("Corr = ", mean_corr_list[[unique_params[7]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[7]]]),
             color = colors[7],size = 3,hjust=0)+
    
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.72, 
             label=paste("Corr = ", mean_corr_list[[unique_params[8]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[8]]]),
             color = colors[8],size = 3,hjust=0)+
    annotate(geom="text", 
             x=0.1*(max(to_plot$true_contribution)- min(to_plot$true_contribution)) + min(to_plot$true_contribution) , 
             y=0.68, 
             label=paste("Corr = ", mean_corr_list[[unique_params[9]]], 
                         ", RMSE = ",mean_RMSE_list[[unique_params[9]]]),
             color = colors[9],size = 3,hjust=0)+
    

    geom_abline(intercept = intercept_list[[unique_params[1]]][1], slope =coefficient_list[[unique_params[1]]][1] , colour=colors[1]) +
    geom_abline(intercept = intercept_list[[unique_params[1]]][2], slope =coefficient_list[[unique_params[1]]][2] , colour=colors[1]) +
    geom_abline(intercept = intercept_list[[unique_params[1]]][3], slope =coefficient_list[[unique_params[1]]][3] , colour=colors[1]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[2]]][1], slope =coefficient_list[[unique_params[2]]][1] , colour=colors[2]) +
    geom_abline(intercept = intercept_list[[unique_params[2]]][2], slope =coefficient_list[[unique_params[2]]][2] , colour=colors[2]) +
    geom_abline(intercept = intercept_list[[unique_params[2]]][3], slope =coefficient_list[[unique_params[2]]][3] , colour=colors[2]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[3]]][1], slope =coefficient_list[[unique_params[3]]][1] , colour=colors[3]) +
    geom_abline(intercept = intercept_list[[unique_params[3]]][2], slope =coefficient_list[[unique_params[3]]][2] , colour=colors[3]) +
    geom_abline(intercept = intercept_list[[unique_params[3]]][3], slope =coefficient_list[[unique_params[3]]][3] , colour=colors[3]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[4]]][1], slope =coefficient_list[[unique_params[4]]][1] , colour=colors[4]) +
    geom_abline(intercept = intercept_list[[unique_params[4]]][2], slope =coefficient_list[[unique_params[4]]][2] , colour=colors[4]) +
    geom_abline(intercept = intercept_list[[unique_params[4]]][3], slope =coefficient_list[[unique_params[4]]][3] , colour=colors[4]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[5]]][1], slope =coefficient_list[[unique_params[5]]][1] , colour=colors[5]) +
    geom_abline(intercept = intercept_list[[unique_params[5]]][2], slope =coefficient_list[[unique_params[5]]][2] , colour=colors[5]) +
    geom_abline(intercept = intercept_list[[unique_params[5]]][3], slope =coefficient_list[[unique_params[5]]][3] , colour=colors[5]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[6]]][1], slope =coefficient_list[[unique_params[6]]][1] , colour=colors[6]) +
    geom_abline(intercept = intercept_list[[unique_params[6]]][2], slope =coefficient_list[[unique_params[6]]][2] , colour=colors[6]) +
    geom_abline(intercept = intercept_list[[unique_params[6]]][3], slope =coefficient_list[[unique_params[6]]][3] , colour=colors[6]) +
    
    
  
    geom_abline(intercept = intercept_list[[unique_params[7]]][1], slope =coefficient_list[[unique_params[7]]][1] , colour=colors[7]) +
    geom_abline(intercept = intercept_list[[unique_params[7]]][2], slope =coefficient_list[[unique_params[7]]][2] , colour=colors[7]) +
    geom_abline(intercept = intercept_list[[unique_params[7]]][3], slope =coefficient_list[[unique_params[7]]][3] , colour=colors[7]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[8]]][1], slope =coefficient_list[[unique_params[8]]][1] , colour=colors[8]) +
    geom_abline(intercept = intercept_list[[unique_params[8]]][2], slope =coefficient_list[[unique_params[8]]][2] , colour=colors[8]) +
    geom_abline(intercept = intercept_list[[unique_params[8]]][3], slope =coefficient_list[[unique_params[8]]][3] , colour=colors[8]) +
    
    
    geom_abline(intercept = intercept_list[[unique_params[9]]][1], slope =coefficient_list[[unique_params[9]]][1] , colour=colors[9]) +
    geom_abline(intercept = intercept_list[[unique_params[9]]][2], slope =coefficient_list[[unique_params[9]]][2] , colour=colors[9]) +
    geom_abline(intercept = intercept_list[[unique_params[9]]][3], slope =coefficient_list[[unique_params[9]]][3] , colour=colors[9]) +
    
    
    scale_color_manual("Subsampling type",values = colors)+
    geom_point(alpha = 0.5, aes(shape = as.factor(unknown_or_not)))+ geom_abline() +
    scale_shape_manual("Source Type",values=c(1,16))+
    ggtitle("Trials") +
    theme_minimal() + xlab(paste0("True Percentage")) +
    ylab(paste0("Estimated Percentage")) + 
    theme(text = element_text(size=13),plot.title = element_text(size=13)) 
  
  
  p
  
  
  
  ggsave(plot=p,paste0(out_dir,"PLOT_Schematic3_lm_", "_",study,"_robustness_",
                       inclusion_criteria,  ".jpeg"),
         width = 7, height = 5) 
  
  
}else{
  
  subset = c(0,1.0)
  full_data_orig = full_data
  
  
  for(complexity in c("low","high")){
    #complexity = "low"
    if(complexity == "low"){
      full_data = full_data_orig %>% filter(true_contribution >= subset[1],true_contribution <= subset[2],
                                            source_count <= 5)

    }else{
      full_data = full_data_orig %>% filter(true_contribution >= subset[1],true_contribution <= subset[2],
                                            source_count > 5)
    }
    mean_coefficient_list = list()
    mean_RMSE_list = list()
    mean_intercept_list = list()
    mean_corr_list = list()
    
    # One list per parameter. Each statistic in that list is for one infant and one trial
    coefficient_list = list()
    RMSE_list = list()
    intercept_list = list()
    corr_list = list()
    corr_pvalue_list = list()
    for(u in unique_params){
      
      coefficient_vec = c()
      RMSE_vec = c()
      intercept_vec = c()
      corr_vec = c()
      corr_p_value_vec = c()
      
      for(trial_t in c(1:3)){
        #head(full_data)
        unique_sinks = unique(full_data$sink)
        for(u_s in unique_sinks){
          
          sub_data = full_data %>% filter( 
            params == u, robust_trial == trial_t,sink == u_s)
          if(nrow(sub_data) == 0){
            next()
          }
          
          model_t = lm(  est_contribution ~ true_contribution ,data = sub_data )
          intercept_vec = c(intercept_vec, model_t$coefficients[1])
          
          coefficient_vec = c(coefficient_vec, model_t$coefficients[2])
          cor_stats = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = "pearson")
          
          corr_vec = c( corr_vec, round(cor_stats$estimate,3) )
          corr_p_value_vec = c( corr_p_value_vec, round(cor_stats$p.value,5) )
          RMSE_vec = c(RMSE_vec, sqrt(sum(( sub_data$est_contribution - sub_data$true_contribution)^2)/nrow(sub_data)))
        }  
      }
      mean_coefficient_list[[u]] = round(mean(coefficient_vec),3)
      mean_RMSE_list[[u]] =  round(mean(RMSE_vec),3)
      mean_intercept_list[[u]] =  round(mean(intercept_vec),3)
      mean_corr_list[[u]] =  round(mean(corr_vec),3)
      
      coefficient_list[[u]] = coefficient_vec
      RMSE_list[[u]] =  RMSE_vec
      intercept_list[[u]] =  intercept_vec
      corr_list[[u]] =  corr_vec
      corr_pvalue_list[[u]] = corr_p_value_vec
    }
    
    
    all_trial_frame = full_data
    head(all_trial_frame)
    #require("wesanderson")
    require(RColorBrewer)
    
    #colors = c(brewer.pal(6,"Blues")[4:6], brewer.pal(6,"Reds")[4:6], brewer.pal(6,"Greens")[4:6])
    
    blues   = brewer.pal(6,"Blues")[4:6]
    reds = brewer.pal(6,"Reds")[4:6]
    greens = brewer.pal(6,"Greens")[4:6]
    colors = c(blues[1], reds[1], greens[1],blues[2], reds[2], greens[2] , blues[3], reds[3], greens[3])
    
    #brewer.pal(9, "Paired")
    #colors = wes_palette("Zissou1" ,n = 9)
    
    
    
    # 
    # 
    # coefficient_list[[u]] = coefficient_vec
    # RMSE_list[[u]] =  RMSE_vec
    # intercept_list[[u]] =  intercept_vec
    # corr_list[[u]] =  corr_vec
    # corr_pvalue_list[[u]] = corr_p_value_vec
    
    dataframe_list = list()
    dataframe_long_list = list()
    for(u in names(coefficient_list)){
      #u = names(coefficient_list)[1]
      info_param = unlist(strsplit(u,"_"))
      
      dataframe_temp = data.frame(prop_species = info_param[2],
                                  prop_snvs = info_param[4],
                                  params = u,
                                  mean_coefficient = mean(coefficient_list[[u]],na.rm=TRUE),
                                  sd_coefficient = sd(coefficient_list[[u]],na.rm=TRUE),
                                  mean_RMSE = mean(RMSE_list[[u]],na.rm=TRUE),
                                  mean_intercept = mean(intercept_list[[u]],na.rm=TRUE),
                                  mean_correlation = mean(corr_list[[u]],na.rm=TRUE),
                                  sd_correlation = sd(corr_list[[u]],na.rm=TRUE),
                                  
                                  min_correlation_p = min(corr_pvalue_list[[u]],na.rm=TRUE),
                                  max_correlation_p = max(corr_pvalue_list[[u]],na.rm=TRUE))
      
      dataframe_long_temp = data.frame(prop_species = info_param[2],
                                  prop_snvs = info_param[4],
                                  params = u,
                                  correlation = corr_list[[u]])
      
      
      dataframe_list[[u]] =dataframe_temp
      dataframe_long_list[[u]] =dataframe_long_temp
      
      
    }
    to_plot =  do.call(rbind,dataframe_list)
    to_plot_long =  do.call(rbind,dataframe_long_list)
    
    
    
    to_plot$params = factor(  to_plot$params ,levels = unique_params)
    to_plot$prop_species<- sapply(  to_plot$prop_species, function(x){
      paste0(as.numeric(x)*100, " %")
    })
    to_plot$prop_snvs<- sapply(  to_plot$prop_snvs, function(x){
      paste0(as.numeric(x)*100, " %")
    })
    to_plot$prop_snvs= factor(  to_plot$prop_snvs ,levels =  c("100 %", "50 %", "10 %"))
    to_plot$prop_species= factor(  to_plot$prop_species, levels =  c("100 %", "50 %", "10 %"))
    
    
    to_plot_long$params = factor(  to_plot_long$params ,levels = unique_params)
    to_plot_long$prop_species<- sapply(  to_plot_long$prop_species, function(x){
      paste0(as.numeric(x)*100, " %")
    })
    to_plot_long$prop_snvs<- sapply(  to_plot_long$prop_snvs, function(x){
      paste0(as.numeric(x)*100, " %")
    })
    to_plot_long$prop_snvs= factor(  to_plot_long$prop_snvs ,levels =  c("100 %", "50 %", "10 %"))
    to_plot_long$prop_species= factor(  to_plot_long$prop_species, levels =  c("100 %", "50 %", "10 %"))
    
    
    #brewer.pal(9, "Paired")
    require(wesanderson)
    colors = wes_palette("Zissou1" ,n = 5)[c(1,3,5)]
    
 
    
    # p <- ggplot(to_plot, aes(fill=prop_snvs, y=mean_correlation, x=prop_species)) + 
    #   geom_bar(position="dodge", stat="identity") + 
    #   geom_errorbar(aes(ymin=mean_correlation - sd_correlation, 
    #                     ymax=mean_correlation + sd_correlation),
    #                 width=.2,                    # Width of the error bars
    #                 position=position_dodge(.9))  + 
    #   scale_color_manual("Subsampling type",values = colors)+
    #   scale_fill_manual("Percentage\nSNVs Used",values = colors)+
    #   ggtitle("Trials") +
    #   theme_minimal() + xlab(paste0("Percentage Species Used")) +
    #   ylab(paste0("Mean Correlation(Truth,Estimate)")) + 
    #   theme(text = element_text(size=13),plot.title = element_text(size=13))
    # 
    # p
    # 
    
    
    
    
    # ggsave(plot=p,paste0(out_dir,"PLOT_Schematic3_lm_", "_",study,"_robustness_",
    #                      inclusion_criteria, "_complexity_", complexity,".jpeg"),
    #        width = 7, height = 5) 
    
    
    # comparisons = list(c("species_1_snvs_1" ,"species_1_snvs_0.5"),
    #                    c("species_1_snvs_1","species_1_snvs_0.1"),
    #                    c("species_1_snvs_1","species_0.5_snvs_1"),
    #                    c("species_1_snvs_1","species_0.1_snvs_1"))
    my_comparisons = list(c("100 %","50 %"),

                          c("50 %","10 %"),
                          c("100 %","10 %"))
    
    #my_comparisons = list(c(1,2))
    
   
    #to_plot_long$prop_species
  
    
    # p <- ggplot(to_plot_long,aes(fill=prop_snvs, y=correlation, x=prop_snvs)) +
    #   geom_bar( position="dodge", stat="summary", fun= "mean") + 
    #   stat_summary(fun.data = mean_se, geom = "errorbar",
    #                width = 0.2,position=position_dodge(.9)) + 
    #   scale_color_manual("Subsampling type",values = colors)+
    #   scale_fill_manual("Percentage\nSNVs Used",values = colors)+
    #   facet_wrap(~prop_species,strip.position="bottom") +
    #   ggtitle(paste0("Community type: ", complexity)) +
    #   theme_minimal() + xlab(paste0("Percentage Species Used")) +
    #   ylab(paste0("Mean Correlation(Truth,Estimate)")) + 
    #   theme(text = element_text(size=13),plot.title = element_text(size=13),
    #         axis.ticks.x = element_blank()) + 
    #   stat_compare_means(comparisons = my_comparisons,
    #                      method = "wilcox.test",label = "p.format",hide.ns = FALSE) 
    # p
    # 
    if(complexity == "low"){
      stars_list = c("**","*","****","*","**","****","**","***","****")
      
    }else{
      stars_list = c("*","**","****","*","***","****","***","***","****")
      
    }
    p <- ggbarplot(to_plot_long,fill="prop_snvs", y="correlation", x="prop_species",
                   add = c("mean_se"),
                   position = position_dodge(width=0.8)) +
      scale_fill_manual("Percentage\nSNVs Used",values = colors) + 
      ggtitle(paste0("Community type: ", complexity)) +
      xlab(paste0("Percentage Species Used")) +
      ylab(paste0("Mean Correlation(Truth,Estimate)")) +
      geom_bracket(
        xmin = 0.75, xmax = 1.72, y.position = 1.16,
        label = stars_list[1],vjust= 0.8,
        tip.length = 0.02,color = colors[1]
      ) + 
      geom_bracket(
        xmin = 1.73, xmax = 2.75, y.position = 1.16,
        label = stars_list[2],vjust= 0.8,
        tip.length = 0.02, color = colors[1]
      ) + 
      geom_bracket(
        xmin = 0.75, xmax = 2.75, y.position = 1.2,
        label = stars_list[3],vjust= 0.8,
        tip.length = 0.02, color = colors[1]
      ) + 
      geom_bracket(
        xmin = 1.0, xmax = 1.98, y.position = 1.07,
        label = stars_list[4],vjust= 0.8,
        tip.length = 0.02,color = "#a78f10"
      ) + 
      geom_bracket(
        xmin = 2.0, xmax = 3.0, y.position = 1.07,
        label = stars_list[5],vjust= 0.8,
        tip.length = 0.02, color = "#a78f10"
      ) + 
      geom_bracket(
        xmin = 1, xmax = 3, y.position = 1.1,
        label = stars_list[6],vjust= 0.8,
        tip.length = 0.02, color = "#a78f10"
      ) + 
      geom_bracket(
        xmin = 1.25, xmax = 2.23, y.position = 0.97,
        label = stars_list[7],vjust= 0.8,
        tip.length = 0.02, color = colors[3]
      ) + 
      geom_bracket(
        xmin = 2.25, xmax = 3.25, y.position = 0.97,
        label = stars_list[8],vjust= 0.8,
        tip.length = 0.02, color = colors[3]
      ) + 
      geom_bracket(
        xmin = 1.25, xmax = 3.25, y.position = 1.01,
        label = stars_list[9],vjust= 0.8,
        tip.length = 0.02, color = colors[3]
      )
    p
    
  
 
    ggsave(plot=p,paste0(out_dir,"PLOT_lm_", "_",study,"_robustness_",
                         inclusion_criteria, "_complexity_", complexity,".jpeg"),
           width = 7, height = 5) 
    
    
      
   
    


  }

  
  
  
}
  

  
  
