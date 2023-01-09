instrain_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/marker_snvs/"
orig_instrain_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/inStrain_compare_tables/"
require(dplyr)
require(ggplot2)
use_all_marker_alleles = TRUE
if(use_all_marker_alleles){
  second_axis = "Nayfach\n% Infant's Marker Alleles Shared with Mother"
  
}else{
  second_axis = "Nayfach\n% Infant's Species With >= 5% Marker Allele Sharing"
  
}

data_list = list()
study  = "SimulationD"
out_dir = paste0('~/Documents/FEASTX/SimulationDFiles/pretty_plots/')
corr_metric = "pearson" #"spearman" 
meta_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/metadata/"
metadata= read.csv(paste0(meta_dir,"metadata_merge_seed21.csv"),stringsAsFactors = FALSE)
metadata2= read.csv(paste0(meta_dir,"metadata_merge_seed24.csv"),stringsAsFactors = FALSE)
metadata = rbind(metadata,metadata2)
temp = read.csv(paste0(instrain_dir ,"nayfach_summary_data.csv"))
#head(temp)
#temp %>% filter(is.na(percentage_infant_species_transmitted_from_mother))

temp$sink = sapply(temp$infant_id,function(x){
  
  unlist(metadata %>% filter(run_accession == x) %>% select(study_id))[1]
})

fix_ = c()
true_props_ = c()
for( r in 1:nrow(temp)){
  
  
  res  = metadata %>% filter(study_id == temp[r,"sink"],run_accession == temp[r,"mother_id"])  
           
  #print(length(res))
  fix_ = c(fix_,unlist(res %>% select(cohort)))
  
  true_props_ = c(true_props_ ,unlist(res %>% select(True_prop))) 
}


temp$sources = fix_
  

temp$true_props_ = true_props_
head(temp)


if(use_all_marker_alleles){
  instrain_data = temp[,c(5,6,7,8)]
}else{
  instrain_data = temp[,c(4,6,7,8)]
}


head(instrain_data)
colnames(instrain_data) = c("est_contribution","sink","sources","true_contribution")
nrow(instrain_data)
head(instrain_data)

instrain_data$datatype = "Nayfach"

temp1 = read.csv(paste0(orig_instrain_dir ,"seed21_data.csv"))
temp2 = read.csv(paste0(orig_instrain_dir ,"seed24_data.csv"))
head(temp1)
temp_feast = rbind(temp1,temp2)
temp_feast = temp_feast %>% filter(sources != "Unknown")

total_data = rbind(instrain_data[,colnames(temp_feast)], temp_feast)




total_data$est_contribution = 100*total_data$est_contribution
total_data$true_contribution = 100*total_data$true_contribution

#total_data %>% filter(is.na(est_contribution))
total_data %>% filter(true_contribution > 80)
total_data = total_data %>% filter(!(sink %in% c("Transmission_24_1_1","Transmission_24_3_1")))



for( percentage_type  in c("small","large","all")){
  #percentage_type = "all"
  if(percentage_type == "small"){
    to_plot = total_data %>% filter(true_contribution < 0.2)
  }else if(percentage_type == "large"){
    to_plot = total_data %>% filter(true_contribution >= 0.2)
  }else{
    to_plot = total_data
  }

  
  
  to_plot$datatype = sapply(to_plot$datatype ,function(x){
    if(x == "otus_only"){
      return("Species_FEAST")
    }else if(x == "snps_only"){
      return("SNV-FEAST")
    }else{
      return("Nayfach")
    }
  })
  cor_stats = list()
  RMSE = list()
  stats_vec = list()
  
  for(u in unique(to_plot$datatype)){
    #head(all_trial_frame)
    #u = "snps_only"
    #u = "Nayfach"
    sub_data  = to_plot %>% filter(datatype ==u, !is.na(est_contribution))
    
    
    cor_stats[[u]] = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = corr_metric)
    RMSE[[u]] = sqrt(sum(( sub_data$est_contribution - sub_data$true_contribution)^2)/nrow(sub_data))
    stats_vec[[u]] = c( round(cor_stats[[u]]$estimate,3),
                        signif(cor_stats[[u]]$p.value, digits=3),round(RMSE[[u]],3))
  }
  #"Species_FEAST",
  to_plot = to_plot %>% filter(datatype %in%  c("SNV-FEAST","Nayfach"))
  
  to_plot$datatype = factor(to_plot$datatype, levels = c("SNV-FEAST","Nayfach"))
  

  require("wesanderson")
  colors = wes_palette("Zissou1" ,n = 4)[c(1,4)]
  
  
  # annotate(geom="text", 
  #          x=0.0*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=1.0,
  #          label=paste("R = ", stats_vec[["Species_FEAST"]][1],", p= ",stats_vec[["Species_FEAST"]][2], 
  #                      ", RMSE = ",stats_vec[["Species_FEAST"]][3]),
  #          color = colors[1],size = 3,hjust=0)+
  
 
  p <- ggplot(to_plot, aes(x = true_contribution, y =est_contribution,
                           color=datatype, group = datatype)) +
    
    annotate(geom="text", 
             x=0.0*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=99, 
             label=paste("R = ", stats_vec[["SNV-FEAST"]][1],", p= ",stats_vec[["SNV-FEAST"]][2], 
                         ", RMSE = ",stats_vec[["SNV-FEAST"]][3]),
             color = colors[1], size=4,hjust=0)+
    annotate(geom="text", 
             x=0.0*(max(to_plot$true_contribution)- min(to_plot$true_contribution)), y=94, 
             label=paste("R = ", stats_vec[["Nayfach"]][1],", p= ",stats_vec[["Nayfach"]][2], 
                         ", RMSE = ",stats_vec[["Nayfach"]][3]),
             color = colors[2], size=4,hjust=0)+
    scale_y_continuous(
      
      # Features of the first axis
      name = "SNV-FEAST\nEstimated % Infant's Microbiome Explained",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~.*1,name=second_axis)
    ) +
    scale_color_manual("Method",values = colors)+
    geom_smooth(method=lm, show.legend = FALSE,linewidth=0.5) +
    geom_point(alpha = 0.9)+ geom_abline() +
    theme_minimal() + xlab(paste0("True % Infant's Microbiome Explained")) +
    ylab(paste0("Estimated Percentage")) + 
    theme(plot.title = element_text(size=13),
          axis.text=element_text(size=9),
          axis.text.y = element_text(colour = colors[1]),
          axis.text.y.right = element_text(color = colors[2]),
          axis.title.y = element_text(colour = colors[1]),
          axis.title.y.right = element_text(colour = colors[2]))
          
          
  
  p
  #p
  
  
  
  
  
  ggsave(plot=p,paste0(out_dir,"PLOT_lm_",corr_metric, "_",study,"_Nayfach_vs_FEAST",
  "_percentage_", percentage_type, "_use_all_marker_alleles_", use_all_marker_alleles,".pdf"),
         width = 6.5, height = 5) 
}

