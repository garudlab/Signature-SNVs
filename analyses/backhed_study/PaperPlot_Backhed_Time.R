forward_time = rep(FALSE,4) #c(rep(TRUE,2),FALSE) #c(FALSE,FALSE,FALSE)
ncm = c(0,0,0)
seed = c(17,18,19)#c(20,21,22)
x_axis = "SNVs" # "SNVs" #  #"SNVs only" # paste0(inclusion_criteria,ncm)
inclusion_criteria = rep("snps_only",3) #rep("snps_only",3) # ## c("otus_only", "combo", "snps_only")# rep("otus_only",3) # #
opt_filter = FALSE
filter =c("MZ","MZ","MZ")
uniqueness = 1
threshold = rep("1.0",3)#,0.8)#c(0.8,0.8) #
min_reads =c(10,10,10) #c(10,20) 
familythreshold = 0.1
#all_types_metrics[[settings[i]]]
start =1 # NA
study = "Backhed"
dropout = 0
use_bin = FALSE
sinks =   c("B","Mo4","Mo12") #c("B","4M","12M")  #  # rep("Infant",3)# #c("Mo4","Mo4","Mo4") #rep("Mo12",3) #c("B","B","B")#,"B")#rep("Mo12",3) #
input_dir = paste0('~/Documents/FEASTX/',study,'Files/')
settings = sinks #c("DZ","MZ")#,"MZ2") #c("otus_only", "combo", "snps_only")#c("MZ","MZ","MZ") #c("Accurate species", "False negative species","False negative species")#paste0(inclusion_criteria,ncm)
time_filter = FALSE
if(time_filter){
  time_stamp = c("week1","week2","week3")
  settings = time_stamp
}
arrangements = c("","","") #rep("NonRandom",3)
all_types = list()
all_types_metrics = list()

#settings = paste0(seed,inclusion_criteria,"_",threshold,"_",min_reads)#paste0(inclusion_criteria,"_",threshold) #c("otus_only", "snps_only", "combo")
#paste0(inclusion_criteria,"_",sinks) #inclusion_criteria #sinks 
require(dplyr)

for(i in 1:length(settings)){
#i=1
  arrangement = arrangements[i]
  #i = 2
  
  
  
  if(min_reads[i] > 0){
    
    file_string = paste0("_seed_",seed[i],"_uniqueness_",uniqueness,"_threshold_",threshold[i],"_familythresh_",
                         familythreshold,"_minreads_",min_reads[i])
  }else{
    file_string = paste0("_seed_",seed[i],"_uniqueness_",uniqueness,"_threshold_",threshold[i],"_familythresh_",
                         familythreshold)
  }
  # if(inclusion_criteria[i] == "otus_only"){
  #   ALLtimes = readRDS(paste0(input_dir,"FEAST_box_input_", inclusion_criteria[i],"_start_", start,file_string,"_bin_",use_bin,".rds"))
  #   
  # }else{
  #   
  if(forward_time[i]){
    if(ncm[i] > 0){
      
      ALLtimes = readRDS(paste0(input_dir,"snv_feast_results/FEAST_", arrangement,"timebox_input_ncm",ncm[i], "_",inclusion_criteria[i],"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[i],".rds"))
      
    }else{
      ALLtimes = readRDS(paste0(input_dir,"snv_feast_results/FEAST_", arrangement,"timebox_input_", inclusion_criteria[i],"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[i],".rds"))
      
    }
    
  }else{
    if(ncm[i] > 0){
      ALLtimes = readRDS(paste0(input_dir,"FEASTX_output/FEAST_", arrangement,"box_input_ncm",ncm[i], "_",inclusion_criteria[i],"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[i],".rds"))
      
    }else{
      #ALLtimes = readRDS(paste0(input_dir,"FEAST_", arrangement,"box_input_", inclusion_criteria[i],"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[i],".rds"))
      ALLtimes_path = paste0(input_dir,"FEASTX_output/FEAST_box_input_", inclusion_criteria[i],"_start_", start,file_string, "_dropout_", dropout,"_sink_",sinks[i],".rds")
      
      #ALLtimes = readRDS(paste0(input_dir,"FEAST_",  arrangement_name,"box_input_", inclusion_criteria,"_start_", start,file_string,"_bin_", use_bin,"sink_",sink,".rds"))
      ALLtimes = readRDS(ALLtimes_path)
      names(ALLtimes) = paste0("Infant_",1:length(ALLtimes))
      
      path1 = paste0(input_dir,"FEASTX_output/FEAST_box_input_snps_only_start_", start,file_string, "_dropout_", dropout,"_sink_",sinks[i],".rds")
      path2 = paste0(input_dir,"FEASTX_output/FEAST_box_input_otus_only_start_", start,file_string, "_dropout_", dropout,"_sink_",sinks[i],".rds")
      data1 = readRDS(path1 )
      names(data1) = paste0("Infant_",1:length(data1))
      data1 = data1[!is.na(data1)]
      infants1 = names(data1)
      
      data2 = readRDS(path2 )
      names(data2) = paste0("Infant_",1:length(data2))
      data2 = data2[!is.na(data2)]
      infants2 = names(data2)
      intersect_infants = intersect(infants1,infants2)
      require(dplyr)
      ALLtimes = ALLtimes[intersect_infants]#ALLtimes %>% filter(host %in% intersect_infants)
      length(ALLtimes)
      
      
    }
    
  }
  
  ALLtimes = ALLtimes[!is.na(ALLtimes)]
  if(study == 'Backhed'){
    fam_unknown = lapply(ALLtimes, function(x){
      #x = ALLtimes[[1]]
      if(length(x) > 0){
        #x =  ALLtimes[[1]]
        fam = data.frame(x %>% filter(variable %in% c("B_B","Mo4_Mo4")))
        if(nrow(fam) > 0){

          
          family = sum(fam %>% select(value))
        }else{
          family = 0
        }
        
        
        mommy = x %>% filter(variable %in% c("M_M"))
        if(nrow(mommy) > 0){
          mother = sum(mommy %>% select(value))
        }else{
          #return(list(NA,NA))
          mother = NA
        }
        
        unknown = sum(x %>% filter(variable %in% c("Unknown")) %>% select(value))
        not_mom = sum(x %>% filter(!(variable %in% c("M_M"))) %>% select(value))
        unrelated_mom = sum(x %>% filter(variable %in% c("M0_M0","M1_M1","M2_M2")) %>% select(value))
        return(list(family,unknown,mother, not_mom,unrelated_mom))
      }
      
    })
    metrics = do.call(rbind,fam_unknown)
    colnames(metrics) = c("Family", "unknown","mother","NotMom","unrelated mother")
    ALLtimes= do.call(rbind,ALLtimes)
    ALLtimes$variable = row.names(ALLtimes)
    
    rename_variable = lapply(ALLtimes$variable,function(x){
      if(grepl("_Mo4",x)){
        return("4M")
      }else if(grepl("_Mo12",x)){
        return("12M")
        
      }else if(grepl("M_M",x)){
        return("M")
        
      }else if(grepl("_M0",x)){
        return("RandomM_1")
      }
      else if(grepl("_M1",x)){
        return("RandomM_2")
      }
      else if(grepl("_M2",x)){
        return("RandomM_3")
      }else if(grepl("_B",x)){
        return("B")
      }else{
        return("Unknown")
      }
    })
    
  }
 
  
  
  
  
  
  ALLtimes$variable = unlist(rename_variable)
  
  all_types[[settings[i]]] = ALLtimes
  all_types_metrics[[settings[i]]] = metrics
}


all_types_metrics_orig = all_types_metrics
all_types_metrics
require(reshape2)

#filter_out_na
for(i in 1:length(settings)){

  if(nrow(all_types_metrics[[settings[i]]]) > 1){
    num_temp = apply(all_types_metrics[[settings[i]]],2,as.numeric)
    valid_rows = !is.na(rowSums(num_temp))
    temp =data.frame(num_temp[valid_rows,])
  }else{
    #temp = all_types_metrics[[settings[i]]]
    test = as.numeric(all_types_metrics[[settings[i]]])
    temp =data.frame(t(test))
    names(temp) = colnames( all_types_metrics[[settings[i]]])
    
    #test = data.frame(all_types_metrics[[settings[i]]])
    
    #test$dtype = "p"
  }
  
  
  temp$dtype = settings[i]
  all_types_metrics[[settings[i]]] = temp
}

fam_unknown_metrics = do.call(rbind,all_types_metrics)
library(ggpubr)
to_plot  = fam_unknown_metrics
my_comparisons <- list(c(settings[1],settings[2]),c(settings[1],settings[3]),c(settings[2],settings[3]) ) #list( c(settings[1],settings[2]))#


if(inclusion_criteria[1]== "otus_only"){
  clean_inclusion_criteria = "Species"
  
  
}else if(inclusion_criteria[1]== "snps_only"){
  clean_inclusion_criteria = "SNVs"
  
  
}else{
  clean_inclusion_criteria = inclusion_criteria[1]
}
if(study == "Backhed"){
  
  for( k in c("Family", "unknown","mother","NotMom","unrelated.mother")){
    #k="Unrelated.Mom"
    head(to_plot)
    
    ylab_ = paste0("% from ", k)
    if(k == "Family"){
      ylab_ = paste0("% from previous")
    }
    ylab_ = gsub("\\."," ",ylab_)
    p <- ggboxplot(to_plot, x = "dtype", y = k,
                   color = "dtype", palette = "jco",add =  "jitter") +xlab("Sink") + ylab(ylab_) +
      theme(text = element_text(size=20),legend.position="none",aspect.ratio=4/3)+ylim(0,1.25) + 
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")
   
    print(p)
  
    
    # x = unlist(to_plot %>% filter(dtype == "B") %>% select(Mother))
    # y = unlist(to_plot %>% filter(dtype == "Mo12") %>% select(Mother))
    # wilcox.test(x,y,paired = FALSE)
    ggsave(plot=p,paste0(input_dir,"FEASTX_plots/FEAST_box_",inclusion_criteria[1],"_",k,"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[1],".pdf")) 
    
    
    unique(to_plot$dtype)
    ## add time element
    time = sapply(to_plot$dtype, function(x){
      if(x == "B"){
        return(0)
      }else if(x == "Mo4"  | x == "4M"){
        return(4)
      }else if(x == "Mo12" | x == "12M"){
        return(12)
      }
      
    })
   to_plot$time = time
   #plot(to_plot$time,to_plot$Mom)
   #install.packages("ggbeeswarm")
   require(ggbeeswarm)
   require(tidyverse)
   # p <- ggplot(to_plot, aes(x = time, y = Mom),
   #                color = "dtype", palette = "jco") + geom_jitter() + xlab("Sink") + 
   #   theme(text = element_text(size=20),legend.position="none")+ylim(0,1.25) +
   #   stat_summary(fun.data= mean_cl_normal) + 
   #   geom_smooth(method='lm')
   # print(p)
   # 
   to_plot[,k] = 100*to_plot[,k]
   fit <- lm(to_plot[,k] ~ to_plot$time)
   summary(fit)

   p <- ggplot(to_plot, aes(x = time, y = to_plot[,k]),
               color = "dtype", palette = "jco") + geom_point(alpha = 0.5, position=position_jitter( width=.5),colour="4F57FF") + 
     ylab(ylab_) +  xlab("Time (months)") +ylim(0,100) + 
     geom_abline(slope = coef(fit)[[2]], intercept = coef(fit)[[1]])+
     labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                        "Intercept =",signif(fit$coef[[1]],5 ),
                        "\nSlope =",signif(fit$coef[[2]], 5),
                        "P =",signif(summary(fit)$coef[2,4], 5))) + 
     theme(axis.text=element_text(size=15),
           axis.title=element_text(size=15,face="bold"),aspect.ratio=0.75)
   print(p)
   ggsave(plot=p,paste0(input_dir,"FEASTX_plots/LM_box_",inclusion_criteria[1],"_",k,"_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[1],".pdf")) 
   
   
   p<-p+geom_point(alpha = 0.01, position=position_jitter(height=.5, width=.5))
  }
  
  
  # p <- ggboxplot(to_plot, x = "dtype", y = "Mom",
  #                color = "dtype", palette = "jco",add = "jitter")+ xlab("Sink") + 
  #   ggtitle(clean_inclusion_criteria) + 
  #   stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif") +
  #   theme(text = element_text(size=20),plot.title = element_text(size=20),
  #         aspect.ratio=1,
  #         legend.position = "none")
  # p
  # ggsave(plot=p,paste0(input_dir,"FEAST_box_start_", start,file_string,"_bin_", use_bin,"sink_",sinks[1],".pdf"),
  #        width=7, height=6) 
  
  
  
  # head(to_plot)
  # to_plot$Mom_Over_Not = log(to_plot$Mom/to_plot$NotMom)
  # p <- ggboxplot(to_plot, x = "dtype", y = "Mom_Over_Not",
  #                color = "dtype", palette = "jco",
  #                add = "jitter") +xlab(x_axis) + 
  #   ggtitle("Proportion that is Mom:Not Mom")+ stat_compare_means(comparisons = my_comparisons,method = "wilcox.test") 
  # p
}

