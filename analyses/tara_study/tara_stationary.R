
require(dplyr)
library(tidyr)
require(ggbeeswarm)
library(ggpubr)
study = "Tara"
type_plot_ = "dots"
plot_always_ = FALSE
code_dir = '~/Documents/FEASTX/VAST/pipeline/scripts_source_tracking/'
input_dir = paste0('~/Documents/FEASTX/',study,'Files/')
out_dir = paste0('~/Documents/FEASTX/',study,'Files/pretty_plots/')
dir.create(out_dir)


source(paste0(code_dir,"plotting_tara_functions.R"))
plotting_dots_main <- function(args, onesource_TF,labels_on_off,plot_always){
  
  all_results = load_data(args, onesource_TF)
  to_plot_data = plotting_dots(args,all_results,input_dir,labels_on_or_off,plot_always)
  return(to_plot_data)
}

plotting_box_main <- function(args1, args2,onesource_TF, labels_on_off,plot_always){
  all_results1 = load_data(args1,onesource_TF)
  all_results2 = load_data(args2,onesource_TF)
  sink1 = args1[3]
  sink2 = args2[3]
  
  
  source1 = args1[4]
  source2 = args2[4]
  
  plotty  <-plotting_box(sink1,sink2,source1, source2, onesource_TF,all_results1, all_results2,
                         input_dir,labels_on_or_off,plot_always)
  
  if(plot_always){
    ggsave(plot=plotty  ,paste0(out_dir,"PLOT_box_onesource",onesource_TF,"_",study,"_sink1_",sink1, "_sink2_", sink2,
                                "_",inclusion_criteria,  ".pdf"),
           height = 5.5, width = 5.5)
  }
  return(plotty)
   
}



########################DOTS ####################
if(type_plot_ == "dots"){
  
  
  
  labels_on_or_off = "off"
  inclusion_criteria = "snps_only"
  onesource_TF_ = TRUE
  distance_metric_ =  "GeosphereDistance" #  "CurrentsTravelTime" #
  analysis_type_ =  "contribution_of_ocean"
  
  # CUSTOMIZATION 
  # c(9,"IO","MS")
  # c(10,"MS","IO")
  # c(11,"NAO","SPO")
  # c(15,"SPO","NAO")
  
  # END CUSTOMIZTION

 # my_args = list()
  #seed_sink_source =c(15,"SPO","NAO")
  #my_args[[1]] = c("contribution_of_ocean",seed_sink_source, "CurrentsTravelTime",inclusion_criteria )
  #my_args[[2]] = c("contribution_of_ocean",seed_sink_source, "GeosphereDistance",inclusion_criteria )
  # 
  # my_args[[3]] = c( "correlation_of_sites",seed_sink_source, "CurrentsTravelTime",inclusion_criteria )
  # my_args[[4]] = c( "correlation_of_sites",seed_sink_source, "GeosphereDistance",inclusion_criteria )
  # 
  # for(a in 1:length(my_args)){
  #   print(a)
  #   plotting_dots_main(my_args[[a]],labels_on_or_off)
  # }
  
  
  
  group_plot = TRUE
  if(group_plot){
    
    my_args = list()
    seeds=c(9, 10, 11, 12, 13, 14, 15)
    sinks=c("IO", "MS", "NAO", "NPO", "RS", "SAO", "SPO")
    onesources=c("IO", "MS", "NAO", "NPO", "RS", "SAO", "SPO")
    counting = 0
    for(s in 1:length(seeds)){
      seed_ = seeds[s]
      sink_ = sinks[s]
      for(onesource_ in onesources){
        
        if(sink_ != onesource_){
          counting = counting + 1
          my_args[[counting]] = c(analysis_type_,c(seed_,sink_,onesource_), distance_metric_,inclusion_criteria )
          
        }
        
      }
    }
    
    all_plots_list = list()
    
    for(a in 1:length(my_args)){
      temp = plotting_dots_main(my_args[[a]],onesource_TF_,labels_on_or_off,plot_always_)
      temp$sink_ocean = my_args[[a]][3]
      temp$source_ocean = my_args[[a]][4]
      all_plots_list[[a]] = temp
    }
    all_plots_data = do.call(rbind,all_plots_list)
    
    to_plot = all_plots_data %>% filter(distance > -1)
    to_plot$distance  = as.numeric(to_plot$distance/(10^6))
    
    y_max = max(to_plot$contribution,na.rm=TRUE)
    
    # Standard LIN REG
    lm = lm( contribution ~ distance,data=to_plot)
    stats_vec = summary(lm)
    
    
    # LM model
    head(to_plot)
    to_plot$sink_source_pair = paste0(to_plot$sink_ocean,"_",to_plot$source_ocean)
    
    # LINEAR MIXED MODEL
    library(nlme)
    model_re <- lme(contribution ~ distance ,random=~+1|sink_ocean,data=to_plot)
    #model_re = lmer(contribution ~ distance + (1|sink_ocean),data = to_plot)
    model_summary = summary(model_re)
    model_anova = anova(model_re)
    model_re_coefficient = model_summary$coefficients$fixed[2]
    model_re_intercept = model_summary$coefficients$fixed[1]
    model_re_p_value = model_anova[2,4]
    
    
    
    if(labels_on_or_off == "on"){
      label_to_plot = paste0(source_station,"->", sink_station)
    }else{
      label_to_plot = NA
    }
    
    #geom_smooth(method = "lm", linetype = "dashed", size = 1,se=F) +
    # annotate(geom="text", 
    #          x=min(to_plot$distance) + 0.1*(max(to_plot$distance)- min(to_plot$distance)), y=105, 
    #          label=paste("Slope = ",signif(stats_vec$coefficients["distance","Estimate"],digits=2),
    #                      ", Adj R2 = ", signif(stats_vec$adj.r.squared,digits=2),
    #                      ", p = ",signif(stats_vec$coefficients["distance",4],digits=2)),
    #          color = "#00A08A",size = 3,hjust=0) + 
    
    
    if(signif(model_re_p_value,digits=2) == 0){
      
      label_line = paste("Slope = ",signif(model_re_coefficient,digits=2),
                         ", p < 1e-16")
    }else{
      label_line = paste("Slope = ",signif(model_re_coefficient,digits=2),
                         ", p = ",signif(model_re_p_value,digits=2))
    }
    ggp <- ggplot(to_plot, aes(x=distance, y= contribution,label = label_to_plot)) + 
      theme_bw() + geom_text(size=2.7,vjust= -1,hjust=0.6) +
      geom_point(alpha = 0.3) + 
      xlab("Distance (thousands of km)")+
      ylab("% Contribution of a Source Ocean") + 
      geom_abline(intercept =  model_re_intercept, slope =model_re_coefficient , 
                  colour="blue") +
      
      
      annotate(geom="text", 
               x=0, y=105, 
               label=label_line,
               color = "#00A08A",size = 5,hjust=0) + 
      ylim(0,105) + ggtitle(paste0("All seas, one source, ",inclusion_criteria, "\n",
                                   "contribution ~ distance + (1| sink_ocean)"))+
      theme(text = element_text(size=15),aspect.ratio=1,plot.title = element_text(size=15)) 
    #ylim(0,1.15*y_max)
    ggp
    
    
    
    
   
    
   
    ggsave(plot=ggp  ,paste0(out_dir,"PLOT_MANY_lm_",analysis_type_ ,"_",distance_metric_, "_",study, "_one_source_",
                             onesource_TF_,"_",inclusion_criteria,  ".pdf"),
           height = 5.5, width = 5.5) 
    
    
    
    
    head(to_plot)
    
  }
}
#################### BOX #########
if(type_plot_ == "box"){
    
  # seed_sink_source1 =c(9,"IO","MS")
  # seed_sink_source2 =c(10,"MS","IO")
  # # 
  #seed_sink_source1 =c(13,"RS","MS")
  #seed_sink_source2 =c(10,"MS","RS")
  # 
  
  labels_on_or_off = "on"
  inclusion_criteria = "snps_only"
  onesource_TF_ = TRUE
  group_plot = TRUE
  analysis_type_ = "contribution_of_ocean"
  if(group_plot){
    
    my_args = list()
    my_args2 = list()
    
    existing_pairings = c()
    seeds=c(9, 10, 11, 12, 13, 14, 15)
    sinks=c("IO", "MS", "NAO", "NPO", "RS", "SAO", "SPO")
    
    counting = 0
    for(s in 1:length(sinks)){
      seed_ = seeds[s]
      sink_ = sinks[s]
      for(os in 1:length(sinks)){
        seed2_ =  seeds[os]
        sink2_  = sinks[os]
        if(sink_ != sink2_){
          
          
          
          new_pairing = paste0(sink_,"_",sink2_)
          if(!(new_pairing %in% existing_pairings)){
            counting = counting + 1
            existing_pairings  = c(existing_pairings,paste0(sink_,"_",sink2_), paste0(sink2_,"_",sink_))
            my_args[[counting]] = c(analysis_type_,c(seed_,sink_,sink2_), NA,inclusion_criteria )
            my_args2[[counting]] = c(analysis_type_,c(seed2_,sink2_,sink_), NA,inclusion_criteria )
            
          }
          
        }
        
      }
    }
    
    all_plots_list = list()
    
    for(a in 1:length(my_args)){
      sorted_order = paste0(sort(c(my_args[[a]][3],my_args2[[a]][3])),collapse = "_")
      given_order = paste0(c(my_args[[a]][3],my_args2[[a]][3]),collapse = "_")
      if(given_order == sorted_order ){
        temp = plotting_box_main(my_args[[a]], my_args2[[a]],onesource_TF_,labels_on_or_off,plot_always_)
        # because sink comes first
        temp$ocean_set = paste0(c(my_args2[[a]][3],my_args[[a]][3]),collapse = " -> ")
      }else{
        temp = plotting_box_main(my_args2[[a]], my_args[[a]],onesource_TF_,labels_on_or_off,plot_always_)
        temp$ocean_set = paste0(c(my_args[[a]][3],my_args2[[a]][3]),collapse = " -> ")
      }
      
      
      
      all_plots_list[[a]] = temp
    }
    
    all_plots_data = do.call(rbind,all_plots_list)
    
    head(all_plots_data)
    
    
    all_plots_data$ocean_set = factor( all_plots_data$ocean_set, levels = sort(as.character(unique(all_plots_data$ocean_set))))
    all_plots_data$Direction = factor(all_plots_data$Direction , levels = c("Forward","Backward"))
    all_plots_data$value = 100*as.numeric(all_plots_data$value)
    
    tangent = TRUE
    if(tangent){
      sub_data1 = all_plots_data %>% filter(ocean_set == "RS -> MS")
      sub_data1$FormalOcean = "Red &\nMediterranean Sea"
      
      sub_data2 = all_plots_data %>% filter(ocean_set == "SPO -> NAO")
      sub_data2$FormalOcean = "South Pacific &\nNorth Atlantic Sea"
      
      
      sub_data12 = rbind(sub_data1,sub_data2)

      
      p <- ggboxplot( sub_data12, x="Sink", y="value", color= "Direction", 
                      add = "jitter", size = 0.5) +
        ggtitle(paste0(inclusion_criteria)) +
        xlab("Ocean Pair") + ylim(0,100) + 
        ylab("% Contribution of a Source Ocean") +
        theme_classic() + 
        theme(axis.title=element_text(size=14),text =element_text(size=12),
              title = element_text(size=12),
              axis.text.x = element_text(angle= 45,hjust=1, color = "black", 
                                         size = 12),
              axis.text.y = element_text(color = "black",size=12),
              legend.position = "none",
              strip.background =element_rect(fill="grey95"),
              plot.margin = margin(0,0.75,0,0.9, "cm"),
              strip.text.x = element_text(size = 12)) +
        facet_wrap(~FormalOcean, scale = "free")+
        geom_bracket(aes(label=p.signif, xmin=1, xmax=2),
                     # this takes some fiddling to look right
                     y.position=max(sub_data12$value),
                     # match this to dodge width
                     bracket.shorten=0.4,
                     # get the p-values from the stats
                     data=compare_means(formula=value~Sink, 
                                        group.by="FormalOcean", data=sub_data12) )
      
      
      ggsave(plot=p  ,paste0(out_dir,"PLOT_FEW_box_",onesource_TF_,"_",study,
                             "_",inclusion_criteria,  ".pdf"),
             height = 5, width = 7)
      
    }
    #list(c("RS->MS","MS->RS"),c("SPO->NAO","NAO->SPO")
    # p2 <- ggplot(all_plots_data, aes(x=Sink, y=value, fill=Direction)) + 
    #   geom_boxplot() +
    #   facet_wrap(~ocean_set, scale="free")
    # 
    # p2
    #"Sink"
    #position=position_dodge(0.75)
    p <- ggboxplot( all_plots_data, x="Direction", y="value", color= "Direction", 
                    add = "jitter", size = 1) +
      ggtitle(paste0(inclusion_criteria)) +
      xlab("Ocean Pair") + ylim(c(0,100)) + 
      ylab("% Contribution of a Source Ocean") +
      theme_classic() + 
      theme(axis.title=element_text(size=18),text =element_text(size=18),
            title = element_text(size=18),
            axis.text.x = element_text(angle= 45,hjust=1, color = "black", size = 13),
            axis.text.y = element_text(color = "black"),
            legend.position = "none",
            strip.background =element_rect(fill="grey95"),
            plot.margin = margin(0,0.75,0,0.9, "cm"),
            strip.text.x = element_text(size = 13)) +
      
      stat_compare_means(comparisons = list(c("Forward","Backward")),method = "wilcox.test",label = "p.signif",hide.ns = FALSE) +
      facet_wrap(~ocean_set, scale = "free",ncol = 7)
    p
    
    ggsave(plot=p  ,paste0(out_dir,"PLOT_MANY_box_",onesource_TF_,"_",study,
                                "_",inclusion_criteria,  ".pdf"),
           height = 12, width = 20)
    
    
    
    p <- ggboxplot( all_plots_data, x="Direction", y="value", color= "Direction", 
                    add = "jitter", size = 1) +
      ggtitle(paste0(inclusion_criteria)) +
      xlab("Ocean Pair") + ylim(c(0,100)) + 
      ylab("% Contribution of a Source Ocean") +
      theme_classic() + 
      theme(axis.title=element_text(size=18),text =element_text(size=18),
            title = element_text(size=18),
            axis.text.x = element_text(angle= 45,hjust=1, color = "black", size = 13),
            axis.text.y = element_text(color = "black"),
            legend.position = "none",
            strip.background =element_rect(fill="grey95"),
            plot.margin = margin(0,0.75,0,0.9, "cm"),
            strip.text.x = element_text(size = 13)) +
      
      stat_compare_means(comparisons = list(c("Forward","Backward")),method = "wilcox.test",label = "p.format",hide.ns = FALSE) +
      facet_wrap(~ocean_set, scale = "free",ncol = 7)
    p
    
    ggsave(plot=p  ,paste0(out_dir,"PLOTformat_MANY_box_",onesource_TF_,"_",study,
                           "_",inclusion_criteria,  ".pdf"),
           height = 12, width = 20)
    
  }
  else{
    
    
    # END CUSTOMIZTION
    my_args1 = list()
    my_args2= list()
    my_args1 = c("contribution_of_ocean",seed_sink_source1, NA,inclusion_criteria )
    my_args2 = c("contribution_of_ocean",seed_sink_source2, NA,inclusion_criteria )
    
    plotting_box_main(my_args1 ,my_args2,onesource_TF_,labels_on_or_off)
  }  
}
  
  
  

