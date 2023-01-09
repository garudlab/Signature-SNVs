require(dplyr)
library(tidyr)
require(ggbeeswarm)
library(ggpubr)
study = "Backhed"
input_dir = paste0('~/Documents/FEASTX/',study,'Files/')
out_dir = paste0(input_dir, "pretty_plots/")


plot_type = "flippy"

seeds = c("103","105","111","112","121","122")
sinks = c("B","M","4M","M","12M","M")
configs = c("M->B","B->M","M->4M","4M->M","M->12M","12M->M")
inclusion_criteria ="snps_only" #"combo" ###"combo"##
one_source = ""
study = "Backhed"
uniqueness = 1
threshold = "1.0"
familythreshold = 0.1
start =1 # NA
dropout= 0
min_reads = 10

#load data
all_data = list()

for(s in 1:length(seeds)){
  seed = seeds[s]
  sink = sinks[s]
  config = configs[s]
  
  file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                       familythreshold,"_minreads_",min_reads)
  
  
  if(one_source != ""){
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_one_source_",one_source,"_dropout_", dropout,"_sink_",sink,".csv")
    
  }else{
    ALLtimes_path = paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".csv")
    
  }
  all_results = read.csv(ALLtimes_path,stringsAsFactors = FALSE,row.names = 1)
  new_rownames = sapply(row.names( all_results ),function(x){unlist(strsplit(x,"_")[[1]])[1]})
  row.names( all_results ) = new_rownames
  
  all_results$sources = row.names(all_results)
  all_results_long = gather(all_results, key = "sink",value = "est_contribution",1:(ncol(all_results)-1))
  all_results_long$config = config
  all_data[[s]] = all_results_long
}
head(all_data)
all_data_df = do.call(rbind,all_data)


if(plot_type == "flippy"){
  configs = list(c("M->B","B->M"),c("M->4M","4M->M"),c("M->12M","12M->M"))
  source_vecs = list(c("M","B"),c("M","Mo4"),c("M","Mo12"))
  
  for(config_num in 1:length(configs)){
    #config_num = 1
    to_plot  = all_data_df %>% filter(config %in% configs[[config_num]], sources %in% source_vecs[[config_num]])
    to_plot$config = factor(to_plot$config, levels = configs[[config_num]])
    to_plot$est_contribution = 100* to_plot$est_contribution
    
    to_plot%>% group_by(config) %>% summarize(
      est_contrib_ = mean(est_contribution, na.rm=TRUE))
    
    p <- ggboxplot(to_plot, x = "config", y =   "est_contribution",
                   color = "config", palette = "jco",add =  "jitter") + ylab("% from source") +
      xlab("source -> sink") +
      theme(text = element_text(size=15),legend.position="none",
            axis.title=element_text(size=15,face="bold"),aspect.ratio=0.75,
            axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0,105) + 
      stat_compare_means(comparisons = list(c(1,2)),method = "wilcox.test",label = "p.signif")
    p
    ggsave(plot=p,paste0(out_dir, paste(c(source_vecs[[config_num]],inclusion_criteria),collapse= "_"),".png"),width = 6, height = 4) 
    
    p <- ggboxplot(to_plot, x = "config", y =   "est_contribution",
                   color = "config", palette = "jco",add =  "jitter") + ylab("% Explained by Source") +
      theme(text = element_text(size=20),legend.position="right",
            aspect.ratio=4/3,
            axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0,105) + 
      stat_compare_means(comparisons = list(c(1,2)),method = "wilcox.test",label = "p.format")
    p
    ggsave(plot=p,paste0(out_dir, paste(c(source_vecs[[config_num]],inclusion_criteria),collapse= "_"),"_format.png"),width = 6.5, height = 5) 
    
  }
  
  
}
  
 
if(plot_type == "time"){
  sources_of_interest =  list(c("M"),c("M0","M1","M2"), "Unknown")
  sources_of_interest_named = c("Mother","Random Mothers", "Unknown")
  for( i in 1:length(sources_of_interest)){
    to_plot_pre =  all_data_df %>% filter(config %in% c("M->B","M->4M","M->12M"), sources %in% sources_of_interest[[i]] )
    to_plot= to_plot_pre  %>% group_by(sink,config) %>% summarize(
      est_contrib_ = sum(est_contribution, na.rm=TRUE))
    source_i = sources_of_interest_named[i]
    
    to_plot$est_contribution = 100*to_plot$est_contrib_
    
    
    time = sapply(to_plot$config, function(x){
      if(x == "M->12M"){
        return(12)
      }else if(x == "M->4M" ){
        return(4)
      }else if(x == "M->B"){
        return(0.13)
      }
      
    })
    to_plot$time = time
    head(to_plot)
    
    fit <- lm(est_contribution ~ time, data = to_plot)
    summary(fit)
    
    library(nlme)
    
    model_re <- lme(est_contribution ~ time, data = to_plot ,random=~+1|sink)
    #model_re = lmer(contribution ~ distance + (1|sink_ocean),data = to_plot)
    model_summary = summary(model_re)
    model_anova = anova(model_re)
    model_re_coefficient = model_summary$coefficients$fixed[2]
    model_re_intercept = model_summary$coefficients$fixed[1]
    model_re_p_value = model_anova[2,4]
    
    
    p <- ggplot(to_plot, aes(x = time, y = est_contribution),
                color = "config", palette = "jco") + geom_point(alpha = 0.5, position=position_jitter( width=.5),colour="4F57FF") + 
      ylab(paste0( "% Infant Explained by ",source_i)) +  xlab("Time (months)") +ylim(0,100) + 
      geom_abline(slope = coef(fit)[[2]], intercept = coef(fit)[[1]])+
      annotate(geom="text", 
               x=5, y=86, 
               label=paste("Slope ==",format(signif(fit$coef[[2]], 3),scientific=TRUE)),color = "black", size=4,hjust=0, parse=TRUE,hjust = 0) + 
      annotate(geom="text", 
               x=5, y=78, 
               label=paste("p ==",format(signif(summary(fit)$coef[2,4], 3), scientific=TRUE)),color = "black", size=4,hjust=0, parse=TRUE,hjust = 0) + 
      annotate(geom="text", 
               x=5, y=70, 
               label=paste("R^{2} == ",format(signif(summary(fit)$adj.r.squared, 3), scientific=TRUE)),color = "black", size=4,hjust=0, parse=TRUE,hjust = 0) + 
      ggtitle(inclusion_criteria)+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold"),aspect.ratio=0.75)
    p
    
    
    ggsave(plot=p,paste0(out_dir,"LM_box_",inclusion_criteria,"_", source_i,".pdf")) 
    
    
    
  }
  
}

