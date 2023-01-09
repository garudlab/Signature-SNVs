calculate_metrics <- function(all_trial_frame,unique_params){
  cor_stats = list()
  RMSE = list()
  stats_vec = list()
  
  
  for(u in unique_params){
    sub_data  = all_trial_frame %>% filter(params == u)
    cor_stats[[u]] = cor.test( sub_data$est_contribution, sub_data$true_contribution,method = "spearman")
    RMSE[[u]] = sqrt(sum(( sub_data$est_contribution - sub_data$true_contribution)^2)/nrow(sub_data))
    stats_vec[[u]] = c( round(cor_stats[[u]]$estimate,3),
                        signif(cor_stats[[u]]$p.value, digits=3),round(RMSE[[u]],3))
  }
  return(stats_vec)
}