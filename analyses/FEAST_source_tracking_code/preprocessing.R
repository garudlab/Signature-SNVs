# Processing scripts for SourceTracking
process_snv_data <- function(snv_matrix = NA,metadata = NA,select_samples = NA,
                             dropout = 0,non_zero_species = NA, robustness_testing = FALSE,
                             proportion_species = NA, proportion_snvs = NA, robustness_testing_seed = NA){
  
  
  #snv_matrix = matrix_dat
  #metadata = fam_ref
  #select_samples = select_samples
  
  
  
  print("original dimensions")
  print(dim(snv_matrix))
  
  keep_name =  snv_matrix[,1]
  #mat_index = which(fam_ref %in% select_samples)
  snv_matrix = snv_matrix[,2:ncol(snv_matrix)]
  snv_matrix = apply(snv_matrix, 2,as.integer)
  
  snv_matrix[is.na(snv_matrix)] = 0
  #dim(snv_matrix)
  colnames(snv_matrix) = select_samples
  row.names(snv_matrix) = unlist(keep_name)
  
  if(dropout < 0){
    rnames = row.names(snv_matrix)
    rnames_spt = strsplit(rnames,"\\|")
    just_species = do.call(rbind,rnames_spt)[,1]
    just_species = sapply(just_species,function(x){
      temp = gsub("Ref_","",gsub("Alt_","",x))
      return(temp)
    })
    chosen_snps = sapply(just_species,function(x){
      if(x %in% non_zero_species){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })
    print(sum(chosen_snps))
    print("dim before dropout")
    print(dim(snv_matrix))
    snv_matrix = snv_matrix[chosen_snps,]
    print("dim after dropout")
    print(dim(snv_matrix))
  } else if (robustness_testing){
    
    rnames = row.names(snv_matrix)
    
    rnames_spt = strsplit(rnames,"\\|")
    just_species = do.call(rbind,rnames_spt)[,1]
    just_species = sapply(just_species,function(x){
      temp = gsub("Ref_","",gsub("Alt_","",x))
      return(temp)
    })
    
    unique_species = sort(unique(just_species))
    
    just_snvs =row.names(snv_matrix)
    just_snvs = gsub("Alt_","",gsub("Ref_","",just_snvs ))
    #num_species_sample = round(proportion_species * length(unique_species))
    #subsample_species = sample(unique_species , num_species_sample)
    species_selexn = as.logical(rbinom(length(unique_species),size = 1,prob = proportion_species))
    subsample_species = unique_species[species_selexn]
    # Subsample species
    
    time_1 = Sys.time()
    
    chosen_snvs = rep(FALSE,length(just_species))
    for(ss in subsample_species){
      #print(ss)
      #ss = subsample_species[1]
      if(ss %in% just_species){
        indexes = which(just_species == ss)
        
        
        
        
        if(proportion_snvs < 1.0){
          
          ss_snv =  just_snvs[indexes]
          unique_ss_snv  = unique(ss_snv)
          
          snvs_selexn = as.logical(rbinom(length(unique_ss_snv),size = 1,prob = proportion_snvs))
          subsample_snvs = unique_ss_snv[snvs_selexn]
          #browser()
          chosen_snvs[which((just_species == ss) & (just_snvs %in% subsample_snvs))] = TRUE
          #print(sum(chosen_snvs))
          #print(table(just_snvs[which((just_snvs == ss) & (just_snvs %in% subsample_snvs))]))
          
          #chosen_snvs[indexes] = as.logical(rbinom(length(indexes),size = 1,prob = proportion_snvs))
        }else{
          chosen_snvs[as.vector(indexes)] = TRUE
          
          #chosen_snvs[indexes[1:3]] = TRUE
          #chosen_snvs[15480]
          
        }
      }
      
        
     
    }
    
  
    print("sum chosen _snvs")
    print(sum(chosen_snvs))
    print("dim before subsample species")
    print(dim(snv_matrix))
    snv_matrix = snv_matrix[chosen_snvs,]
    print(dim(snv_matrix))
    
   

  }
  return(snv_matrix)
}
unique_snvs_filter <- function(snv_matrix = NA){
  print("CHECKING FOR UNIQUE SNV")
  print(dim(snv_matrix))
 
  vector = c(1,rep(0,ncol(snv_matrix)-1))
  bool_matrix = t(t(snv_matrix) == vector)
  unique_check = bool_matrix[,1]
  
  # get row names of "unique" to remove alt and ref
  unique_snv = row.names(snv_matrix)[unlist(unique_check)]
  unique_snv = gsub("Alt_","",gsub("Ref_","",unique_snv))
  all_snv = gsub("Alt_","",gsub("Ref_","",row.names(snv_matrix)))
  bool_vector = !(all_snv %in% unique_snv)
  snv_matrix = snv_matrix[bool_vector,]
  print(dim(snv_matrix))
  return(snv_matrix)
}

customize_data_one_source <- function(one_source,snv_matrix, feast_metadata){
  
  keep_index = which((feast_metadata$Env == one_source))
  if(length(keep_index) != 1){
    keep_index = which(grepl(one_source,feast_metadata$Env))
  }
  feast_metadata = feast_metadata[c(1,keep_index),]
  
  snv_matrix = snv_matrix[,c(1,keep_index)]
  
  if( length(keep_index) == 1){
    feast_metadata = rbind(feast_metadata,c("Dummy","Source",1))
    row.names(feast_metadata)[3] = "Dummy"
    
    snv_matrix = cbind(snv_matrix,0)
    head(snv_matrix)
    colnames(snv_matrix)[3] = "Dummy"
    
    snv_matrix = rbind(snv_matrix,c(0,rep(1,(ncol(snv_matrix) -1))))
    
  }
  return(list(snv_matrix, feast_metadata))
}

make_feast_metadata <- function(fam_id,select_samples, sink_name){
  col_id = rep(fam_id,length(select_samples))
  col_env = names(select_samples)
  col_sourcesink = rep("Source",length(select_samples))
  col_sourcesink[which(col_env== sink_name)] = "Sink"
  
  custom_metadata = data.frame('Env' = col_env, 'SourceSink' = col_sourcesink, 'id' = col_id)
  row.names(custom_metadata) = col_env
  
  return(custom_metadata)
}