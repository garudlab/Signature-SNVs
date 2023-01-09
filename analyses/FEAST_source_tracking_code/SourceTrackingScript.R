

args = commandArgs(trailingOnly=TRUE)
#args =  c("snps_only",1, 22, "1.0", "B", 10, "SimulationD" ,0, "", 0.5,0.5, 3)   #, "M" )
#args =  c("snps_only",1, 22, "1.0", "B", 10, "SimulationD" ,0, "")# , 0.5,0.5, 3)   #, "M" )

#-arg snps_only -arg 1 -arg $seed -arg $thresh -arg 0 -arg $tp -arg $minreads -arg $study -arg $arran -arg 0

#args =  c("otus_only",1, 9, "1.0", "IO", 5, "Tara" ,0, "MS" )
#args =  c("otus_only",1, 10, "1.0", "MS", 5, "Tara" ,0, "IO" )
#args =  c("snps_only",1, 9, "1.0", "IO", 5, "Tara" ,0, "MS" )
#args =  c("snps_only",1, 111, "1.0", "4M", 10, "Backhed" ,0) #, "M" )
#args =  c("snps_only",1, 111, "1.0", "M", 10, "Backhed" ,0, "Mo4" )

#args =  c("snps_only",1, 103, "1.0", "B", 10, "Backhed" ,0, "M" )

#args =  c("snps_only",1, 105, "1.0", "M", 10, "Backhed" ,0, "B" )
# args =  c("snps_only",1, 10, "1.0", "Infant", 10, "Brooks" ,0 )

#args =  c("snps_only",1, 10, "1.0", "Infant", 5, "Brooks" ,0)
local_bool = FALSE
# plot_local =TRUE
#local_bool = TRUE # 
#local_bool = FALSE
plot_local =FALSE # plot_local =TRUE
snv_count = c()

debug_mode = FALSE
filter_unique_SNV = TRUE
#
inclusion_criteria =args[1] #"combo"#"otus_only" #"combo"#snps_only"#"otus_only" # 
start=as.integer(args[2])
print(args)
seed=as.integer(args[3])
uniqueness=1
threshold=args[4]
sink = args[5]
min_reads = args[6]
study = args[7]
dropout = as.integer(args[8])
if(length(args) > 8){
  if(args[9] == "o"){
    one_source = ""
  }else{
    one_source = args[9]
  }
  
}else{
  one_source = ""
}

if(length(args) > 9){
  robustness_testing = TRUE
  proportion_species = as.numeric(args[10])
  proportion_snvs = as.numeric(args[11])
  robustness_testing_seed =  as.numeric(args[12])
}else{
  robustness_testing = FALSE
}

print("Robust ness testing")
print(robustness_testing)
familythreshold=0.1
select_genus = FALSE
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
lapply(Packages, library, character.only = TRUE)

sink_alt = sink
if(sink == "4M"){
  sink_alt = "X4M"
  
}else if (sink == "12M"){
  sink_alt = "X12M"
}



#Rscript SourceTrackingScript.R otus_only 1 10 1.0 B 10 SimulationD -1 
require(data.table)
## END packages

if(local_bool){
  input_dir =  paste0('~/Documents/FEASTX/',study,'Files/')
  code_dir = '/Users/leahbriscoe/Documents/FEASTX/VAST/pipeline/'
  overall_dir = '~/Documents/FEASTX/VAST/pipeline/'
}else if(grepl("Suez",study)){
  input_dir =  paste0('/u/scratch/b/briscoel/',study,'Files/')
  code_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/VAST/pipeline/scripts_snv_selection/'
  overall_dir = '/u/scratch/b/briscoel/'
}else{
  input_dir =  paste0('/u/home/b/briscoel/project-ngarud/FEASTX/',study,'Files/')
  code_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/VAST/pipeline/'
  overall_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/'
}
if(((seed >= 10) & (seed <= 13)) & study == "SimulationD"){
  seed_data = 7
}else if((seed %in% c(115,215,109,209)) & study == "Tara"){
  seed_data = (seed  %% 100)
}else if((seed==22 | seed == 23) & study == "SimulationD"){
  seed_data = 21
}else if((seed==25 | seed == 26) & study == "SimulationD"){
  seed_data = 24
}else{
  seed_data = seed
}

source(paste0(code_dir ,"scripts_source_tracking/FEAST_source_tracking/preprocessing.R"))
dir.create(paste0(input_dir,"snv_feast_results/"))


##### PREPARING METADATA
# Family refers to the groups
fam_reference = read.csv(paste0(input_dir,"sink_source_config/Pipe_ChosenFamily_seed",seed,".csv"),stringsAsFactors = FALSE)
if(study == "SimulationD"){
  meta_merge = read.csv(paste0(input_dir,"metadata/metadata_merge_seed",seed,".csv"),stringsAsFactors = FALSE)
  
}
fam_ids = sort(unique(fam_reference$family_id))

####### Input Private SNVs FILE NAMING ######

file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                     familythreshold,"_minreads_",min_reads)
if(one_source != ""){
  file_string = paste0(file_string,"_one_source_",one_source)
}
if(robustness_testing){
  file_string = paste0(file_string,"_robustnesstest", robustness_testing_seed, "_species_",proportion_species,"_snv_",proportion_snvs,"_")
}

inputfile_string = paste0("_seed_",seed_data,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                          familythreshold,"_minreads_",min_reads)

###### SIMULATION METADATA
if (dropout < 0){
  transmission_0_1 = read.csv(paste0(input_dir,"species_lists/species_presence_absence_seed",seed,".csv"),
                              row.names = 1,stringsAsFactors = FALSE)
}




### PREPARING COUNT DATA
if(inclusion_criteria == "otus_only"){
  if(study ==  "SimulationD"){
    counts_ab = fread(paste0(input_dir,"species/count_reads_", seed_data,".txt"),sep="\t",stringsAsFactors = FALSE)
    
  }else{
    counts_ab = fread(paste0(input_dir,"species/count_reads.txt"),sep="\t",stringsAsFactors = FALSE)
    
  }
  save_rows = counts_ab$species_id
  counts_ab = as.matrix(counts_ab[,2:ncol(counts_ab)])
  row.names(counts_ab) = save_rows
  print("original counts ab dim")
  print(dim(counts_ab))
  counts_ab_orig = counts_ab

}
  

set.seed(0)

if(robustness_testing){
  set.seed(robustness_testing_seed)
}
ALLtimes = list()

for(fam_id_i in c(start:length(fam_ids))){
  
  
  print(paste("fam number",fam_id_i))
  fam_id = fam_ids[fam_id_i]
  print(paste("Fam id",fam_id))
  fam_ref = fam_reference %>% filter(family_id == fam_id)
  
  ##################  SELECT SAMPLES: based on arrangement ###########
  select_samples = fam_ref[2:length(fam_ref)]
  print(paste("number sources",length(select_samples)))
  ##################  END SELET       ################## 
  
  
  ################### FOR DROPOUT SIMULATION: FILTER OUT SPECIES ##########
  if (dropout < 0){
    non_zero_species = names(transmission_0_1[fam_id,transmission_0_1[fam_id,] != 0])
    print("non zero species")
    print(non_zero_species)
  }else{
    non_zero_species = NA
  }
    
  ################### LOAD: SNPs ##################
  if(inclusion_criteria == "snps_only"){
    path_check = paste0(input_dir,"private_snvs_merged/CATTED_fam_",fam_id,inputfile_string,"_counts.bz2")
    
    if(study == "Tara" | (  (study ==  "SimulationD") & (seed < 24 )  )){
      path_check = paste0(input_dir,"private_snvs_merged/CATTED_fam_",fam_id,inputfile_string,"_counts.csv.bz2")
      
    }
    
    if(((seed == 12) & study == "SimulationD") | ((seed %in% c(109,115)) & study == "Tara")){
      path_check = paste0(input_dir,prefix,"1_CATTED_fam_",fam_id,inputfile_string,".csv.bz2")
    }
    if(((seed == 13) & study == "SimulationD") | ((seed %in% c(209,215)) & study == "Tara")){
      path_check = paste0(input_dir,prefix,"4_CATTED_fam_",fam_id,inputfile_string,".csv.bz2")
      
    }
    print(path_check)
    if(!file.exists(file.path(path_check))){
      
      print(paste0("Family file ", fam_id, " absent"))
      next
    }

    matrix_dat <- tryCatch(fread(path_check,stringsAsFactors = FALSE,header = TRUE,
                                 quote ="",na.strings=c("NA","NaN", " ","nan","")), error=function(e){
                                   print("not file")
                                   return(e)
                                 })
    
    if(inherits(matrix_dat, "error")){
      ALLtimes[[paste0("fam_" , fam_id)]] = NA
      next
    }
    else if (nrow(matrix_dat) == 0){
      ALLtimes[[paste0("fam_" , fam_id)]] = NA
      next
    }else{
      matrix_dat =as.matrix(matrix_dat)
    }
    
    matrix_dat = process_snv_data(snv_matrix = matrix_dat,metadata = fam_ref,
                     select_samples = select_samples,
                     dropout  = dropout,non_zero_species = non_zero_species,
                     robustness_testing = robustness_testing,
                     proportion_species = proportion_species , 
                     proportion_snvs = proportion_snvs , robustness_testing_seed = robustness_testing_seed)
  } else if(inclusion_criteria == "otus_only"){
    intersectselect_samples = intersect(colnames(counts_ab),as.character(select_samples))
    matrix_dat =  counts_ab[,as.character(select_samples)]
    if (dropout < 0){
      print("non zero species")
      print(length(non_zero_species))
      intersect_transmission = intersect(row.names(counts_ab_orig),non_zero_species)
      print(" intersect_transmission")
      print( intersect_transmission)
      counts_ab = counts_ab_orig[intersect_transmission,,drop=FALSE]
      print("Dim(counts_ab")
      print(dim(counts_ab))
    }
  }
  
  #  MAKE METADATA FOR SOURCE TRACK
  
  t1 = Sys.time()  
  meta_dat= make_feast_metadata(fam_id = fam_id,
                                 select_samples=select_samples, 
                                 sink_name = sink_alt) 
  colnames(matrix_dat) = row.names(meta_dat)
  
  
  #CODE SKIP THIS FAMILY ID:  If only one column in private snv matrix with non zero total counts, skip
  if(sum(colSums(matrix_dat)!=0) == 1){
    next
  }
  #CODE SKIP THIS FAMILY ID: If no SNVs left, skip
  if(dim(matrix_dat)[1] == 0){
    next
  }
  
  if(study == "SimulationD"){
    meta_fam = meta_merge %>% filter(study_id == fam_id)
    print("meta fam")
    print(meta_fam)
  }
  if(debug_mode){
    overall_dir = '~/Documents/FEASTX'
    source(paste0(overall_dir,"FEAST/R/FEAST.R"))
    FEAST_output <- tryCatch(FEAST(C = t(matrix_dat), metadata =meta_dat, different_sources_flag = 0, dir_path =input_dir,
                                   outfile="demo",overall_dir =paste0(overall_dir,"FEAST/R/")), error=function(e){
                                     print(e)
                                     print("feast not work")
                                     return(e)})
  }else{
    library(FEAST)
    
    if(filter_unique_SNV){
      matrix_dat = unique_snvs_filter(snv_matrix = matrix_dat)
    }
    
    
    colsum_sources = colSums(matrix_dat[,2:ncol(matrix_dat)])
    csum = colSums(matrix_dat)
    coverage_min =min(csum[csum!=0])
   
    if(sum(colsum_sources > 1) == 1){
      matrix_dat = rbind(matrix_dat,c(0,rep(1,(ncol(matrix_dat) -1))))
    }
    
   
  

    if(one_source != ""){
      customized_data = customize_data_one_source(one_source = one_source,
                                snv_matrix = matrix_dat, 
                                feast_metadata = meta_dat)
      matrix_dat = customized_data[[1]]
      meta_dat = customized_data[[2]]
    }
    
    #?FEAST
    # write.table(meta_dat,"~/Documents/FEASTX/Signature_SNVs/example_2/FEAST_metadata.csv",sep=",",
    #             quote = FALSE)
   
    FEAST_output <- tryCatch(FEAST(C = t(matrix_dat), metadata =meta_dat, different_sources_flag = 0, dir_path =input_dir,
                                   outfile="demo",COVERAGE =coverage_min), error=function(e){
                                     print(e)
                                     print("feast not work")
                                     return(e)})

  }
  
  head(matrix_dat)

  if(inherits(FEAST_output, "error")){
    ALLtimes[[paste0("fam_" , fam_id)]] = NA
    #next
  }else{
    print("FEAST TEST")
    FEAST_output = FEAST_output[c(1:(nrow(meta_dat)))]
    orig_names_output = names(FEAST_output)
    FEAST_output = do.call(rbind,FEAST_output)
    print(Sys.time() - t1)
    all_times = data.frame(variable =  orig_names_output , value = FEAST_output[,1])
    print(all_times)
    
    ALLtimes[[paste0("fam_" , fam_id)]] = all_times
    # for writing to table

    write.table(all_times,paste0(input_dir,"snv_feast_results/FEAST_box_input_", 
                                 inclusion_criteria,"_start_", start,file_string, 
                                 "_dropout_", dropout,"_sink_",sink,"_fam_",fam_id,".csv"),sep = ",",
                col.names = FALSE,
                row.names= FALSE,quote = FALSE)
    
    snv_count = c(snv_count,dim(matrix_dat)[1])
  }
  
  saveRDS(ALLtimes,paste0(input_dir,"snv_feast_results/FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".rds"))
    

  names_ALLtimes = paste0("Family_num_",start:length(ALLtimes))
  
  ALLtimes_format = lapply(ALLtimes[!is.na(ALLtimes) ], function(x){
    y = x[,"value",drop=FALSE]
    return(y)
    
  })
 
  ALLtimes_table = do.call(cbind,ALLtimes_format)
  colnames(ALLtimes_table ) = names(ALLtimes)[!is.na(ALLtimes)]
  
  write.table(ALLtimes_table,paste0(input_dir,"snv_feast_results/FEAST_box_input_", 
                               inclusion_criteria,"_start_", start,file_string, 
                               "_dropout_", dropout,"_sink_",sink,".csv"),sep = ",",
              col.names = TRUE,
              row.names=TRUE,quote = FALSE)
  
  write.table(snv_count,paste0(input_dir,"snv_feast_results/FEAST_box_input_", 
                                    inclusion_criteria,"_start_", start,file_string, 
                                    "_dropout_", dropout,"_sink_",sink,"_SNVCOUNT.csv"),sep = ",",
              col.names = TRUE,
              row.names=TRUE,quote = FALSE)
  
}


