# For each experiment, make a merged private snvs file and place in directory called merged_snvs , at the same level as your
# sink_source.csv. e.g. example_2/merged_snvs/Signature_SNVs_merged_fam_Transmission_21_3_1_counts.csv.bz2

### USER INPUT.  Modify the "sink", and "input_dir" or try using the example data in example_2 (as shown below)
# If you want to remove SNVs that only appear in the sink and no in any source, enter TRUE for filter_unique_SNV
# this can reduce over-estimation of the unknown 
# 
sink = "B"
input_dir = "/Users/leahbriscoe/Documents/FEASTX/Signature-SNVs/example_2"
input_dir = paste0(input_dir,"/")
snvs_dir = paste0(input_dir,"merged_snvs/")
filter_unique_SNV = FALSE

### if needed, run these commands to install FEAST
#install.packages("devtools")
# Install feast using following line
#devtools::install_github("cozygene/FEAST")

Packages <- c("R.utils","Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", 
              "gridExtra", "data.table","FEAST","here","rstudioapi")


lapply(Packages, library, character.only = TRUE)


# TO RUN FROM COMMAND LINE, UNCOMMMENT AND MODIFY THE CODE BELOW
# Rscript SourceTrackingScript.R --sink MS --input_dir ~/Documents/FEASTX/TaraFiles/
# args <- cmdArgs()
# sink <- cmdArg("sink", 42L)
# input_dir <- cmdArg("input_dir", 42L)

current_path = dirname(getSourceEditorContext()$path  )
source(paste0(current_path,"/preprocessing.R"),chdir = TRUE)


### CREATE OUTPUT DIRECTORIES
snv_count = c()
dir.create(paste0(input_dir,"snv_feast_results/"))


##### PREPARING METADATA
# Family refers to the groups
fam_reference = read.csv(paste0(input_dir,"sink_source.csv"),stringsAsFactors = FALSE)
fam_ids = sort(unique(fam_reference$family_id))

####### Input Private SNVs FILE NAMING ######

file_string = paste0("_sink_",sink)
set.seed(0)

ALLtimes = list()
for(fam_id_i in c(1:length(fam_ids))){
  #fam_id_i = 1
  print(paste("fam number",fam_id_i))
  fam_id = fam_ids[fam_id_i]
  print(paste("Fam id",fam_id))
  fam_ref = fam_reference %>% filter(family_id == fam_id)
  
  ##################  SELECT SAMPLES: based on arrangement ###########
  select_samples = fam_ref[2:length(fam_ref)]
  print(paste("number sources",length(select_samples)))
  
  ################### LOAD: SNPs ##################
  path_check = paste0(snvs_dir,"Signature_SNVs_merged_fam_",fam_id, "_counts.csv.bz2")
  if(!file.exists(file.path(path_check))){
    print(paste0("Family file ", fam_id, " absent"))
    next
  }
  
  matrix_dat <- tryCatch(fread(path_check,stringsAsFactors = FALSE,header = FALSE,
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
  
  # FORMAT THE SNV DATA 
  matrix_dat = process_snv_data(snv_matrix = matrix_dat,metadata = fam_ref,
                                select_samples = select_samples)
  
  
  #  MAKE METADATA FOR SOURCE TRACKING
  
  t1 = Sys.time()  
  meta_dat= make_feast_metadata(fam_id = fam_id,
                                select_samples=select_samples, 
                                sink_name = sink) 
  colnames(matrix_dat) = row.names(meta_dat)
  
  
  #CODE SKIP THIS FAMILY ID:  If only one column in private snv matrix with non zero total counts, skip
  if(sum(colSums(matrix_dat)!=0) == 1){
    next
  }
  
  if(dim(matrix_dat)[1] == 0){
    next
  }
  
  
  if(filter_unique_SNV){
    matrix_dat = unique_snvs_filter(snv_matrix = matrix_dat)
  }
  
  
  colsum_sources = colSums(matrix_dat[,2:ncol(matrix_dat)])
  csum = colSums(matrix_dat)
  coverage_min =min(csum[csum!=0])
  
  if(sum(colsum_sources > 1) == 1){
    matrix_dat = rbind(matrix_dat,c(0,rep(1,(ncol(matrix_dat) -1))))
  }
  
  FEAST_output <- tryCatch(FEAST(C = t(matrix_dat), metadata =meta_dat, different_sources_flag = 0, dir_path =input_dir,
                                 outfile="demo",COVERAGE =coverage_min), error=function(e){
                                   print(e)
                                   print("error from FEAST")
                                   return(e)})
  

  
  head(matrix_dat)

  if(inherits(FEAST_output, "error")){
    ALLtimes[[paste0("fam_" , fam_id)]] = NA
    #next
  }else{
    print("FEAST RESULTS")
    FEAST_output = FEAST_output[c(1:(nrow(meta_dat)))]
    orig_names_output = names(FEAST_output)
    FEAST_output = do.call(rbind,FEAST_output)
    print(Sys.time() - t1)
    all_times = data.frame(variable =  orig_names_output , value = FEAST_output[,1])
    print(all_times)
    
    ALLtimes[[paste0("fam_" , fam_id)]] = all_times
    # for writing to table
    
    # write.table(all_times,paste0(input_dir,"snv_feast_results/FEAST_input", 
    #                              file_string,"_fam_",fam_id,".csv"),sep = ",",
    #             col.names = FALSE,
    #             row.names= FALSE,quote = FALSE)

    
    
    snv_count = c(snv_count,dim(matrix_dat)[1])
  }
}

  
saveRDS(ALLtimes,paste0(input_dir,"snv_feast_results/FEAST_output", file_string, ".rds"))


names_ALLtimes = paste0("Family_num_",1:length(ALLtimes))

ALLtimes_format = lapply(ALLtimes[!is.na(ALLtimes) ], function(x){
  y = x[,"value",drop=FALSE]
  return(y)
  
})

ALLtimes_table = do.call(cbind,ALLtimes_format)
colnames(ALLtimes_table ) = names(ALLtimes)[!is.na(ALLtimes)]

write.table(ALLtimes_table,paste0(input_dir,"snv_feast_results/FEAST_box_input", 
                                  file_string, ".csv"),sep = ",",
            col.names = TRUE,
            row.names=TRUE,quote = FALSE)

write.table(snv_count,paste0(input_dir,"snv_feast_results/FEAST_box_input", 
                             file_string, "_SNVCOUNT.csv"),sep = ",",
            col.names = TRUE,
            row.names=TRUE,quote = FALSE)




