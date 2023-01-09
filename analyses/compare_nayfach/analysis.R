marker_snvs_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/marker_snvs/"
meta_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/metadata/"
metadata= read.csv(paste0(meta_dir,"metadata_merge_seed27.csv"),stringsAsFactors = FALSE)
all_species = list.files(marker_snvs_dir)
infant_marker_alleles = list()
mother_marker_alleles = list()
mother_infant_marker_alleles = list()
mother_infant_marker_alleles[["Escherichia_coli_53492"]]
for(s in all_species){
  print(s)
  #s = all_species[1]
  
  infant_marker_alleles[[s]] = list()
  mother_marker_alleles[[s]] = list()
  mother_infant_marker_alleles[[s]] = list()
  
    
  filenames = c("infant",  "mother","mother_infant")
  for(f in filenames){
    file_path = paste0(marker_snvs_dir,s,"/", f, "_marker_allele_strain__",s,".csv")
    if(file.size(file_path) > 1){
      data_s = read.csv(file_path)
    }else{
      print(paste("skip", file_path))
      next
    }
    
    
    
    for(cl in colnames(data_s)){
      #cl = colnames(data_s)[3]
      
      if( f == "infant"){
        infant_marker_alleles[[s]][[cl]] = data_s[,cl][(data_s[,cl] != "")]
      }else if (f == "mother"){
        mother_marker_alleles[[s]][[cl]] = data_s[,cl][(data_s[,cl] != "")]
        
      }else{
        mother_infant_marker_alleles[[s]][[cl]] = data_s[,cl][(data_s[,cl] != "")]
      }
    
    }
  }
  
  
  
}


head(metadata)
unique_study_id = unique(metadata$study_id)



## 1 line per infant_mother_species triplet
#data1_total_sharing_percentage_per_species = c() ## shared/ shared + infant only + mother only
#data1_sharing_percentage_infant = c() ## pair /pair + infant only
data1 = list()

temp_data = c(NA,NA,NA,NA, NA, NA, NA, NA, NA)
names(temp_data) = c("infant_id", "mother_id", "infant_presence","mother_presence","total_sharing","percentage_of_infant", 
                     "num_infant_snvs", "num_mother_snvs","num_pair_snvs")
## 1 line per infant_mother
data2_percentage_species_transmitted = c() # using 5% cutoff

for(studyid in unique_study_id){
  print(studyid)
  #studyid = unique_study_id[1]
  infant_id = unlist(metadata %>% filter(cohort == "B", study_id == studyid) %>% select(run_accession))
  study_mothers = unlist(metadata %>% filter(grepl("M",cohort), study_id == studyid) %>% select(run_accession))
  for( mother_id in study_mothers){
    #mother_id = study_mothers[1]
    pairname = paste0(infant_id,"_",mother_id)
    

    for( s in all_species){
      #s = "Escherichia_coli_53492"
      data1[[pairname]][[s]] = temp_data
      data1[[pairname]][[s]]["infant_id"] = infant_id
      data1[[pairname]][[s]]["mother_id"] = mother_id
      
      
      
      if(length(infant_marker_alleles[[s]][[infant_id]]) == 0){
        num_infant_marker_alleles = 0
        data1[[pairname]][[s]]["infant_presence"]  = "absent"
      }else{
        
        num_infant_marker_alleles  = length(unique(infant_marker_alleles[[s]][[infant_id]]))
        data1[[pairname]][[s]]["infant_presence"]  = "present"
        
      }
      
      if(length(mother_marker_alleles[[s]][[mother_id]]) == 0){
        num_mother_marker_alleles = 0
        data1[[pairname]][[s]]["mother_presence"]  =  "absent"
      }else{
        num_mother_marker_alleles  = length(unique(mother_marker_alleles[[s]][[mother_id]]))
        data1[[pairname]][[s]]["mother_presence"]   = "present"
      }
      
      if(length(mother_infant_marker_alleles[[s]][[  pairname ]]) == 0){
        num_mother_infant_marker_alleles = 0
      }else{
        num_mother_infant_marker_alleles = length(unique(mother_infant_marker_alleles[[s]][[pairname]]))
      }
      
      
      ## sanity check - removing any overlap between pair and mother alone since I didn't use infant in background of mom
      if( length(mother_infant_marker_alleles[[s]][[  pairname ]]) > 0   & 
          length(mother_marker_alleles[[s]][[mother_id]]) > 0){
        overlap_pair_mom = intersect(mother_infant_marker_alleles[[s]][[  pairname ]],
                                     mother_marker_alleles[[s]][[mother_id]])
        print(length(overlap_pair_mom))
        
        if (length(overlap_pair_mom) > 0){
          temp_mother_marker_alleles = mother_marker_alleles[[s]][[mother_id]]
          temp_pair_alleles = mother_infant_marker_alleles[[s]][[pairname]]
          temp_mother_marker_alleles  = temp_mother_marker_alleles [!(temp_mother_marker_alleles %in% temp_pair_alleles)]
          print(paste("old length",num_mother_marker_alleles))
          num_mother_marker_alleles  = length(temp_mother_marker_alleles)
          print(paste("new length",num_mother_marker_alleles))
        }
      }
      
      data1[[pairname]][[s]]["total_sharing"] = num_mother_infant_marker_alleles/(num_mother_infant_marker_alleles + num_infant_marker_alleles + num_mother_marker_alleles)
      ## pair/ pair + infant only + mother only
      data1[[pairname]][[s]]["percentage_of_infant"] = num_mother_infant_marker_alleles/(num_mother_infant_marker_alleles + num_infant_marker_alleles )
    
      data1[[pairname]][[s]]["num_infant_snvs"] = num_mother_marker_alleles
      data1[[pairname]][[s]]["num_mother_snvs"] =  num_infant_marker_alleles
      data1[[pairname]][[s]]["num_pair_snvs"] = num_mother_infant_marker_alleles 
   
    }
   
    
  }
  
}

pairname_df_list = list()
summary_data_list = list()
for(pairname in names(data1)){
  #pairname = names(data1)[1]
  pairname_df = data.frame(do.call(rbind,data1[[pairname]]))
  pairname_df$species = row.names(pairname_df)
  as.numeric(pairname_df$percentage_of_infant)
  pairname_df_list [[pairname]] = pairname_df
  
  present_in_infant = pairname_df %>% filter(infant_presence == "present")
  percent_present = sum(as.numeric(present_in_infant$total_sharing) > 0.05)/nrow(present_in_infant )
  percent_snvs_present = sum(as.numeric(present_in_infant$num_pair_snvs))/sum(c(as.numeric(present_in_infant$num_pair_snvs), as.numeric(present_in_infant$num_infant_snvs)))
  
  temp_summary= c(pairname_df[1,1:2 ],nrow(present_in_infant),percent_present, percent_snvs_present )
  names(temp_summary) = c("infant_id" , "mother_id", "num_species_in_infant", "percentage_infant_species_transmitted_from_mother", "percentage_infant_marker_snvs_shared")
  
  
  
  summary_data_list[[pairname]] = temp_summary
  
  
}
rbind_summary_data = do.call(rbind, summary_data_list)
rbind_pairname_df = do.call(rbind, pairname_df_list)
  



write.table(rbind_summary_data,paste0(marker_snvs_dir ,"nayfach_summary_data.csv"),sep = ",",
            col.names = TRUE,
            row.names= TRUE,quote = FALSE)

