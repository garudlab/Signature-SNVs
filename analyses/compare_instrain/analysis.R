instrain_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/inStrain_compare_tables/"
meta_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/metadata/"
metadata= read.csv(paste0(meta_dir,"metadata_merge_seed21.csv"),stringsAsFactors = FALSE)
metadata2= read.csv(paste0(meta_dir,"metadata_merge_seed24.csv"),stringsAsFactors = FALSE)
metadata = rbind(metadata,metadata2)
scaffold_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/scaffold_mappings/"
scaffold_info = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/scaffold_info/"
# scaffold_mapping = read.csv(paste0(scaffold_dir,"all_scaffold_ids.csv"),stringsAsFactors = FALSE,header = FALSE)
# colnames(scaffold_mapping) = c("species","scaffold")
# nrow(scaffold_mapping)
# scaffold_mapping= scaffold_mapping %>% distinct(species,scaffold, .keep_all= TRUE)
# nrow(scaffold_mapping)
# head(scaffold_mapping)

scaffold_mapping = readRDS(paste0(scaffold_dir,"scaffold_map.rds"))

all_files = list.files(instrain_dir)

require(dplyr)
threshold = 0.99999

prop_above_threshold = c()
true_props = c()
acc1_list = c()
acc2_list = c()
prop_species_covered =c()
prop_of_infant_scaffolds_list = c()

for(f in 1:length(all_files)){
  #f = 1
  file = all_files[f]
  #file = "CompareCV5_ERR0S21T14N1_ERR525987_comparisonsTable.tsv"
  print(file)
  acc1 = unlist(strsplit(file,"_"))[2]
  acc2 = unlist(strsplit(file,"_"))[3]
  
  data = read.csv(paste0(instrain_dir,file),stringsAsFactors = FALSE,sep ="\t")
  scaffold_data1 = read.csv(paste0(scaffold_info ,"instrain_",acc1, "_scaffold_info.tsv"),stringsAsFactors = FALSE,sep ="\t")
  #scaffold_data2 = read.csv(paste0(scaffold_info ,"instrain_",acc2, "_scaffold_info.tsv"),stringsAsFactors = FALSE,sep ="\t")
  #scaffold_data1 = scaffold_data1 %>% filter(coverage >= 5) ## for the infant
   #nrow
  #head(scaffold_data1)
  # scaffold_mapping1 = read.csv(paste0(scaffold_dir,"scaffold_mapping_",acc1,".txt"),stringsAsFactors = FALSE,header = FALSE)
  # scaffold_mapping2 = read.csv(paste0(scaffold_dir,"scaffold_mapping_",acc2,".txt"),stringsAsFactors = FALSE,header = FALSE)
  # scaffold_mapping = rbind(scaffold_mapping1, scaffold_mapping2)
  # colnames(scaffold_mapping) = c("species","scaffold")
  # nrow(scaffold_mapping)

  # scaffold_mapping= scaffold_mapping %>% distinct(species,scaffold, .keep_all= TRUE)
  
  t1 = Sys.time()
  scaffold_species = sapply(data$scaffold,function(x){
    scaffold_mapping[[x]]
  })
  print(Sys.time() - t1)
  
  data$species = scaffold_species
  
  
  scaffold_info_species = unlist(sapply(scaffold_data1$scaffold,function(x){
    scaffold_mapping[[x]]
  }))
  
  scaffold_data1$species = scaffold_info_species
  
  
  print("number of contigs with no species identified")
  print(nrow(data) - nrow(data %>% filter(!is.na(species))))
  print("of")
  print(nrow(data))
  
  # found_species = data %>% filter(!is.na(species))
  # unfound_species = data %>% filter(is.na(species))
  # test = sapply(unfound_species$scaffold,function(x){
  #   
  #   unique_species = unique(scaffold_mapping[grepl(substr(x,1,10),scaffold_mapping$scaffold),"species"])
  #   if(length(unique_species) > 1){
  #     print(unique_species)
  #     print(x)
  #   }
  #   
  # })
  # scaffold_mapping %>% filter(species %in% c("Faecalibacterium_prausnitzii_61481", "Thalassiobium_sp_61960"))
  # 
  # 
  # test =  data %>% filter(is.na(species))
  # head(test)
  # nrow(scaffold_mapping)
  # intersect("AHZD01000806",scaffold_mapping$scaffold)
  # test = scaffold_mapping[grepl("AHZD010008", scaffold_mapping$scaffold ),]
  # nrow(test)
  
  
  #print("Size unknown scaffolds")
  
  prop_species_covered = c(prop_species_covered ,nrow(data %>% filter(!is.na(species)))/nrow(data))
  
  
  
  summary_data = data %>% group_by(species) %>% filter(!is.na(species)) %>% summarize(
    popANI_max = max(popANI,na.rm=TRUE),
    popANI_min = min(popANI, na.rm=TRUE),
    popANI_avg = mean(popANI, na.rm=TRUE))
  total_popANI_valid = nrow(summary_data %>% filter(!is.na(popANI_avg)))
  
  scaffold_info_summary = scaffold_data1 %>% group_by(species) %>% filter(!is.na(species)) %>% summarize(
    popANI_reference_avg = mean(popANI_reference,na.rm=TRUE),
    coverage_avg = mean(coverage, na.rm=TRUE))
  scaffold_infant_sufficient_coverage = nrow(scaffold_info_summary %>% filter(coverage_avg >= 5))
  
  
   
  total_passed = nrow(summary_data %>% filter(popANI_avg > threshold))
  
 
  
  prop_above_thresh = total_passed/total_popANI_valid
  prop_above_threshold = c(prop_above_threshold,prop_above_thresh)
  
  
  prop_of_infant_scaffolds = total_passed/scaffold_infant_sufficient_coverage
  prop_of_infant_scaffolds_list = c(prop_of_infant_scaffolds_list, prop_of_infant_scaffolds)
  
  family_id = unlist(metadata %>% filter(run_accession == acc1) %>% select(study_id))[1]
  true_prop = unlist(metadata %>% filter(study_id == family_id, run_accession == acc2) %>% select(True_prop))[1]

  true_props = c(true_props,true_prop)
  acc1_list = c(acc1_list,acc1)
  acc2_list = c(acc2_list,acc2)
  
}
length(true_props)
length(prop_above_threshold)
length(prop_species_covered)
length(acc2_list)
all_trial_frame = data.frame(true_props,prop_above_threshold,prop_species_covered,acc1_list,acc2_list)

#all_trial2 <- readRDS(paste0("/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/inStrain_compare_tables24/","seed24_data.rds"))
#all_trial_frame_ <- rbind(all_trial_frame, all_trial2)
#nrow(all_trial_frame_)
require(ggplot2)
p <- ggplot(all_trial_frame, aes(x = true_props, y =prop_above_threshold)) +

  geom_smooth(method=lm, fill="grey",show.legend = FALSE,size=0.5) +
  geom_point()+ geom_abline() +
  ggtitle("Trials") +
  theme_minimal() + xlab(paste0("True Percentage")) +
  ylab(paste0("Estimated Percentage")) + 
  theme(text = element_text(size=13),plot.title = element_text(size=13)) 


p
write.table(all_trial_frame,paste0(instrain_dir ,"instrain_percent_popANI_above_99999.csv"),sep = ",",
            col.names = TRUE,
            row.names= TRUE,quote = FALSE)

