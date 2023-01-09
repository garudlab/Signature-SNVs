scaffold_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/scaffold_mappings/"

scaffold_mapping = read.csv(paste0(scaffold_dir,"all_scaffold_ids.csv"),stringsAsFactors = FALSE,header = FALSE)
colnames(scaffold_mapping) = c("species","scaffold")
nrow(scaffold_mapping)
scaffold_mapping= scaffold_mapping %>% distinct(species,scaffold, .keep_all= TRUE)
head(scaffold_mapping)
scaffold_map = list()
for(r in 1:nrow(scaffold_mapping)){
  if(r %% 100 == 0){
    print(r)
  }
  scaffold_map[[scaffold_mapping[r,"scaffold"]]] = scaffold_mapping[r,"species"]
}
saveRDS(object = scaffold_map, paste0(scaffold_dir, "scaffold_map.rds"))


