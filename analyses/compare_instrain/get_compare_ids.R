data_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/metadata/"
out_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/accessions_list/"
accession_reference_path = paste0(data_dir,"metadata_merge_seed24.csv")
accession_reference = read.csv(accession_reference_path,stringsAsFactors = FALSE)
head(accession_reference)
require(dplyr)
babies_reference = accession_reference %>% filter(cohort == "B")
acc1 = c()
acc2 = c()
for(b in 1:nrow(babies_reference)){
  baby_acc = babies_reference$run_accession[b]
  baby_study_id = babies_reference$study_id[b]
 
  mom_data = accession_reference %>% filter(study_id == baby_study_id,!(cohort %in% c("UNK1","B")))
  for(m in 1:nrow(mom_data)){
    acc1 <- c(acc1,baby_acc)
    acc2<- c(acc2,mom_data$run_accession[m])
  }
  
}
length(acc1)
acc1 <- c(acc1,"ERR0S21T14N1")
acc2 <- c(acc2,"ERR525987")

write.table(acc1,paste0(out_dir,"instrain_compare1.txt"),sep = ",",
            col.names = FALSE,
            row.names= FALSE,quote = FALSE)

write.table(acc2,paste0(out_dir,"instrain_compare2.txt"),sep = ",",
            col.names = FALSE,
            row.names= FALSE,quote = FALSE)
length(unique(accession_reference$run_accession))
length(unique(acc1))
length(unique(acc2))
