require(dplyr)
#install.packages("dplyr")
simulation_d <- read.csv(paste0("/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/metadata/",
                                "metadata_merge_seed27.csv"))
head(simulation_d)
unique_moms <- unique(unlist(simulation_d %>% filter(cohort != "B") %>% select(run_accession)))

backhed <- read.csv(paste0("/Users/leahbriscoe/Documents/FEASTX/BackhedFiles/metadata_merge.csv"))
head(backhed)
family_of_unique_moms <- unlist(backhed %>% filter(run_accession %in% unique_moms) %>% select(study_id))

all_children_that_can_include <- unlist(backhed %>% filter(!(study_id %in% family_of_unique_moms), cohort != "M") %>% select(run_accession))
all_mothers <- unlist(backhed %>% filter(cohort  == "M") %>% select(run_accession))
nrow(backhed)

backhed_inclusion = c(all_mothers,all_children_that_can_include )
length(backhed_inclusion)


all_hmp_us = read.csv(paste0("/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/Nayfach_analysis/HMP_ids_order.txt"), 
sep = "\t")
head(all_hmp_us)
unique(all_hmp_us$country)
unique(all_hmp_us$subject_id)
hmp_inclusion = unique(all_hmp_us$sample_id)
length(hmp_inclusion)
length(c(backhed_inclusion,hmp_inclusion))
write.table(c(backhed_inclusion,hmp_inclusion), "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/Nayfach_analysis/all_inclusion.csv",
      sep = ",",quote=FALSE)
