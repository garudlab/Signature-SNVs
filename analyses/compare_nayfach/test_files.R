input_dir <- "/Users/leahbriscoe/Documents/FEASTX/BackhedFiles/snps/Blargasaurus_rex/"
inclusion_list <- "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/Nayfach_analysis/all_inclusion.csv"
inclusion_samples = read.csv(paste0(inclusion_list))[,1]

snp_depth  =read.csv(paste0(input_dir, "snps_depth.txt"),sep = "\t")
inclusion_samples = intersect(inclusion_samples,colnames(snp_depth))
write.table(snp_depth[,inclusion_samples],paste0(input_dir, "snps_depth.txt"),
            sep = "\t",quote=FALSE)
snp_freq  =read.csv(paste0(input_dir, "snps_ref_freq.txt"),sep = "\t")

write.table(snp_freq[,inclusion_samples],paste0(input_dir, "snps_ref_freq.txt"),
            sep = "\t",quote=FALSE)
