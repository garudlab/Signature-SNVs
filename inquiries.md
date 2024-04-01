bzcat private_snvs_merged/CATTED_fam_SRF_Above20_IO_TARA_B000000609_seed_9_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_counts.csv.bz2  | awk '$1 ~/^Alt_Alpha/'


bzcat private_snvs_merged/CATTED_fam_SRF_Above20_IO_TARA_B000000609_seed_9_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_counts.csv.bz2  | awk '$1 ~/^Alt_Alpha/'


bzcat private_snvs_merged/CATTED_fam_SRF_Above20_IO_TARA_B000000609_seed_9_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_counts.csv.bz2 | head -n 3 | cut -d "," -f1-4

bzcat private_snvs_merged/CATTED_fam_SRF_Above20_IO_TARA_B000000609_seed_9_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_counts.csv.bz2  | awk '$1 ~/^Ref_Alpha_proteobacterium_62227[|]CP003801[|]16/' | head -n 3 | cut -d "," -f1-4


bzcat private_snvs_merged/CATTED_fam_SRF_Above20_IO_TARA_B000000609_seed_9_uniqueness_1_threshold_1.0_familythresh_0.1_minreads_5_counts.csv.bz2  | awk '$1 ~/^Ref_Alpha_proteobacterium_62227[|]CP003801[|]2233/'  | head -n 3 | cut -d "," -f1-4

cat Alpha_proteobacterium_62227/snps_depth.txt | awk '$1 ~/^CP003801[|]1663[|]T/' | cut -d$'\t' -f1-3
cat Alpha_proteobacterium_62227/snps_depth.txt | awk '$1 ~/^CP003801[|]1696[|]A/' | cut -d$'\t' -f1-3
cat Alpha_proteobacterium_62227/snps_depth.txt | awk '$1 ~/^CP003801[|]2233[|]T/' | cut -d$'\t' -f1-3

cat Alpha_proteobacterium_62227/snps_ref_freq.txt | awk '$1 ~/^CP003801[|]1663[|]T/' | cut -d$'\t' -f1-3
cat Alpha_proteobacterium_62227/snps_ref_freq.txt | awk '$1 ~/^CP003801[|]1696[|]A/' | cut -d$'\t' -f1-3
cat Alpha_proteobacterium_62227/snps_ref_freq.txt | awk '$1 ~/^CP003801[|]2233[|]T/' | cut -d$'\t' -f1-3


IO_TARA_B000000609_seed_9_

SRF_Above20_IO_TARA_B000000609,ERR599057,ERR598993,ERR315862,ERR315857,ERR594288,ERR599123,ERR599029,ERR599136,ERR598955,ERR599142,ERR598989,ERR599030,ERR599162,ERR599116,ERR599049,ERR599106,ERR598991,ERR594328,ERR598984,ERR599126,ERR594286,ERR594310,ERR599038,ERR599069,ERR599160,ERR598954,ERR598992,ERR599119,ERR594339,ERR594344,ERR594296,ERR594311,ERR594287,ERR594347,ERR594326,ERR594292,ERR594306,ERR594307


first 3 columns: ERR599057,ERR598993,ERR315862
 in snps_depth.txt 2,nan,ERR315862
Alt_Alpha_proteobacterium_62227|CP003801|1663|T,3.0,nan,71.0
Alt_Alpha_proteobacterium_62227|CP003801|1696|A,1.0,nan,87.0
Alt_Alpha_proteobacterium_62227|CP003801|2233|T,1.0,nan,0.0
Ref_Alpha_proteobacterium_62227|CP003801|1663|T,17.0,nan,0.0
Ref_Alpha_proteobacterium_62227|CP003801|1696|A,27.0,nan,1.0
Ref_Alpha_proteobacterium_62227|CP003801|2233|T,11.0,nan,121.0


the snps_depth.txt

CP003801|1663|T	14	71
CP003801|1696|A	19	88
CP003801|2233|T	92	121


CP003801|1663|T	0.5714285714285714	0.0
CP003801|1696|A	0.42105263157894735	0.011363636363636364
CP003801|2233|T	1.0	1.0





snps_ref_freq.txt


|CP003801|16/



tHe file name iRef_Alpha_proteobacterium_62227|CP003801|1663|T,17.0,nan,0.0
Ref_Alpha_proteobacterium_62227|CP003801|1696|A,27.0,nan,1.0



python <site-packages_directory>/signature_snvs/signature_snvs_cli.py --species Bacteroides_uniformis_57318_short --min_reads 5 --start_index 1 --end_index 200 --config_file_path config.yaml
