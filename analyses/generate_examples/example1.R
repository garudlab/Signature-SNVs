# Frequencies

ref_freq_data = read.csv("/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/midas_output/snps/Bacteroides_uniformis_57318/snps_ref_freq_pre.txt", 
         sep = "\t",header=TRUE)
head(ref_freq_data)
colnames(ref_freq_data) = c("site_id","Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")
# bakcgorund: almost all 1
# snps 1 to 20 are common (some low coverage though, for later) - they have lower than usual ref freq
# snps 21 to 30 are high in mother 1 and baby 1, rare in rest
# snps 31 to 40 are high in mother 2 and baby 2, rare in rest
# snps 41 to 50 are high in mother 3 and baby 3, rare in rest
# snps 51 to 60 are high in mother 4 and baby 4, rare in rest
# snps 61 to 30 are high in mother 1 only, rare in rest
# snps 71 to 40 are high in mother 2 only, rare in rest
# snps 81 to 50 are high in mother 3 only, rare in rest
# snps 91 to 100 are high in mother 4 only, rare in rest

# bakcgorund: almost all 0
# same pattern

# Set background
set.seed(123)
for(i in 1:100){
  ref_freq_data[i,c("Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")] = rbeta(8,shape1=1, shape2 = 0.01)
}

for(i in 101:200){
  ref_freq_data[i,c("Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")] = round(1-rbeta(8,shape1=1, shape2 = 0.01),7)
}

for(i in 21:30){
  ref_freq_data[i,c("Baby1","Mother1")] = runif(2,0.6,0.8)
}

for(i in 31:40){
  ref_freq_data[i,c("Baby2","Mother2")] = runif(2,0.6,0.8)
}

for(i in 41:50){
  ref_freq_data[i,c("Baby3","Mother3")] = runif(2,0.6,0.8)
}

for(i in 51:60){
  ref_freq_data[i,c("Baby4","Mother4")] = runif(2,0.6,0.8)
}

for(i in 61:70){
  ref_freq_data[i,c("Mother1")] = runif(1,0.6,0.8)
}

for(i in 71:80){
  ref_freq_data[i,c("Mother2")] = runif(1,0.6,0.8)
}

for(i in 81:90){
  ref_freq_data[i,c("Mother3")] = runif(1,0.6,0.8)
}

for(i in 91:100){
  ref_freq_data[i,c("Mother4")] = runif(1,0.2,0.8)
}



for(i in 121:130){
  ref_freq_data[i,c("Baby1","Mother1")] = runif(2,0.2,0.6)
}

for(i in 131:140){
  ref_freq_data[i,c("Baby2","Mother2")] = runif(2,0.2,0.6)
}

for(i in 141:150){
  ref_freq_data[i,c("Baby3","Mother3")] = runif(2,0.2,0.6)
}

for(i in 151:160){
  ref_freq_data[i,c("Baby4","Mother4")] = runif(2,0.2,0.6)
}

for(i in 161:170){
  ref_freq_data[i,c("Mother1")] = runif(1,0.2,0.6)
}

for(i in 171:180){
  ref_freq_data[i,c("Mother2")] = runif(1,0.2,0.6)
}

for(i in 181:190){
  ref_freq_data[i,c("Mother3")] = runif(1,0.2,0.6)
}

for(i in 191:200){
  ref_freq_data[i,c("Mother4")] = runif(1,0.2,0.8)
}



write.table(ref_freq_data,"/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/midas_output/snps/Bacteroides_uniformis_57318/snps_ref_freq.txt", 
            sep = "\t",quote= FALSE, row.names = FALSE)



# bakcgorund: either almost all 0 or almost all 1
# snps 1 to 40 are common (some low coverage though, for later)
# snps 41 to 50 are high in mother 1 and baby 1, rare in rest
# snps 51 to 60 are high in mother 2 and baby 2, rare in rest
# snps 61 to 70 are high in mother 3 and baby 3, rare in rest
# snps 71 to 80 are high in mother 4 and baby 4, rare in rest
# snps 81 to 85 are high in only mother 1 and baby 1
# snps 86 to 90 are high in only mother 2 and baby 2
# snps 91 to 95 are high in only mother 3 and baby 3
# snps 96 to 100 are high in only mother 4 and baby 4

snp_depth = read.csv("/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/midas_output/snps/Bacteroides_uniformis_57318/snps_depth_pre.txt", 
                         sep = "\t",header=TRUE)

colnames(snp_depth) = c("site_id","Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")
# bakcgorund: almost all 1
# snps 1 to 20 are common (some low coverage though, for later) - they have lower than usual ref freq
# snps 21 to 30 are high in mother 1 and baby 1, rare in rest
# snps 31 to 40 are high in mother 2 and baby 2, rare in rest
# snps 41 to 50 are high in mother 3 and baby 3, rare in rest
# snps 51 to 60 are high in mother 4 and baby 4, rare in rest
# snps 61 to 30 are high in mother 1 only, rare in rest
# snps 71 to 40 are high in mother 2 only, rare in rest
# snps 81 to 50 are high in mother 3 only, rare in rest
# snps 91 to 100 are high in mother 4 only, rare in rest

# bakcgorund: almost all 0
# same pattern

# Set background
set.seed(123)
for(i in 1:200){
  snp_depth[i,c("Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")] = rpois(n=8, lambda=10)
}
#pick out random low depth sites in 20% of sites

low_coverage_sites= sample(c(1:200),40 )
for(s in low_coverage_sites){
  snp_depth[s,c("Baby1","Baby2","Baby3","Baby4","Mother1","Mother2","Mother3","Mother4")] = rpois(n=8, lambda=2)
}

write.table(snp_depth,"/Users/leahbriscoe/Documents/FEASTX/snv-feast/example_1/midas_output/snps/Bacteroides_uniformis_57318/snps_depth.txt", 
            sep = "\t",quote= FALSE, row.names = FALSE)






