
# ------------------------------------ #
# ----------- SUPPLEMENT 02  --------- #
# ------------------------------------ #

### HWE TESTS
### Author: Annabel Smith; Di Binyin

# Binyin: Load and tidy workspace and remove everything except necessary objects:
load("D:/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# Annabel: Load and tidy workspace and remove everything except necessary objects:
load("binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# This is the full data set:
ghead(snp_onerow); dim(snp_onerow)

# This is the partially filtered data from STEP 02:
ghead(filtered_data); dim(filtered_data)

## --- *** Run HWE tests *** --- ##

#----run sample data----#

# randomly sample 20 loci for testing:

samp_cols<-which(colnames(filtered_data)[3:ncol(filtered_data)] %in% sample(colnames(filtered_data)[3:ncol(filtered_data)],20))
rand_snp<-filtered_data[,c(1,2,samp_cols)]
rand_snp<-tidy.df(rand_snp)
ghead(rand_snp); dim(rand_snp) 

# this function does the exact tests:
hwe.res<-hwe_exact(rand_snp)
head(hwe.res); dim(hwe.res)

# when you have it working on the test data, do it for all loci

# write.table(hwe.res, file="HWE_test.txt", quote=F, sep="\t", row.names=F)


#----fun full dataset----#

ghead(filtered_data); dim(filtered_data)


hwe.res<-hwe_exact(filtered_data)
head(hwe.res); dim(hwe.res)

write.table(hwe.res, file="HWE_test.txt", quote=F, sep="\t", row.names=F)

#Results

#> # This is the partially filtered data from STEP 02:
#  > ghead(filtered_data); dim(filtered_data)
#site     ind L6 L7 L8 L9 L11 L12 L17 L25
#1 X01b X01b_01  2  1  2  0   0   0   0   0
#2 X01b X01b_02  2  1  2  0   0   0   0   0
#3 X01b X01b_03  2  1  2  0   0   0   0   0
#4 X01b X01b_04  2  1  2  0   0   0   0   0
#5 X01u X01u_01  1  2  0  0   0   0   0   0
#6 X01u X01u_02  2  1  2  0   0   0  NA   0
#[1]    93 40713
#> hwe.res<-hwe_exact(filtered_data)
#> head(hwe.res); dim(hwe.res)
#locus            p        p.adj
#1    L6 1.000000e+00 1.000000e+00
#2    L7 2.314964e-01 1.000000e+00
#3    L8 6.324177e-02 1.000000e+00
#4    L9 1.423518e-23 5.729233e-19
#5   L11 1.369863e-02 1.000000e+00
#6   L12 6.206069e-02 1.000000e+00
#[1] 40711     3




