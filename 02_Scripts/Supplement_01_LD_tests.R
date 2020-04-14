
# ------------------------------------ #
# ----------- SUPPLEMENT 01  --------- #
# ------------------------------------ #

### LD TESTS
### Author: Annabel Smith & Di Binyin

# Binyin: Load and tidy workspace and remove everything except necessary objects:
load("D:/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# Annabel: Load and tidy workspace and remove everything except necessary objects:
load("binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# This is the full data set:
ghead(snp_onerow); dim(snp_onerow)

# Perhaps we should use the partially filtered data from STEP 01, which has 29007 loci? This will be quicker than doing the full data set and we are not interested in those loci which have already been filtered out. 
ghead(filtered_data); dim(filtered_data)

# To get the script working, you could run through a smaller data set, that won't take as long to run. Once you have it working you can work on the full filtered_data data set:

# randomly sample 20 loci for testing:
samp_cols<-which(colnames(filtered_data)[3:ncol(filtered_data)] %in% sample(colnames(filtered_data)[3:ncol(filtered_data)],20))
rand_snp<-filtered_data[,c(1,2,samp_cols)]
rand_snp<-tidy.df(rand_snp)
ghead(rand_snp); dim(rand_snp)

# format genotype file for the calc_LD function from evachang and run test:

sites_to_test<-levels(rand_snp$site)
dat_test<-rand_snp


param_file<-paste("RESULTS/LD_results/parameters","_LD.txt",sep="")

# if using the partially filtered data set, we've already removed monomorphic loci and loci with high levels of missing data. So we can move directly to the linkage tests. In this data set, all sites are in the same region, so we can look at the population as a whole. There is no need to divide the data into separate sites. 

# First remove the site and individual data:
dat_test<-dat_test[,3:ncol(dat_test)]

# geno needs to be m x n where m is the number of markers and n is the number of individuals:
dat_test<-t(dat_test)

ghead(dat_test); dim(dat_test)

ld.thisrun<-calc_LD(dat_test,inds=1:nrow(dat_test), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)

loc_combn<-combn(rownames(dat_test),2)

df_test<-data.frame(loc1=loc_combn[1,],loc2=loc_combn[2,],r2=ld.thisrun$rsq[lower.tri(ld.thisrun$rsq)])

# for 20 loci, there are 190 pairwise comparisons of linkage disequilibrium:
head(df_test); dim(df_test)

hist(df_test$r2)

# there will be MANY, MANY more for the full data set. Once you have it working on the test data set, you can proceed to work on the full data set. 

# we need to decide what a reasonable cut-off is for discarding linked loci. In the PNAS paper, we used 0.75

# the STEP_02 script pulls in a file of linked loci. If the file is very large, it might be good to reduce it to only those locus pairs above the cut-off:
df_test<-df_test[which(df_test$r2>0.75),]
df_test<-tidy.df(df_test)

# write.table(df_test, file="LD_r75_LOCI_FOR_REMOVAL", quote=F, sep="\t", row.names=T)










