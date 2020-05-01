
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

#Running LD Scripts

# Use  partially filtered data from STEP 01 - this will be quicker than doing the full data set and we are not interested in those loci which have already been filtered out. 

#BD note: Use the partially fitered data from STEP 02, after malim ==0.05 process
ghead(filtered_data); dim(filtered_data)

### TEST SCRIPT on sample data:

#Warning: Sample dataset were used in following scripts

# To get the script working, run through a smaller data set, that won't take as long to run. Once you have it working you can work on the full filtered_data data set:

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

### LD ANALYSIS on FULL dataset:

#rand_snp replaced by filtered_data for full data run. 
sites_to_test<-levels(filtered_data$site)
dat_test<-filtered_data
param_file<-paste("RESULTS/LD_results/parameters","_LD.txt",sep="")
dat_test<-dat_test[,3:ncol(dat_test)]

# geno needs to be m x n where m is the number of markers and n is the number of individuals:
dat_test<-t(dat_test)

ghead(dat_test); dim(dat_test)

ld.thisrun<-calc_LD(dat_test,inds=1:nrow(dat_test), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)

loc_combn<-combn(rownames(dat_test),2)

df_test<-data.frame(loc1=loc_combn[1,],loc2=loc_combn[2,],r2=ld.thisrun$rsq[lower.tri(ld.thisrun$rsq)])

# for all loci, there are 182739403 pairwise comparisons of linkage disequilibrium:
head(df_test); dim(df_test)
hist(df_test$r2)

# the file is very large, it might be good to reduce it to only those locus pairs above the cut-off:
df_test<-df_test[which(df_test$r2>0.5),]
df_test<-tidy.df(df_test)

# This step only tells us which loci are linked, not which loci we should remove (the next script does that). 

# Write a file of all the linked loci, reduced to a smaller file size. Previously we were calling this "LD_r50_LOCI_FOR_REMOVAL", but really it should be "LD_r50_linked_loci", because we haven't yet determined which to remove. 

# write.table(df_test, file="LD_r50_LOCI_FOR_REMOVAL", quote=F, sep="\t", row.names=T)

### LD ANALYSIS: Determine which of the linked pairs we should remove:

##BD's Script
LD_dir<-"D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
LD_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
dir(LD_dir)
ld_loc<-read.table(paste(LD_dir, "LD_r75_filtered_data",sep="/"),header=T)

# AS LD: 
LD_dir_AS<-"../Offline_Results/LD_results"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r70_LOCI_FOR_REMOVAL",sep="/"),header=T)
head(ld_loc); dim(ld_loc)

# Decide  a reasonable cut-off is for  linked loci. In the PNAS paper, we used 0.75
# Set cutoff at 0.75:
ld_loc<-ld_loc[which(ld_loc$r2>0.75),]
ld_loc<-tidy.df(ld_loc)
head(ld_loc); dim(ld_loc)

# The determine which loci would need to be removed to ensure no linked loci would occur together:

# Total number of loci in the linked data set:
all_ldloc<-c(as.character(ld_loc$loc1),as.character(ld_loc$loc2))
length(unique(all_ldloc))
head(all_ldloc)

# Frequency with which each of the linked loci occurs:
freq_loc<-data.frame(locus=names(table(all_ldloc)),no_times_total=as.numeric(table(all_ldloc)))
head(freq_loc); dim(freq_loc)
head(ld_loc); dim(ld_loc)

# The following script removes a locus from a linked pair, prioritising the locus which is more frequently linked to other loci:

# start 1000h, finish 2000h - next day!! 24+10 = 34 hr! And this is on the big mac. Something is wrong... Might be that it's online constantly connecting to Github, I need to do some tests to figure out what's going on. 

save.image("../Offline_Results/LD_selection.RData")

removed.loci<-list()

for (i in 1:nrow(ld_loc)){
  
  print(paste("Starting test", i, "of", nrow(ld_loc), "tests", sep=" "))
  
  line.thisrun<-ld_loc[i,]
  l1.thisrun<-as.character(line.thisrun$loc1)
  l2.thisrun<-as.character(line.thisrun$loc2)
  
  freq_l1<-freq_loc$no_times_total[which(freq_loc$locus==l1.thisrun)]
  freq_l2<-freq_loc$no_times_total[which(freq_loc$locus==l2.thisrun)]
  
  # If one of the pair has already been assigned to the "remove" pile, then the pair is OK and can skip to the next test
  if(length(which(c(l1.thisrun,l2.thisrun) %in% unlist(removed.loci)==T))>0) next
  
  # If they have the same frequency in the linkage summary, remove l2:
  if (freq_l1==freq_l2) {
    removed.loci[[i]]<-l2.thisrun
    next
  } # close if same freq
  
  # If they have a different frequency, remove the one with the higher frequency:
  if(freq_l1!=freq_l2) {
    removed.loci[[i]]<-c(l1.thisrun,l2.thisrun)[which(c(freq_l1,freq_l2)==max(freq_l1,freq_l2))]
  } # close different freq
  
} # close for

Sys.time()
save.image("../Offline_Results/LD_selection.RData")

# Summarise results:
rm_loci<-data.frame(locus=unlist(removed.loci))
rm_loci$for_removal<-1
loci_toremove<-as.character(rm_loci$locus)
head(loci_toremove)
length(loci_toremove)
length(removed.loci)
head(rm_loci)
tail(rm_loci)

# write.table(rm_loci,"LD_r75_loci_to_remove.txt",sep="\t",row.names=F,quote=F)

# Check:
ck_rm<-merge(ld_loc,rm_loci,by.x="loc1",by.y="locus",all.x=T, all.y=F)
colnames(ck_rm)[length(colnames(ck_rm))]<-"l1_removal"
ck_rm$l1_removal[which(is.na(ck_rm$l1_removal))]<-0
ck_rm<-merge(ck_rm,rm_loci,by.x="loc2",by.y="locus",all.x=T, all.y=F)
colnames(ck_rm)[length(colnames(ck_rm))]<-"l2_removal"
ck_rm$l2_removal[which(is.na(ck_rm$l2_removal))]<-0

# Sometimes both will be removed, because they both occur somewhere else
check.rows(ck_rm)
head(ck_rm)

# But every line should add to at least 1:
range(rowSums(ck_rm[,which(colnames(ck_rm) %in% c("l1_removal","l2_removal"))]))

# This shows that at least one of each pair will be removed and that no linked pair will be included in the final data set. 

save.image("../Offline_Results/LD_selection.RData")






