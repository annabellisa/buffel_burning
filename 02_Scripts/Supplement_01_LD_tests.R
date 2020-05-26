
# ------------------------------------ #
# ----------- SUPPLEMENT 01  --------- #
# ------------------------------------ #

### LD TESTS
### Author: Annabel Smith & Di Binyin

# Binyin: Load and tidy workspace and remove everything except necessary objects:
load("D:/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# Annabel: Load and tidy workspace and remove everything except necessary objects:
load("../Offline_Results/LD_40711_loci/LD_40711_loci.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat", "filtered_data")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))
# library(parallel)

# This is the full data set:
ghead(snp_onerow); dim(snp_onerow)

# This is the filtered data set on which to do the LD tests (40711 loci)
ghead(filtered_data); dim(filtered_data)

### ^^^^ ____________
# Data set-up:
### ^^^^ ____________

sites_to_test<-levels(filtered_data$site)
dat_ld<-filtered_data
dat_ld<-dat_ld[,3:ncol(dat_ld)]
ghead(dat_ld); dim(dat_ld)

# Make a parameter file so we know what filters were applied before the tests: test change
# param_file<-paste("../Offline_Results/LD_40711_loci/LD_40711_parameters.txt",sep="")
# write("LD tests, round 2",file=param_file,sep="")
# write(as.character(print(Sys.time())),file=param_file,append = T)
# write(paste("number of sites = ",length(sites_to_test),sep=""),file=param_file,append = T)
# write(paste("sites:", sep=""),file=param_file,append = T)
# write(paste(sites_to_test,collapse=", "),file=param_file,append = T)
# write(paste("number of loci = ",ncol(dat_ld),sep=""),file=param_file,append = T)
# write(paste("filters applied = dups removed, mono removed, max missing= 50%"),file=param_file,append = T)

# geno needs to be m x n where m is the number of markers and n is the number of individuals:
dat_ld<-t(dat_ld)
ghead(dat_ld); dim(dat_ld)

### ^^^^ ____________
# LD TESTS
### ^^^^ ____________

# these are the actual tests:
# 2hr 40min on big mac (40711 loci)
print("LD test START:")
print(Sys.time())
ld.thisrun<-calc_LD(dat_ld,inds=1:nrow(dat_ld), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)
print("LD test END")
print(Sys.time())

# Create pairwise locus column names:
# 12 min big mac
print(Sys.time())
loc_combn<-combn(rownames(dat_ld),2)
print(Sys.time())

# Summarise results:
ld_df<-data.frame(loc1=loc_combn[1,],loc2=loc_combn[2,],r2=ld.thisrun$rsq[lower.tri(ld.thisrun$rsq)])

# save.image("../Offline_Results/LD_40711_loci/LD_40711_loci.RData")

# for all loci, there are 828672405 pairwise comparisons of linkage disequilibrium for 40711 loci (182739403 for ~19,000 loci):
head(ld_df); dim(ld_df)
# hist(ld_df$r2)
# table(ld_df$r2>0.5)
 
# the file is very large (the workspace is 15GB), reduce:
ld_df<-ld_df[which(ld_df$r2>0.5),]
ld_df<-tidy.df(ld_df)

# save.image("../Offline_Results/LD_40711_loci/LD_40711_loci.RData")
# reduced size is 9GB - I've deleted this since the results were saved in the text file

# This step only tells us which loci are linked, not which loci we should remove (the next script does that). 

# Write a file of all the linked loci, reduced to a smaller file size. Previously we were calling this "LD_r50_LOCI_FOR_REMOVAL", but really it should be "LD_r50_linked_loci", because we haven't yet determined which to remove. 

# write.table(ld_df, file="../Offline_Results/LD_40711_loci/LD_r50_linked_loci", quote=F, sep="\t", row.names=T)

### ^^^^ ____________
# LD locus selection
### ^^^^ ____________

# AS LD: 
LD_dir_AS<-"../Offline_Results/LD_40711_loci"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r75_linked_loci.txt",sep="/"),header=T)
head(ld_loc); dim(ld_loc)

# I saved 0.75 as a separate file, to reduce the amount of data stored in the workspace when running. Re-load from the first lines in this section. 
# Decide  a reasonable cut-off is for  linked loci. In the PNAS paper, we used 0.75
# Set cutoff at 0.75:
# ld_loc<-ld_loc[which(ld_loc$r2>0.75),]
# ld_loc<-tidy.df(ld_loc)
# head(ld_loc); dim(ld_loc)
# write.table(ld_loc, file="../Offline_Results/LD_40711_loci/LD_r75_linked_loci.txt", quote=F, sep="\t", row.names=T)

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

# First run from limited (19,000) data set:
# start 1000h, finish 2000h - next day!! 24+10 = 34 hr! Have done speed tests and nothing is wrong - the only way to speed it up is to learn to use parallel()

# Starting now (start time was approx 1345 21 May 2020):
# This run didn't work. Computer froze and nothing was saved.  

# Starting again at 26 May 2020, 0910:

Sys.time()
save.image("../Offline_Results/LD_40711_loci/LD_selection.RData")

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

ld_select<-function(pw_data, freq_data, line_to_test){

  # pw_data = the r2 values for all pairs of loci
  # freq_data = the frequency of occurrence for each locus
  # lines_to_test = a sequence (e.g. 1:5000) of lines in pw_data to test

line.thisrun<-pw_data[line_to_test,]
l1.thisrun<-as.character(line.thisrun$loc1)
l2.thisrun<-as.character(line.thisrun$loc2)

freq_l1<-freq_data$no_times_total[which(freq_data$locus==l1.thisrun)]
freq_l2<-freq_data$no_times_total[which(freq_data$locus==l2.thisrun)]

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

return(removed.loci)
  
} # close function

head(freq_loc); dim(freq_loc)
head(ld_loc); dim(ld_loc)

removed.loci<-list()

for (i in 1:5000){

line_to_test<-i
xx<-ld_select(ld_loc, freq_loc, line_to_test)

} # close for loop

str(xx)
head(xx)

rm_loci<-data.frame(locus=unlist(xx))
rm_loci$for_removal<-1
loci_toremove<-as.character(rm_loci$locus)
head(loci_toremove)
length(loci_toremove)
head(rm_loci); dim(rm_loci)
tail(rm_loci)

head(xx,20)
head(xx[unlist(which(lapply(xx,is.null)==T))])

tail(xx)
sample(xx, 6 )


foreach (i= 1:10000, .combine=cbind) %dopar% {
  
ld_select(ld_loc, freq_loc)
  
} # close foreach

# Finish time:
Sys.time()
save.image("../Offline_Results/LD_40711_loci/LD_selection.RData")




# the original script (it works!) but it's too slow:

removed.loci<-list()

head(ld_loc)

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

# Finish time:
Sys.time()
save.image("../Offline_Results/LD_40711_loci/LD_selection.RData")


# up to here:
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

save.image("../Offline_Results/LD_40711_loci/LD_selection.RData")











### **************** 
# BELOW is SCRIPT DEVELOPMENT ONLY.
# We did not use the full panel of loci on this analysis. Even the "full data set" was only 19,000 loci, so we are re-doing the tests on the complete data set. 
### ****************

#Running LD Scripts

# Use  partially filtered data from STEP 01 - this will be quicker than doing the full data set and we are not interested in those loci which have already been filtered out. 

#BD note: Use the partially fitered data from STEP 02, after malim ==0.05 process
ghead(filtered_data); dim(filtered_data)

### TEST SCRIPT on sample data:

#Warning: Sample dataset were used in following scripts

# To get the script working, run through a smaller data set, that won't take as long to run. Once you have it working you can work on the full filtered_data data set:

# randomly sample 20-2000 loci for testing:
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

print("LD test START:")
print(Sys.time())

ld.thisrun<-calc_LD(dat_test,inds=1:nrow(dat_test), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)

print("LD test END")
print(Sys.time())

# 2000 loci = 2 min

# trying to get this to work with mclapply, but not sure how it works yet:

f <- function(i) {
  nrow(dat_test)
}

f()
x1 <- mclapply(1:nrow(dat_test), f)

save1 <- mclapply(1:100, f)

print("LD test START:")
print(Sys.time())

ld.thisrun<-calc_LD(dat_test,inds=1:nrow(dat_test), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)

print("LD test END")
print(Sys.time())


# this goes back to the normal script:

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

#save the file by above the cutoff
ldf<-0.75
#Method 1.#
install.packages("dplyr")
library(dplyr)
ldfilt75<-ld_loc %>%
  filter(r2 > ldf)
head(ldfilt75)
write.table(ldfilt75, file = "LD_r75_filtered_data", quote = F, sep = "\t", row.names = T)

#Method 2.#
ldlift75<-ld_loc[which(ld_loc$r2>0.75),]
ldlift75<-tidy.df(ldlift75)
write.table(ldlift75, file="LD_r75_another_method", quote=F, sep="\t", row.names=T)


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






