
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

# This is the full data set:
ghead(snp_onerow); dim(snp_onerow)

# This is the filtered data set on which to do the LD tests (40711 loci)
ghead(filtered_data); dim(filtered_data)

### ^^^^ ____________
# Data set-up:
### ^^^^ ____________

# Data set-up:
sites_to_test<-levels(filtered_data$site)
dat_ld<-filtered_data
dat_ld<-dat_ld[,3:ncol(dat_ld)]
ghead(dat_ld); dim(dat_ld)

# Make a parameter file so we know what filters were applied before the tests:
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
# these are the actual tests, but the mac is slow and we need to figure out why first:
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

# for all loci, there are 828672405 pairwise comparisons of linkage disequilibrium for 40711 loci (there were previously 182739403 for ~19,000 loci):
head(ld_df); dim(ld_df)
# hist(ld_df$r2)
# table(ld_df$r2>0.5)
 
# the file is very large (the workspace is 15GB), reduce:
ld_df<-ld_df[which(ld_df$r2>0.5),]
ld_df<-tidy.df(ld_df)

# save.image("../Offline_Results/LD_40711_loci/LD_40711_loci.RData")

# reduced size is 9GB - I've deleted this since the results were saved in the text file
head(ld_df); dim(ld_df)
hist(ld_df$r2)

# the file is very large, reduce it to only those locus pairs above the cut-off:
ld_df<-ld_df[which(ld_df$r2>0.5),]
ld_df<-tidy.df(ld_df)

# This step only tells us which loci are linked, not which loci we should remove (the next section does that). 

# Write a file of all the linked loci, reduced to a smaller file size. Previously we were calling this "LD_r50_LOCI_FOR_REMOVAL", but really it should be "LD_r50_linked_loci", because we haven't yet determined which to remove. 

# write.table(ld_df, file="../Offline_Results/LD_40711_loci/LD_r50_linked_loci", quote=F, sep="\t", row.names=T)

### ^^^^ ____________
# LD locus selection
### ^^^^ ____________

# Decide  a reasonable cut-off is for  linked loci. In the PNAS paper, we used 0.75

dir("../Offline_Results/LD_40711_loci")

# no need to run this section if loading the LD_selection.RData workspace below
LD_dir_AS<-"../Offline_Results/LD_40711_loci"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r75_linked_loci.txt",sep="/"),header=T)
head(ld_loc); dim(ld_loc)

# this workspace has the full data set of 0.75 pw comparisons (16235 individual loci, 15040523 pw comparisons)
load("../Offline_Results/LD_40711_loci/LD_selection.RData")

# The determine which loci would need to be removed to ensure no linked loci would occur together:

# the new "quick method":

# Method: We first removed any locus that was linked more than 100 times with another locus (12385 out of 16235 loci). For the remaining 3850 loci, we systematically evaluated each pair of linked loci, and removed the locus with the highest frequency of linkage to other loci in the data set. This process identified 2722 loci that were  unlinked to any other locus after the more highly linked loci were removed. The remaining 13513 loci were removed from all analyses that assumed independence of loci. 

# Total number of loci in the full data set:
all_ldloc<-c(as.character(ld_loc$loc1),as.character(ld_loc$loc2))
length(unique(all_ldloc))
head(all_ldloc)

# Frequency with which each of the linked loci occurs:
freq_loc<-data.frame(locus=names(table(all_ldloc)),no_times_total=as.numeric(table(all_ldloc)))
head(freq_loc); dim(freq_loc)
head(ld_loc); dim(ld_loc)

# order the loci by the frequency with which they are correlated with others:
f2loc<-freq_loc[order(freq_loc$no_times_total,decreasing=T),]
f2loc<-tidy.df(f2loc)
head(f2loc); dim(f2loc)

# plot(f2loc$no_times_total)
# hist(f2loc$no_times_total)

# remove any locus that is highly linked (e.g. with 100 other loci). Removing loci that are very highly linked reduced the data set to a more computationally managable size:
freq_cutoff<-100
highly_linked<-as.character(f2loc$locus[f2loc$no_times_total>freq_cutoff])

# Of the 16235 linked loci, 12385 of them are highly linked. 
head(highly_linked); length(highly_linked)
head(which(duplicated(highly_linked)))

# This leaves 3850 loci with low levels of linkage (low_linked):
f3loc<-f2loc[f2loc$no_times_total<=freq_cutoff,]
f3loc<-tidy.df(f3loc)
# plot(f3loc$no_times_total)
# hist(f3loc$no_times_total)
head(f3loc); dim(f3loc)

low_linked<-as.character(f3loc$locus)
head(low_linked); length(low_linked)

# update the pw comparison data, so it only includes the relevant loci:
head(ld_loc); dim(ld_loc)

# if we include ALL of the low_linked markers here (e.g. c(which(ld_loc$loc1 %in% f3loc$locus),which(ld_loc$loc2 %in% f3loc$locus))), we end up with low_linked markers paired with high_linked markers, thus, more rows than we need. So we need to REMOVE highly_linked from each row in the ld_loc data set:
ld_loc<-ld_loc[-c(which(ld_loc$loc1 %in% highly_linked),which(ld_loc$loc2 %in% highly_linked)),]
ld_loc<-tidy.df(ld_loc)
head(ld_loc); dim(ld_loc)

# with 100 as the cutoff, this reduces from 15040523 pw comparisons to 19250. 

# a problem here is that this removes some of the low_linked loci where all of their pw comparisons were with high_linked loci. So, after updating the frequencies, we need to add them back in as zeros.

# update freq_loc so it only includes the new data set:

# Total number of loci in the new data set:
all_ldloc<-c(as.character(ld_loc$loc1),as.character(ld_loc$loc2))
length(unique(all_ldloc))
head(all_ldloc)

# Frequency with which each of the linked loci occurs:
freq_loc<-data.frame(locus=names(table(all_ldloc)),no_times_total=as.numeric(table(all_ldloc)))
head(freq_loc); dim(freq_loc)

# create df of the low linked loci that do NOT appear in the new freq_loc data set; these are the ones that should be zero - i.e. all of their pairs were with high_linked loci, that have already been removed. 

# ll2 includes all loci which should be zero
ll2<-low_linked[which(!low_linked %in% freq_loc$locus)]
head(ll2); length(ll2)

#  and its length, plus the length of the new frequency dataset give the total number of low_linked loci: 2078+1772==3850; this should be TRUE:
length(ll2)+nrow(freq_loc)==nrow(f3loc)

ll_df<-data.frame(locus=ll2,no_times_total=0)
head(ll_df); dim(ll_df)

# rbind both datasets:
freq_loc<-rbind(freq_loc, ll_df)
range(freq_loc$no_times_total)
# hist(freq_loc$no_times_total)

head(ld_loc); dim(ld_loc)
head(freq_loc); dim(freq_loc)

# then run the main script (quick method only takes 1 min)

# The following script removes a locus from a linked pair, prioritising the locus which is more frequently linked to other loci:

save.image("../Offline_Results/LD_40711_loci/LD_quick_method.RData")

# add removed loci to the pile:

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
save.image("../Offline_Results/LD_40711_loci/LD_quick_method.RData")

head(removed.loci); length(removed.loci)

# Summarise results:
rm_loci<-data.frame(locus=unlist(removed.loci))
rm_loci$for_removal<-1
loci_toremove<-as.character(rm_loci$locus)
head(loci_toremove); length(loci_toremove)

# Add frequencies to removed loci to check where they sat on the spectrum of linkage level (there is quite a range):
rm_loci<-merge(rm_loci, freq_loc, by="locus", all.x=T, all.y=F)
head(rm_loci); dim(rm_loci)
hist(rm_loci$no_times_total)

# this process identified 1128 loci to be removed (in addition to the "highly_linked" markers identified previously); the process thus "saved" 2722 loci from removal (nrow(f3loc)-nrow(rm_loci)) from a total of 16235 that were used in the LD selection process
saved_loci<-data.frame(locus=as.character(freq_loc$locus[which(!freq_loc$locus %in% rm_loci$locus)]))
saved_loci<-merge(saved_loci, freq_loc, by="locus", all.x=T, all.y=F)

# apart from one or two exceptions, most of the "saved" loci had very low levels of linkage with other loci
# this essentially "validates" the process of removing the highly_linked markers to begin with, without checking if any of them could be saved. None of the saved markers were linked with > 100 other markers and only one or two had high-ish levels of linkage (50, 76, etc). We might have lost one or two from the mass removal of highly_linked, but unlikely to make much difference, and a reasonable trade-off given the savings in computational time. 
range(saved_loci$no_times_total)
mean(saved_loci$no_times_total)
hist(saved_loci$no_times_total)
head(saved_loci); dim(saved_loci)

# the final step is to add all highly linked loci to the rm_loci df:
hgll<-data.frame(locus=highly_linked, for_removal=1)
head(hgll); dim(hgll)

# we've overwritten the frequencies for highly_linked, but since we no longer need this column (it was just for checking), we can remove it from the main df and rbind them together:
rm_loci$no_times_total<-NULL
rm_loci<-rbind(rm_loci, hgll)
head(rm_loci); dim(rm_loci)

# the final result is 13513 loci that will need to be removed from the neutral analyses

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

save.image("../Offline_Results/LD_40711_loci/LD_quick_method.RData")




