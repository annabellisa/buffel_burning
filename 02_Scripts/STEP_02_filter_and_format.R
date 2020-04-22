
# ------------------------------------ #
# ------------- STEP 02  ------------- #
# ------------------------------------ #

### filter SNP loci & format for different software
### Author: Annabel Smith & Di Binyin

#Set Working dir
setwd("D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter")
setwd("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter")

# Binyin: Load and tidy workspace and remove everything except necessary objects:
load("D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# Annabel: Load and tidy workspace and remove everything except necessary objects:
load("binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

#########################################
####	     FULL DATA SET:    		 ####
#########################################

###-->> Set data:
data_name<-"snp_onerow"
filtered_data<-get(data_name)

# --- *** Discard duplicated *** --- #
dup_loc<-as.character(linf$locus[linf$duplicate==1])
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% dup_loc)]
filtered_data<-tidy.df(filtered_data)
ghead(filtered_data); dim(filtered_data)

####	     FILTER LOCI:    		 ####

# --- *** Filter monomorphic loci *** --- #
filtered_data<-mono_loci(filtered_data,3)
ghead(filtered_data); dim(filtered_data)

# --- *** DartSeq Quality Control (QC) filters *** --- #

# Filter loci with high missing data rate (see remarks in missing_data function):
###-->> Set maximum missing data:
##missing_data==1-CallRate
cr<-0.2
missing_sum<-missing_data(filtered_data,3,cr)
m_summary<-missing_sum$miss_sum
#hist(m_summary$missing_data)
filtered_data<-missing_sum$filt_dat

# Filter loci with low reproducibility:
###-->> Set RepAvg:
ra<-0.95
repavg95<-linf$locus[which(linf$RepAvg<ra)]
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% repavg95)]
filtered_data<-tidy.df(filtered_data)
ghead(filtered_data); dim(filtered_data)



# --- *** Minor Allele Frequency (MAF) filters *** --- #
maf_sum<-maf_summary(filtered_data)
head(maf_sum)
#hist(maf_sum$maf[maf_sum$maf<0.1])

# Filter loci with extreme maf:
###-->> Set maf limit:
malim<-0.05
filtered_data<-maf_filter(maf_sum,filtered_data,malim)
ghead(filtered_data); dim(filtered_data)

save.image("Partially Filtered Data.RData")

# --- ***Linkage disequilibrium (LD) filters *** --- #

----#Supplement_01_LD_test.R#----

##BD's Script
LD_dir<-"D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
LD_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
dir(LD_dir)
ld_loc<-read.table(paste(LD_dir, "LD_r50_LOCI_FOR_REMOVAL",sep="/"),header=T)

# AS LD: 
LD_dir_AS<-"../Offline_Results/LD_results"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r50_LOCI_FOR_REMOVAL",sep="/"),header=T)

head(ld_loc)
View(ld_loc)

#Steps
#1. Read in the file you saved with the LD_results as a data frame
#2. Set the LD cutoff - 0.7 or 0.75, or whatever you decide
#3. Identify the loci from the results which occur above this cutoff
#4. Subset the SNP data set to exclude the loci in disequilibrium

ldf<-0.75
#Method 1.#
install.packages("tidyverse")
library(tidyverse)
library(tibble)
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

# Hi Binyin - what you're doing above is saving the LD results to a new file, which is good. But the following lines are needed to do the actual filtering. I.e. now that we've got the LD results saved, we can use them to filter out the loci that are too correlated, according to our cut-off. 


----#Changes#----
##BD's Script
LD_dir<-"D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
LD_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
dir(LD_dir)
ld_loc<-read.table(paste(LD_dir, "LD_r75_filtered_data",sep="/"),header=T)

# AS LD: 
LD_dir_AS<-"../Offline_Results/LD_results"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r75_filtered_data",sep="/"),header=T)

ldfilt<-as.character(ld_loc$loc2)
print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-tidy.df(filtered_data)
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% ldfilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data)
save.image("binyin_winter.RData")

# --- *** HWE filters *** --- #

#Supplement_02_HWE_test.R#

#Annabel
hwe_dir<-"../../ANALYSIS_RESULTS/ALL_pops_HWE_test"
dir(hwe_dir)
hwe_res<-read.table(paste(hwe_dir,"hwe_all_pops.txt",sep="/"),header=T)
head(hwe_res)

#Binyin
hwe_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/HWE_results"
dir(hwe_dir)

#See-sup#2
ghead(filtered_data); dim(filtered_data)
hwe.res<-hwe_exact(filtered_data)
head(hwe.res); dim(hwe.res)
write.table(hwe.res, file="HWE_test1.txt", quote=F, sep="\t", row.names=F)

# we need to make some changes here. The previous project I did, we looked for loci that were consistently out of HWE in > 5 populations. For ours however. we're treating them as a single population so I think we can just go with the p value. 

# I've made some preliminary changes here, but we need to look at the results, check the script and make some decisions:

# Filter loci with HWD:
hwe_flag<-T
hwe_cutoff<-0.1 # we need to decide on the cutoff
hwefilt<-as.character(hwe.res$locus[hwe.res$p.adj<hwe_cutoff])
print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% hwefilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

save.image("HWE_filter_data.Rdata")

# --- *** NEUTRALITY filter *** --- #

ghead(filtered_data); dim(filtered_data)

# Directory with results:
sel_dir<-"../../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION"
dir(sel_dir)

# Outlier loci:
res<-read.table(paste(sel_dir,"outliers_all.txt",sep="/"),header=T)
head(res)

# Filter outlier loci from BayeScan, PCAdapt and LFMM (not bayenv):
outl_loci<-as.character(res$locus[c(which(res$bs_outl==1),which(res$pca_outl==1),which(res$lfmm_outl==1))])
outl_loci<-outl_loci[-which(duplicated(outl_loci))]
head(outl_loci)
length(outl_loci)

# Filter outlier loci:
neutral_flag<-T
print(paste("no loci before neutral filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% outl_loci)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after neutral filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# Create data set with ONLY non-neutral loci:
neutral_flag<-F
print(paste("no loci before neutral filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,c(1,2,which(colnames(filtered_data) %in% outl_loci))]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after neutral filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# close filter loci ----

####   	 	 FORMAT GENEPOP:    	 ####

# Single row data:
# 0 = Reference allele homozygote (0101)
# 1 = SNP allele homozygote (0202)
# 2 = heterozygote (0102)

data<-filtered_data

# List parameters:
headline<-"enter_data_name_here"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-T
param.callrate<-T
param.MAF<-T
param.LD<-T
param.HWE<-T
param.neu<-T
param.dup<-T

ghead(data); dim(data)

# This makes three files: the genepop file, the parameter file and the locus file:

# Takes < 1.5 hr for full DPlan18
# < 1 min for 2500 loci
# 5 min for 53 x 18321
format_genepop(data,headline)

# close format genepop ----

####  FORMAT PLINK (for STRUCTURE): ####

# Use this for PLINK analyses and for STRUCTURE

###-->> Set data:
data_to_plink<-filtered_data

# ~ 5-10 MINS
formatted_ped<-format_plink_ped(snp_data=data_to_plink,locus_data=linf,remove_og=NULL,remove_cultivar=NULL)
ghead(formatted_ped)

check_plink_ped(orig_data=data_to_plink,plink_data=formatted_ped)

# write.table(formatted_ped,"ntham_filt1.ped",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ****** .map file ****** ~~~~ ##
formatted_map<-format_plink_map(ped_file=formatted_ped,locus_data=linf)
head(formatted_map)

# write.table(formatted_map,"ntham_filt1.map",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ***** locus info file ***** ~~~~ ##
plink_locus_info<-data.frame(lind=1:length(colnames(data_to_plink)[3:ncol(data_to_plink)]),locus=colnames(data_to_plink)[3:ncol(data_to_plink)])
head(plink_locus_info)

# write.table(plink_locus_info,"ntham_filt1_loci.txt",sep="\t",row.names=F,quote=F)

# Save parameters to file:
# List parameters:
data<-filtered_data
headline<-"ntham_filt1"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-T
param.callrate<-T
param.MAF<-T
param.LD<-T
param.HWE<-T
param.neu<-T
param.dup<-T

# The original plink parameter file was write_parameters() in the format_plink.R library but the genepop one is working better

gp_param(data,headline)

# close format plink ----

####  	 		FORMAT LFMM:	   	 ####

# Single row data:
# 0 = Reference allele homozygote (0101)
# 1 = SNP allele homozygote (0202)
# 2 = heterozygote (0102)

# For lfmm 0,1,2 represents the number of alleles. 
# From http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/note.pdf: The number of alleles can be the number of reference alleles or the number of derived alleles as long as a same choice is made for an entire column

ghead(filtered_data)
dim(filtered_data)

# List parameters:
data<-filtered_data
headline<-"lfmm_filt2"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-T
param.callrate<-T
param.MAF<-T
param.LD<-T
param.HWE<-T
param.neu<-F
param.dup<-T

ghead(data)
dim(data)

# Use genepop and structure converter:
data_lfmm<-data[,3:length(data)]
ghead(data_lfmm)

# Re-code genotypes:
# This only works if the is.na goes in the first place:
data_lfmm<-apply(
data_lfmm,2,function(x)
ifelse(is.na(x),"9",
ifelse(x=="2","1",
ifelse(x=="1","0",
ifelse(x=="0","2",x
))))
)
ghead(data_lfmm)

# check:
loci.now<-sample(colnames(data),3)
inds.now<-sample(1:nrow(data),3)
data[inds.now,colnames(data) %in% loci.now]
data_lfmm[inds.now,colnames(data_lfmm) %in% loci.now]

# Coded so the number is the number of reference alleles:
# 0 = 2
# 1 = 0
# 2 = 1

ghead(data_lfmm)

# Write data:
write.table(data_lfmm,file="lfmm.txt",row.names=F,col.names=F,quote=F,sep=" ")

# From format_genepop library:
gp_param(data,headline)

# Output locus info index:
lfmm_loci_filt2<-data.frame(lind=1:length(colnames(data)[3:ncol(data)]),locus=colnames(data)[3:ncol(data)])
head(lfmm_loci_filt2)

# write.table(lfmm_loci_filt2,"lfmm_loci_filt2.txt",row.names=F,quote=F,sep="\t")

ghead(data)
head(data[,1:2])
# Output site and individual data to merge with enviro data:

write.table(data[,1:2],"lfmm_site.txt",row.names=F,quote=F,sep="\t")

# close format lfmm ----

















