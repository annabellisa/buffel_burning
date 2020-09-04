
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

# Hi Binyin, I've restored the original data. Maybe we should save "binyin_winter.RData" for the full data set (69,799 loci) and when we do any filtering - we save those to different RData files? I think you've already started doing this in your "Datasets" folder which is great. 

# save.image("binyin_winter.RData")

#########################################
####	     FULL DATA SET:    		 ####
#########################################

###-->> Set data:
data_name<-"snp_onerow"
filtered_data<-get(data_name)

# --- *** Discard duplicated *** --- #

# 22923 out of 69799 loci (33%) were identified as duplicate sequences by BLAST. This seemed very high, so I did some random manual checks of the blast results and the results are correct: we do have many SNPs occurring on the same sequence. 
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
#Orginally
missing_cutoff<-0.5
missing_sum<-missing_data(filtered_data,3,missing_cutoff)
m_summary<-missing_sum$miss_sum
head(m_summary); dim(m_summary)
range(m_summary$missing)
filtered_data<-missing_sum$filt_dat
ghead(filtered_data); dim(filtered_data)

# save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/0.5 Datasets_filter_and_format/Max Cutoff.RData")

#hist(m_summary$missing_data)
#hist(1-m_summary$missing_data)
#hist(linf$CallRate)

#Remarks:
#when missing_cutoff<-0.5
#[1] "Original data: 93 individuals; 46875 loci"
#[1] "Filtered data: 93 individuals; 40711 loci"
#[1] "6164 loci with more than 50 % missing data removed"
#missing_cutoff<-0.2
#[1] "Original data: 93 individuals; 46875 loci"
#[1] "Filtered data: 93 individuals; 24579 loci"
#[1] "22296 loci with more than 20 % missing data removed"

# Jumpt to the test, following scripts are the filters


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

# I usually hashtag the write.table() and save.image() functions so we don't acidentally overwrite saved files

# write.table(filtered_data, file = "Partially_filtered_data_after_MAF", quote = F, sep = "\t", row.names = T)

# save.image("Partially Filtered Data After MAF.RData")

# --- ***Linkage disequilibrium (LD) filters *** --- #

ldf<-0.75

#Supplement_01_LD_test.R#
##BD's Script#
LD_dir<-"D:/OneDrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
LD_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LD_results"
dir(LD_dir)
ld_loc<-read.table(paste(LD_dir, "LD_r75_loci_to_remove.txt",sep="/"),header=T)

# AS LD: 
LD_dir_AS<-"RESULTS/LD_results"
dir(LD_dir_AS)
ld_loc<-read.table(paste(LD_dir_AS, "LD_r75_loci_to_remove.txt",sep="/"),header=T)
head(ld_loc); dim(ld_loc)

ldfilt<-as.character(ld_loc$locus)
head(ldfilt); length(ldfilt)

print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% ldfilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# --- *** HWE filters *** --- #

#Supplement_02_HWE_test.R#

#Annabel
hwe_dir<-"RESULTS/HWE_results"
dir(hwe_dir)
hwe_res<-read.table(paste(hwe_dir,"HWE_test.txt",sep="/"),header=T)
head(hwe_res)

#Binyin
hwe_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/HWE_results"
dir(hwe_dir)
hwe_res<-read.table(paste(hwe_dir,"HWE_test.txt",sep="/"),header=T)
head(hwe_res)

# Filter loci with HWD:
hwe_flag<-T
hwe_cutoff<-0.1 # we need to decide on the cutoff
hwefilt<-as.character(hwe_res$locus[hwe_res$p<hwe_cutoff])

# with a p cutoff of 0.1, there are 27655 (p) and 20644 (p.adj) loci identified as out of HWE. This is pretty much all of our loci, after applying the other filters. It leaves us with only 2000 or so loci for analysis. I would assume that the HWE results were influenced by the structure in the data and are thus not reliable. 

# we have two options: (1) skip this analysis as it's possibly inappropriate for our data. We're mainly doing this to look for genotyping issues, rather than biological issues. And, with so many other quality filters applied to the data, it's possibly not necessary; possibly a hang-over from my microsat days where genotyping errors and quality were harder to detect. (2) re-do the HWE analysis on a subset of the data, which form a single genetic cluster, e.g. the 19 individuals in buf01 and buf02. 

# I just checked 4 x recent papers that used DartSeq SNPs and none of them tested for HWE. This could be the end of my HWE testing days. 

head(hwefilt)
length(hwefilt)
dim(filtered_data)

print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% hwefilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# save.image("HWE_filter_data.Rdata")

# --- *** NEUTRALITY filter *** --- #

# 1. BAYESCAN - command line program - AS
# 2. LFMM - R - BD
# 3. PCAdapt - R - BD

# for all neutrality tests, remove: duplicated and monomorphic loci, but leave EVERYTHING ELSE (including low call rate - I changed my mind on this. 
# If I'm not wrong, PCAdapt uses bed, bim, fam files (i.e. PLINK files) and LFMM a special format below. The next step is to format the files for these programs. 


# Directory with results:
sel_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/00_Data"
sel_dir<-"D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/00_Data"


sel_dir<-"../../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION"
dir(sel_dir)


# filtered_data<-read.table(paste(sel_dir,"Partially_filtered_data",sep="/"),header=T)

ghead(filtered_data); dim(filtered_data)


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
headline<-"Genpop_Diversity_Original"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.dup<-T
param.mono<-T
param.callrate<-T
param.repavg<-T
param.MAF<-T
param.LD<-T
param.HWE<-F
param.neu<-T



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

# ~ 3 MINS on bigmac; 40711 loci Cenchrus
formatted_ped<-format_plink_ped(snp_data=data_to_plink,locus_data=linf,remove_og=NULL,remove_cultivar=NULL)
ghead(formatted_ped); dim(formatted_ped)

check_plink_ped(orig_data=data_to_plink,plink_data=formatted_ped)

# write.table(formatted_ped,"Cenchrus_filt1.ped",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ****** .map file ****** ~~~~ ##
formatted_map<-format_plink_map(ped_file=formatted_ped,locus_data=linf)
head(formatted_map)

# write.table(formatted_map,"Cenchrus_filt1.map",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ***** locus info file ***** ~~~~ ##
plink_locus_info<-data.frame(lind=1:length(colnames(data_to_plink)[3:ncol(data_to_plink)]),locus=colnames(data_to_plink)[3:ncol(data_to_plink)])
head(plink_locus_info)

# write.table(plink_locus_info,"Cenchrus_filt1_loci.txt",sep="\t",row.names=F,quote=F)

# Save parameters to file:
# List parameters:
data<-filtered_data
headline<-"Cenchrus_filt1"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-F
param.callrate<-T
param.MAF<-F
param.LD<-F
param.HWE<-F
param.neu<-F
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
param.repavg<-F
param.callrate<-T
param.MAF<-F
param.LD<-F
param.HWE<-F
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

####   	 	 CHANGE SITE-names (BayeScan):	   	 ####

ghead(filtered_data); dim(filtered_data)
old.site.name<-filtered_data$site
new.site.name<-filtered_data$site
# this is a hack that uses the first three characters and replaces them directly; not generalisable to other data sets and use with caution
new.site.name<-substr(new.site.name,start=1,stop=3)
filtered_data$site<-as.factor(new.site.name)

# close change site names ----

####   	 	 FORMAT STRUCTURE:	   	 ####
  
  # Single row data:
  # 0 = Reference allele homozygote (0101)
  # 1 = SNP allele homozygote (0202)
  # 2 = heterozygote (0102)
  
ghead(filtered_data); dim(filtered_data)

  # List parameters:
  data<-filtered_data
  headline<-"Cenchrus_filt2"
  param.sites<-levels(data$site)
  param.nosites<-length(param.sites)
  param.noloci<-ncol(data)-2
  param.noindiv<-nrow(data)
  param.mono<-T
  param.repavg<-F
  param.callrate<-T
  param.MAF<-F
  param.LD<-F
  param.HWE<-F
  param.neu<-F
  param.dup<-T
  
  ghead(data)
  dim(data)
  
  # This is MUCH faster than format_genepop (only a few mins)
  format_structure(data,headline)
  
  # Output locus info index:
  bs_loci_filt2<-data.frame(lind=1:length(colnames(data)[3:ncol(data)]),locus=colnames(data)[3:ncol(data)])
  head(bs_loci_filt2)
  
  # write.table(bs_loci_filt2,"bs_loci_filt2.txt",row.names=F,quote=F,sep="\t")
  
  # close format strucutre ----















