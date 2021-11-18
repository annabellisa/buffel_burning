
# ------------------------------------ #
# ------------- STEP 02  ------------- #
# ------------------------------------ #

### filter SNP loci & format for different software
### Updated after first peer-review

### Author: Annabel Smith & Di Binyin

# Load and tidy workspace and remove everything except necessary objects:
load("buffel_burning.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# "buffel_burning.RData" is for the full data set (69,799 loci); filtered data sets should be saved to different RData files.

# save.image("buffel_burning.RData")

# FULL DATA SET:  

###-->> Set data:
data_name<-"snp_onerow"
filtered_data<-get(data_name)
ghead(filtered_data); dim(filtered_data)

###-->> FILTER # 1: duplicated, monomorphic, missing data

# --- *** Discard duplicated *** --- #

# 22921 out of 69799 loci (33%) were identified as duplicate sequences by BLAST. This seemed very high, so I did some random manual checks of the blast results and the results are correct: we do have many SNPs occurring on the same sequence. We are removing these to reduce the chance of physical linkage. For a given sequence, we are keeping the SNP with the highest call rate (i.e. those assigned as 'duplicates' on the locus info file have lower call rates, see STEP_01 script for details).  
dup_loc<-as.character(linf$locus[linf$duplicate==1])
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% dup_loc)]
filtered_data<-tidy.df(filtered_data)
ghead(filtered_data); dim(filtered_data)

# FILTER LOCI:    	# ----

# --- *** Filter monomorphic loci *** --- #
filtered_data<-mono_loci(filtered_data,3)
ghead(filtered_data); dim(filtered_data)

# --- *** DartSeq Quality Control (QC) filters *** --- #

# Filter loci with high missing data rate:

###-->> Set maximum missing data:
## missing_data==1-CallRate
missing_cutoff<-0.5
missing_sum<-missing_data(filtered_data,3,missing_cutoff)
m_summary<-missing_sum$miss_sum
head(m_summary); dim(m_summary)
range(m_summary$missing)
filtered_data<-missing_sum$filt_dat
ghead(filtered_data); dim(filtered_data)

###-->> FILTER # 2: reproducibility, MAF, LD, paralogs

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

# Filter loci with extreme maf:
###-->> Set maf limit:
malim<-0.05
filtered_data<-maf_filter(maf_sum,filtered_data,malim)
ghead(filtered_data); dim(filtered_data)

# --- ***Linkage disequilibrium (LD) filters *** --- #

# Set cut-off
ldf<-0.75

# Run Supplement_01_LD_tests.R #

# AS LD: 
LD_dir_AS<-"04_RESULTS/Filtering/LD_results"
ld_loc<-read.table(paste(LD_dir_AS, "LD_r75_loci_to_remove.txt",sep="/"),header=T)
head(ld_loc); dim(ld_loc)

ldfilt<-as.character(ld_loc$locus)
head(ldfilt); length(ldfilt)

print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% ldfilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# --- *** PARALOG filter *** --- #

# Filter putative paralogs by K clusters from structure because these might correspond to different cytotypes. 

# Make genind object to get Hobs and FIS for each K cluster separately. 
library(adegenet); library(hierfstat)

# Make genind object (takes time so have included it in a separate workspace; filtered_data in this workspace has all of the above filters already applied)
load("03_Workspaces/paralog_filter.RData")
genind_all

# gp_dir<-"00_Data/Genepop_Files"
# dir(gp_dir)
# genind_all<-read.genepop(file=paste(gp_dir,"Genepop_all_loci_byK3.gen",sep="/"), ncode=2L,quiet=FALSE)
# save.image("03_Workspaces/paralog_filter.RData")

# The genind object names populations using the last individual in the list for each K; I updated this before formatting the genepop file so they would relate to the K3 clusters here. 
genind_all@pop

# Get basic stats from hierfstat (3 mins, saved in paralog_filter.RData):
# gendiv_all <- basic.stats(genind_all, diploid = TRUE, digits = 2)
head(gendiv_all$Fis)
head(gendiv_all$Ho)

fis_all<-as.data.frame(gendiv_all$Fis)
Hobs_all<-as.data.frame(gendiv_all$Ho)
head(fis_all)
head(Hobs_all)
table(rownames(fis_all)==rownames(Hobs_all))

# Get list of loci that have Hobs > 80% in any cluster:
Hobs80<-unique(c(rownames(Hobs_all[which(Hobs_all$K1>0.8),]),rownames(Hobs_all[which(Hobs_all$K2>0.8),]),rownames(Hobs_all[which(Hobs_all$K3>0.8),])))
head(Hobs80); length(Hobs80)
# 14718 loci fit this category but most (11662) have already been filtered out with the above quality control

# save.image("03_Workspaces/paralog_filter.RData")

# Filter loci with extreme observed heterozygosity (80%, following Reynes et al. 2021 MER):

length(which(colnames(filtered_data) %in% Hobs80))

filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% Hobs80)]
filtered_data<-tidy.df(filtered_data)
ghead(filtered_data); dim(filtered_data)

###-->> FILTER # 3: neutrality

# --- *** NEUTRALITY filter *** --- #

# Neutrality tests were run in two different programs (see STEP_04_loci_under_selection.R): LFMM and PCAdapt

# Neutrality tests were run using the Filter # 2 output

# PCAdapt uses bed, bim, fam files (i.e. PLINK files) and LFMM a special format. See the format scripts below for these programs.

# Directory with results:
sel_dir<-"00_Data"
dir(sel_dir)

ghead(filtered_data); dim(filtered_data)

# Outlier loci (from PCAdapt and LFMM):
res<-read.table(paste(sel_dir,"outliers_all.txt",sep="/"),header=T)
head(res,3)

# Names and length of outliers:
outl_loci<-as.character(res$locus[c(which(res$pca_outl==1),which(res$lfmm_outl==1))])
outl_loci<-outl_loci[-which(duplicated(outl_loci))]
head(outl_loci)
length(outl_loci)

# Filter outlier loci (for neutral analyses):
neutral_flag<-T
print(paste("no loci before neutral filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% outl_loci)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after neutral filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# Create data set with ONLY non-neutral loci (for adaptive diversity analyses):
neutral_flag<-F
print(paste("no loci before neutral filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,c(1,2,which(colnames(filtered_data) %in% outl_loci))]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after neutral filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

# close filter loci ----

# FORMAT DARTSEQ:    	# ----

# For analyses that require DartSeq format (e.g. our genetic diversity analysis), the data can be written directly, without any further processing:
write.table(filtered_data, "dartseqNN_filt4.txt", quote=F, row.names=F, sep="\t")

# close format DartSeq ----

# FORMAT genepop by K3 (all sites):    	# ----

# Format by K3 with all individuals in the same file:

dir("00_Data")
kclust<-read.table("00_Data/K_genetic_clusters_Cenchrus_filt4.txt", header=T)
kclust<-kclust[,c("indiv", "K3")]
head(kclust)

# These should all be true:
table(kclust$indiv %in% filtered_data$ind)
ghead(filtered_data); dim(filtered_data)
head(kclust); dim(kclust)
filtered_data<-merge(filtered_data, kclust, by.x="ind", by.y="indiv", all.x=T, all.y=F)
filtered_data$site<-NULL
filtered_data<-filtered_data[,c(which(colnames(filtered_data) %in% c("K3", "ind")),grep("L",colnames(filtered_data)))]
filtered_data<-filtered_data[,c(2,1,grep("L", colnames(filtered_data)))]
colnames(filtered_data)[which(colnames(filtered_data)=="K3")]<-"site"
filtered_data$site<-as.factor(filtered_data$site)

# Only run the this section if you need the population names to be meaningful (i.e. for the paralog filter):
ghead(filtered_data); dim(filtered_data)
filtered_data$ind<-as.character(filtered_data$ind)
filtered_data$ind[filtered_data$site==1]<-"K1"
filtered_data$ind[filtered_data$site==2]<-"K2"
filtered_data$ind[filtered_data$site==3]<-"K3"
filtered_data$ind<-as.factor(filtered_data$ind)

# close format genepop K3 ----

# Format Genepop with K in separate files, individuals grouped into sites ---- 

dir("00_Data")
kclust<-read.table("00_Data/K_genetic_clusters_Cenchrus_filt4.txt", header=T)
kclust<-kclust[,c("indiv", "K3")]
head(kclust)

# Format separate files for each K:

# These should all be true:
table(kclust$indiv %in% filtered_data$ind)
ghead(filtered_data); dim(filtered_data)
head(kclust); dim(kclust)

filtered_data<-merge(filtered_data, kclust, by.x="ind", by.y="indiv", all.x=T, all.y=F)
ghead(filtered_data); dim(filtered_data)

fd1<-filtered_data[which(filtered_data$K3==1),]
fd2<-filtered_data[which(filtered_data$K3==2),]
fd3<-filtered_data[which(filtered_data$K3==3),]

fd1<-tidy.df(fd1)
fd2<-tidy.df(fd2)
fd3<-tidy.df(fd3)

# Remove the K column:
fd1$K3<-NULL
fd2$K3<-NULL
fd3$K3<-NULL

# Remove sites with < 4 individuals
# There are too few individuals in K2 to do the IBD analysis, so do this on K1 and K3 only
fd1<-fd1[-which(fd1$site %in% names(which(table(fd1$site)<4))),]
fd2<-fd2[-which(fd2$site %in% names(which(table(fd2$site)<4))),]
fd3<-fd3[-which(fd3$site %in% names(which(table(fd3$site)<4))),]
fd1<-tidy.df(fd1)
fd2<-tidy.df(fd2)
fd3<-tidy.df(fd3)

# Check it:
table(fd1$site)
table(fd2$site)
table(fd3$site)

# Flip ind and site cols:
fd1<-fd1[,c(2,1,3:ncol(fd1))]
fd3<-fd3[,c(2,1,3:ncol(fd3))]

ghead(fd1); dim(fd1)
ghead(fd3); dim(fd3)

# --- *** Filter monomorphic loci *** --- #
fd1<-mono_loci(fd1,3)
ghead(fd1); dim(fd1)

fd3<-mono_loci(fd3,3)
ghead(fd3); dim(fd3)

# Remove all loci which are all NAs (only applies to fd3):

fd3<-fd3[,-which(colnames(fd3) %in% names(unlist(apply(fd3[,3:ncol(fd3)], 2, function(x) which(length(which(is.na(x)))==nrow(fd3))))))]

# Check it:
table(fd1$site)
table(fd3$site)

ghead(fd1); dim(fd1)
ghead(fd3); dim(fd3)

# close format genepop separate K ----

####   	 	 FORMAT GENEPOP:    	 ####

# Single row data:
# 0 = Reference allele homozygote (0101)
# 1 = SNP allele homozygote (0202)
# 2 = heterozygote (0102)

data<-filtered_data

# parameter flags for param file:
headline<-"Genepop_filt4_non_neutral"
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
param.plog<-T
param.HWE<-F
param.neu<-F

ghead(data); dim(data)

# format_genepop makes three files: the genepop file, the parameter file and the locus file:

format_genepop(data,headline)

# close format genepop ----

####  FORMAT PLINK (for STRUCTURE): ####

# Use for PLINK analyses and for STRUCTURE

###-->> Set data:
data_to_plink<-filtered_data

# Cenchrus on big mac: ~ 3 MINS 40711 loci; 1 min 20159 loci 
formatted_ped<-format_plink_ped(snp_data=data_to_plink,locus_data=linf,remove_og=NULL,remove_cultivar=NULL)
ghead(formatted_ped); dim(formatted_ped)

check_plink_ped(orig_data=data_to_plink,plink_data=formatted_ped)

# write.table(formatted_ped,"Cenchrus_filt4.ped",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ****** .map file ****** ~~~~ ##
formatted_map<-format_plink_map(ped_file=formatted_ped,locus_data=linf)
head(formatted_map)

# write.table(formatted_map,"Cenchrus_filt4.map",quote=F,row.names=F,col.names=F,sep=" ")

## ~~~~ ***** locus info file ***** ~~~~ ##
plink_locus_info<-data.frame(lind=1:length(colnames(data_to_plink)[3:ncol(data_to_plink)]),locus=colnames(data_to_plink)[3:ncol(data_to_plink)])
head(plink_locus_info); dim(plink_locus_info)

# write.table(plink_locus_info,"Cenchrus_filt4_loci.txt",sep="\t",row.names=F,quote=F)

# parameter flags for param file:
data<-filtered_data
headline<-"Cenchrus_filt4"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-T
param.callrate<-T
param.MAF<-T
param.LD<-T
param.plog<-T
param.HWE<-F
param.neu<-T
param.dup<-T

# my original plink parameter file was write_parameters() in the format_plink.R library but the genepop one is working better

gp_param(data,headline)

# close format plink ----

####  	 		FORMAT LFMM:	   	 ####

# Single row data:
# 0 = Reference allele homozygote (0101)
# 1 = SNP allele homozygote (0202)
# 2 = heterozygote (0102)

# For lfmm 0,1,2 represents the number of alleles. 
# From http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/note.pdf: The number of alleles can be the number of reference alleles or the number of derived alleles as long as a same choice is made for an entire column

ghead(filtered_data); dim(filtered_data)

# parameter flags for param file:
data<-filtered_data
headline<-"lfmm_filt3"
param.sites<-levels(data$site)
param.nosites<-length(param.sites)
param.noloci<-ncol(data)-2
param.noindiv<-nrow(data)
param.mono<-T
param.repavg<-T
param.callrate<-T
param.MAF<-T
param.LD<-T
param.plog<-T
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
lfmm_loci_filt3<-data.frame(lind=1:length(colnames(data)[3:ncol(data)]),locus=colnames(data)[3:ncol(data)])
head(lfmm_loci_filt3); dim(lfmm_loci_filt3)

# write.table(lfmm_loci_filt3,"lfmm_loci_filt3.txt",row.names=F,quote=F,sep="\t")

ghead(data)
head(data[,1:2])
# Output site and individual data to merge with enviro data:

write.table(data[,1:2],"lfmm_site.txt",row.names=F,quote=F,sep="\t")

# close format lfmm ----

####   	 	 FORMAT STRUCTURE:	   	 ####
  
  # Single row data:
  # 0 = Reference allele homozygote (0101)
  # 1 = SNP allele homozygote (0202)
  # 2 = heterozygote (0102)
  
ghead(filtered_data); dim(filtered_data)

# parameter flags for param file:
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















