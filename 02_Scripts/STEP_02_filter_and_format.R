
# ------------------------------------ #
# ------------- STEP 02  ------------- #
# ------------------------------------ #

### filter SNP loci & format for different software
### Author: Annabel Smith & Di Binyin

#Set Working dir
setwd("D:/GitHub/Binyin_Winter")
setwd("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter")

# Binyin: Load and tidy workspace and remove everything except necessary objects:
load("D:/GitHub/Binyin_Winter/binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))
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
cr<-0.5
missing_sum<-missing_data(filtered_data,3,cr)
m_summary<-missing_sum$miss_sum
# hist(m_summary$missing_data)
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
# hist(maf_sum$maf[maf_sum$maf<0.1])

# Filter loci with extreme maf:
###-->> Set maf limit:
malim<-0.05
filtered_data<-maf_filter(maf_sum,filtered_data,malim)
ghead(filtered_data); dim(filtered_data)

# save.image("binyin_winter.RData")

# --- ***Linkage disequilibrium (LD) filters *** --- #

# Filter loci with were in LD in > 5 populations with a correlation of 0.75 (see Supplement_01_LD_tests.R for details):


LD_dir<-"../../ANALYSIS_RESULTS/LINKAGE_DISEQUILIBRIUM/LD_parameters"
dir(LD_dir)

ld_loc<-read.table(paste(LD_dir, "LD_r75_over5pop_LOCI_FOR_REMOVAL.txt",sep="/"),header=T)
head(ld_loc)

# Filter loci with LD in > 5 pops:
###-->> Set ld limit:
ldf<-0.75
ldfilt<-as.character(ld_loc$locus)
print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% ldfilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data)

# --- *** HWE filters *** --- #

hwe_dir<-"../../ANALYSIS_RESULTS/ALL_pops_HWE_test"
dir(hwe_dir)

hwe_res<-read.table(paste(hwe_dir,"hwe_all_pops.txt",sep="/"),header=T)
head(hwe_res)

# we need to make some changes here. The previous project I did, we looked for consistent loci that were consistently out of HWE in > 5 populations. For ours however, I think we can just go with the p value. 

# Summarise hwe tests:
loci_outin5<-data.frame(locus=names(table(hwe_res$loci_out)[which(table(hwe_res$loci_out)>5)]),no_pops_out=as.numeric(table(hwe_res$loci_out)[which(table(hwe_res$loci_out)>5)]))
head(loci_outin5)
dim(loci_outin5)

# Filter loci with HWED in > 5 pops:
hwe_flag<-T
hwefilt<-as.character(loci_outin5$locus)
print(paste("no loci before ld filt = ",dim(filtered_data)[2],sep=""))
filtered_data<-filtered_data[,-which(colnames(filtered_data) %in% hwefilt)]
filtered_data<-tidy.df(filtered_data)
print(paste("no loci after ld filt = ",dim(filtered_data)[2],sep=""))
ghead(filtered_data); dim(filtered_data)

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

















