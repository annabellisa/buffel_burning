
## -- ** INDIVIDUAL HETEROZYGOSITY:

# load libraries
library(lme4)
library(tidyverse)



# Calculate individual heterozygosity
# Use dartseq format files, stored in the Genepop folder:

gp_dir<-"D:\\Onedrive\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\RESULTS\\Diversity_and_Distance\\NF_Format"
gp_dir<-"D:\\Onedrive\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\RESULTS\\Diversity_and_Distance\\N_NF_Format"



dir(gp_dir)


# NF
ds_filt2<-read.csv(paste(gp_dir,"dartseq_format_NF.txt",sep="/"),header=T)

# NNF
ds_filt2<-read.csv(paste(gp_dir,"dartseq_format_N_NF.txt",sep="/"),header=T)



ghead(ds_filt2) # check the file format, csv works 
dim(ds_filt2)

# Codes for onerow formatted data:
# 0 = Reference allele homozygote
# 1 = SNP allele homozygote
# 2 = heterozygote

# 3 MINS on laptop

ih_out<-data.frame(ds_filt2[,1:2],ind_het=NA)
head(ih_out,25)

for(i in 1:nrow(ih_out)){
  ind.cons<-as.character(ds_filt2[i,3:length(ds_filt2)])
  t.cons<-data.frame(ind.cons=as.numeric(names(table(ind.cons))),count=as.numeric(table(ind.cons)),stringsAsFactors=F)
  
  if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] 
  else stop("no zero category")
  
  if(length(which(t.cons$ind.cons==2))==0) ih_out[i,3]<-"no_heterozygotes" 
  else ih_out[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)


} # close for i, mins

# need further investigation 
class(ih_NF$ind_het)


ih_NNF$ind_het<-as.numeric(ih_NNF$ind_het, na.rm = TRUE)
ggplot(data = ih_NNF, aes(x = site, y = ind_het)) +
  geom_boxplot() +
theme_classic()





save.image("C:\\Users\\s4467005\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\03_Workspaces\\STEP04_divdist_wksp.RData")
head(ih_out)

# write.table(ih_out,"ih.txt",quote=F,row.names=F,sep="\t")

# close genetic diversity


# Models

load("C:\\Users\\s4467005\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\03_Workspaces\\STEP04_divdist_wksp.RData")

# ih_out dataset integration (same way as AR) full join (trt,longitude)

mod1<-lmer(ih_het~treatment + (1|site) , data = ih_out)

# CV coef of variance on site level (group_by site)



