
## -- ** INDIVIDUAL HETEROZYGOSITY:

# Calculate individual heterozygosity
# Use dartseq format files, stored in the Genepop folder:
dir(gp_dir)

# NF
ds_filt2<-read.table(paste(gp_dir,"dartseq_format_NF.txt",sep="/"),header=T)

# NNF
ds_filt2<-read.table(paste(gp_dir,"dartseq_format_N_NF.txt",sep="/"),header=T)



# ghead(ds_filt2) # Error
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

if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] else stop("no zero category")

if(length(which(t.cons$ind.cons==2))==0) ih_out[i,3]<-"no_heterozygotes" else ih_out[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)

} # close for i
save.image("../04_workspaces/STEP04_divdist_wksp")
head(ih_out)

# write.table(ih_out,"ih.txt",quote=F,row.names=F,sep="\t")

# close genetic diversity



