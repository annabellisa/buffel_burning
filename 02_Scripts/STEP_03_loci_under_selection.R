
# ------------------------------------ #
# ------------- STEP 03  ------------- #
# ------------------------------------ #

### Post-process BayeScan (run in command line) and run PCAdapt and LFMM to detect loci under selection

### Author: Annabel Smith, Binyin Di

# AS Load and tidy workspace and remove everything except necessary objects:
load("binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# BD load R.Data and Functions
load("D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/binyin_winter.RData")

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

load("03_Workspaces/STEP03_loci_under_selection.RData")

#  POST-PROCESS BAYESCAN results:    	# ----

# *** ANALYSE OUTLIER LOCI:

# Rename data:
gt_data<-snp_onerow

# Directory with bayescan results:

bs_dir<-"D:\\Onedrive\\OneDrive - The University of Queensland\\Offline Winter Project\\Cenchrus_filt2_BS_po400_RESULTS" 

bs_dir<-"/Users/annabelsmith/OneDrive - The University of Queensland/02_TEACHING/STUDENTS/Binyin_Winter/Analysis/Offline_Results/BayeScan/BayeScan_ANALYSIS/Cenchrus_filt2"

dir(bs_dir)

# Define bayescan outliers from FST Outlier file:
bs_fst<-paste(bs_dir,"Cenchrus_filt2_BS_po400_RESULTS/Cenchrus_filt2_BS_po400_fst.txt",sep="/") 
bsres<-plot_bayescan(res=bs_fst,FDR=0.05) 
bs_outl<-bsres$outliers
bs_n_outl<-bsres$nb_outliers

# Make MARKER file to separate non-neutral loci:

# Get locus index (this is the locus index file that is made in format_structure() function for bayescan):

bslinf<-read.table("D:\\Onedrive\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\RESULTS\\STRUCTURE\\STRUCTURE_DIR\\Cenchrus_filt2\\Cenchrus_filt2_loci.txt", header = TRUE)

bslinf<-read.table(paste(bs_dir,"bs_loci_filt2.txt",sep="/"), header = TRUE) 
head(bslinf); dim(bslinf)

bsoutl<-data.frame(lind=bs_outl,outl=1)
head(bsoutl); dim(bsoutl)

# Outlier loci:
bslinf<-merge(bslinf,bsoutl,by="lind", all.x = T, all.y = F)
bslinf$outl[which(is.na(bslinf$outl))]<-0
head(bslinf); dim(bslinf)
check.rows(bslinf)

colnames(bslinf)[which(colnames(bslinf)=="outl")]<-"bs_outl"
table(bslinf$bs_outl)

# write.table(bslinf,file="bayescan_outliers.txt",quote=F,row.names=F,sep="\t")

# *** PLOT DIAGNOSTICS:

bs_sel<-paste(bs_dir,"Cenchrus_filt2_BS_po400_RESULTS/Cenchrus_filt2_BS_po400.sel",sep="/")
seldat<-read.table(bs_sel,colClasses="numeric")
# save.image("03_Workspaces/STEP03_loci_under_selection.RData")

# Plot log likelihood:
dev.new(width=5,height=5, noRStudioGD = T,dpi=100) 
par(mfrow=c(1,1),mar=c(4,4,2,1),mgp=c(2,0.5,0))
parameter<-"logL" # a few minutes 
plot(density(seldat[[parameter]]),xlab="log Likelihood",main="BayeScan posterior distribution", font.main=1)

# Plot FST:
dev.new(width=5,height=8, noRStudioGD = T,dpi=80) 
par(mfrow=c(3,3),mar=c(2,2,0.2,0.2),mgp=c(2,0.5,0))
for (i in grep("Fst",colnames(seldat))){
par.thisrun<-colnames(seldat)[i]
plot(density(seldat[[par.thisrun]]),xlab=par.thisrun,main="",cex.axis=0.75)
legend("bottom",legend=par.thisrun,cex=1,bty="n")
}

dim(seldat)

# Plot alpha for a random selection of loci:
alphacols<-colnames(seldat)[grep("alpha",colnames(seldat))]
outl_alpha<-paste("alpha",bs_outl,sep="")
head(outl_alpha)
table(outl_alpha %in% alphacols)

dev.new(width=12,height=8, noRStudioGD = T,dpi=80) 
par(mfrow=c(5,10),mar=c(2,2,0.2,0.2),mgp=c(2,0.5,0))
for (i in 1:50){
col.thisrun<-sample(outl_alpha,1)
par.thisrun<-colnames(seldat)[which(colnames(seldat)==col.thisrun)]
plot(density(seldat[[par.thisrun]]),xlab=par.thisrun,main="",cex.axis=0.75)
legend("bottom",legend=par.thisrun,cex=1,bty="n")
}

# *** PLOT OUTLIER LOCI:

loc.toanalyse<-as.character(bslinf$locus[which(bslinf$bs_outl==1)]) 
length(loc.toanalyse)
head(loc.toanalyse)

# Make sure site in genotype and site data match:
site_data<-sdat
site_data$site<-sub("buf","X",site_data$site)
ghead(gt_data)
head(site_data,3)
table(gt_data$site %in% site_data$site)

# PLOT GENOTYPE FREQUENCIES BY LOCATION:

loc.toanalyse # for Cenchrus there's only 5 outliers from BS:
length(loc.toanalyse)

# Location for heatmap files:	
out.dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/BayeScan/Cenchrus_filt2_BS_po400"
dir(out.dir)

# this plots by longitude (see plot_freq_long folder for results); use "jpg" for windows or "pdf" for mac
out.dir<-"RESULTS/BayeScan/plot_freq_long"
plot_freq_long(loci=loc.toanalyse,genotype_data = gt_data,site_data = site_data,out.dir = out.dir,number_to_plot = 5, out_file_type="pdf") 

out.dir<-"RESULTS/BayeScan/plot_freq_burn"
plot_freq_long(loci=loc.toanalyse,genotype_data = gt_data,site_data = site_data,out.dir = out.dir,number_to_plot = 5, out_file_type="pdf") 

# close BayeScan ----

#  PCADAPT:    	# ----

# tutorial:
# browseVignettes("pcadapt")
# or: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html  

library(pcadapt)
library(qvalue)
library(grDevices)

# to install qvalue:
# https://www.bioconductor.org/packages/release/bioc/html/qvalue.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

path_to_file <- "RESULTS/PCAdapt/PCAdapt_files"

path_to_file<- "D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/PCAdapt/PCAdapt_files"

dir(path_to_file)

filename <- read.pcadapt(paste(path_to_file,"Cenchrus_filt1.bed",sep="/"), type = "bed")

### --- *** SET INITIAL K *** --- ###

K <- 19 # Begin with the number of sample sites
x <- pcadapt(input = filename, K = 19) 

# CHOOSE K FROM SCREE PLOT:

dev.new(width=4,height=4,dpi=160,pointsize=12, noRStudioGD = T)
par(mar=c(4,4,2,0.5))
plot(1:K,x$singular.values^2,xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n", main = "PCA scree plot", font.main=1)
title(xlab="Principal component",mgp=c(2.5,1,0))
grid()
lines(1:K,x$singular.values^2)
points(1:K,x$singular.values^2,pch=20)

# CHOOSE K FROM PCs:

# Get population names from the fam file:
famf<-read.table(paste(path_to_file,"Cenchrus_filt1.fam",sep="/"),header=F) # strange fam file
head(famf)
poplist.names <- as.character(famf$V1)
head(poplist.names)

# Reproduce inbuilt PC plot: 

dev.new(width=10,height=4,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(2,5),mar=c(4,4,1,1),mgp=c(2.5,1,0))

for (i in seq(1,ncol(x$scores),2)){
  pc1<-x$scores[,i]
  pc2<-x$scores[,i+1]
  plot(pc1,pc2,col=rainbow(length(unique(poplist.names))),xlab="",ylab="",pch=20)
  title(xlab=paste("PC",i,sep=""),cex.lab=1)
  title(ylab=paste("PC",i+1,sep=""),cex.lab=1)
}
par(xpd=NA)
legend(0.6, 0.55,legend=unique(poplist.names),col=rainbow(length(unique(poplist.names))),pch=20,ncol=3)
par(xpd=F)

### --- *** SET K *** --- ###

K<-2 # return 
x <- pcadapt(input = filename, K = K) 
summary(x)

### --- *** CHECK LOADINGS FOR LD *** --- ###

dev.new(width=10,height=4,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (i in 1:ncol(x$loadings)){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
  }
pc1<-x$scores[,1]
pc2<-x$scores[,2]
pc3<-x$scores[,3]

# the loadings for PC3 suggest LD is an issue which is not surprising because we have not filtered LD loci for this initial stage; however, when I tried re-running the PCA with thinning, the problem on the loadings was not resolved. 

# Get outliers based on q values:

loc_dat<-read.table(paste(path_to_file,"Cenchrus_filt1_loci.txt",sep="/"),header=T) 
loc_dat$pvalue<-x$pvalues

# there are plenty of NAs:
length(which(is.na(loc_dat$pvalue)))

# Get qvalues
loc_dat$qvalue<-qvalue(loc_dat$pvalue)$qvalue
length(which(is.na(loc_dat$qvalue)))

# Get outliers
alpha <- 0.05
loc_dat$outlier <- 0
loc_dat$outlier[which(loc_dat$qvalue < alpha)]<-1
table(loc_dat$outlier)

# Associate outliers with PCs: library(pcadapt)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)

snp_pc <- get.pc(x,outliers)
head(snp_pc); dim(snp_pc)
range(snp_pc$PC)
range(snp_pc$SNP)

# Add them to the main data frame:
loc_dat<-merge(loc_dat, snp_pc, by.x="lind", by.y="SNP", all.x=T, all.y=F)
loc_dat$PC[which(is.na(loc_dat$PC))]<-0
check.rows(loc_dat)

# most are outlying on PC2:
table(loc_dat$PC)

# write.table(loc_dat,"outliers_PCAdapt_from_snp_pc.txt",sep="\t",row.names=F,quote=F)

save.image("03_Workspaces/STEP03_loci_under_selection.RData")

# Reproduce manhattan plot:
head(loc_dat); dim(loc_dat)

dev.new(width=10,height=4,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,1),mar=c(3,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(loc_dat),-log(loc_dat$pvalue,10),type="n", xlab="", ylab="-log10(pvalue)")
points(1:nrow(loc_dat),-log(loc_dat$pvalue,10),pch=20, cex=0.5,col=ifelse(as.factor(loc_dat$qvalue<0.05)==T,"red","black"))

# close PCAdapt ----

#  LFMM:    	# ----

library(lfmm) # https://bcm-uga.github.io/lfmm/articles/lfmm
library(qvalue)

lfmm_dir<-"RESULTS/LFMM/LFMM_files"
lfmm_dir<-"D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LFMM/lfmm_files"

dir(lfmm_dir)

# Genotype data:
lfdat<-as.matrix(read.table(paste(lfmm_dir,"lfmm.txt",sep="/"),header=F)) # a few seconds 
colnames(lfdat)<-NULL
ghead(lfdat)

# Site data for each individual in the genotype file:
lfsite<-read.table(paste(lfmm_dir,"lfmm_site.txt",sep="/"),header=T)
head(lfsite)

# Enviro data for each site:

# Add site info
sdat$site_code <- sub("buf","X",sdat$site)
head(sdat) 

# just need site_code, burnt_unburnt, and long:
lfsd<-sdat[which(sdat$site_code %in% unique(lfsite$site)),c(which(colnames(sdat)=="site"),which(colnames(sdat)=="burn_unburnt"),which(colnames(sdat)=="long"):which(colnames(sdat)=="site_code"))] 
lfsd<-lfsd[order(lfsd$site),]
lfsd<-tidy.df(lfsd)
head(lfsd)

# Make sure they're in the same order, these should all be T:
lfsd$site_code==unique(lfsite$site)

# This is for the biplot labels:
rownames(lfsd)<-lfsd[,"site_code"] 

# Make PCs:

# burn_unburn to binary code: burnt is 1, otherwise, 0/
lfsd$burn_unburnt<-ifelse(lfsd$burn_unburnt =="b",1,0)
lf_pc<-princomp(lfsd[,c(2,3)],cor=T) 
summary(lf_pc)
# biplot(lf_pc,xlab="PC1",ylab="PC2",cex=0.7)
lf_pc$loadings
head(lf_pc$scores)

# Add PCs to individual level data:
lfpcs<-data.frame(site=rownames(lfsd),lf_pc$scores[,1])
lfpcs<-tidy.df(lfpcs)
head(lfpcs)

lfsd2<-lfsite
lfsd3<-merge(lfsd2, lfsd, by.y = "site_code", by.x = "site", all.x = TRUE, all.y = FALSE)

#check:
lfsite$ind == lfsd3$ind

# Reduce columns to burnt_unburnt and long - these are the environmental variables we're using in the model:
lfs<-lfsd3[,c(4,5)]
head(lfs)
head(lfs); dim(lfs)

# The environmental matrix (X):
lfs<-as.matrix(lfs[,1:2])
head(lfs)

# The genotype matrix (Y):
ghead(lfdat)

# choose K:
pc <- prcomp(lfdat)
head(pc$sdev)
length(pc$sdev)
summary(pc)$importance[,1:5]
str(pc)
head(summary(pc)$importance[2,])

# Plot prop var explained for comparison with pcadapt (instead of stdev^2 used in example):
dev.new("",4,4,dpi=160,pointsize=12)
par(mar=c(4,4,0.5,0.5))
plot(1:40,summary(pc)$importance[2,][1:40], main = "LFMM scree plot",xlab="Principal component ",ylab="Proportion variance explained",pch=20,las=1,type="n")
title(xlab="PC",mgp=c(2.5,1,0))
grid()
lines(1:40,summary(pc)$importance[2,][1:40])
points(1:40,summary(pc)$importance[2,][1:40],pch=20)

# The screes from the PCAdapt and LFMM are similar; K between 3 and 5

# Ridge estimates and GWAS tests
head(lfs) # The environmental matrix (X):
ghead(lfdat) # The genotype matrix (Y):

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = lfdat, X = lfs, K = 5)
summary(mod.lfmm)
head(mod.lfmm$B)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = lfdat, X = lfs,  lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue 
head(pvalues)
length(pvalues)

# Direct effect sizes estimated from latent factor models:
efs.lfmm1<-effect_size(Y = lfdat, X = as.matrix(lfs[,1]), lfmm = mod.lfmm) # a few minutes
efs.lfmm2<-effect_size(Y = lfdat, X = as.matrix(lfs[,2]), lfmm = mod.lfmm)

head(efs.lfmm1)
head(efs.lfmm2)

lfmm_loci<-read.table(paste(lfmm_dir,"lfmm_loci_filt2.txt",sep="/"),header=T)
head(lfmm_loci)
dim(lfmm_loci)

# Add p and q values to locus info:
lfres<-data.frame(locus=lfmm_loci$locus,pvalues)
colnames(lfres)<-paste(colnames(lfres),"_p",sep="")

lfqs<-apply(as.matrix(lfres[,2:ncol(lfres)]),2,function(x) qvalue(x)$qvalues)
colnames(lfqs)<-paste(substr(colnames(lfqs),1,2),"_q",sep="")
lfres<-cbind(lfres,lfqs)
head(lfres)

# Get outliers based on q values:
alpha <- 0.05
lf_outl<-apply(as.matrix(lfres[,2:ncol(lfres)]),2,function(x) which(x<alpha))

lfres$bu_outl<-ifelse(rownames(lfres) %in% lf_outl$bu_q,1,0)
lfres$lo_outl<-ifelse(rownames(lfres) %in% lf_outl$lo_q,1,0)
head(lfres)

# no loci were outliers along the burn gradient with a q threshold of 0.05
# 35 loci were outliers along the longitudinal gradient with a q threshold of 0.05

table(lfres$bu_outl)
table(lfres$lo_outl)

# Reproduce manhattan plot:
lfres$never_outl<-rowSums(lfres[,which(colnames(lfres)=="bu_outl"):ncol(lfres)])
lfres$never_outl<-rowSums(lfres[,8,drop = FALSE]) 
head(lfres)

# BD ggplot:
lfres_df_lo<-lfres %>% subset(lo_q<alpha)
lfres_df_lo$rowname<-lfres %>% subset(lo_q<alpha) %>% row.names()

lfres_df_bu<-lfres %>% subset(bu_q<alpha)
lfres_df_bu$rowname<-lfres %>% subset(bu_q<alpha) %>% row.names()


ggplot(data = lfres,mapping = aes(x= 1:nrow(lfres), y = -log(long_p,10))) +
  geom_point(aes(x= 1:nrow(lfres),-log(long_p,10))) +
  geom_point(data = lfres_df_lo, aes(x = as.numeric(rowname),-log(long_p,10)), colour = "red") +
 # geom_point(data = lfres_df_bu, aes(x = as.numeric(rowname),-log(burn_unburnt,10), colour = "green"))+
  labs(y ="-log10(p value)",
       x ="locus") +
  theme_article() 
  


# AS manhattan (30 Oct 2020):
quartz("",8,4,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),-log(lfres$long_p,10),type="n")
plot(1:nrow(lfres),-log(lfres$long_p,10),type="n",ylab="-log10(p value)",las=1,xlab="Locus",ylim= c(0,max(-log(lfres[,2:4],10))), main = "LFMM outlier loci for two PCs")

points(1:nrow(lfres),-log(lfres$long_p,10),pch=20,col=as.factor(lfres$lo_q<0.05))

points(rownames(lfres)[lfres$bu_q<alpha],-log(lfres$burn_unburnt_p[lfres$bu_q<alpha],10),pch=20,col="green")
points(rownames(lfres)[lfres$lo_q<alpha],-log(lfres$long_p[lfres$lo_q<alpha],10),pch=20,col="blue")
legend("topleft",legend=c("PC longtitude"),pch=20,col=c("blue"))
# legend("topleft",legend=c("PC1 burnt or unburnt","PC2 longtitude"),pch=20,col=c("green","blue"))


# Plot effect sizes:
lfres$efs1<-efs.lfmm1
lfres$efs2<-efs.lfmm2
head(lfres)


# BD ggplot
lfres$rowname<-lfres %>% row.names()
lfres_df_efs<-lfres %>% subset(never_outl == 0)

ggplot(data = lfres_df_efs,mapping = aes(x= rowname)) +
  geom_point(aes(x=rowname, y=efs1)) +
  geom_point(aes(x=rowname, y=efs2)) +   
  geom_hline(yintercept = 0,colour = "purple") +
  geom_point(data = lfres_df_lo, aes(x = as.numeric(rowname),efs2), colour = "red") +
  geom_point(data = lfres_df_bu, aes(x = as.numeric(rowname),efs1), colour = "green")+
  labs(y ="Effect Size",
       x ="Locus") +
  theme_article()



quartz("",8,4,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),lfres$efs1,type="n",ylab="Effect size",las=1,xlab="Locus",ylim=c(min(c(lfres$efs1,lfres$efs2,lfres$efs3))+-0.5,max(c(lfres$efs1,lfres$efs2,lfres$efs3))+0.5))

points(rownames(lfres)[lfres$never_outl==0],lfres$efs1[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],lfres$efs2[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))

points(rownames(lfres)[lfres$bu_q<alpha],lfres$efs1[lfres$bu_q<alpha],pch=20,col="red")
points(rownames(lfres)[lfres$lo_q<alpha],lfres$efs2[lfres$lo_q<alpha],pch=20,col="green")

abline(0,0,col="purple")

# xofs<-rep(c(700,-700),length(which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1)))[1:length(which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1))]
# yofs<-rep(c(0.04,-0.04),length(which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1)))[1:length(which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1))]

text(which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1)+xofs,lfres$efs1[which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1)]+yofs,labels=lfres$locus_p[which(abs(lfres$efs1)>0.5 & lfres$bu_outl==1)])

## QQ plot:
qqplot(rexp(length(pvalues), rate = log(10)),   -log10(pvalues), xlab = "Expected quantile",  pch = 19, cex = .4)
abline(0,1)

# Extra outlier loci plot from plot_freq_location()
lfmm_loc.toanalyse<-as.character(lfres$locus_p[which(lfres$lo_outl==1)]) 
length(lfmm_loc.toanalyse)
head(lfmm_loc.toanalyse)

# Location for files:	
out.dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LFMM"
out.dir

head(lfmm_loc.toanalyse)
ghead(gt_data)
head(sdat,3)
dir(out.dir)

# close LFMM ----

#  COMBINE RESULTS:    	# ----

head(bslinf)
head(pca_loc)
head(lfr)

dim(bslinf)
dim(pca_loc)
dim(lfr)

outl_all<-merge(bslinf,pca_loc,by=c("lind","locus"),all.x=T,all.y=F)
outl_all<-merge(outl_all,lfr,by=c("locus"),all.x=T,all.y=F)
outl_all<-outl_all[order(outl_all$lind),]
outl_all<-tidy.df(outl_all)
head(outl_all)

outl_all$lfmm_outl<-ifelse(rowSums(outl_all[,c("lfmm_PC1","lfmm_PC2")])>0,1,0)
outl_all[outl_all$lfmm_PC1==1,]
check.rows(outl_all)

# There is very little overlap among three methods:
table(rowSums(outl_all[,c("bs_outl","pca_outl","lfmm_outl")]))

# write.table(outl_all,"outliers_all.txt",sep="\t",row.names=F,quote=F)

# close combine ----





