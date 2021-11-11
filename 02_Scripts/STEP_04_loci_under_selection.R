
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### PCAdapt and LFMM to detect loci under selection
### Updated after first peer-review

### Author: Annabel Smith, Binyin Di

# AS Load and tidy workspace and remove everything except necessary objects:
load("buffel_burning.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

load("03_Workspaces/STEP04_loci_under_selection.RData")

#  PCADAPT:    	# ----

# tutorial:
# browseVignettes("pcadapt")
# or: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html  

library(pcadapt)
library(qvalue)

# to install qvalue:
# https://www.bioconductor.org/packages/release/bioc/html/qvalue.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

path_to_file <- "04_RESULTS/Loci_under_selection/PCAdapt_revised_Nov2021/PCAdapt_files"

dir(path_to_file)

filename <- read.pcadapt(paste(path_to_file,"Cenchrus_filt3.bed",sep="/"), type = "bed")

### --- *** SET INITIAL K *** --- ###

K <- 19 # Begin with the number of sample sites
x <- pcadapt(input = filename, K = 19) 
summary(x)

# This is the proportion variance explained
# It's described as a vector containing the K ordered square root of the proportion of variance explained by each PC (https://bcm-uga.github.io/pcadapt/articles/pcadapt.html.
# But it's the square root of the following that's written in the function (https://github.com/bcm-uga/pcadapt/blob/master/R/pcadapt.R)
pvarexp<-x$singular.values[1:K]^2 / sum(x$singular.values^2)

# CHOOSE K FROM SCREE PLOT:
dev.new(width=4,height=4,dpi=160,pointsize=12, noRStudioGD = T)
par(mar=c(4,4,1,0.5))
plot(1:K,pvarexp,xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n", main = "", font.main=1)
title(xlab="Principal component",mgp=c(2.5,1,0))
grid()
lines(1:K,pvarexp)
points(1:K,pvarexp,pch=20)

# CHOOSE K FROM PCs:

# Get population names from the fam file:
famf<-read.table(paste(path_to_file,"Cenchrus_filt3.fam",sep="/"),header=F) # strange fam file
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

K<-3 # return 
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

# Get outliers based on q values:

loc_dat<-read.table(paste(path_to_file,"Cenchrus_filt3_loci.txt",sep="/"),header=T) 
loc_dat$pvalue<-x$pvalues
head(loc_dat)

# Get outliers based on qvalues
loc_dat$qvalue <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(loc_dat$qvalue < alpha)
loc_dat$outlier <- 0
loc_dat$outlier[which(loc_dat$qvalue < alpha)]<-1
length(outliers)
table(loc_dat$outlier)
head(loc_dat); range(loc_dat$outlier)

# Associate outliers with PCs: library(pcadapt)
snp_pc <- get.pc(x,outliers)
head(snp_pc); dim(snp_pc)
range(snp_pc$PC)
range(snp_pc$SNP)

# Add them to the main data frame:
loc_dat<-merge(loc_dat, snp_pc, by.x="lind", by.y="SNP", all.x=T, all.y=F)
loc_dat$PC[which(is.na(loc_dat$PC))]<-0
check.rows(loc_dat)

# most are outlying on PC2, none are outlying on PC3:
table(loc_dat$PC)

# write.table(loc_dat,"outliers_PCAdapt_from_snp_pc.txt",sep="\t",row.names=F,quote=F)

# save.image("03_Workspaces/STEP04_loci_under_selection.RData")

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

lfmm_dir<-"04_RESULTS/Loci_under_selection/LFMM_revised_Nov2021/LFMM_files"
dir(lfmm_dir)

# Genotype data:
lfdat<-as.matrix(read.table(paste(lfmm_dir,"lfmm.txt",sep="/"),header=F)) # a few seconds 
colnames(lfdat)<-NULL
ghead(lfdat); dim(lfdat)

# Site data for each individual in the genotype file:
lfsite<-read.table(paste(lfmm_dir,"lfmm_site.txt",sep="/"),header=T)
head(lfsite)

# Enviro data for each site:

# Add site info
sdat$site_code <- sub("buf","X",sdat$site)
head(sdat,2); dim(sdat) 

# enviro data to include: site_code, burnt_unburnt, and long:
lfsd<-sdat[which(sdat$site_code %in% unique(lfsite$site)),c(which(colnames(sdat)=="site"),which(colnames(sdat)=="burn_unburnt"),which(colnames(sdat)=="long"):which(colnames(sdat)=="site_code"))] 
lfsd<-lfsd[order(lfsd$site),]
lfsd<-tidy.df(lfsd)
lfsd$site<-NULL
head(lfsd); dim(lfsd)

# Make sure they're in the same order, these should all be T:
lfsd$site_code==unique(lfsite$site)

# This is for the biplot labels:
rownames(lfsd)<-lfsd[,"site_code"] 

# burn_unburn to binary code: burnt is 1, otherwise, 0
lfsd$burn_unburnt<-ifelse(lfsd$burn_unburnt =="b",1,0)

# merge individual level (nrow==93) data with site data (nrow==19):
lfsd2<-lfsite
head(lfsd); dim(lfsd)
head(lfsd2); dim(lfsd2)

lfsd3<-merge(lfsd2, lfsd, by.x = "site", by.y = "site_code",  all.x = TRUE, all.y = FALSE)

#check:
table(lfsite$ind == lfsd3$ind) # should all be T
head(lfsd3); dim(lfsd3)

# Add K cluster:
kdat<-read.table("00_Data/K_genetic_clusters_Cenchrus_filt2.txt",header=T)
# Check individual names correspond:
table(kdat$indiv %in% lfsd3$ind) # should all be T
table(lfsd3$ind %in% kdat$indiv ) # should all be T
kdat<-kdat[,c("indiv","K3")]
head(kdat,3); dim(kdat)
# Merge with lfsd3

lfsd3<-merge(lfsd3, kdat, by.x="ind", by.y="indiv", all.x=T, all.y=F)
head(lfsd3); dim(lfsd3)

# Reduce columns to burnt_unburnt, long and K - these are the environmental variables we're using in the model:
lfs<-lfsd3[,c("burn_unburnt","long", "K3")]
head(lfs); dim(lfs)

# The environmental matrix (X):
lfs<-as.matrix(lfs[,1:3])
head(lfs); dim(lfs)

# The genotype matrix (Y):
ghead(lfdat); dim(lfdat)

# CHOOSE K FROM SCREE PLOT:
pc <- prcomp(lfdat)
head(pc$sdev); length(pc$sdev)
summary(pc)$importance[,1:5]
str(pc)

dev.new(width=4,height=4,dpi=160,pointsize=12, noRStudioGD = T)
par(mar=c(4,4,1,0.5))
plot(1:20,summary(pc)$importance[2,][1:20],xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n", main = "", font.main=1)
title(xlab="Principal component",mgp=c(2.5,1,0))
grid()
lines(1:20,summary(pc)$importance[2,][1:20])
points(1:20,summary(pc)$importance[2,][1:20],pch=20)

# The scree from PCAdapt suggests K==3, while the scree from LFMM suggests K==5

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
efs.lfmm3<-effect_size(Y = lfdat, X = as.matrix(lfs[,3]), lfmm = mod.lfmm)
# save.image("03_Workspaces/STEP04_loci_under_selection.RData")

head(efs.lfmm1)
head(efs.lfmm2)
head(efs.lfmm2)

lfmm_loci<-read.table(paste(lfmm_dir,"lfmm_loci_filt3.txt",sep="/"),header=T)
head(lfmm_loci)
dim(lfmm_loci)

# save.image("03_Workspaces/STEP04_loci_under_selection.RData")

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
lfres$K_outl<-ifelse(rownames(lfres) %in% lf_outl$K3_q,1,0)
head(lfres)

table(lfres$bu_outl)
table(lfres$lo_outl)
table(lfres$K_outl)

# no loci were outliers along the burn gradient with a q threshold of 0.05
# 12 loci were outliers along the longitudinal gradient with a q threshold of 0.05
# 22 loci were outliers along the K gradient with a q threshold of 0.05

# Reproduce manhattan plot:
lfres$never_outl<-rowSums(lfres[,which(colnames(lfres)=="bu_outl"):ncol(lfres)])
lfres$never_outl<-rowSums(lfres[,8,drop = FALSE]) 
head(lfres)

# AS manhattan:
dev.new(width=10,height=4,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),-log(lfres$long_p,10),type="n")
plot(1:nrow(lfres),-log(lfres$long_p,10),type="n",ylab="-log10(p value)",las=1,xlab="Locus",ylim= c(0,max(-log(lfres[,2:4],10))), main = "")

points(1:nrow(lfres),-log(lfres$long_p,10),pch=20,col=as.factor(lfres$lo_q<0.05))

# none were outlying on the burn gradient, so plot longitude and K outliers only:
points(rownames(lfres)[lfres$bu_q<alpha],-log(lfres$burn_unburnt_p[lfres$bu_q<alpha],10),pch=20,col="blue")
points(rownames(lfres)[lfres$lo_q<alpha],-log(lfres$long_p[lfres$lo_q<alpha],10),pch=20,col="green")
points(rownames(lfres)[lfres$K3_q<alpha],-log(lfres$K3_p[lfres$K3_q<alpha],10),pch=20,col="cornflowerblue")
legend("topleft",legend=c("longtitude", "K"),pch=20,col=c("green","cornflowerblue"))

# Plot effect sizes:
lfres$efs1<-efs.lfmm1
lfres$efs2<-efs.lfmm2
lfres$efs3<-efs.lfmm3
head(lfres,3); dim(lfres)

# write.table(lfres,"lfmm_outliers.txt",sep="\t",quote=F, row.names=F)

dev.new(width=8,height=4, noRStudioGD = T,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),lfres$efs1,type="n",ylab="Effect size",las=1,xlab="Locus",ylim=c(min(c(lfres$efs1,lfres$efs2,lfres$efs3))+-0.5,max(c(lfres$efs1,lfres$efs2,lfres$efs3))+0.5))

points(rownames(lfres)[lfres$never_outl==0],lfres$efs1[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],lfres$efs2[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],lfres$efs3[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))

points(rownames(lfres)[lfres$bu_q<alpha],lfres$efs1[lfres$bu_q<alpha],pch=20,col="red")
points(rownames(lfres)[lfres$lo_q<alpha],lfres$efs2[lfres$lo_q<alpha],pch=20,col="green")
points(rownames(lfres)[lfres$K3_q<alpha],lfres$efs3[lfres$K3_q<alpha],pch=20,col="cornflowerblue")

abline(0,0,col="purple")

legend("topleft",legend=c("longtitude","K"),pch=20,col=c("green","cornflowerblue"))

## QQ plot:
dev.new(width=8,height=4, noRStudioGD = T,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
qqplot(rexp(length(pvalues), rate = log(10)),   -log10(pvalues), xlab = "Expected quantile",  pch = 19, cex = .4)
abline(0,1)

# save.image("03_Workspaces/STEP04_loci_under_selection.RData")

# close LFMM ----

#  COMBINE RESULTS:    	# ----

pcadap_linf<-loc_dat[,c("lind","locus","outlier")]
colnames(pcadap_linf)[which(colnames(pcadap_linf)=="outlier")]<-"pca_outl"
lfmm_linf<-lfres[,c("locus_p","bu_outl","lo_outl","K_outl")]
lfmm_linf$lfmm_outl<-rowSums(lfmm_linf[,c("bu_outl","lo_outl","K_outl")])
table(lfmm_linf$lfmm_outl==lfmm_linf$bu_outl)
table(lfmm_linf$lfmm_outl==lfmm_linf$K_outl)
table(lfmm_linf$lfmm_outl==lfmm_linf$lo_outl)

lfmm_linf<-lfmm_linf[,c("locus_p","lfmm_outl")]

head(pcadap_linf,3); dim(pcadap_linf) # PCAdapt
head(lfmm_linf,3); dim(lfmm_linf) # LFMM

table(pcadap_linf$pca_outl)
table(lfmm_linf$lfmm_outl)

outl_all<-pcadap_linf

outl_all<-merge(outl_all,lfmm_linf,by.x="locus",by.y="locus_p",all.x=T,all.y=F)
outl_all<-outl_all[order(outl_all$lind),]
outl_all<-tidy.df(outl_all)
head(outl_all); dim(outl_all)
check.rows(outl_all)

# LFMM = 34 outliers
# PCAdapt = 5005 outliers
# Overlap between methods = 9 outliers
table(rowSums(outl_all[,c("pca_outl","lfmm_outl")]))

# write.table(outl_all,"outliers_all.txt",sep="\t",row.names=F,quote=F)

# save.image("03_Workspaces/STEP04_loci_under_selection.RData")

# VENN diagram:

library(VennDiagram)

# Summarise overlap:
res2<-outl_all[,c("pca_outl","lfmm_outl")]
res2$sum1<-rowSums(cbind(res2[,1],res2[,2]))

res2$sum1<-ifelse(res2$sum1==1,0,res2$sum1)
res2$sum1<-ifelse(res2$sum1==2,1,res2$sum1)
head(res2); table(res2$sum1)

res3<-colSums(res2)
res3

dev.new(width=4,height=4,dpi=160,pointsize=12, noRStudioGD = T)
par(mar=c(4,4,1,0.5))
venn.plot<-draw.pairwise.venn(res3[1],res3[2],res3[3], category=c("PCAdapt","LFMM"),scaled=F,fill=rgb(0,0,0,0.5),fontfamily="sans",cat.fontfamily="sans",cex=1, cat.pos=c(2,10),lwd=1)

pdf(file="venn.pdf",width=4,height=4,pointsize=12)
par(mar=c(4,4,1,1))
grid.draw(venn.plot)
dev.off()

# close combine ----





