
# ------------------------------------ #
# ------------- STEP 03  ------------- #
# ------------------------------------ #

### Post-process results from analyses to detect loci under selection
### Author: Annabel Smith, Binyin Di

# Load workspace (not updated for Cenchrus):
# load("03_Workspaces/STEP03_sel_wksp")
# load("03_Workspaces/STEP03_sel_wksp_BD.RData")

# AS Load and tidy workspace and remove everything except necessary objects:
load("binyin_winter.RData"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# BD load R.Data and Functions
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/binyin_winter.RData")

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

#  ANALYSE BAYESCAN:    	# ----

# *** ANALYSE OUTLIER LOCI:

# Rename data:
gt_data<-snp_onerow

# Directory with bayescan results:
bs_dir<-"../Offline_Analysis/BayeScan/Cenchrus_filt1/Cenchrus_BS_po400_RESULTS"
bs_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/Offline Winter Project/BayeScan/Cenchrus_BS_po400_RESULTS"
dir(bs_dir)

bs_fst<-paste(bs_dir,"Cenchrus_BS_po400_fst.txt",sep="/")

# FST Outlier loci:
bsres<-plot_bayescan(bs_fst,FDR=0.05)
bs_outl<-bsres$outliers
bs_n_outl<-bsres$nb_outliers

# MARKER file to separate non-neutral loci:

# Get locus index (this is the locus index file that is made in format_structure() function for bayescan):
bslinf<-read.table("../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION/BayeScan/bayescan_filt3/bs_loci_filt3.txt",header=T)
bslinf<-read.table("C:/Users/s4467005/OneDrive - The University of Queensland/Offline Winter Project/BayeScan/bs_loci_filt1.txt", header = TRUE) # bs_loci_filt1.txt
head(bslinf)

bsoutl<-data.frame(lind=bs_outl,outl=1)
head(bsoutl)

# Outlier loci:
bslinf<-merge(bslinf,bsoutl,by="lind",all.x=T,all.y=F)
bslinf$outl[which(is.na(bslinf$outl))]<-0
head(bslinf)
table(bslinf$outl)
check.rows(bslinf)
colnames(bslinf)[which(colnames(bslinf)=="outl")]<-"bs_outl"

# write.table(bslinf,file="bayescan_outliers.txt",quote=F,row.names=F,sep="\t")

# *** PLOT DIAGNOSTICS:

bs_sel<-paste(bs_dir,"Cenchrus_BS_po400.sel",sep="/")
seldat<-read.table(bs_sel,colClasses="numeric")

# Plot log likelihood:
parameter<-"logL" # a few minutes
plot(density(seldat[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution PO=400"))

# Plot FST:
quartz("",5,8,dpi=100) # Error
par(mfrow=c(10,6),mar=c(2,2,0.2,0.2),mgp=c(2,0.5,0))
for (i in grep("Fst",colnames(seldat))){
par.thisrun<-colnames(seldat)[i]
plot(density(seldat[[par.thisrun]]),xlab=par.thisrun,main="",cex.axis=0.75)
legend("bottom",legend=par.thisrun,cex=1,bty="n")
}

dim(seldat)

alphacols<-colnames(seldat)[grep("alpha",colnames(seldat))]
head(alphacols)
outl_alpha<-paste("alpha",bs_outl,sep="")
head(outl_alpha)
table(outl_alpha %in% alphacols)

# TRUE 
# 2561

# Plot alpha for a random selection of loci:
quartz("",12,8,dpi=80)
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

# Put files here:	
out.dir<-"../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION/BayeScan/bayescan_filt3/bayescan_filt3_50_random_heatmaps"
out.dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/BayeScan"

out.dir
## Make sure they're all in site data

# sdat$site_code <- sub("buf","X",sdat$site)
# unique(gt_data$site) %in% sdat$site_code # add a new column

# PLOT GENOTYPE FREQUENCIES BY LOCATION:

# Takes about 15 seconds for 100

head(loc.toanalyse)
ghead(gt_data)
head(sdat,3)
dir(out.dir)


plot_freq_location(loc.toanalyse,gt_data,sdat,out.dir,50) # Error: requirement for site_data [c("site_code","native","country","region","latitude")] needs to be changed.

# From the manual:
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
dir(bs_dir)

# you can plot population specific Fst coefficient by setting
parameter<-"Fst1"

# if you test for selection, you can plot the posterior for alpha coefficient for selection:
parameter<-"alpha41121"

# you also have access to the likelihood with:
parameter<-"logL"

# if you have the package "boa" installed, you can very easily obtain Highest Probability 
library(boa)

# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
boa.hpd(seldat[[parameter]],0.05)

# Lower Bound Upper Bound 
# -995514.0   -994271.8 


# close BayeScan ----



#########################################
####	    	PCADAPT:	    	 ####
#########################################
{

  
# Tutorial: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html  
  

install.packages("pcadapt")  
install.packages("https://cran.r-project.org/src/contrib/Archive/qvalue/qvalue_1.26.0.tar.gz", repos = NULL, type = "source")  

library(pcadapt)
library(qvalue)

# tutorial:
# browseVignettes("pcadapt")

path_to_file <- "../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION/PCAdapt/pcadapt_filt2/pcadapt_files"
path_to_file<- "C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1"
path_to_file<- "D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1"

dir(path_to_file)

filename <- read.pcadapt(paste(path_to_file,"Cenchrus_filt1.bed",sep="/"), type = "bed")

### --- *** CHOOSE K *** --- ###
K <- 22 # Capitalised K
  
x <- pcadapt(input = filename, K = 22) 

# CHOOSE K FROM SCREE PLOT: reproduce built in plot: plot(x, option = "screeplot")

# plot(x, option = "screeplot")

# install.packages("grDevices")
# library(grDevices)

quartz("",4,4,dpi=160,pointsize=12) # Error
par(mar=c(4,4,0.5,0.5))
plot(1:K,x$singular.values^2,xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
title(xlab="PC",mgp=c(2.5,1,0))
grid()
lines(1:K,x$singular.values^2)
points(1:K,x$singular.values^2,pch=20)

# CHOOSE K FROM PCs:

# Get population names from the fam file:
famf<-read.table(paste(path_to_file,"Cenchrus_filt1.fam",sep="/"),header=F)
head(famf)
poplist.names <- as.character(famf$V1)
head(poplist.names)

# Reproduce inbuilt PC plot: plot(x, option = "scores", i = 1, j = 2, pop = poplist.names)

quartz("",10,10,dpi=80,pointsize=14) # Error
par(mfrow=c(5,5),mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (i in seq(1,ncol(x$scores),2)){
pc1<-x$scores[,i]
pc2<-x$scores[,i+1]
plot(pc1,pc2,col=rainbow(length(unique(poplist.names))),xlab="",ylab="",pch=20)
title(xlab=paste("PC",i,sep=""),cex.lab=1)
title(ylab=paste("PC",i+1,sep=""),cex.lab=1)
}
par(xpd=NA)
legend(-5, -2,legend=unique(poplist.names),col=rainbow(length(unique(poplist.names))),pch=20,ncol=9)

# 10-25 pcs needed to explain variation

# Score plot - Example
# poplist.int<-c(rep(1,50), rep(2,50), rep(3,50))
# print(poplist.int)

# poplist.names<-c(rep("POP1", 50), rep("POP2",50), rep("POP3",50))
# print(poplist.names)

# plot(x, option = "scores", pop = poplist.int)
# plot(x, option = "scores", pop = poplist.names)




### --- *** SET K *** --- ###

library(pcadapt)

K<-3 # return # 182
x <- pcadapt(input = filename, K = K) 
summary(x)


# > summary(x)
# Length Class  Mode   
# scores             930 -none- numeric
# singular.values     10 -none- numeric
# loadings        407110 -none- numeric
# zscores         407110 -none- numeric
# af               40711 -none- numeric
# maf              40711 -none- numeric
# chi2.stat        40711 -none- numeric
# stat             40711 -none- numeric
# gif                  1 -none- numeric
# pvalues          40711 -none- numeric
# pass             34752 -none- numeric



# Check for LD 
dev.new("",10,4,dpi=80,pointsize=14) # Error
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (i in 1:ncol(x$loadings))
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

# Reproduce PC plot:

head(x$scores)

dev.new("",10,4) # Error
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5,1,0))
pc1<-x$scores[,1]
pc2<-x$scores[,2]
pc3<-x$scores[,3]

colour_codes<-poplist.names
# colour_codes<-unique(poplist.names) # no need
colour_codes[grep("b",unique(poplist.names))]<-"red"
colour_codes[grep("u",unique(poplist.names))]<-"blue"


# PC should be a square shape when exporting the plots

plot(pc1,pc2,col=colour_codes,xlab="pc1",ylab="pc2",pch=20, cex=3)
text(pc1,pc2, labels = poplist.names)

plot(pc1,pc3,col=colour_codes,xlab="pc1",ylab="pc3",pch=20, cex=3)
text(pc1,pc3, labels = poplist.names)


plot(pc2,pc3,col=colour_codes,xlab="pc2",ylab="pc3",pch=20, cex=3)
text(pc2,pc3, labels = poplist.names)

# title(xlab=paste("PC",i,sep=""),cex.lab=1)
# title(ylab=paste("PC",i+1,sep=""),cex.lab=1)

par(xpd=NA)
legend(-5,-1,legend=unique(poplist.names),col=rainbow(length(unique(poplist.names))),pch=20,ncol=9)

# lattice methods 

install.packages("lattice", dependencies = TRUE)
library(lattice)
show.settings()

install.packages("gridExtra")
library(gridExtra)
require(gridExtra) 
require(lattice)

plot1<-xyplot(pc2~pc1, scales=list(cex=1, col="red"),
              col=colour_codes, 
              xlab="pc1", ylab="pc2",
              main="pc1 v pc2",
              panel.text(pc1, pc2, labels = poplist.names))

  
plot2<-xyplot(pc2~pc3, scales=list(cex=1, col="red"),
              col=colour_codes,
              xlab="pc3", ylab="pc2",
              main="pc2 v pc3")  
  
plot3<-xyplot(pc1~pc3, scales=list(cex=1, col="red"),
              col=colour_codes,
              xlab="pc3", ylab="pc1",
              main="pc1 v pc3")

grid.arrange(plot1,plot2,plot3, ncol=3) 

text(pc1~pc3, labels = unique(poplist.names))

# alternative, multi-panel by ggpairs() 

library("tidyr")
pc_long <- x$scores %>%
  data.frame() %>%
  rename(PC1 = 1, PC2 = 2, PC3 = 3) %>%
  pivot_longer(cols = 1:3, 
               names_to = "PC",
               values_to = "Scores")
# class(x$scores)
  
pc_mat<-x$scores
names(pc_mat)<-c("PC1", "PC2", "PC3")
pairs(x$scores, upper.panel = NULL)

# GGally::ggpairs()
pc_df<-data.frame(pc_mat)
install.packages("GGally")
library(GGally)
ggpairs(pc_df,  aes(colour = colour_codes), columns = 1:3) +
  theme_bw() +
  scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme(text = element_text(size = 25)) 


# Tidyverse, ggplot
install.packages("directlabels", repos = "http://r-forge.r-project.org", dependencies = TRUE)
install.packages("grid")
install.packages("quadprog")
library(directlabels)



library(tidyverse)
library(ggrepel)
install.packages("egg")
library(egg)


# colour-trtments
# label-site_name

site_name<-str_extract(poplist.names,"[0-9]+")
colour_codes[grep("b",poplist.names)]<-"red"
colour_codes[grep("u",poplist.names)]<-"blue"
plotdata <- data.frame(pc1,pc2,pc3,colour_codes,site_name)

plot1<-ggplot(plotdata,mapping = aes(x = pc1, y = pc2, colour = colour_codes),
       xlab="pc1", ylab="pc2",
       main="pc1 v pc2") +
  geom_point(size = 2,alpha = 0.25) +  scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme_bw() +
  theme(text = element_text(size = 25)) 
  
# theme(axis.text = element_text(size = 20))               # Axis text size
# For more: https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend

# label method 1
library(gghighlight)

plot1_1<- plot1 +
  #geom_label_repel(aes(label = site_name),nudge_y = .05) +
  gghighlight() +
  facet_wrap(~site_name)


# label method 2
plot1 + 
geom_label_repel(aes(label = site_name),
                 box.padding   = 0.01, 
                 point.padding = 0.025) +
  theme_article()

# segment.color = 'grey50'
# label method 3
plot1 +
  theme_bw()+
  geom_text(aes(label=ifelse(pc2<-0.1,as.character(poplist.names),'')),hjust=0,vjust=0)


plot2<-ggplot(plotdata,mapping = aes(x = pc1, y = pc3, colour = colour_codes),
              xlab="pc1", ylab="pc3",
              main="pc1 v pc3") +
  geom_point(size = 2,alpha = 0.25) + scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme_bw() +
  theme(text = element_text(size = 25)) 

plot3<-ggplot(plotdata,mapping = aes(x = pc2, y = pc3, colour = colour_codes),
              xlab="pc2", ylab="pc3",
              main="pc2 v pc3") +
  geom_point(size = 2,alpha = 0.25) + scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme_bw() +
  theme(text = element_text(size = 25)) 

# Original plots 

ggplot(mapping = aes(x = pc1, y = pc2, colour = poplist.names, shape = poplist.names),
              xlab="pc1", ylab="pc2",
              main="pc1 v pc2") +
  geom_point() +
  geom_text(aes(label = poplist.names))

# note: y~x | x,y
# library(Rmisc)
# multiplot(plot1_1, plot2_1, plot3_1, cols=3)


# For future, 3D images: 
install.packages("https://cran.r-project.org/src/contrib/Archive/plot3Drgl/plot3Drgl_1.0.tar.gz", repos = NULL, type = "source", dependencies = TRUE) 
install.packages("plot3Drgl", dependencies = TRUE)
install.packages("rgl")
library(rgl)
library(plot3Drgl)

# plotrgl()


# Get outliers based on q values:

## be careful with NA, qvalue: no NA, p or pval has NA, and is.na_remove is pval without NAs
loci <- read.pcadapt(paste(path_to_file,"Cenchrus_filt1.bed",sep="/"), type = "bed")

pval <-x$pvalues

length(pval)
# write.table(pval,file="pval.txt",quote=F,row.names=F,sep="\t")
is.na_remove<-x$pvalues[!is.na(x$pvalues)]
qval <- qvalue(is.na_remove)$qvalues

length(qval)

# write.table(qval,file="qval.txt",quote=F,row.names=F,sep="\t")

# --------------------------------------------
# https://statisticsglobe.com/r-is-na-function/
# Detect if there are any NAs
any(is.na(pval))

# Locate NA in data set via which()
which(is.na(pval))


# if game
for (i in 1:length(pval)) {
 if(is.na(pval[i])) {
   print("Damn, it's an NA")
 }  
  else {
    print("Wow, that's awesome")
  }
}

ifelse(is.na(pval), "Damn, it's an NA", "WOW, that's awesome")

# ---------------------------



alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
# [1st run]4070 [2nd run]7040

# Associate outliers with PCs: library(pcadapt)
snp_pc <- get.pc(x, outliers)
head(snp_pc)

# write.table(snp_pc,"outliers_PCAdapt_from_snp_pc.txt",sep="\t",row.names=F,quote=F)

# Reproduce manhattan plot:

qp<-data.frame(q = qval,p = is.na_remove)
head(qp)
plot(1:nrow(qp),-log(qp$p,10),type="n")
points(1:nrow(qp),-log(qp$p,10),pch=20,col=as.factor(qp$q<0.05))
# plot(1:nrow(qp),qp$q,pch=20,col=as.factor(qp$q<0.05))

# get locus indices for outliers:
pca_loc<-read.table("../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION/PCAdapt/pcadapt_filt2/pcadapt_filt2_loci.txt",header=T)
pca_loc<-read.table("D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1/Cenchrus_filt1_loci.txt", header = T)
pca_loc<-read.table("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1/Cenchrus_filt1_loci.txt", header = T)
pca_loc$pca_outl<-ifelse(pca_loc$lind %in% outliers,1,0)
table(pca_loc$pca_outl)


maindir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1"
orig_loci<-read.table(paste(maindir, "Cenchrus_filt1_loci.txt", sep="/"), header=T)
head(orig_loci)
head(pca_loc)

table(orig_loci$locus == pca_loc$locus)

dir()


# [1st run]

# > table(pca_loc$pca_outl)

# 0     1 
# 32653  8058 


# [2nd run]
# > table(pca_loc$pca_outl)

# 0     1 
# 33671  7040 



} # close pcadapt

#########################################
####	    	 LFMM:	    	 	 ####
#########################################
{
library(lfmm)
library(qvalue)

lfmm_dir<-"../ANALYSIS_RESULTS/LOCI_UNDER_SELECTION/LFMM/lfmm_filt2"
lfmm_dir<-"D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LFMM"
lfmm_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/LFMM"
dir(lfmm_dir)

# Genotype data:
lfdat<-as.matrix(read.table(paste(lfmm_dir,"lfmm.txt",sep="/"),header=F))
colnames(lfdat)<-NULL
ghead(lfdat)

# Site data for each individual in the genotype file:
lfsite<-read.table(paste(lfmm_dir,"lfmm_site.txt",sep="/"),header=T)
head(lfsite)

# Enviro data for each site:

# Add site info
# sdat$site_code <- sub("buf","X",sdat$site)
# head(sdat$site_code) 


# lfsd<-sdat[which(sdat$site_code %in% unique(lfsite$site)),c(which(colnames(sdat)=="site_code"),which(colnames(sdat)=="elevation"):which(colnames(sdat)=="sm"))] # What are they? 
lfsd<-sdat[which(sdat$site_code %in% unique(lfsite$site)),c(which(colnames(sdat)=="site_code"))]

# lfsd<-lfsd[order(lfsd$site_code),]

lfsd<-tidy.df(lfsd)
head(lfsd)

# Make sure they're in the same order, these should all be T:
lfsd$site_code==unique(lfsite$site)

# This is for the biplot labels:
rownames(lfsd)<-lfsd[,"site_code"] 

# Make PCs:
lf_pc<-princomp(lfsd[,2:ncol(lfsd)],cor=T)
summary(lf_pc)
# biplot(lf_pc,xlab="PC1",ylab="PC2",cex=0.7)
lf_pc$loadings
head(lf_pc$scores)

# Add PCs to individual level data:
lfpcs<-data.frame(site=rownames(lfsd),lf_pc$scores[,1:3])
lfpcs<-tidy.df(lfpcs)
head(lfpcs)

lfs<-data.frame(site=lfsite[,1])
head(lfs)
lfs<-merge(lfs,lfpcs,by="site",all.x=T,all.y=F)
head(lfs); dim(lfs)

# The environmental matrix (X):
lfs<-as.matrix(lfs[,2:ncol(lfs)])
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
quartz("",4,4,dpi=160,pointsize=12)
par(mar=c(4,4,0.5,0.5))
plot(1:40,summary(pc)$importance[2,][1:40],xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
title(xlab="PC",mgp=c(2.5,1,0))
grid()
lines(1:40,summary(pc)$importance[2,][1:40])
points(1:40,summary(pc)$importance[2,][1:40],pch=20)

# The screes from the two methods are remarkably similar and K is similarly hard to choose, being somewhere between 10 and 25

# Ridge estimates and GWAS tests
head(lfs) # The environmental matrix (X):
ghead(lfdat) # The genotype matrix (Y):

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = lfdat, X = lfs, K = 10)
summary(mod.lfmm)
head(mod.lfmm$B)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = lfdat, X = lfs,  lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue 
head(pvalues)
length(pvalues)

# Direct effect sizes estimated from latent factor models:
efs.lfmm1<-effect_size(Y = lfdat, X = as.matrix(lfs[,1]),  lfmm = mod.lfmm)
efs.lfmm2<-effect_size(Y = lfdat, X = as.matrix(lfs[,2]),  lfmm = mod.lfmm)
efs.lfmm3<-effect_size(Y = lfdat, X = as.matrix(lfs[,3]),  lfmm = mod.lfmm)

head(efs.lfmm1)
length(efs.lfmm3)

lfmm_loci<-read.table(paste(lfmm_dir,"lfmm_loci_filt2.txt",sep="/"),header=T)
head(lfmm_loci)
dim(lfmm_loci)

# Add p and q values to locus info:
lfres<-data.frame(locus=lfmm_loci$locus,pvalues)
colnames(lfres)<-paste(colnames(lfres),"_p",sep="")

lfqs<-apply(lfres[,2:ncol(lfres)],2,function(x)qvalue(x)$qvalues)
colnames(lfqs)<-paste(substr(colnames(lfqs),1,6),"_q",sep="")
lfres<-cbind(lfres,lfqs)
head(lfres)

# Get outliers based on q values:
alpha <- 0.1
lf_outl<-apply(lfres[,5:ncol(lfres)],2,function(x) which(x<alpha))

lfres$Comp.1_outl<-ifelse(rownames(lfres) %in% lf_outl$Comp.1_q,1,0)
lfres$Comp.2_outl<-ifelse(rownames(lfres) %in% lf_outl$Comp.2_q,1,0)
lfres$Comp.3_outl<-ifelse(rownames(lfres) %in% lf_outl$Comp.3_q,1,0)
head(lfres)

table(lfres$Comp.1_outl)
table(lfres$Comp.2_outl)
table(lfres$Comp.3_outl)

# Reproduce manhattan plot:
lfres$never_outl<-rowSums(lfres[,which(colnames(lfres)=="Comp.1_outl"):ncol(lfres)])
head(lfres)

quartz("",8,4,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),-log(lfres$Comp.1_p,10),type="n",ylab="-log10(p value)",las=1,xlab="locus",ylim=c(0,max(-log(lfres[,2:4],10))))

points(rownames(lfres)[lfres$never_outl==0],-log(lfres$Comp.1_p[lfres$never_outl==0],10),pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],-log(lfres$Comp.2_p[lfres$never_outl==0],10),pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],-log(lfres$Comp.3_p[lfres$never_outl==0],10),pch=20,cex=1,col=rgb(0,0,0,0.08))

points(rownames(lfres)[lfres$Comp.1_q<alpha],-log(lfres$Comp.1_p[lfres$Comp.1_q<alpha],10),pch=20,col="red")
points(rownames(lfres)[lfres$Comp.2_q<alpha],-log(lfres$Comp.2_p[lfres$Comp.2_q<alpha],10),pch=20,col="green")
points(rownames(lfres)[lfres$Comp.3_q<alpha],-log(lfres$Comp.3_p[lfres$Comp.3_q<alpha],10),pch=20,col="blue")

legend("topleft",legend=c("PC1 mt,sp,mm,sm","PC2 st,ap","PC3 Elevation"),pch=20,col=c("red","green","blue"))

# Plot effect sizes:
lfres$efs1<-efs.lfmm1
lfres$efs2<-efs.lfmm2
lfres$efs3<-efs.lfmm3
head(lfres)

quartz("",8,4,dpi=100)
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(1:nrow(lfres),lfres$efs1,type="n",ylab="Effect size",las=1,xlab="Locus",ylim=c(min(c(lfres$efs1,lfres$efs2,lfres$efs3))+-0.5,max(c(lfres$efs1,lfres$efs2,lfres$efs3))+0.5))

points(rownames(lfres)[lfres$never_outl==0],lfres$efs1[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],lfres$efs2[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))
points(rownames(lfres)[lfres$never_outl==0],lfres$efs3[lfres$never_outl==0],pch=20,cex=1,col=rgb(0,0,0,0.08))

points(rownames(lfres)[lfres$Comp.1_q<alpha],lfres$efs1[lfres$Comp.1_q<alpha],pch=20,col="red")
points(rownames(lfres)[lfres$Comp.2_q<alpha],lfres$efs2[lfres$Comp.2_q<alpha],pch=20,col="green")
points(rownames(lfres)[lfres$Comp.3_q<alpha],lfres$efs3[lfres$Comp.3_q<alpha],pch=20,col="blue")

abline(0,0,col="purple")

xofs<-rep(c(700,-700),length(which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1)))[1:length(which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1))]
yofs<-rep(c(0.04,-0.04),length(which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1)))[1:length(which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1))]

text(which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1)+xofs,lfres$efs1[which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1)]+yofs,labels=lfres$locus_p[which(abs(lfres$efs1)>0.5 & lfres$Comp.1_outl==1)])

## QQ plot:
qqplot(rexp(length(pvalues), rate = log(10)),   -log10(pvalues), xlab = "Expected quantile",  pch = 19, cex = .4)
abline(0,1)

# Make results summary
# Old version was "lfr", saved in wksp but not defined in script
# lfres however, has been updated:

lfr_old<-lfr
lfr_new<-lfres[,c("locus_p","Comp.1_outl","Comp.2_outl","Comp.3_outl")]
colnames(lfr_new)<-c("locus","lfmm_PC1","lfmm_PC2","lfmm_PC3")
head(lfr_new)

# the loci detected with the new enviro data are DIFFERENT!
head(lfr_old); dim(lfr_old)
head(lfr_new); dim(lfr_new)
check.rows(cbind(lfr_old,lfr_new))
cbind(lfr_old,lfr_new)[which(lfr_old$lfmm_PC1==1),]
cbind(lfr_old,lfr_new)[which(lfr_old$lfmm_PC2==1),]
cbind(lfr_old,lfr_new)[which(lfr_old$lfmm_PC3==1),]

head(lfres)

lfr<-lfr_new

} # close lfmm

#########################################
####    	 COMBINE RESULTS:	  	 ####
#########################################
{

head(bslinf)
head(pca_loc)
head(lfr)

outl_all<-merge(bslinf,pca_loc,by=c("lind","locus"),all.x=T,all.y=F)
outl_all<-merge(outl_all,lfr,by=c("locus"),all.x=T,all.y=F)
outl_all<-outl_all[order(outl_all$lind),]
outl_all<-tidy.df(outl_all)
head(outl_all)

outl_all$lfmm_outl<-ifelse(rowSums(outl_all[,c("lfmm_PC1","lfmm_PC2","lfmm_PC3")])>0,1,0)
outl_all[outl_all$lfmm_PC1==1,]
check.rows(outl_all)

# There is very little overlap among three methods:
table(rowSums(outl_all[,c("bs_outl","pca_outl","lfmm_outl")]))

# write.table(outl_all,"outliers_all.txt",sep="\t",row.names=F,quote=F)


} # close combine





