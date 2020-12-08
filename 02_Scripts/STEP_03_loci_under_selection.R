
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

load("D:\\Onedrive\\OneDrive - The University of Queensland\\GitHub\\Binyin_Winter\\03_Workspaces\\divdist_ALL.RData")



#  POST-PROCESS BAYESCAN results:    	# ----

# *** ANALYSE OUTLIER LOCI:

# Rename data:
gt_data<-snp_onerow

# Directory with bayescan results:
bs_dir<-"D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/Offline_Results/Old Files//BayeScan/Cenchrus_filt2_BS_po400_RESULTS" # check fst file?
bs_dir<-"D:\\Onedrive\\OneDrive - The University of Queensland\\Offline Winter Project\\Cenchrus_BS_po400_RESULTS" 


bs_dir<-"../Offline_Results/BayeScan/Cenchrus_filt2_BS_po400_RESULTS"
dir(bs_dir)

bs_fst<-paste(bs_dir,"Cenchrus_filt2_BS_po400_fst.txt",sep="/") #check fst file

# FST Outlier loci:
# dev.new() 
bsres<-plot_bayescan(res=bs_fst,FDR=0.05) # error 
bs_outl<-bsres$outliers
bs_n_outl<-bsres$nb_outliers

# MARKER file to separate non-neutral loci:

# Get locus index (this is the locus index file that is made in format_structure() function for bayescan):

bslinf<-read.table("D:/Onedrive/OneDrive - The University of Queensland/Offline Winter Project/Old Files/BayeScan/bs_loci_filt1.txt", header = TRUE)
bslinf<-read.table("../Offline_Results/BayeScan/BayeScan_ANALYSIS/Cenchrus_filt2/bs_loci_filt2.txt", header = TRUE) 

head(bslinf); dim(bslinf)
bsoutl<-data.frame(lind=bs_outl,outl=1)
head(bsoutl); dim(bsoutl)

# Outlier loci:
bslinf<-merge(bslinf,bsoutl,by="lind", all.x = T, all.y = F)
bslinf$outl[which(is.na(bslinf$outl))]<-0
head(bslinf)
table(bslinf$outl)
check.rows(bslinf)

colnames(bslinf)[which(colnames(bslinf)=="outl")]<-"bs_outl"

# write.table(bslinf,file="bayescan_outliers.txt",quote=F,row.names=F,sep="\t")

# *** PLOT DIAGNOSTICS:

bs_sel<-paste(bs_dir,"Cenchrus_filt2_BS_po400.sel",sep="/")
seldat<-read.table(bs_sel,colClasses="numeric")

# Plot log likelihood:
parameter<-"logL" # a few minutes
plot(density(seldat[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution PO=400"))

# Plot FST:
quartz("",5,8,dpi=100) # Error
par(mfrow=c(3,3),mar=c(2,2,0.2,0.2),mgp=c(2,0.5,0))
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
# 5

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

# Location for heatmap files:	
out.dir<-"../Offline_Results/RESULTS/BayeScan/new_heatmaps"
out.dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/BayeScan/Cenchrus_filt2_BS_po400"
out.dir

## Make sure they're all in site data

# unique(gt_data$site) %in% sdat$site_code # add a new column

# PLOT GENOTYPE FREQUENCIES BY LOCATION:

# Takes about 15 seconds for 100

head(loc.toanalyse)
ghead(gt_data)
head(sdat,3)
dir(out.dir)
length(loc.toanalyse)

sdat$site_code <- sub("buf","X",sdat$site)
site_data<-sdat

# this plots by longitude (see plot_freq_long folder for results)
plot_freq_long(loc.toanalyse,gt_data,sdat,out.dir,5) 

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

# close BayeScan ----

#  PCADAPT:    	# ----

# tutorial:
# browseVignettes("pcadapt")

# Tutorial: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html  

library(pcadapt)
library(qvalue)
library(grDevices)

path_to_file <- "RESULTS/PCAdapt/PCAdapt_files"

path_to_file<- "D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/PCAdapt/PCAdapt_files"

dir(path_to_file)

filename <- read.pcadapt(paste(path_to_file,"Cenchrus_filt1.bed",sep="/"), type = "bed")

### --- *** CHOOSE K *** --- ###
K <- 22 # Capitalised K
  
x <- pcadapt(input = filename, K = 22) 

# CHOOSE K FROM SCREE PLOT: reproduce built in plot: plot(x, option = "screeplot")

quartz("",4,4,dpi=160,pointsize=12) # Error
par(mar=c(4,4,0.5,0.5))
plot(1:K,x$singular.values^2,xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
title(xlab="PC",mgp=c(2.5,1,0))
grid()
lines(1:K,x$singular.values^2)
points(1:K,x$singular.values^2,pch=20)

# CHOOSE K FROM PCs:

# Get population names from the fam file:
famf<-read.table(paste(path_to_file,"Cenchrus_filt1.fam",sep="/"),header=F) # strange fam file
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
legend(0, -2,legend=unique(poplist.names),col=rainbow(length(unique(poplist.names))),pch=20,ncol=9)

### --- *** SET K *** --- ###

K<-3 # return 
x <- pcadapt(input = filename, K = K) 
summary(x)



# Check for LD 
dev.new("",10,4,dpi=80,pointsize=14) # Error
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5,1,0))
for (i in 1:ncol(x$loadings)){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
  }


pc1<-x$scores[,1]
pc2<-x$scores[,2]
pc3<-x$scores[,3]



# colour-trtments
# label-site_name
library(stringr)
site_name<-str_extract(poplist.names,"[0-9]+")
colour_codes<-NA
colour_codes[grep("b",poplist.names)]<-"red"
colour_codes[grep("u",poplist.names)]<-"blue"
treatments<-NA
treatments[grep("b",poplist.names)]<-"Burnt"
treatments[grep("u",poplist.names)]<-"Unburnt"

plotdata <- data.frame(pc1,pc2,pc3,colour_codes,site_name,treatments)




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

library(tidyr)
library(tidyverse)
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
# install.packages("GGally")
library(GGally)


pc_df<-pc_df %>%
  pivot_longer(cols = 1:3, 
               names_to = "colour_codes",
               values_to = "Scores")


ggpairs(pc_df,  aes(colour = colour_codes)) +
  theme_bw() +
  scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme(text = element_text(size = 25)) 

# Tidyverse, ggplot
install.packages("quadprog")
install.packages("directlabels", repos = "http://r-forge.r-project.org", dependencies = TRUE)
install.packages("grid")
install.packages("egg")

library(directlabels)
library(tidyverse)
library(ggrepel)
library(egg)



plot1<-ggplot(plotdata,mapping = aes(x = pc1, y = pc2, colour = treatments),
       xlab="pc1", ylab="pc2",
       main="pc1 v pc2") +
  geom_point(size = 2,alpha = 0.25) +  scale_color_manual(values=c("Burnt"="red","Unburnt"="blue")) +
  theme_article() +
  theme(legend.position="bottom", text = element_text(size = 25)) +
  labs(x = expression(italic("PC 1")),
       y = expression(italic("PC 2"))) +
  theme_article()


plot2<-ggplot(plotdata,mapping = aes(x = pc1, y = pc3, colour = treatments),
              xlab="pc1", ylab="pc3",
              main="pc1 v pc3") +
  geom_point(size = 2,alpha = 0.25) +  scale_color_manual(values=c("Burnt"="red","Unburnt"="blue")) +
  theme_article() +
  theme(legend.position="bottom", text = element_text(size = 25)) +
  labs(x = expression(italic("PC 1")),
       y = expression(italic("PC 3"))) +
  theme_article()

plot3<-ggplot(plotdata,mapping = aes(x = pc2, y = pc3, colour = treatments),
              xlab="pc2", ylab="pc3",
              main="pc2 v pc3") +
  geom_point(size = 2,alpha = 0.25) +  scale_color_manual(values=c("Burnt"="red","Unburnt"="blue")) +
  theme_article() +
  theme(legend.position="bottom", text = element_text(size = 25)) +
  labs(x = expression(italic("PC 2")),
       y = expression(italic("PC 3"))) +
  theme_article()


ggarrange(plot1,plot2,plot3)
  
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
  theme_article()+
  geom_text(aes(label=ifelse(pc2<-0.1,as.character(poplist.names),'')),hjust=0,vjust=0)

plot2<-ggplot(plotdata,mapping = aes(x = pc1, y = pc3, colour = colour_codes),
              xlab="pc1", ylab="pc3",
              main="pc1 v pc3") +
  geom_point(size = 2,alpha = 0.25) + scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme_article() +
  theme(text = element_text(size = 25)) 

plot3<-ggplot(plotdata,mapping = aes(x = pc2, y = pc3, colour = colour_codes),
              xlab="pc2", ylab="pc3",
              main="pc2 v pc3") +
  geom_point(size = 2,alpha = 0.25) + scale_color_manual(values=c("red"="red","blue"="blue")) +
  theme_article() +
  theme(text = element_text(size = 25)) 

# Original plots 

ggplot(mapping = aes(x = pc1, y = pc2, colour = site_name, shape = treatments),
              xlab="pc1", ylab="pc2",
              main="pc1 v pc2") +
  geom_point(size = 10) +
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
snp_pc <- get.pc(x, outliers) # 3 PCs as K=3?
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
pca_loc<-read.table("D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/PCAdapt/PCAdapt_files/Cenchrus_filt1_loci.txt", header = T)

pca_loc$pca_outl<-ifelse(pca_loc$lind %in% outliers,1,0)
table(pca_loc$pca_outl)

maindir<-"D:/Onedrive/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/PCAdapt/PCAdapt_files/"
orig_loci<-read.table(paste(maindir, "Cenchrus_filt1_loci.txt", sep="/"), header=T)
head(orig_loci)
head(pca_loc)

table(orig_loci$locus == pca_loc$locus)

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
plot(1:40,summary(pc)$importance[2,][1:40],xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
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
plot(1:nrow(lfres),-log(lfres$long_p,10),type="n",ylab="-log10(p value)",las=1,xlab="locus",ylim= c(0,max(-log(lfres[,2:4],10))))

points(1:nrow(lfres),-log(lfres$long_p,10),pch=20,col=as.factor(lfres$lo_q<0.05))

points(rownames(lfres)[lfres$bu_q<alpha],-log(lfres$burn_unburnt_p[lfres$bu_q<alpha],10),pch=20,col="green")
points(rownames(lfres)[lfres$lo_q<alpha],-log(lfres$long_p[lfres$lo_q<alpha],10),pch=20,col="blue")
legend("topleft",legend=c("PC1 burnt or unburnt","PC2 longtitude"),pch=20,col=c("green","blue"))

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





