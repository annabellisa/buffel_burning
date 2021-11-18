
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### Diversity & distances
### Author: Annabel Smith & Binyin Di

# AS Load workspace:
load("03_workspaces/STEP04_divdist_ALL.RData")

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# Load libraries:
library(diveRsity);library(hierfstat);library(adegenet); library(ecodist); library(AICcmodavg);library(geosphere)

#  GENIND object & site data:    	# ----

# Load Genepop files
gp_dir<-"00_Data/Genepop_Files/Genepop_gen_files"
dir(gp_dir)

# Make genind objects:
# Updated Nov 2021 after peer review

# ~~ Neutral 
genind_neutral<-read.genepop(file=paste(gp_dir,"Genepop_filt4.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_neutral

# ~~ Non-Neutral (~ 10 seconds for 5030 loci)
genind_nonneutral<-read.genepop(file=paste(gp_dir,"Genepop_filt5_nonneutral.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_nonneutral

# Load site data:
sdat<-read.table(paste("00_data/Cenchrus_site_data.txt",sep=""),header=T)
sdat<-sdat[sdat$sequenced==1,]
sdat<-tidy.df(sdat)
head(sdat,2); dim(sdat)

# Format to match genetic data:
sdt<-sdat[,c("block","site","year","burn_unburnt","time_since_burn","no_samples","lat","long")]
colnames(sdt)<-c("block","site","year","burn","TSF","no_samples","lat","long")
sdt$burn<-factor(sdt$burn, levels=c("u","b"))
sdt$burn2<-as.character(sdt$burn)
sdt$burn2[which(sdt$block=="buf11")]<-"b2"
sdt$burn2<-factor(sdt$burn2, levels=c("u","b","b2"))
head(sdt, 3); dim(sdt)

# Load cluster assignment data:
kdat<-read.table("00_data/K_genetic_clusters_Cenchrus_filt4.txt", header=T)
kd<-kdat[,c("indiv", "K3", "K4")]
head(kd,3); dim(kd)

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# close genind object ----

# RESULT:
# See parameter files in gp_dir for filters
genind_neutral # all filters + neutral markers only (15965 loci)
genind_nonneutral # all filters + non-neutral markers only  (5030 loci)
head(sdat,3); dim(sdat) # site data

### -- *** 
# Question 1: Does fire affect spatial genetic structure?

#  HIERARCHICAL CLUSTERING:    	# ----

# see code in Supplement_02_plot_structure.R where this is re-run and plotted over structure results (not run here)

# hclust (in base R)
# Consensus UPGMA dendrogram (see Acquadro et al. 2017)
# https://popgen.nescent.org/2015-05-18-Dist-SNP.html
# https://adegenet.r-forge.r-project.org/files/Glasgow2015/practical-introphylo.1.0.pdf
# https://dyerlab.github.io/applied_population_genetics/genetic-distances.html

# re-do distance matrix on raw data:
ddir<-"00_Data/Filtered_DartSeq_format"
ddat<-read.table(paste(ddir, "dartseq_filt4.txt", sep="/"), header=T)
ghead(ddat)

# data
clust_dat<-ddat
rownames(clust_dat)<-clust_dat$ind
clust_dat<-clust_dat[,3:length(clust_dat)]
ghead(clust_dat)

# distance matrix and hierarchical clustering:
hc_dist<-dist(x = clust_dat, method="euclidean")
euc_clust<-hclust(hc_dist)
str(euc_clust)

# See Supplement_02_plot_structure.R for plotting the hclust object

# close hclust ----

# FST & Mantel tests (SITE LEVEL ANALYSIS):    	# ----

### -- *** CALCULATE FST:

# Load Genepop files
gp_dir<-"00_Data/Genepop_Files/Genepop_gen_files"
dir(gp_dir)

# Get FST by site ALL individuals ALL K:
print(Sys.time())
fst<-diffCalc(paste(gp_dir,"Genepop_filt4.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())
# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# Get FST for each K separately (K1 and K3 only):
print(Sys.time())
fst_K1<-diffCalc(paste(gp_dir,"Genepop_filt4_K1.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())

print(Sys.time())
fst_K3<-diffCalc(paste(gp_dir,"Genepop_filt4_K3.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# Get FST by K3 clusters, all individuals:
print(Sys.time())
fstK3clust<-diffCalc(paste(gp_dir,"Genepop_filt4_byK3.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())
# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

head(fst$pairwise$Fst)
head(fst_K1$pairwise$Fst)
head(fst_K3$pairwise$Fst)

# FST between clusters:
head(fstK3clust$pairwise$Fst)

fst_clust_name<-rownames(fstK3clust$pairwise$Fst)
fst_clust_name<-substr(fst_clust_name,1,nchar(fst_clust_name)-1)

# Get K data for interpretation:
dir("00_Data")
kclust<-read.table("00_Data/K_genetic_clusters_Cenchrus_filt4.txt", header=T)
kclust<-kclust[,c("indiv", "K3")]
head(kclust)

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

kclust_3<-kclust[which(kclust$indiv %in% fst_clust_name),]
kclust_3<-tidy.df(kclust_3)
kclust_3

fstK3clust2<-fstK3clust$pairwise$Fst
rownames(fstK3clust2)<-substr(rownames(fstK3clust2),1,nchar(rownames(fstK3clust2))-1)
colnames(fstK3clust2)<-substr(colnames(fstK3clust2),1,nchar(colnames(fstK3clust2))-1)

rownames(fstK3clust2)<-paste("K",kclust_3$K3[unlist(lapply(rownames(fstK3clust2) , FUN=function(x) which(kclust_3$indiv %in% x)))],sep="")
colnames(fstK3clust2)<-paste("K",kclust_3$K3[unlist(lapply(colnames(fstK3clust2) , FUN=function(x) which(kclust_3$indiv %in% x)))],sep="")
fstnK3clust<-rownames(fstK3clust2)

fst_K3clust_df<-data.frame(pop1=combn(fstnK3clust,2)[1,],pop2=combn(fstnK3clust,2)[2,],fst=fstK3clust2[lower.tri(fstK3clust2)])

# COLOURS and K name UPDATE Nov 2021
# K1 == Site 1, 2, 6u, 7u, 8u, 10 == PURPLE (#BEAED4)
# K2 == Site 7b, 8b == ORANGE (#FDC086)
# K3 == Site 3, 5, 6b, 11b1 == GREEN (#7FC97F)

# c("#7FC97F","#BEAED4","#FDC086")
# this is green purple orange

# they're plotted in this order on the structure plot: K3 (green), K1 (purple), K2 (orange)
# confusion arises because in the K clusters file that we're using for all K subset analyses, they're listed in order from site 1 onwards. On the structure plot, they're ordered from west to east, so site 11 comes up first (thus K3, green gets plotted first)
# but the colours and Ks as listed here are correct
fst_K3clust_df
dev.new(8,8,pointsize=16, dpi=80,noRStudioGD = T)
plot(1:3, 1:3, type="p", pch=20, cex=4,col=c("#7FC97F","#BEAED4","#FDC086"))

# Convert square matrix to column matrix:
fstn<-rownames(fst$pairwise$Fst)
fstn<-substr(fstn,1,(nchar(fstn)-2))
fstn1<-rownames(fst_K1$pairwise$Fst)
fstn1<-substr(fstn1,1,(nchar(fstn1)-2))
fstn3<-rownames(fst_K3$pairwise$Fst)
fstn3<-substr(fstn3,1,(nchar(fstn3)-2))

# All sites
ind_endnum<-which(unlist(gregexpr("[0-9]",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endnum<-fstn[ind_endnum]
fstn[ind_endnum]<-substr(with_endnum,1,nchar(with_endnum)-1)

ind_endund<-which(unlist(gregexpr("_",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endund<-fstn[ind_endund]
fstn[ind_endund]<-substr(with_endund,1,nchar(with_endund)-1)

# K1
ind_endnum1<-which(unlist(gregexpr("[0-9]",substr(fstn1,nchar(fstn1),nchar(fstn1))))>0)
with_endnum1<-fstn1[ind_endnum1]
fstn1[ind_endnum1]<-substr(with_endnum1,1,nchar(with_endnum1)-1)

ind_endund1<-which(unlist(gregexpr("_",substr(fstn1,nchar(fstn1),nchar(fstn1))))>0)
with_endund1<-fstn1[ind_endund1]
fstn1[ind_endund1]<-substr(with_endund1,1,nchar(with_endund1)-1)

# K3
ind_endnum3<-which(unlist(gregexpr("[0-9]",substr(fstn3,nchar(fstn3),nchar(fstn3))))>0)
with_endnum3<-fstn3[ind_endnum3]
fstn3[ind_endnum3]<-substr(with_endnum3,1,nchar(with_endnum3)-1)

ind_endund3<-which(unlist(gregexpr("_",substr(fstn3,nchar(fstn3),nchar(fstn3))))>0)
with_endund3<-fstn3[ind_endund3]
fstn3[ind_endund3]<-substr(with_endund3,1,nchar(with_endund3)-1)

# Combine into data frames

# All sites:
fst_df<-data.frame(pop1=combn(fstn,2)[1,],pop2=combn(fstn,2)[2,],fst=fst$pairwise$Fst[lower.tri(fst$pairwise$Fst)],gst=fst$pairwise$gst[lower.tri(fst$pairwise$gst)],Gst=fst$pairwise$Gst[lower.tri(fst$pairwise$Gst)],GGst=fst$pairwise$GGst[lower.tri(fst$pairwise$GGst)],D=fst$pairwise$D[lower.tri(fst$pairwise$D)])
head(fst_df)

# K1:
fst_df1<-data.frame(pop1=combn(fstn1,2)[1,],pop2=combn(fstn1,2)[2,],fst=fst_K1$pairwise$Fst[lower.tri(fst_K1$pairwise$Fst)],gst=fst_K1$pairwise$gst[lower.tri(fst_K1$pairwise$gst)],Gst=fst_K1$pairwise$Gst[lower.tri(fst_K1$pairwise$Gst)],GGst=fst_K1$pairwise$GGst[lower.tri(fst_K1$pairwise$GGst)],D=fst_K1$pairwise$D[lower.tri(fst_K1$pairwise$D)])
head(fst_df1)

# K3:
fst_df3<-data.frame(pop1=combn(fstn3,2)[1,],pop2=combn(fstn3,2)[2,],fst=fst_K3$pairwise$Fst[lower.tri(fst_K3$pairwise$Fst)],gst=fst_K3$pairwise$gst[lower.tri(fst_K3$pairwise$gst)],Gst=fst_K3$pairwise$Gst[lower.tri(fst_K3$pairwise$Gst)],GGst=fst_K3$pairwise$GGst[lower.tri(fst_K3$pairwise$GGst)],D=fst_K3$pairwise$D[lower.tri(fst_K3$pairwise$D)])
head(fst_df3)

mean(fst_df$fst)
range(fst_df$fst,na.rm=T)
mean(fst_df1$fst)
range(fst_df1$fst,na.rm=T)
mean(fst_df3$fst)
range(fst_df3$fst,na.rm=T)

# write.table(fst_df3,"fst_K3.txt",row.names=F,quote=F,sep="\t")

### -- *** ADD GEOGRAPHIC DISTANCE:

# The following formatting of FST (adding geog dist, etc) has been run and saved in fst_and_distances_all_sites.txt (until "Analyse FST")

dat_dir<-"04_RESULTS/Diversity_and_Distance/FST"
dir(dat_dir)

# load fst data:
pwpop<-read.table(paste(dat_dir,"fst_all_sites.txt",sep="/"),header=T)
pwpop1<-read.table(paste(dat_dir,"fst_K1.txt",sep="/"),header=T)
pwpop3<-read.table(paste(dat_dir,"fst_K3.txt",sep="/"),header=T)
head(pwpop)
head(pwpop1)
head(pwpop3)

# load site data:

# add site code to match fst data:
sdat$pop<-sdat$site
sdat$pop<-paste("X",substr(x=sdat$pop,start = 4,stop = nchar(as.character(sdat$pop))), sep="")
head(sdat,3)

# Check all sites in data have site data (should all be T):
table(unique(c(pwpop$pop1,pwpop$pop2)) %in% sdat$pop)
table(unique(c(pwpop1$pop1,pwpop1$pop2)) %in% sdat$pop)
table(unique(c(pwpop3$pop1,pwpop3$pop2)) %in% sdat$pop)

sll<-sdat[,c("pop","lat","long")]
head(sll); dim(sll)
head(pwpop,2); dim(pwpop)

# All SITES:
m1<-merge(pwpop,sll,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("lat","long")]<-c("lat1","lon1")
m1<-merge(m1,sll,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("lat","long")]<-c("lat2","lon2")
m1<-m1[order(m1$pop1,m1$pop2),]
m1<-tidy.df(m1)
m1<-m1[,c(2,1,3:length(m1))]
m1$geog_dist<-distGeo(m1[,c("lon1","lat1")],m1[,c("lon2","lat2")])
# lat dist is redundant here. For Cenchrus, long dist would make more sense, but that will be captured in pure geographic distance:
m1$lat_dist<-abs(m1$lat1)-abs(m1$lat2)
head(m1,2)
head(sdat,4)

# Add same site distance:
ndis<-sdat[,c("pop","block")]
ndis$block<-paste("X",substr(x=sdat$block,start = 4,stop = nchar(as.character(sdat$block))), sep="")
head(ndis)
m2<-merge(m1,ndis,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("block")]<-c("block1")
m2<-merge(m2,ndis,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("block")]<-c("block2")
m2<-m2[order(m2$pop1,m2$pop2),]
m2<-tidy.df(m2)
m2<-m2[,c(2,1,3:length(m2))]
head(m2,3)

m2$same_block<-m2$block1==m2$block2
m2$same_block<-ifelse(m2$same_block==T,1,0)
check.rows(m2)

# K1:
head(pwpop1); dim(pwpop1)
m1K1<-merge(pwpop1,sll,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m1K1)[colnames(m1K1) %in% c("lat","long")]<-c("lat1","lon1")
m1K1<-merge(m1K1,sll,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m1K1)[colnames(m1K1) %in% c("lat","long")]<-c("lat2","lon2")
m1K1<-m1K1[order(m1K1$pop1,m1K1$pop2),]
m1K1<-tidy.df(m1K1)
head(m1K1)

m1K1<-m1K1[,c(2,1,3:length(m1K1))]
m1K1$geog_dist<-distGeo(m1K1[,c("lon1","lat1")],m1K1[,c("lon2","lat2")])
# lat dist is redundant here. For Cenchrus, long dist would make more sense, but that will be captured in pure geographic distance:
m1K1$lat_dist<-abs(m1K1$lat1)-abs(m1K1$lat2)
head(m1K1,2); dim(m1K1)
head(sdat,4)

# Add same site distance:
head(ndis); dim(ndis)
m2K1<-merge(m1K1,ndis,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m2K1)[colnames(m2K1) %in% c("block")]<-c("block1")
m2K1<-merge(m2K1,ndis,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m2K1)[colnames(m2K1) %in% c("block")]<-c("block2")
m2K1<-m2K1[order(m2K1$pop1,m2K1$pop2),]
m2K1<-tidy.df(m2K1)
m2K1<-m2K1[,c(2,1,3:length(m2K1))]
head(m2K1,3)

m2K1$same_block<-m2K1$block1==m2K1$block2
m2K1$same_block<-ifelse(m2K1$same_block==T,1,0)
head(m2K1,3); dim(m2K1)

# K3:
head(pwpop3); dim(pwpop3)
m1K3<-merge(pwpop3,sll,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m1K3)[colnames(m1K3) %in% c("lat","long")]<-c("lat1","lon1")
m1K3<-merge(m1K3,sll,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m1K3)[colnames(m1K3) %in% c("lat","long")]<-c("lat2","lon2")
m1K3<-m1K3[order(m1K3$pop1,m1K3$pop2),]
m1K3<-tidy.df(m1K3)
head(m1K3); dim(m1K3)

m1K3<-m1K3[,c(2,1,3:length(m1K3))]
m1K3$geog_dist<-distGeo(m1K3[,c("lon1","lat1")],m1K3[,c("lon2","lat2")])
# lat dist is redundant here. For Cenchrus, long dist would make more sense, but that will be captured in pure geographic distance:
m1K3$lat_dist<-abs(m1K3$lat1)-abs(m1K3$lat2)
head(m1K3,2); dim(m1K3)
head(sdat,4)

# Add same site distance:
head(ndis); dim(ndis)
m2K3<-merge(m1K3,ndis,by.x="pop1",by.y="pop",all.x=T,all.y=F)
colnames(m2K3)[colnames(m2K3) %in% c("block")]<-c("block1")
m2K3<-merge(m2K3,ndis,by.x="pop2",by.y="pop",all.x=T,all.y=F)
colnames(m2K3)[colnames(m2K3) %in% c("block")]<-c("block2")
m2K3<-m2K3[order(m2K3$pop1,m2K3$pop2),]
m2K3<-tidy.df(m2K3)
m2K3<-m2K3[,c(2,1,3:length(m2K3))]
head(m2K3,3)

m2K3$same_block<-m2K3$block1==m2K3$block2
m2K3$same_block<-ifelse(m2K3$same_block==T,1,0)
head(m2K3,3); dim(m2K3)

head(m2,3); dim(m2)
head(m2K1,3); dim(m2K1)
head(m2K3,3); dim(m2K3)

# write.table(m2K3,"m2K3.txt",row.names=F,quote=F,sep="\t")
# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# Analyse FST:

dir(dat_dir)
# load fst data:
fst_all<-read.table(paste(dat_dir,"fst_and_distances_all_sites.txt",sep="/"),header=T)
fst_K1<-read.table(paste(dat_dir,"fst_and_distances_K1.txt",sep="/"),header=T)
fst_K3<-read.table(paste(dat_dir,"fst_and_distances_K3.txt",sep="/"),header=T)
head(fst_all,2); dim(fst_all)
head(fst_K1,2); dim(fst_K1)
head(fst_K3,2); dim(fst_K3)

summary(fst_all$fst)
summary(fst_K1$fst)
summary(fst_K3$fst)

# sites with the very large FSTs correspond to the structure clusters... they're probably different lines
fst_all[,1:3]

# mantel test:
mant_all<-mantel(formula = fst~geog_dist, data = fst_all)
mant_all
mant_K1<-mantel(formula = fst~geog_dist, data = fst_K1)
mant_K1
mant_K3<-mantel(formula = fst~geog_dist, data = fst_K3)
mant_K3

# plot FST:
# plot FST:
dev.new(6,8,dpi=90, pointsize=16,noRStudioGD = T)
par(mfrow=c(3,1),mar=c(4,4,0.5,0.5), mgp=c(2.7,1,0), oma=c(0,0,1.5,0))
plot(fst_all$geog_dist/1000, fst_all$fst, pch=20, xlab="", xlim=c(0,100),ylab=expression("Genetic distance ("*italic("F")[ST]*")"), las=1)
title(xlab="Geographic distance (km)",mgp=c(2.2,1,0))
mtext(text="(a) All sites",side=3, at=-10, adj=0,line=0.5, cex=0.7)

plot(fst_K1$geog_dist/1000, fst_K1$fst, pch=20, xlab="", xlim=c(0,100),ylab=expression("Genetic distance ("*italic("F")[ST]*")"),ylim=c(-0.2,1), las=1, col="#BEAED4")
title(xlab="Geographic distance (km)",mgp=c(2.2,1,0))
mtext(text="(b) Genetic cluster K1",side=3, at=-10,adj=0, line=0.5, cex=0.7)

plot(fst_K3$geog_dist/1000, fst_K3$fst, pch=20, xlab="", xlim=c(0,100),ylab=expression("Genetic distance ("*italic("F")[ST]*")"),ylim=c(-0.2,1), las=1, col="#7FC97F")
title(xlab="Geographic distance (km)",mgp=c(2.2,1,0))
mtext(text="(c) Genetic cluster K3",side=3,at=-10,adj=0, line=0.5, cex=0.7)

# close FST ----

# PCA:    	# ----

  # PCA adegenet:
  K_test<-22 # set big K
  X_test <- scaleGen(genind_neutral, scale = FALSE, NA.method = "mean")
  pcaX_test <- dudi.pca(X_test, cen = FALSE, scale = FALSE, scannf = FALSE, nf = K_test)
  summary(pcaX_test)
  
  # Calculate variance explained
  
  # This tutorial is from one of the authors of ade4:
  # https://pbil.univ-lyon1.fr/R/pdf/course2.pdf
  var_expl_test <- pcaX_test$eig/sum(pcaX_test$eig)
  head(cumsum(var_expl_test),20)
  
  ### --- *** CHOOSE K *** --- ###
  
  # CHOOSE K FROM SCREE PLOT (they are very consistent with PCAs from outlier tests):
  
  # Of eigen values:
  dev.new(height=4,width=4,dpi=160,noRStudioGD = T,pointsize=12)
  par(mar=c(4,4,0.5,0.5))
  plot(1:K_test,pcaX_test$eig[1:22],xlab="",ylab="Eigenvalue",pch=20,las=1,type="n")
  title(xlab="PC",mgp=c(2.5,1,0))
  grid()
  lines(1:K_test,pcaX_test$eig[1:22])
  points(1:K_test,pcaX_test$eig[1:22],pch=20)
  
  # Or variance explained:
  dev.new(height=4,width=4,dpi=160,noRStudioGD = T,pointsize=12)
  par(mar=c(4,4,0.5,0.5))
  plot(1:K_test,var_expl_test[1:22],xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
  title(xlab="PC",mgp=c(2.5,1,0))
  grid()
  lines(1:K_test,var_expl_test[1:22])
  points(1:K_test,var_expl_test[1:22],pch=20)
  
  # K=3
  
  ### --- *** SET K & re-run *** --- ###
  
  K<-3
  
  X <- scaleGen(genind_neutral, scale = FALSE, NA.method = "mean")
  pcaX <- dudi.pca(X, cen = FALSE, scale = FALSE, scannf = FALSE, nf = K)
  summary(pcaX)
  
  # Calculate variance explained
  var_expl <- pcaX$eig/sum(pcaX$eig)
  head(cumsum(var_expl),10)
  
  ### *** plot PCA
  
  head(pcaX$li)
  rownames(pcaX$li)
  head(sdat,2)
  
  # The genind object stores both individual names from the genpop file and population names, which it appears to take from the last individual in each population (this could be fixed by updating the pop names in the genepop format script but, for now, I'm just dealing with it here)
  
  ldat<-data.frame(rowpca=rownames(pcaX$li),genind_pop=pop(genind_neutral))
  ldat$site_code<-sapply(ldat$genind_pop,function(x) substr(x,1,gregexpr("_",x)[[1]][1]-1))
  ldat$site_code<-as.character(ldat$site_code)
  ldat$ind<-1:nrow(ldat)
  ldat<-tidy.df(ldat)
  
  sdt2<-sdt[,c("site","year","burn","burn2","lat","long")]

  head(ldat,3); dim(ldat)
  head(sdt2,2); dim(sdt2)
  
  # Align site names in two data sets:
  sdt2$site<-paste("X",substr(sdt2$site, 4, nchar(as.character(sdt2$site))),sep="")
  
  # should all be TRUE:
  table(ldat$site_code %in% sdt2$site); table(sdt2$site %in% ldat$site_code )
  
  # Add country code, re-merge and re-order:
  ldat<-merge(ldat,sdt2,by.x="site_code",by.y="site",all.x=T,all.y=F)
  ldat<-ldat[order(ldat$ind),]
  ldat<-tidy.df(ldat)
  
  # Should all be TRUE:
  table(rownames(pcaX$li) == ldat$rowpca)
  
  # Add PCs to site data:
  ldat<-cbind(ldat, pcaX$li)
  head(ldat,3); dim(ldat)
  
  dev.new(height=6, width=6, noRStudioGD = T, dpi=100, pointsize=12)
  par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5,1,0), oma=c(0,0,1.5,0))

    plot(1:K_test,var_expl_test[1:22],xlab="",ylab="Proportion variance explained",pch=20,las=1,type="n")
  title(xlab="Principal Component",mgp=c(2.3,1,0))
  grid()
  lines(1:K_test,var_expl_test[1:22])
  points(1:K_test,var_expl_test[1:22],pch=20)
  mtext("(a)", side = 3, line=0.5, cex = 1, adj=0)
  
  blankplot()
  
  par(xpd=NA)
  legend(x=5.4, y=0.8, legend = c("road unburnt","road burnt","burnt site 11"), pt.cex=1.5,col=c("black","red","orange"),pch=20)
  par(xpd=F)
  
  # colour order is 1,2,3==black,red,green,
  # burnt(2)==green,unburnt(1)==black,burn2(3)==green,
  # need it to be black, red, orange
  
  c.now<-data.frame(burn2=ldat$burn2,burn_order=as.numeric(ldat$burn2))
  c.now$new.col<-c.now$burn_order
  c.now$new.col<-ifelse(ifelse(ifelse(c.now$new.col==1,"black",c.now$new.col)==2,"red",ifelse(c.now$new.col==1,"black",c.now$new.col))==3,"orange",ifelse(ifelse(c.now$new.col==1,"black",c.now$new.col)==2,"red",ifelse(c.now$new.col==1,"black",c.now$new.col)))
  head(c.now)
  
  plot(ldat$Axis1, ldat$Axis2, col=c.now$new.col, pch=20, xlab="PC 1", ylab="PC 2",cex=1.5, las=1, mgp=c(2.3,1,0))
  mtext("(b)", side = 3, line=0.5, cex = 1, adj=0)
  
  plot(ldat$Axis1, ldat$Axis3, col=c.now$new.col, pch=20, xlab="PC 1", ylab="PC 3" ,cex=1.5, las=1, mgp=c(2.3,1,0))
  mtext("(c)", side = 3, line=0.5, cex = 1, adj=0)
  
  # save.image("03_Workspaces/STEP04_divdist_ALL.RData")
  
# close PCA ----

# Question 2: Does fire affect genetic diversity?
  
 #  CALCULATE SITE-LEVEL Genetic diversity: 	# ----
  
  # The following section has been run and saved in "Genetic_Diversity_ALL.txt"
  
  load("03_workspaces/STEP04_divdist_ALL.RData")
  
  genind_neutral # all filters + neutral markers only (15965 loci)
  genind_nonneutral # all filters + non-neutral markers only  (5030 loci)
  head(sdat,3); dim(sdat)
  
  # Calculate genetic diversity per population in hierfstat:
  
  # neutral:
  gendiv_neutral <- basic.stats(genind_neutral, diploid = TRUE, digits = 2)
  str(gendiv_neutral) # 1 min
  head(gendiv_neutral$Fis); dim(gendiv_neutral$Fis)
  head(gendiv_neutral$Ho)
  tail(gendiv_neutral$Ho)
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  gd_neutral<-data.frame(site=names(apply(gendiv_neutral$Ho,2,mean,na.rm=T)),max_n=apply(gendiv_neutral$n.ind.samp,2,max,na.rm=T),Ho=apply(gendiv_neutral$Ho,2,mean,na.rm=T),He=apply(gendiv_neutral$Hs,2,mean,na.rm=T),Fis=apply(gendiv_neutral$Fis,2,mean,na.rm=T))
  gd_neutral<-tidy.df(gd_neutral)
  
  gd_neutral$site<-substr(gd_neutral$site,1,nchar(as.character(gd_neutral$site))-3)
  head(gd_neutral); dim(gd_neutral)
  
  # plot(gd_neutral$Ho, gd_neutral$He)
  
  # non-neutral
  gendiv_nonneutral <- basic.stats(genind_nonneutral, diploid = TRUE, digits = 2)
  
  str(gendiv_nonneutral) # a few seconds
  head(gendiv_nonneutral$Ho)
  tail(gendiv_nonneutral$Ho)
  
  gd_nonneutral<-data.frame(site=names(apply(gendiv_nonneutral$Ho,2,mean,na.rm=T)),max_n=apply(gendiv_nonneutral$n.ind.samp,2,max,na.rm=T),Ho=apply(gendiv_nonneutral$Ho,2,mean,na.rm=T),He=apply(gendiv_nonneutral$Hs,2,mean,na.rm=T),Fis=apply(gendiv_nonneutral$Fis,2,mean,na.rm=T))
  gd_nonneutral<-tidy.df(gd_nonneutral)
  head(gd_nonneutral,3); dim(gd_nonneutral)
  
  gd_nonneutral$site<-substr(gd_nonneutral$site,1,nchar(as.character(gd_nonneutral$site))-3)
  head(gd_nonneutral)
  
  # save.image("03_workspaces/STEP_04_kinship.RData")
  
  ## -- ** ALLELIC RICHNESS:
  
  head(sdat,3); dim(sdat)
  genind_neutral
  genind_nonneutral
  
  print(Sys.time())
  ar_neutral<- allelic.richness(genind_neutral, diploid = TRUE) # 1 min
  print(Sys.time())
  
  print(Sys.time())
  ar_nonneutral<- allelic.richness(genind_nonneutral, diploid = TRUE) # a few seconds
  print(Sys.time())
  
  str(ar_neutral)
  head(ar_neutral$Ar,2)
  
  str(ar_nonneutral)
  head(ar_nonneutral$Ar,2)
  
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  # Summarise
  
  # Neutral:
  head(ar_neutral$Ar,2)
  ar_neutral$min.all
  ar_neutral_res<-data.frame(site=levels(genind_neutral@pop),ar_neutral=apply(ar_neutral$Ar,2,mean,na.rm=T))
  head(ar_neutral_res); dim(ar_neutral_res)
  
  # Non-Neutral:
  head(ar_nonneutral$Ar,2)
  ar_nonneutral$min.all
  ar_nonneutral_res<-data.frame(site=levels(genind_nonneutral@pop),ar_nonneutral=apply(ar_nonneutral$Ar,2,mean,na.rm=T))
  head(ar_nonneutral_res); dim(ar_nonneutral_res)
  
  # combine allelic richness into table with other genetic diversity data:
  
  # fix up the names on AR data frames:
  ar_neutral_res$site<-substr(ar_neutral_res$site,1,nchar(as.character(ar_neutral_res$site))-3)
  head(ar_neutral_res)
  ar_nonneutral_res$site<-substr(ar_nonneutral_res$site,1,nchar(as.character(ar_nonneutral_res$site))-3)
  head(ar_nonneutral_res)
  
  # check (all four should be all T):
  gd_neutral$site %in% ar_neutral_res$site
  ar_neutral_res$site %in% gd_neutral$site
  
  gd_nonneutral$site %in% ar_nonneutral_res$site
  ar_nonneutral_res$site %in% gd_nonneutral$site
  
  # merge
  head(gd_neutral,3); dim(gd_neutral)
  head(ar_neutral_res,3); dim(ar_neutral_res)
  gd_neutral<-merge(gd_neutral, ar_neutral_res, by="site", all.x=T, all.y=F)
  
  head(gd_nonneutral,3); dim(gd_nonneutral)
  head(ar_nonneutral_res,3); dim(ar_nonneutral_res)
  gd_nonneutral<-merge(gd_nonneutral, ar_nonneutral_res, by="site", all.x=T, all.y=F)
  
  # combine neutral and non-neutral:
  
  colnames(gd_neutral)[colnames(gd_neutral) %in% c("Ho","He","Fis")]<-paste(colnames(gd_neutral)[colnames(gd_neutral) %in% c("Ho","He","Fis")],"_neutral",sep="")
  colnames(gd_nonneutral)[colnames(gd_nonneutral) %in% c("Ho","He","Fis")]<-paste(colnames(gd_nonneutral)[colnames(gd_nonneutral) %in% c("Ho","He","Fis")],"_nonneutral",sep="")
  
  head(gd_neutral,3); dim(gd_neutral)
  head(gd_nonneutral,3); dim(gd_nonneutral)
  
  gd_all<-gd_neutral
  gd_all$max_n<-NULL
  gd_all<-merge(gd_all,gd_nonneutral, by="site", all.x=T, all.y=F)
  gd_all<-gd_all[,c(c(which(colnames(gd_all) %in% c("site","max_n"))),c(which(!colnames(gd_all) %in% c("site","max_n"))))]
  head(gd_all,3); dim(gd_all)
  
  # write.table(gd_all, "gd_all.txt", row.names=F, quote=F, sep="\t")
  # save.image("03_workspaces/STEP_04_kinship.RData")
  
  # Format data for analysis:
  
  # Site-level data:
  gd_all<-read.table("RESULTS/Diversity_and_distance/Genetic_Diversity/Genetic_Diversity_ALL.txt", header=T)
  head(gd_all,3); dim(gd_all)
  
  # simplify site data
  head(sdt4)
  
  # merge with genetic diversity data:
  table(sdt4$site %in% gd_all$site)
  table(gd_all$site %in% sdt4$site)
  gd_all<-merge(gd_all, sdt4, by="site", all.x=T, all.y=F)
  
  # add second treatment variable:
  gd_all$treatment<-gd_all$burn
  gd_all$trt<-as.character(gd_all$treatment)
  gd_all$trt[grep("b",gd_all$trt)]<-"b1"
  gd_all$trt[grep("X11",gd_all$site)]<-"b2"
  gd_all$trt<-factor(gd_all$trt, levels=c("u","b1","b2"))
  head(gd_all,3); dim(gd_all)
  
  #### Accounting for population genetic structure
  
  # Given the strong background genetic structure in the data, we must control for phylogeny. 
  # We modelled background genetic structure using two methods: the mean site-level diversity score, calculated from individual level diversity scores and the most common K per site
  # Regardless of the method used, there was no effect of background genetic structure on site-level allelic richness. Both methods showed the treatment only model to be superior, although it was not better than a null model. 
  # For the update during the revision, I have removed the diversity score, just focussing on the K method. 
  
  # Assignment probabilities from STRUCTURE:
  site_assig<-read.table("00_Data/assig_prob_K3_Cenchrus_filt4.txt", header=T)
  # site_assig<-cbind(site_assig, matrix(data=c(1,2,3),ncol=3,nrow=nrow(site_assig),byrow=T))
  # site_assig$genstr<-apply(site_assig[,3:8],1,function(x) weighted.mean(x[4:6], x[1:3]))
  head(site_assig); dim(site_assig)
  
  ## -- ** METHOD 2: K
  
  # ih_dat has K per individual; it was merged with individual data in the Individual genetic diversity section below; they are the correct Ks, as in, they reflect the K with the highest assignment probability per individual; check site 7b which has a mixture of K=1 and K=3:
  head(ih_dat,3); dim(ih_dat)
  head(gd_all,2); dim(gd_all)
  
  # Use which.max to get the mode - the most common K for each site:
  Ksite<-data.frame(site=names(tapply(ih_dat$K3,ih_dat$site, FUN = function(x) names(which.max(summary(x))))),most_likelyK=tapply(ih_dat$K3,ih_dat$site, FUN = function(x) names(which.max(summary(x)))))
  Ksite<-tidy.df(Ksite)
  head(Ksite); dim(Ksite)
  head(gd_all)
  
  gd_all<-merge(gd_all, Ksite, by="site", all.x=T, all.y=F)
  # Check it against the structure plot
  check.rows(gd_all[,c("site","most_likelyK")])
  
  # close calc genetic diversity ----
  
#  CALCULATE individual genetic diversity	# ----
  
  # Calculate ind het on dart seq format:
  ddir<-"00_Data/Filtered_DartSeq_format"
  dir(ddir)
  ddat<-read.table(paste(ddir, "dartseq_filt4_neutral.txt", sep="/"), header=T)

  # Neutral data set:
  ghead(ddat); dim(ddat)
  
  # Non-neutral data set:
  nndat<-read.table(paste(ddir,"dartseq_filt4_non_neutral.txt", sep="/"), header=T)
  ghead(nndat); dim(nndat)
  
  # Codes for onerow formatted data:
  # 0 = Reference allele homozygote
  # 1 = SNP allele homozygote
  # 2 = heterozygote
  
  # Neutral
  # 3 MINS on laptop 
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  ih_out<-data.frame(ddat[,1:2],ind_het=NA)
  head(ih_out,25)
  
  for(i in 1:nrow(ih_out)){
    ind.cons<-as.character(ddat[i,3:length(ddat)])
    t.cons<-data.frame(ind.cons=as.numeric(names(table(ind.cons))),count=as.numeric(table(ind.cons)),stringsAsFactors=F)
    
    if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] else stop("no zero category")
    
    if(length(which(t.cons$ind.cons==2))==0) ih_out[i,3]<-"no_heterozygotes" else ih_out[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)
    
  } # close for i
  
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  # Non-neutral
  
  ih_out_nn<-data.frame(nndat[,1:2],ind_het=NA)
  head(ih_out_nn,15); dim(ih_out_nn)
  ghead(nndat); dim(nndat)
  
  for(i in 1:nrow(ih_out_nn)){
    ind.cons<-as.character(nndat[i,3:length(nndat)])
    t.cons<-data.frame(ind.cons=as.numeric(names(table(ind.cons))),count=as.numeric(table(ind.cons)),stringsAsFactors=F)
    
    if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] else stop("no zero category")
    
    if(length(which(t.cons$ind.cons==2))==0) ih_out_nn[i,3]<-"no_heterozygotes" else ih_out_nn[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)
    
  } # close for i
  
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  ih_out_nn$ind_het<-as.numeric(ih_out_nn$ind_het)
  
  ghead(ddat); dim(ddat)
  ghead(nndat); dim(nndat)
  
  head(ih_out); dim(ih_out)
  head(ih_out_nn); dim(ih_out_nn)
  
  # add site data and cluster assignment:
  # Format to match ind het data:
  sdt4<-sdt
  sdt4$site<-paste("X",substr(sdt4$site,4,nchar(sdt4$site)),sep="")
  head(sdt4, 3); dim(sdt4)
  
  # Check: should all be T
  table(sdt4$site %in% ih_out$site)
  table(ih_out$site %in% sdt4$site)
  
  kd$site<-substr(kd$indiv,1,(unlist(gregexpr("_",kd$indiv))-1))
  head(kd,10); dim(kd)
  
  # Check: should all be T
  table(kd$site %in% ih_out$site)
  table(ih_out$site %in% kd$site)
  
  head(ih_out,6); dim(ih_out)
  head(ih_out_nn,6); dim(ih_out_nn)
  
  # Before merging, check they're in the same order:
  table(ih_out$ind==ih_out_nn$ind)
  
  ih_dat<-data.frame(ih_out,ind_het_nonneut=ih_out_nn$ind_het_nonneut)
  
  # Check all sites are present in new df:
  table(ih_dat$site %in% sdt4$site)
  
  # Merge with site data:
  ih_dat<-merge(ih_dat, sdt4, by="site", all.x=T, all.y=F)
  head(ih_dat,3); dim(ih_dat)
  
  # merge K on individual:
  # kd has the correct individual level K values:
   head(kd,3); dim(kd)
  table(kd$indiv %in% ih_dat$ind)
  table(ih_dat$ind %in% kd$indiv)
  
  ih_dat<-merge(ih_dat, kd, by.x="ind", by.y="indiv", all.x=T, all.y=T)
  ih_dat$K3<-as.factor(ih_dat$K3)
  ih_dat$K4<-as.factor(ih_dat$K4)
  
  # Remove the extra site:
  table(ih_dat$site.x==ih_dat$site.y)
  ih_dat$site.y<-NULL
  colnames(ih_dat)[which(colnames(ih_dat)=="site.x")]<-"site"
  head(ih_dat,3); dim(ih_dat)
  
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
# close individual diversity ----

# Question 3: Does fire influence the mode of reproduction?

#  Individual genetic distance & relatedness	# ----

genind_neutral

ddir<-"00_Data"
dir(ddir)

# SNPRelate (uses plink format - can use STRUCTURE files - Cenchrus_filt2)
# http://corearray.sourceforge.net/tutorials/SNPRelate/

library("SNPRelate")
# To install this package, start R (version "4.0") and enter:
# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("SNPRelate")

str.dir<-"/Users/annabelsmith/Documents/00_UQ_offline/buffel_burning/04_RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt4"
dir(str.dir)

# Read structure files:
bed.fn <- paste(str.dir,"Cenchrus_filt4.bed",sep="/")
fam.fn <- paste(str.dir,"Cenchrus_filt4.fam",sep="/")
bim.fn <- paste(str.dir,"Cenchrus_filt4.bim",sep="/")

# Make gds format
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "Cenchrus_filt4.gds")
snpgdsSummary(paste(ddir,"Cenchrus_filt4.gds",sep="/"))

# Import gds format :
genofile <- snpgdsOpen(paste(ddir,"Cenchrus_filt4.gds",sep="/"))
summary(genofile)

# Estimating IBD Using PLINK method of moments (MoM)
ibd <- snpgdsIBDMoM(genofile, maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

# Estimating IBD Using Maximum Likelihood Estimation (MLE)
set.seed(100)
snp.id <- sample(ibd$snp.id, 50)  # random 1500 SNPs
ibd_mle <- snpgdsIBDMLE(genofile, maf=0.05, missing.rate=0.05, num.thread=2)
ibd_mle.coeff <- snpgdsIBDSelection(ibd_mle)
head(ibd_mle.coeff); dim(ibd_mle.coeff)

# plot(ibd.coeff$kinship,ibd_mle.coeff$kinship)

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# combine:
neuc_df$kinship_mle<-ibd_mle.coeff$kinship
neuc_df$kinship_mom<-ibd.coeff$kinship
head(neuc_df)

table(neuc_df$kinship_mle>0.48)
table(neuc_df$kinship_mom>0.45)

# Are all samples WITHIN sites clones?
# Kinship coefficient == 1/2 * relatedness coefficient, such that:
# clone = 0.5
# full sib = 0.25
# parent offspring = 0.25

# Add site1 and site2
neuc_df$site1<-substr(neuc_df$ind1,1,unlist(gregexpr("_",neuc_df$ind1))-1)
neuc_df$site2<-substr(neuc_df$ind2,1,unlist(gregexpr("_",neuc_df$ind2))-1)
check.rows(neuc_df[,c(1,2,length(neuc_df)-1,length(neuc_df))])

# Add block1 and block2
neuc_df$block1<-substr(neuc_df$site1,1,3)
neuc_df$block2<-substr(neuc_df$site2,1,3)
check.rows(neuc_df[,c(1,2,length(neuc_df)-1,length(neuc_df))])

# Add same flags:
neuc_df$same_site<-ifelse(neuc_df$site1==neuc_df$site2,1,0)
neuc_df$same_block<-ifelse(neuc_df$block1==neuc_df$block2,1,0)
check.rows(neuc_df,3)

# what proportion of samples are clones, if we infer that kinship > 0.45 is a clone?
table(neuc_df$kinship_mle>0.45)
table(neuc_df$kinship_mle[neuc_df$same_block==1]>0.45)
table(neuc_df$kinship_mle[neuc_df$same_site==1]>0.45)

quartz("",6,6,dpi=100)
par(mfrow=c(2,2),mar=c(4,4,2,1), mgp=c(2.5,1,0))
hist(neuc_df$kinship_mle, main="", font.main=1, xlab="", ylab="",las=1)
arrows(0.45,0,0.45,2500,length=0, col="red", lwd=1.5)
mtext(text = "(a) all samples", side=3, line=0.5, adj=0)

text(0.43, 2000,paste("proportion of pairs\n> 0.45 = ",round(table(neuc_df$kinship_mle>0.45)[2]/sum(table(neuc_df$kinship_mle>0.45)), 2),sep=""),col="red", adj=1)

hist(neuc_df$kinship_mle[neuc_df$same_block==1], main="", font.main=1, xlab="", ylab="",las=1, ylim=c(0,200))
arrows(0.45,0,0.45,200,length=0, col="red", lwd=1.5)
mtext(text = "(b) within location", side=3, line=0.5, adj=0)

text(0.43, 160,paste("proportion of pairs\n> 0.45 = ",round(table(neuc_df$kinship_mle[neuc_df$same_block==1]>0.45)[2]/sum(table(neuc_df$kinship_mle[neuc_df$same_block==1]>0.45)), 2),sep=""),col="red", adj=1)

hist(neuc_df$kinship_mle[neuc_df$same_site==1], main="", font.main=1, xlab="", ylab="",las=1, ylim=c(0,100))
arrows(0.45,0,0.45,100,length=0, col="red", lwd=1.5)
mtext(text = "(c) within site", side=3, line=0.5, adj=0)

text(0.43, 80,paste("proportion of pairs\n> 0.45 = ",round(table(neuc_df$kinship_mle[neuc_df$same_site==1]>0.45)[2]/sum(table(neuc_df$kinship_mle[neuc_df$same_site==1]>0.45)), 2),sep=""),col="red", adj=1)

# Fewer than half of the samples within the same site are less than full sibs:
table(neuc_df$kinship_mle[neuc_df$same_site==1]<0.25)
table(neuc_df$kinship_mle[neuc_df$same_block==1]<0.25)

# although most clones are within the same block or site, there are still 800 pairs of plants that are clones from different blocks, i.e several km apart:
table(neuc_df$kinship_mle[neuc_df$same_block==0]>0.45)
head(neuc_df,3); dim(neuc_df)

# save.image("03_Workspaces/STEP_04_kinship.RData")

# Prepare data for analysis:

# is there a difference in the probability of being a clone (> 0.45) between the burnt and unburnt sites?

head(neuc_df,3); dim(neuc_df)

kinsh<-neuc_df

# set clone cutoff
kinsh$clone<-ifelse(kinsh$kinship_mle>0.45,1,0)

table(sdt4$site %in% kinsh$site1)
table(sdt4$site %in% kinsh$site2)

# Add burn categories to distance matrix:
head(sdt4,3); dim(sdt4)

bdat<-sdt4[,c("site","burn","burn2")]
colnames(bdat)[2:3]<-c("b2L","b3L")
head(bdat,3); dim(bdat)

kinsh<-merge(kinsh, bdat, by.x="site1", by.y="site", all.x=T, all.y=F)
colnames(kinsh)[which(colnames(kinsh)==c("b2L","b3L"))]<-c("b2L_s1","b3L_s1")
kinsh<-merge(kinsh, bdat, by.x="site2", by.y="site", all.x=T, all.y=F)
colnames(kinsh)[which(colnames(kinsh)==c("b2L","b3L"))]<-c("b2L_s2","b3L_s2")

head(kinsh,3); dim(kinsh)

# Add K genetic clusters:
head(ih_dat,3); dim(ih_dat)
indK<-ih_dat[,c("ind","K3")]
head(indK,3); dim(indK)

table(indK$ind %in% unique(c(kinsh$ind1, kinsh$ind2)))
table(unique(c(kinsh$ind1, kinsh$ind2)) %in% indK$ind)

kinsh<-merge(kinsh, indK, by.x="ind1", by.y="ind", all.x=T, all.y=F)
colnames(kinsh)[which(colnames(kinsh)==c("K3"))]<-c("K3_s1")
kinsh<-merge(kinsh, indK, by.x="ind2", by.y="ind", all.x=T, all.y=F)
colnames(kinsh)[which(colnames(kinsh)==c("K3"))]<-c("K3_s2")

# add same flag for K:
kinsh$same_K<-ifelse(kinsh$K3_s1==kinsh$K3_s2,1,0)
head(kinsh,3); dim(kinsh)

# remove pw comparisons from different blocks:
kinsh2<-kinsh[-which(kinsh$same_block==0),]

# add same flags for burn category:
kinsh2$same_b2L<-ifelse(kinsh2$b2L_s1==kinsh2$b2L_s2,1,0)
kinsh2$same_b3L<-ifelse(kinsh2$b3L_s1==kinsh2$b3L_s2,1,0)

# remove pw comparisons from different sites:
kinsh3<-kinsh2[-which(kinsh2$same_site==0),]

# remove pw comparisons from different K:
kinsh4<-kinsh[-which(kinsh$same_K==0),]

head(kinsh2,3); dim(kinsh2) # same block
head(kinsh3,3); dim(kinsh3) # same site
head(kinsh4,3); dim(kinsh4) # same K

# Focussing on within site comparisons (i.e. kinsh3), for each site calculate the proportion of pairwise distances that indicate clonality, as a measure of the rate of clonality:

prop_clone<-aggregate(clone~site1, data=kinsh3, FUN=function(x) length(x[x==1])/length(x))
colnames(prop_clone)<-c("site","prop_clone")
head(prop_clone); dim(prop_clone)

dir()
# save.image("03_Workspaces/STEP_04_kinship.RData")

# check:
site.now<-sample(kinsh3$site1,1)
d.now<-kinsh3[kinsh3$site1==site.now,c("site1","site2","ind1","ind2","kinship_mle","clone")]
d.now
prop_clone[prop_clone$site==site.now,]
length(d.now$clone[d.now$clone==1])/length(d.now$clone)

# add site data to clone proportions:
head(sdt4)
prop_clone<-merge(prop_clone, sdt4, by="site", all.x=T, all.y=F)
head(prop_clone); dim(prop_clone)

# Merge with the diversity score and K to account for background genetic structure:

head(gd_all,2); dim(gd_all)
gd_ds<-gd_all[,c("site","most_likelyK")]
prop_clone<-merge(prop_clone, gd_ds, by="site", all.x=T, all.y=F)
head(prop_clone); dim(prop_clone)

# save.image("03_Workspaces/STEP_04_kinship.RData")

# close individual distances ----

# Visualise distribution of FIS	(i.e. level of clonality) # ----

# The paralog filter workspace has the full set of FIS values for all loci (40711)
# Analyse distribution of FIS, following Reynes et al. 2021 MER):
load("03_Workspaces/paralog_filter.RData")


# Plot 

summary(fis_all$K1)
summary(fis_all$K2)
summary(fis_all$K3)

sd(fis_all$K1, na.rm=T)/mean(fis_all$K1, na.rm=T)
sd(fis_all$K2, na.rm=T)/mean(fis_all$K1, na.rm=T)
sd(fis_all$K3, na.rm=T)/mean(fis_all$K1, na.rm=T)

head(fis_all,3); dim(fis_all)
fisK1<-fis_all$K1[-which(is.na(fis_all$K1))]
fisK2<-fis_all$K2[-which(is.na(fis_all$K2))]
fisK3<-fis_all$K3[-which(is.na(fis_all$K3))]

length(fisK1); length(fisK2); length(fisK3)

dev.new(width=8, height=8, dpi=80, pointsize=20, noRStudioGD = T)
par(mfrow=c(2,2), mar=c(4,4,2,1), mgp=c(2.5,1,0))
h1<-hist(fisK1, xlab=expression(""*italic("F")[IS]*""), main="",las=1,yaxt="n", ylab="Number of loci (/1000)")
title(main="(a) K1", font.main=1, adj=0)
axis(side=2, at=pretty(x=h1$counts,n=3), labels=pretty(x=h1$counts,n=3)/1000, las=1)

h2<-hist(fisK2, xlab=expression(""*italic("F")[IS]*""), main="",las=1,yaxt="n", ylab="Number of loci (/1000)")
title(main="(b) K2", font.main=1, adj=0)
axis(side=2, at=pretty(x=h2$counts,n=5), labels=pretty(x=h2$counts,n=5)/1000, las=1)

h3<-hist(fisK3, xlab=expression(""*italic("F")[IS]*""), main="",las=1,las=1,yaxt="n", ylab="Number of loci (/1000)")
title(main="(c) K3", font.main=1, adj=0)
axis(side=2, at=pretty(x=h3$counts,n=5), labels=pretty(x=h3$counts,n=5)/1000, las=1)


# close FIS distribution ----

#  SUMMARY PLOTS: genetic diversity & reproductive mode	# ----

# save.image("03_Workspaces/STEP_04_kinship.RData")

# Site level genetic diversity, including most likely K
head(gd_all,3); dim(gd_all)

# The proportion of clones per site, including most likely K
head(prop_clone); dim(prop_clone)

# Individual heterozygosity, including actual K:
head(ih_dat,3); dim(ih_dat)

# Pairwise kinship coefficients:
head(neuc_df,3); dim(neuc_df)

# COLOURS and K name UPDATE Nov 2021
# K1 == Site 1, 2, 6u, 7u, 8u, 10 == PURPLE (#BEAED4)
# K2 == Site 7b, 8b == ORANGE (#FDC086)
# K3 == Site 3, 5, 6b, 11b1 == GREEN (#7FC97F)

# this is purple, orange, green (i.e. K1, K2, K3)
col.order<-c("#BEAED4","#FDC086","#7FC97F")

# PLOT GENETIC DIVERSITY: neutral markers:

dev.new(width=7.2,height=9.4,pointsize=16,dpi=80,noRStudioGD = T)
par(mfrow=c(4,3), mar=c(4,4,1,1), oma=c(0,0,1,1))
boxplot(ar_neutral~burn2, data=gd_all, las=1, ylab="allelic richness", xaxt="n", xlab="")
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(0.7,2,3.3), labels = c("unburnt","burnt","site 11"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Fire category", mgp=c(2,1,0))
mtext("(a)",side=3, at=0.5, line=0.3, cex=0.7)

boxplot(ar_neutral~most_likelyK, data=gd_all, las=1, ylab="allelic richness", xaxt="n", xlab="",col=col.order)
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(1,2,3), labels = c("1","2","3"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
mtext("(b)",side=3, at=0.5, line=0.3, cex=0.7)

full.int<-levels(interaction(gd_all$burn2,gd_all$most_likelyK))
ex.int<-levels(interaction(gd_all$burn2,gd_all$most_likelyK,drop=T))
col.fullint<-rep(col.order,rep(3,length(col.order)))
col.ind<-which(full.int %in% ex.int)

b1<-boxplot(ar_neutral~burn2+most_likelyK, data=gd_all, las=2, ylab="allelic richness", cex.axis=1,col=col.fullint[col.ind], xaxt="n", xlab="", drop=T)
axis(side=1, at=1:length(ex.int), labels = rep("",length(ex.int)), cex.axis=0.7)
axis(side=1, at=1:length(ex.int), labels = c("U","B","B","U","B","11"), tick = F, cex.axis=1, mgp=c(3,0.7,0), gap.axis = 0.2)
arrows(c(2.5,3.5),0,c(2.5,3.5),2,code=0,lty=2,col="grey50")
title(xlab="Fire category x K", mgp=c(2,1,0))
mtext("(c)",side=3, at=0.5, line=0.3, cex=0.7)
mtext(c("K1","K2","K3"), side=3, at=c(1.5,3,5), line=0.2, cex=0.7)

# individual heterozygosity
head(ih_dat,3)
ih_dat[c("ind","K3")]

boxplot(ind_het_neut~burn2, data=ih_dat, las=1, ylab="individual heterozygosity", xaxt="n", xlab="")
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(0.7,2,3.3), labels = c("unburnt","burnt","site 11"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Fire category", mgp=c(2,1,0))
mtext("(d)",side=3, at=0.5, line=0.3, cex=0.7)

b2<-boxplot(ind_het_neut~K3, data=ih_dat, las=1, ylab="individual heterozygosity", xaxt="n", xlab="",col=col.order)
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(1,2,3), labels = c("1","2","3"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
mtext("(e)",side=3, at=0.5, line=0.3, cex=0.7)

full.int2<-levels(interaction(ih_dat$burn2,ih_dat$K3))
ex.int2<-levels(interaction(ih_dat$burn2,ih_dat$K3,drop=T))
col.ind2<-which(full.int2 %in% ex.int2)

b3<-boxplot(ind_het_neut~burn2+K3, data=ih_dat, las=2, ylab="individual heterozygosity", cex.axis=1,col=col.fullint[col.ind2], xaxt="n", xlab="", drop=T)
axis(side=1, at=1:length(full.int), labels = rep("",length(full.int)), cex.axis=0.7)
axis(side=1, at=1:length(col.ind2), labels = rep(c("U","B","11"),3)[col.ind2], tick = F, cex.axis=0.8, mgp=c(3,0.7,0),gap.axis=0.1)
arrows(c(3.5,5.5),0,c(3.5,5.5),2,code=0,lty=2,col="grey50")
title(xlab="Fire category x K", mgp=c(2,1,0))
mtext("(f)",side=3, at=0.5, line=0.3, cex=0.7)
mtext(c("K1","K2","K3"), side=3, at=c(2,4.5,7), line=0.2, cex=0.7)

# Proportion clones
head(prop_clone); dim(prop_clone)

boxplot(prop_clone~burn2, data=prop_clone, las=1, ylab="proportion clones", xaxt="n", xlab="")
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(0.7,2,3.3), labels = c("unburnt","burnt","site 11"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Fire category", mgp=c(2,1,0))
mtext("(g)",side=3, at=0.5, line=0.3, cex=0.7)

b4<-boxplot(prop_clone~most_likelyK, data=prop_clone, las=1, ylab="proportion clones", xaxt="n", xlab="",col=col.order)
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(1,2,3), labels = c("1","2","3"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
mtext("(h)",side=3, at=0.5, line=0.3, cex=0.7)

full.int3<-levels(interaction(prop_clone$burn2,prop_clone$most_likelyK))
ex.int3<-levels(interaction(prop_clone$burn2,prop_clone$most_likelyK,drop=T))
col.ind3<-which(full.int3 %in% ex.int3)

b5<-boxplot(prop_clone~burn2+most_likelyK, data=prop_clone, las=2, ylab="proportion clones", cex.axis=1,col=col.fullint[col.ind3], xaxt="n", xlab="", drop=T)
axis(side=1, at=1:length(full.int3), labels = rep("",length(full.int3)), cex.axis=0.7)
axis(side=1, at=1:length(col.ind3), labels = rep(c("U","B","11"),3)[col.ind3], tick = F, cex.axis=1, mgp=c(3,0.7,0),gap.axis=0.1)
arrows(c(2.5,3.5),0,c(2.5,3.5),2,code=0,lty=2,col="grey50")
title(xlab="Fire category x K", mgp=c(2,1,0))
mtext("(i)",side=3, at=0.5, line=0.3, cex=0.7)
mtext(c("K1","K2","K3"), side=3, at=c(1.5,3,5), line=0.2, cex=0.7)

# Pairwise kinship coefficients
head(neuc_df,3); dim(neuc_df) # all pairwise coefs
head(kinsh3,2); dim(kinsh3) # same burn
head(kinsh4,3); dim(kinsh4) # same K

boxplot(kinship_mle~b3L_s1, data=kinsh3, las=1, ylab="kinship coefficient", xaxt="n", xlab="")
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(0.7,2,3.3), labels = c("unburnt","burnt","site 11"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Fire category", mgp=c(2,1,0))
mtext("(j)",side=3, at=0.5, line=0.3, cex=0.7)

b6<-boxplot(kinship_mle~K3_s1, data=kinsh4, las=1, ylab="kinship coefficient", xaxt="n", xlab="",col=col.order)
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(1,2,3), labels = c("1","2","3"), tick = F, cex.axis=1, mgp=c(3,0.7,0))
title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
mtext("(k)",side=3, at=0.5, line=0.3, cex=0.7)

full.int4<-levels(interaction(kinsh4$b3L_s1,kinsh4$K3_s1))
ex.int4<-levels(interaction(kinsh4$b3L_s1,kinsh4$K3_s1,drop=T))
col.ind4<-which(full.int4 %in% ex.int4)

b5<-boxplot(kinship_mle~b3L_s1+K3_s1, data=kinsh4, las=2, ylab="kinship coefficient", cex.axis=1,col=col.fullint[col.ind4], xaxt="n", xlab="", drop=T)
axis(side=1, at=1:length(full.int4), labels = rep("",length(full.int4)), cex.axis=0.7)
axis(side=1, at=1:length(col.ind4), labels = rep(c("U","B","11"),3)[col.ind4], tick = F, cex.axis=1, mgp=c(3,0.7,0),gap.axis=0.1)
arrows(c(3.5,5.5),0,c(3.5,5.5),2,code=0,lty=2,col="grey50")
title(xlab="Fire category x K", mgp=c(2,1,0))
mtext("(l)",side=3, at=0.5, line=0.3, cex=0.7)
mtext(c("K1","K2","K3"), side=3, at=c(2,4.5,7), line=0.2, cex=0.7)

# close summary plot ----

#  ANALYSE Genetic diversity & reproductive mode:	# ----

# save.image("03_Workspaces/STEP_04_analysis.RData")

## ANALYSIS:

# Site level genetic diversity, including most likely K
head(gd_all,3); dim(gd_all)

# The proportion of clones per site, including most likely K
head(prop_clone); dim(prop_clone)

# Individual heterozygosity, including actual K:
head(ih_dat,3); dim(ih_dat)

# Model ALLELIC RICHNESS:

# Make K a factor
gd_all$most_likelyK<-as.factor(gd_all$most_likelyK)
head(gd_all,3); dim(gd_all)

ar_neut.1<-lm(ar_neutral~treatment*most_likelyK, data=gd_all)
summary(ar_neut.1)
anova(ar_neut.1)

ar_nonneut.1<-lm(ar_nonneutral~treatment*most_likelyK, data=gd_all)
summary(ar_nonneut.1)
anova(ar_nonneut.1)

# Model INDIVIDUAL HETEROZYGOSITY:

head(ih_dat,3); dim(ih_dat)
table(ih_dat$K3, ih_dat$burn2) # rank deficient data structure
summary(ih_dat$ind_het_neut)

library(lme4)
library(lmerTest)

ih_neut.1<-lmer(ind_het_neut~burn*K3+(1|site), REML=F,data=ih_dat)
summary(ih_neut.1)
anova(ih_neut.1)

ih_nonneut.1<-lmer(ind_het_nonneut~burn*K3+(1|site), REML=F,data=ih_dat)
summary(ih_nonneut.1)
anova(ih_nonneut.1)

# Model PROPORTION CLONES:
head(prop_clone); dim(prop_clone)
library(mgcv)

pclone.neut1<-gam(prop_clone~burn*most_likelyK, data=prop_clone, family=betar)
summary(pclone.neut1)
anova(pclone.neut1)

# Predict from model:
nd_ar<-data.frame(treatment=levels(gd_all$treatment), most_likelyK=as.factor(rep(c(1,2,3), rep(2,3))))

nd_pc<-data.frame(burn=levels(prop_clone$burn), most_likelyK=as.factor(rep(c(1,2,3), rep(2,3))))

nd_ih<-data.frame(burn=as.factor(levels(ih_dat$burn)), K3=as.factor(rep(c(1,2,3), rep(2,3))))

# For site-level analysis (ar and pclone) there is no unburnt in K2 (rank deficient data structure)
table(gd_all$most_likelyK, gd_all$treatment)
table(prop_clone$most_likelyK, prop_clone$burn)
table(ih_dat$K3, ih_dat$burn) 

# Generate predictions for these levels but overwrite them after, to help with lining up the plots:

head(nd_ar)
head(nd_ih)
head(nd_pc)

ar_neut.pred<-predCI(ar_neut.1,nd_ar)
ar_nonneut.pred<-predCI(ar_nonneut.1,nd_ar)

ih_neut.pr<-predictSE(ih_neut.1, newdata = nd_ih, se.fit=T)
ih_neut.pr<-data.frame(nd_ih, fit=ih_neut.pr$fit, se=ih_neut.pr$se.fit)
ih_neut.pr$lci<-ih_neut.pr$fit-(1.96*ih_neut.pr$se)
ih_neut.pr$uci<-ih_neut.pr$fit+(1.96*ih_neut.pr$se)

ih_nonneut.pr<-predictSE(ih_nonneut.1, newdata = nd_ih, se.fit=T)
ih_nonneut.pr<-data.frame(nd_ih, fit=ih_nonneut.pr$fit, se=ih_nonneut.pr$se.fit)
ih_nonneut.pr$lci<-ih_nonneut.pr$fit-(1.96*ih_nonneut.pr$se)
ih_nonneut.pr$uci<-ih_nonneut.pr$fit+(1.96*ih_nonneut.pr$se)

pclone.pred<-predict(pclone.neut1, newdata = nd_pc, se.fit=T, type="response")
pclone.pred<-data.frame(nd_pc, fit=pclone.pred$fit, se=pclone.pred$se.fit)
pclone.pred$lci<-pclone.pred$fit-(1.96*pclone.pred$se)
pclone.pred$uci<-pclone.pred$fit+(1.96*pclone.pred$se)

ar_neut.pred
ar_nonneut.pred
pclone.pred
ih_neut.pr
ih_nonneut.pr

# Extract model coefficients for paper:

ar_neut.coef<-round(summary(ar_neut.1)$coefficient[,c(1)],3)
ar_nonneut.coef<-round(summary(ar_nonneut.1)$coefficient[,c(1)],3)
ih_neut.coef<-round(summary(ih_neut.1)$coefficient[,c(1)],3)
ih_nonneut.coef<-round(summary(ih_nonneut.1)$coefficient[,c(1)],3)
pclone.coef<-round(summary(pclone.neut1)$p.coeff,3)

ar_neut.se<-round(summary(ar_neut.1)$coefficient[,c(2)],3)
ar_nonneut.se<-round(summary(ar_nonneut.1)$coefficient[,c(2)],3)
ih_neut.se<-round(summary(ih_neut.1)$coefficient[,c(2)],3)
ih_nonneut.se<-round(summary(ih_nonneut.1)$coefficient[,c(2)],3)
pclone.se<-round(summary(pclone.neut1)$p.coeff,3)

# For site-level analysis (ar and pclone) there is no unburnt in K2 (rank deficient data structure); the p clone (gam) model returns this as zero; but the lm models leave it out (so need to add back in)

rownames(summary(ar_neut.1)$coefficient)

ar_neut.coef<-c(ar_neut.coef[1:4],NA,ar_neut.coef[5])
ar_nonneut.coef<-c(ar_nonneut.coef[1:4],NA,ar_nonneut.coef[5])
ar_neut.se<-c(ar_neut.se[1:4],NA,ar_neut.se[5])
ar_nonneut.se<-c(ar_nonneut.se[1:4],NA,ar_nonneut.se[5])

coef.df<-data.frame(term=c("Intercept","Burnt","K2","K3","Burnt K2","Burnt K3"),ar_neut.coef,ar_neut.se,ar_nonneut.coef,ar_nonneut.se,ih_neut.coef, ih_neut.se,ih_nonneut.coef,ih_nonneut.se,pclone.coef, pclone.se)
coef.df<-tidy.df(coef.df)
coef.df

p.df<-data.frame(term=c("Fire category","K","Fire category x K"),ar_neut.p=round(anova(ar_neut.1)$"Pr(>F)"[1:3],3),ar_nonneut.p=round(anova(ar_nonneut.1)$"Pr(>F)"[1:3],3),ih_neut.p=round(anova(ih_neut.1)$"Pr(>F)",3),ih_nonneut.p=round(anova(ih_nonneut.1)$"Pr(>F)",3),pclone.p=summary(pclone.neut1)$pTerms.table[,3])

# write.table(coef.df,"coef.txt", row.names = F, quote=F, sep="\t")
# write.table(p.df,"pdf.txt", row.names = F, quote=F, sep="\t")

# save.image("03_Workspaces/STEP_04_analysis.RData")

# PLOT MODEL ESTIMATES:

ar.fullint<-unique(interaction(gd_all$treatment,gd_all$most_likelyK))
ar.exint<-unique(interaction(gd_all$treatment,gd_all$most_likelyK,drop=T))
col.ind<-which(full.int %in% ex.int)

dev.new(width=7,height=7,pointsize=16,dpi=80,noRStudioGD = T)
par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,1,1))

plot(1:nrow(p_neut),p_neut$fit, ylim=c(min(p_neut$lci),max(p_neut$uci)), xlim=c(0.5,nrow(p_neut)+0.5), las=2, ylab="estimated allelic richness", cex.axis=1,col=col.fullint[col.ind], xaxt="n", xlab="", pch=20, cex=2)
axis(side=1, at=1:length(ar.exint), labels = rep("",length(ar.exint)), cex.axis=0.7)
axis(side=1, at=1:nrow(p_neut), labels = c("U","B","B","U","B"), tick = F, cex.axis=1, mgp=c(3,0.7,0), gap.axis = 0.2)
arrows(1:nrow(p_neut),p_neut$lci,1:nrow(p_neut),p_neut$uci,code=3, length=0.2, angle=90)
points(1:nrow(p_neut),p_neut$fit, pch=20, cex=2, col=col.fullint[col.ind])
arrows(c(2.5,3.5),0,c(2.5,3.5),2,code=0,lty=2,col="grey50")
title(xlab="Fire category x K", mgp=c(2,1,0))
mtext("(c)",side=3, at=0.5, line=0.3, cex=0.7)
mtext(c("K1","K2","K3"), side=3, at=c(1.5,3,5), line=0.2, cex=0.7)



# close analyse genetic diversity ----



