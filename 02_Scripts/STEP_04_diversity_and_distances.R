
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
  
  # non-neutral
  gendiv_nonneutral <- basic.stats(genind_nonneutral, diploid = TRUE, digits = 2)
  
  str(gendiv_nonneutral) # a few seconds
  head(gendiv_nonneutral$Ho)
  tail(gendiv_nonneutral$Ho)
  
  gd_nonneutral<-data.frame(site=names(apply(gendiv_nonneutral$Ho,2,mean,na.rm=T)),max_n=apply(gendiv_nonneutral$n.ind.samp,2,max,na.rm=T),Ho=apply(gendiv_nonneutral$Ho,2,mean,na.rm=T),He=apply(gendiv_nonneutral$Hs,2,mean,na.rm=T),Fis=apply(gendiv_nonneutral$Fis,2,mean,na.rm=T))
  gd_nonneutral<-tidy.df(gd_nonneutral)
  head(gd_nonneutral); dim(gd_nonneutral)
  
  gd_nonneutral$site<-substr(gd_nonneutral$site,1,nchar(as.character(gd_nonneutral$site))-3)
  head(gd_nonneutral)
  
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
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
  
  # check:
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
  # save.image("03_workspaces/STEP04_divdist_ALL.RData")
  
  # close calc genetic diversity ----
  
#  ANALYSE SITE-LEVEL Genetic diversity:	# ----
  
  # Analyse the influence of fire treatment on site-level genetic diversity. 
  
  # Site-level data:
  gd_all<-read.table("RESULTS/Diversity_and_distance/Genetic_Diversity/Genetic_Diversity_ALL.txt", header=T)
  head(gd_all,3)
  
  # simplify site data
  sdt3<-sdat[,c("block","site","pop","burn_unburnt","no_samples","lat","long")]
  sdt3$site<-paste("X",substr(sdt3$site, 4, nchar(as.character(sdt3$site))),sep="")
  
  # merge with genetic diversity data:
  sdt3$site %in% gd_all$site
  gd_all$site %in% sdt3$site
  gd_all<-merge(gd_all, sdt3, by="site", all.x=T, all.y=F)
  
  # add second treatment variable:
  gd_all$treatment<-gd_all$burn_unburnt
  gd_all$burn_unburnt<-NULL
  gd_all$trt<-as.character(gd_all$treatment)
  gd_all$trt[grep("b",gd_all$trt)]<-"b1"
  gd_all$trt[grep("X11",gd_all$site)]<-"b2"
  gd_all$trt<-factor(gd_all$trt, levels=c("u","b1","b2"))
  head(gd_all,3); dim(gd_all)
  
  #### Accounting for population genetic structure
  
  # Given the strong background genetic structure in the data, we must control for phylogeny. 
  # We modelled background genetic structure using two methods: the mean site-level diversity score, calculated from individual level diversity scores and the most common K per site
  # Regardless of the method used, there was no effect of background genetic structure on site-level allelic richness. Both methods showed the treatment only model to be superior, although it was not better than a null model. 

  # Assignment probabilities from STRUCTURE:
  site_assig<-read.table("00_Data/assig_prob_K3_Cenchrus_filt2.txt", header=T)
  site_assig<-cbind(site_assig, matrix(data=c(1,2,3),ncol=3,nrow=nrow(site_assig),byrow=T))
  site_assig$genstr<-apply(site_assig[,3:8],1,function(x) weighted.mean(x[4:6], x[1:3]))
  head(site_assig)
  range(site_assig$genstr)

  ## -- ** METHOD 1: Admixture Diversity Score:
  
  ### ---- provide calculations on the level of admixture per population
  
  # Working from equations in Harismendy et al. 2019, Journal of the American Medical Informatics Association, 26(5), 2019, 457–46
  
  # 93 individuals 
  head(site_assig); dim(site_assig)

  # Associated site data (19 sites):
  head(sdt2,2); dim(sdt2)
  head(gd_all,3); dim(gd_all)
  
  # Use individual level data for the analysis:
  head(site_assig); dim(site_assig)

  # This is the Diversity Score from Harismendy. We don't need the first eqn (the cumulative admixture fraction), because ours already sum to 1. 
  
  # This diversity score is the same as the shannon diversity index but with a scaling factor to account for the number of clusters (Hmax). That is, it makes the diversity score relative to complete evenness. 
  
  DS<-function(x,K) {
    Hmax<-K*((1/K)*log(1/K))
    (-sum(x*log(x)))/-Hmax
  }
  
  # SET K:
  K<-3
  
  # Calculate admixture diversity for all individuals:
  ds_div<-apply(site_assig[3:5], 1, DS, K=3)
  range(ds_div)
  site_assig$ds_div<-ds_div
  head(site_assig); dim(site_assig)
  
  # Summarise site-level admixture diversity:
  addiv_site<-aggregate(ds_div~site,FUN=mean,data=site_assig)
  gd_all<-merge(gd_all, addiv_site, by="site", all.x=T, all.y=F)
  head(gd_all,2); dim(gd_all)
  
  ## -- ** METHOD 2: K
  
  # ih_dat has K per individual; it was merged with individual data in the Individual genetic diversity section below; they are the correct Ks, as in, they reflect the K with the highest assignment probability per individual; check site 7b which has a mixture of K=1 and K=3:
  head(ih_dat,3); dim(ih_dat)
  head(gd_all,2); dim(gd_all)
  
  # Use which.max to get the mode - the most common K for each site:
  Ksite<-data.frame(site=names(tapply(ih_dat$K3,ih_dat$site, FUN = function(x) names(which.max(summary(x))))),K3=tapply(ih_dat$K3,ih_dat$site, FUN = function(x) names(which.max(summary(x)))))
  Ksite<-tidy.df(Ksite)
  
  gd_all<-merge(gd_all, Ksite, by="site", all.x=T, all.y=F)
  # Check it against the structure plot
  check.rows(gd_all[,c("site","K3")])
  
  ## ANALYSIS:
  head(gd_all,2); dim(gd_all)

  ### With K:
  {
  # neutral models:
  mod7.a<-lm(ar_neutral~1, data=gd_all)
  mod7.b<-lm(ar_neutral~trt, data=gd_all)
  mod7.c<-lm(ar_neutral~K3, data=gd_all)
  mod7.d<-lm(ar_neutral~trt+K3, data=gd_all)
  mod7.e<-lm(ar_neutral~trt*K3, data=gd_all)
  AICc(mod7.a); AICc(mod7.b); AICc(mod7.c); AICc(mod7.d); AICc(mod7.e)
  
  summary(mod7.b); anova(mod7.b)

  nd_neut<-data.frame(trt=levels(gd_all$trt))
  p_neut<-predict(mod7.b, newdata = nd_neut, se.fit=T)
  p_neut<-data.frame(nd_neut, fit=p_neut$fit, se=p_neut$se.fit)
  p_neut$lci<-p_neut$fit-(1.96*p_neut$se)
  p_neut$uci<-p_neut$fit+(1.96*p_neut$se)
  p_neut
  
  # non-neutral models:
  mod8.a<-lm(ar_nonneutral~1, data=gd_all)
  mod8.b<-lm(ar_nonneutral~trt, data=gd_all)
  mod8.c<-lm(ar_nonneutral~K3, data=gd_all)
  mod8.d<-lm(ar_nonneutral~trt+K3, data=gd_all)
  mod8.e<-lm(ar_nonneutral~trt*K3, data=gd_all)
  AICc(mod8.a); AICc(mod8.b); AICc(mod8.c); AICc(mod8.d); AICc(mod8.e)
  
  summary(mod8.b); anova(mod8.b)
  
  nnd_neut<-data.frame(trt=levels(gd_all$trt))
  p_nneut<-predict(mod8.b, newdata = nnd_neut, se.fit=T)
  p_nneut<-data.frame(nnd_neut, fit=p_nneut$fit, se=p_nneut$se.fit)
  p_nneut$lci<-p_nneut$fit-(1.96*p_nneut$se)
  p_nneut$uci<-p_nneut$fit+(1.96*p_nneut$se)
  p_nneut
  
  ## PLOT neutral and non-neutral together
  
  dev.new(height=4,width=9,dpi=100, noRStudioGD = T,pointsize=16)
  par(mfrow=c(1,2),mar=c(3,4,1.5,0.5), mgp=c(2.8,0.8,0))
  plot(1:3, p_neut$fit, ylim=c(min(p_neut$lci), max(p_neut$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Neutral allelic richness", xlab="", xaxt="n", col="black")
  arrows(1:3, p_neut$lci,1:3, p_neut$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
  text(1,1.3, labels=paste("P = ",round(anova(mod1.b)$"Pr(>F)"[1],2),sep=""), adj=0)
  mtext("(a)", side=3, line=0.5, adj=0)
  
  plot(1:3, p_nneut$fit, ylim=c(min(p_nneut$lci), max(p_nneut$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Non-neutral allelic richness", xlab="", xaxt="n", col="black")
  arrows(1:3, p_nneut$lci,1:3, p_nneut$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
  text(1,1.3, labels=paste("P = ",round(anova(mod2.b)$"Pr(>F)"[1],2),sep=""), adj=0)
  mtext("(b)", side=3, line=0.5, adj=0)
  
  } # close K
  
  # With diversity score:
  {
  # neutral models:
  mod1.a<-lm(ar_neutral~1, data=gd_all)
  mod1.b<-lm(ar_neutral~trt, data=gd_all)
  mod1.b2<-lm(ar_neutral~ds_div, data=gd_all)
  mod1.c<-lm(ar_neutral~trt+ds_div, data=gd_all)
  mod1.d<-lm(ar_neutral~trt*ds_div, data=gd_all)
  AICc(mod1.a); AICc(mod1.b); AICc(mod1.b2); AICc(mod1.c); AICc(mod1.d)
  
  summary(mod1.b); anova(mod1.b)
  summary(mod1.c); anova(mod1.c)
  
  nd_neut<-data.frame(trt=levels(gd_all$trt))
  p_neut<-predict(mod1.b, newdata = nd_neut, se.fit=T)
  p_neut<-data.frame(nd_neut, fit=p_neut$fit, se=p_neut$se.fit)
  p_neut$lci<-p_neut$fit-(1.96*p_neut$se)
  p_neut$uci<-p_neut$fit+(1.96*p_neut$se)
  p_neut
  
  # non-neutral models:
  mod2.a<-lm(ar_nonneutral~1, data=gd_all)
  mod2.b<-lm(ar_nonneutral~trt, data=gd_all)
  mod2.b2<-lm(ar_nonneutral~ds_div, data=gd_all)
  mod2.c<-lm(ar_nonneutral~trt+ds_div, data=gd_all)
  mod2.d<-lm(ar_nonneutral~trt*ds_div, data=gd_all)
  AICc(mod2.a); AICc(mod2.b); AICc(mod2.b2); AICc(mod2.c); AICc(mod2.d)
  
  summary(mod2.b); anova(mod2.b)
  
  nnd_neut<-data.frame(trt=levels(gd_all$trt))
  p_nneut<-predict(mod2.b, newdata = nnd_neut, se.fit=T)
  p_nneut<-data.frame(nnd_neut, fit=p_nneut$fit, se=p_nneut$se.fit)
  p_nneut$lci<-p_nneut$fit-(1.96*p_nneut$se)
  p_nneut$uci<-p_nneut$fit+(1.96*p_nneut$se)
  p_nneut
  
  ## PLOT neutral and non-neutral together
  
  dev.new(height=4,width=9,dpi=100, noRStudioGD = T,pointsize=16)
  par(mfrow=c(1,2),mar=c(3,4,1.5,0.5), mgp=c(2.8,0.8,0))
  plot(1:3, p_neut$fit, ylim=c(min(p_neut$lci), max(p_neut$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Neutral allelic richness", xlab="", xaxt="n", col="black")
  arrows(1:3, p_neut$lci,1:3, p_neut$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
  text(1,1.3, labels=paste("P = ",round(anova(mod1.b)$"Pr(>F)"[1],2),sep=""), adj=0)
  mtext("(a)", side=3, line=0.5, adj=0)
  
  plot(1:3, p_nneut$fit, ylim=c(min(p_nneut$lci), max(p_nneut$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Non-neutral allelic richness", xlab="", xaxt="n", col="black")
  arrows(1:3, p_nneut$lci,1:3, p_nneut$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
  text(1,1.3, labels=paste("P = ",round(anova(mod2.b)$"Pr(>F)"[1],2),sep=""), adj=0)
  mtext("(b)", side=3, line=0.5, adj=0)
  
  } # close diversity score
  
  # close analyse genetic diversity ----
  
#  Individual genetic diversity	# ----
  
  # Neutral data set:
  ghead(ddat); dim(ddat)
  
  # Non-neutral data set:
  nndat<-read.table(paste("00_Data/Filtered_DartSeq_format","dartseq_filt3_adaptive.txt", sep="/"), header=T)
  ghead(nndat); dim(nndat)
  
  # Codes for onerow formatted data:
  # 0 = Reference allele homozygote
  # 1 = SNP allele homozygote
  # 2 = heterozygote
  
  # Neutral
  # 3 MINS on laptop 
  
  ih_out<-data.frame(ddat[,1:2],ind_het=NA)
  head(ih_out,25)
  
  for(i in 1:nrow(ih_out)){
    ind.cons<-as.character(ddat[i,3:length(ddat)])
    t.cons<-data.frame(ind.cons=as.numeric(names(table(ind.cons))),count=as.numeric(table(ind.cons)),stringsAsFactors=F)
    
    if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] else stop("no zero category")
    
    if(length(which(t.cons$ind.cons==2))==0) ih_out[i,3]<-"no_heterozygotes" else ih_out[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)
    
  } # close for i
  
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
  
  ih_out_nn$ind_het<-as.numeric(ih_out_nn$ind_het)
  
  ghead(ddat); dim(ddat)
  ghead(nndat); dim(nndat)
  
  # add site data and cluster assignment:
  
  # sdt[,which(is.na(colnames(sdt)))]<-NULL
  head(sdt, 3); dim(sdt)
  head(ih_out,6); dim(ih_out)
  head(ih_out_nn,6); dim(ih_out_nn)
  
  colnames(ih_out)[which(colnames(ih_out)=="ind_het")]<-"ind_het_neut"
  colnames(ih_out_nn)[which(colnames(ih_out_nn)=="ind_het")]<-"ind_het_nonneut"
  
  ih_dat<-data.frame(ih_out,ind_het_nonneut=ih_out_nn[,3])
  ih_dat<-merge(ih_dat, sdt, by="site", all.x=T, all.y=F)
  head(ih_dat,3); dim(ih_dat)
  
  # merge K on individual:
  # kd has the correct individual level K values:
   head(kd,3); dim(kd)
  table(kd$indiv %in% ih_dat$ind)
  table(ih_dat$ind %in% kd$indiv)
  
  ih_dat<-merge(ih_dat, kd, by.x="ind", by.y="indiv", all.x=T, all.y=T)
  ih_dat$K3<-as.factor(ih_dat$K3)
  ih_dat$K4<-as.factor(ih_dat$K4)
  head(ih_dat,3); dim(ih_dat)
  
  # Add admixture diversity score:
  head(site_assig,2); dim(site_assig)
  sassig<-site_assig[,c("indiv","ds_div")]
  head(sassig); dim(sassig)
  table(sassig$indiv %in% ih_dat$ind)
  
  ih_dat<-merge(ih_dat, sassig, by.x="ind", by.y="indiv", all.x=T, all.y=F)
  
  # MODEL individual heterozygosity ~ burn category:
  
  # NEUTRAL
  head(ih_dat,3); dim(ih_dat)
  table(ih_dat$K3, ih_dat$burn2) # rank deficient data structure
  summary(ih_dat$ind_het_neut)
  
  mod3.a<-lmer(ind_het_neut~1+(1|site), REML=F, data=ih_dat)
  mod3.b<-lmer(ind_het_neut~burn2+(1|site),REML=F, data=ih_dat)
  mod3.c<-lmer(ind_het_neut~K3+(1|site),REML=F, data=ih_dat)
  mod3.d<-lmer(ind_het_neut~burn2+K3+(1|site), REML=F,data=ih_dat)
  mod3.e<-lmer(ind_het_neut~burn2*K3+(1|site), REML=F,data=ih_dat)
  
  AICc(mod3.a); AICc(mod3.b); AICc(mod3.c); AICc(mod3.d); AICc(mod3.e)

  summary(mod3.c); anova(mod3.c)
  summary(mod3.d); anova(mod3.d)
  
  # plot fire + K for fun (mod3.d)
  nd_fk<-data.frame(K3=rep(as.factor(c(1,2,3)),rep(3,3)), burn2=factor(c("u","b","b2"), levels=c("u","b","b2")))
  p_fk<-predictSE(mod3.d, newdata = nd_fk, se.fit=T)
  p_fk<-data.frame(nd_fk, fit=p_fk$fit, se=p_fk$se.fit)
  p_fk$lci<-p_fk$fit-(1.96*p_fk$se)
  p_fk$uci<-p_fk$fit+(1.96*p_fk$se)
  p_fk
  
  dev.new(height=4,width=4,noRStudioGD = T,dpi=100, pointsize=20)
  xofs.now<-0.1
  par(mfrow=c(1,1),mar=c(3,3.5,0.5,0.1), mgp=c(2.6,0.8,0))
  plot(1:3-xofs.now, p_fk$fit[p_fk$K3==1], ylim=c(min(p_fk$lci), max(p_fk$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Individual heterozygosity", xlab="", xaxt="n", col=switch.col2[1])
  arrows(1:3-xofs.now, p_fk$lci[p_fk$K3==1],1:3-xofs.now, p_fk$uci[p_fk$K3==1], code=3, angle=90, length=0.05,lwd=1.5)
  points(1:3-xofs.now, p_fk$fit[p_fk$K3==1], pch=20, col=switch.col2[1])
  arrows(1:3, p_fk$lci[p_fk$K3==2],1:3, p_fk$uci[p_fk$K3==2], code=3, angle=90, length=0.05,lwd=1.5)
  points(1:3, p_fk$fit[p_fk$K3==2], pch=20, col=switch.col2[2])
  arrows(1:3+xofs.now, p_fk$lci[p_fk$K3==3],1:3+xofs.now, p_fk$uci[p_fk$K3==3], code=3, angle=90, length=0.05,lwd=1.5)
  points(1:3+xofs.now, p_fk$fit[p_fk$K3==3], pch=20, col=switch.col2[3])
  axis(side=1, at=c(1:3), labels = c("u","b","b2"))
  title(xlab="Fire category", mgp=c(2,1,0))

  # Adding the diversity score does not improve the model fit:
  mod3.f<-lmer(ind_het~K3+ds_div+(1|site),REML=F, data=ih_dat)
  summary(mod3.f)
  AICc(mod3.c); AICc(mod3.f)
  
  # NON-NEUTRAL
  head(ih_dat,3); dim(ih_dat)

  mod6.a<-lmer(ind_het_nonneut~1+(1|site), REML=F, data=ih_dat)
  mod6.b<-lmer(ind_het_nonneut~burn2+(1|site),REML=F, data=ih_dat)
  mod6.c<-lmer(ind_het_nonneut~K3+(1|site),REML=F, data=ih_dat)
  mod6.d<-lmer(ind_het_nonneut~burn2+K3+(1|site), REML=F,data=ih_dat)
  mod6.e<-lmer(ind_het_nonneut~burn2*K3+(1|site), REML=F,data=ih_dat)
  
  AICc(mod6.a); AICc(mod6.b); AICc(mod6.c); AICc(mod6.d); AICc(mod6.e)
  
  summary(mod6.c); anova(mod6.c)
  summary(mod6.d); anova(mod6.d)
  
  # MODEL ESTIMATES:
  
  # NEUTRAL
  nd_ih<-data.frame(K3=as.factor(c(1,2,3)))
  p_ih<-predictSE(mod3.c, newdata = nd_ih, se.fit=T)
  p_ih<-data.frame(nd_ih, fit=p_ih$fit, se=p_ih$se.fit)
  p_ih$lci<-p_ih$fit-(1.96*p_ih$se)
  p_ih$uci<-p_ih$fit+(1.96*p_ih$se)
  p_ih
  
  # NON-NEUTRAL
  nd_ih_nn<-data.frame(K3=as.factor(c(1,2,3)))
  p_ih_nn<-predictSE(mod6.c, newdata = nd_ih_nn, se.fit=T)
  p_ih_nn<-data.frame(nd_ih_nn, fit=p_ih_nn$fit, se=p_ih_nn$se.fit)
  p_ih_nn$lci<-p_ih_nn$fit-(1.96*p_ih_nn$se)
  p_ih_nn$uci<-p_ih_nn$fit+(1.96*p_ih_nn$se)
  p_ih_nn
  
  # PLOT TOP MODEL: genetic diversity ~ genetic cluster:
  
  # COLOURS:
  # The main structure plot used the colour order in switch.col and the structure plot in the dendro used switch.col2 - rearranged so they would match
  # switch.col is green, purple, orange
  switch.col<-c("#7FC97F", "#BEAED4", "#FDC086")
  
  # switch.col2 is orange, green, purple
  switch.col2<-c("#FDC086", "#7FC97F", "#BEAED4")
  
  # switch.col2 is in the correct order for plotting K1:K3
  ih_dat[which(ih_dat$ind=="X11b1_02"),] # K2 == green
  ih_dat[which(ih_dat$ind=="X08b_03"),] # K1 == orange
  ih_dat[which(ih_dat$ind=="X02b_03"),] # K3 == purple
  # plot(1:3, 1:3, col=switch.col2, cex=2, pch=16)
  
  head(ih_dat,3); dim(ih_dat)
  
  # NEUTRAL
  dev.new(height=4,width=4,noRStudioGD = T,dpi=100, pointsize=20)
  par(mfrow=c(1,1),mar=c(3,3.5,0.5,0.1), mgp=c(2.6,0.8,0))
  plot(1:3, p_ih$fit, ylim=c(min(p_ih$lci), max(p_ih$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Individual heterozygosity", xlab="", xaxt="n", col=switch.col2)
  arrows(1:3, p_ih$lci,1:3, p_ih$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c(1,2,3))
  title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
  text(2,0.12, labels="P < 0.001", adj=0)
  points(1:3, p_ih$fit, col=switch.col2, cex=2)
  points(1:3, p_ih$fit, col=switch.col2, cex=2,pch=20)
  
  # NON-NEUTRAL
  dev.new(height=4,width=4,noRStudioGD = T,dpi=100, pointsize=20)
  par(mfrow=c(1,1),mar=c(3,3.5,0.5,0.1), mgp=c(2.6,0.8,0))
  plot(1:3, p_ih_nn$fit, ylim=c(min(p_ih_nn$lci), max(p_ih_nn$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Individual heterozygosity", xlab="", xaxt="n", col=switch.col2)
  arrows(1:3, p_ih_nn$lci,1:3, p_ih_nn$uci, code=3, angle=90, length=0.05,lwd=1.5)
  axis(side=1, at=c(1:3), labels = c(1,2,3))
  title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
  text(2,0.12, labels="P < 0.001", adj=0)
  points(1:3, p_ih_nn$fit, col=switch.col2, cex=2)
  points(1:3, p_ih_nn$fit, col=switch.col2, cex=2,pch=20)
  
  # save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# close individual diversity ----
  
# Question 3: Does fire influence the mode of reproduction?
  
#  Individual genetic distance & relatedness	# ----

genind_neutral

# Euclidean distance (from adegenet, works on genind object)
# See Shirk et al. 2017 - Euc dist performs as well as other distance measures:
neu_euc_dist<-dist(genind_neutral, method="euclidean")

# get names BEFORE as.matrix:
neu_euc_names<-combn(attr(neu_euc_dist, "Labels"),2)
neu_euc<-as.matrix(neu_euc_dist)
neuc_df<-data.frame(ind1=neu_euc_names[1,], ind2=neu_euc_names[2,], dist_euc=neu_euc[lower.tri(neu_euc)])

# check:
head(neuc_df)
neu_euc[1:5,1:5]
neu_euc[90:93,1:5]
neuc_df[90:93,]

# Bray-Curtis (from vegan, works on data.frame):

ddir<-"/Users/annabelsmith/Documents/00_UQ_offline/Binyin_Winter/00_Data/Filtered_DartSeq_format"
ddat<-read.table(paste(ddir, "dartseq_filt2.txt", sep="/"), header=T)
ghead(ddat); dim(ddat)
neu_bray<-vegdist(x = ddat[,3:length(ddat)], method="bray",na.rm=T)
neuc_df$dist_bray<-neu_bray
head(neuc_df); dim(neuc_df)
# these are very highy correlated:
# plot(neuc_df$dist_euc, neuc_df$dist_bray)
# cor.test(neuc_df$dist_euc, neuc_df$dist_bray)

# SNPRelate (uses plink format - can use STRUCTURE files - Cenchrus_filt2)
# http://corearray.sourceforge.net/tutorials/SNPRelate/

library("SNPRelate")
# To install this package, start R (version "4.0") and enter:
# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("SNPRelate")

str.dir<-"/Users/annabelsmith/Documents/00_UQ_offline/Binyin_Winter/RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt2"
dir(str.dir)

# Read structure files:
bed.fn <- paste(str.dir,"Cenchrus_filt2.bed",sep="/")
fam.fn <- paste(str.dir,"Cenchrus_filt2.fam",sep="/")
bim.fn <- paste(str.dir,"Cenchrus_filt2.bim",sep="/")

# Make gds format
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "Cenchrus_filt2.gds")
snpgdsSummary(paste(ddir,"Cenchrus_filt2.gds",sep="/"))

# Import gds format :
genofile <- snpgdsOpen(paste(ddir,"Cenchrus_filt2.gds",sep="/"))
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

# Add proportion matching genotypes

head(neuc_df,3); dim(neuc_df)
ghead(ddat); dim(ddat)

neuc_df$prop_match<-NA

### prop_match is the proportion of NON-NA genotypes that match

## WARNING: takes approx 50 min

for (i in 1:nrow(neuc_df)){
  
  pair.thisrun<-neuc_df[i,c(1,2)]
  pair1.thisrun<-ddat[which(ddat$ind==pair.thisrun$ind1),3:length(ddat)]
  pair2.thisrun<-ddat[which(ddat$ind==pair.thisrun$ind2),3:length(ddat)]
  comp.thisrun<-pair1.thisrun==pair2.thisrun
  t1<-table(comp.thisrun)
  neuc_df$prop_match[i]<-t1[2]/sum(t1)
  
} # close i

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# compare methods:

quartz("",6,6,dpi=100)
par(mfrow=c(2,2),mar=c(4,4,2,1), mgp=c(2.5,1,0))
plot(neuc_df$dist_euc, neuc_df$dist_bray, pch=20, xlab="Euclidean", ylab="Bray-Curtis")
plot(neuc_df$dist_euc, neuc_df$kinship_mle, pch=20, xlab="Euclidean", ylab="Kinship MLE")
plot(neuc_df$dist_bray, neuc_df$kinship_mle, pch=20, xlab="Bray-Curtis", ylab="Kinship MLE")
plot(neuc_df$kinship_mle, neuc_df$prop_match, pch=20, xlab="Kinship MLE", ylab="Proportion matching")

plot(neuc_df$kinship_mom, neuc_df$kinship_mle,  pch=20, xlab="Kinship MoM", ylab="Kinship MLE")

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

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

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# close individual distances ----

#  Analyse Individaul distance & relatedness	# ----

# is there a difference in the probability of being a clone (> 0.45) between the burnt and unburnt sites?

head(neuc_df,3); dim(neuc_df)

kinsh<-neuc_df

# set clone cutoff
kinsh$clone<-ifelse(kinsh$kinship_mle>0.45,1,0)

sdt$site %in% kinsh$site1
sdt$site %in% kinsh$site2

# Add burn categories to distance matrix:
bdat<-sdt[,c("site","burn","burn2")]
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

indK$ind %in% unique(c(kinsh$ind1, kinsh$ind2))
unique(c(kinsh$ind1, kinsh$ind2)) %in% indK$ind

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

# check:
site.now<-sample(kinsh3$site1,1)
d.now<-kinsh3[kinsh3$site1==site.now,c("site1","site2","ind1","ind2","kinship_mle","prop_match","clone")]
d.now
prop_clone[prop_clone$site==site.now,]
length(d.now$clone[d.now$clone==1])/length(d.now$clone)

# add site data to clone proportions:
head(sdt2)
prop_clone<-merge(prop_clone, sdt2, by="site", all.x=T, all.y=F)
# plot(prop_clone$burn2, prop_clone$prop_clone)

# Merge with the diversity score and K to account for background genetic structure:

head(gd_all,2); dim(gd_all)
gd_ds<-gd_all[,c("site","ds_div","K3")]
prop_clone<-merge(prop_clone, gd_ds, by="site", all.x=T, all.y=F)
head(prop_clone); dim(prop_clone)

# What is the effect of burn category on the probability of being a clone (site level analysis with binomial models)?

mod5.bNULL<-glm(prop_clone~1, data=prop_clone, family="binomial")
mod5.b1<-glm(prop_clone~burn2, data=prop_clone, family="binomial")
mod5.b2<-glm(prop_clone~K3, data=prop_clone, family="binomial")
mod5.b3<-glm(prop_clone~burn2+K3, data=prop_clone, family="binomial")
mod5.b4<-glm(prop_clone~burn2*K3, data=prop_clone, family="binomial")
summary(mod5.b1)$coefficients
anova(mod5.b1)
AICc(mod5.bNULL); AICc(mod5.b1); AICc(mod5.b2); AICc(mod5.b3); AICc(mod5.b4)

nd.b1<-data.frame(burn2=factor(c("u","b","b2"),level=c("u","b","b2")))
pr.b1<-predict(mod5.b1, newdata = nd.b1, type="response", se.fit=T)
pr.b1<-data.frame(nd.b1, fit=pr.b1$fit, se=pr.b1$se.fit)
pr.b1$lci<-pr.b1$fit-(1.96*pr.b1$se)
pr.b1$uci<-pr.b1$fit+(1.96*pr.b1$se)
# pr.b1$lci[which(pr.b1$lci<0)]<-0

# get effects and 95% CI on link scale
pr.b1_EFS<-predict(mod5.b1, newdata = nd.b1, type="link", se.fit=T)
pr.b1_EFS<-data.frame(nd.b1, fit=pr.b1_EFS$fit, se=pr.b1_EFS$se.fit)
pr.b1_EFS$lci<-pr.b1_EFS$fit-(1.96*pr.b1_EFS$se)
pr.b1_EFS$uci<-pr.b1_EFS$fit+(1.96*pr.b1_EFS$se)

# PLOT ESTIMATED site-level proportion CLONE
dev.new(height=4,width=7,dpi=100, noRStudioGD = T,pointsize=15)
par(mfrow=c(1,2),mar=c(3,4,2,1), mgp=c(2.5,0.8,0))
plot(1:3, pr.b1$fit, ylim=c(min(pr.b1$lci), max(pr.b1$uci)), pch=20, xlim=c(0.75, 3.25), las=1, ylab="Proportion asexual individuals", xlab="", xaxt="n", col="black")
arrows(1:3, pr.b1$lci,1:3, pr.b1$uci, code=3, angle=90, length=0.05,lwd=1.5)
axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
mtext("(a) response scale",side=3, line=0.5, adj=0)

# plot effect sizes:
par(mgp=c(2.3,0.8,0))
plot(1:3,pr.b1_EFS$fit,  ylim=c(min(pr.b1_EFS$lci), max(pr.b1_EFS$uci)), pch=20, las=1, ylab="Proportion asexual individuals", xlab="", col="black",xaxt="n",xlim=c(0.75, 3.25))
arrows(1:3,pr.b1_EFS$lci,1:3, pr.b1_EFS$uci, code=3, angle=90, length=0.05,lwd=1.5)
arrows(0,0, 4,0, code=3, angle=90, length=0,lwd=1.5)
axis(side=1, at=c(1:3), labels = c("unburnt","burnt","site 11"))
mtext("(b) logit scale",side=3, line=0.5, adj=0)

# Focussing on within K comparisons (i.e. kinsh4), for each K, calculate the proportion of pairwise distances that indicate clonality, as a measure of the rate of clonality:

plot(kinsh3$b3L_s1, kinsh3$kinship_mle)
head(kinsh3,2); dim(kinsh3) # same burn
head(kinsh4,3); dim(kinsh4) # same K

prop_cloneK<-aggregate(clone~K3_s1, data=kinsh4, FUN=function(x) length(x[x==1])/length(x))
colnames(prop_cloneK)<-c("K","prop_clone")
prop_cloneK

# PLOT pair-wise kinship coefficients within sites across fire cateogries and K to visualise variation in the data that might not have been captured by the site-level proportion variable:

dev.new(height=4,width=8,noRStudioGD = T,dpi=100, pointsize=20)
par(mfrow=c(1,2),mar=c(3,3.5,1.5,0.5), mgp=c(2.3,0.8,0))

plot(kinsh3$b3L_s1, kinsh3$kinship_mle, las=1, ylab="Kinship coefficient", xlab="", xaxt="n")
axis(side=1, at=c(1:3), labels = c("","",""), cex.axis=0.7)
axis(side=1, at=c(1,2,3), labels = c("unburnt","burnt","site 11"), tick = F, cex.axis=0.8)
title(xlab="Fire regime", mgp=c(2,1,0))
mtext("(a)",side=3, at=0.5, line=0.3)

plot(as.factor(kinsh4$K3_s1), kinsh4$kinship_mle, pch=20, las=1, ylab="Kinship coefficient", xlab="", xaxt="n")
axis(side=1, at=c(1:3), labels = c(1,2,3))
title(xlab="Genetic cluster (K)", mgp=c(2,1,0))
mtext("(b)",side=3, at=0.5, line=0.3)

# save.image("03_Workspaces/STEP04_divdist_ALL.RData")

# close analyse distance ----

# Analyse distribution of FIS	# ----

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
