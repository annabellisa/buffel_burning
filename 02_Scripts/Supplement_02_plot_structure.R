
# ------------------------------------ #
# ----------- SUPPLEMENT 03  --------- #
# ------------------------------------ #

# Plot structure results:
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

library("RColorBrewer")
# library(rworldmap)
# library("plotrix")

# --- *** DATA SET UP *** --- #

# SET WORKING DIRECTORIES, FILES AND TIDY DATA:

# The project dir is the location of the structure input files:
proj_dir<-"RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt2"
dir(proj_dir)

# The results dir contains the structure results, usually within proj_dir:
res_dir<-paste(proj_dir,"Cenchrus_filt2_results",sep="/")
dir(res_dir)

file_name<-"Cenchrus_filt2"

# this is usually .fam but can also be .str:
str_file<-paste(proj_dir,dir(proj_dir)[grep(".fam",dir(proj_dir))],sep="/")

# str_sites is site data that corresponds to genetic assignment probabilities:
str_sites<-read.table(str_file,colClasses=c(rep("character",2),rep("NULL",dim(read.table(str_file))[2]-2)),header=F)

dat_dir<-"00_Data"

# read main site data file:
sdat<-read.table(paste(dat_dir,"Cenchrus_site_data.txt",sep="/"),header=T)
sdat<-sdat[,c("block","site","year","burn_unburnt","time_since_burn","no_samples","lat","long")]
sdat<-tidy.df(sdat)
sdat<-sdat[order(sdat$block,sdat$site),]
sdat<-tidy.df(sdat)
head(sdat,2); dim(sdat)

#########################################
#  SET K AND ASSIGNMENT PROBABILITIES:  #
#########################################

# Set K and run this whole section

# SET K:
K<-3

# Get assignment probabilities:
# use out_dir or res_dir, depending on where results are:
assig<-read.table(paste(res_dir,dir(res_dir)[grep(paste("\\.",K,".meanQ",sep=""),dir(res_dir))],sep="/"),header=F)

# Combine site_data and assigment probs:
site_assig<-cbind(str_sites,assig)
head(site_assig)
head(sdat,2)

# this is the order for .fam files
colnames(site_assig)<-c("site","indiv",paste("assig",1:K,sep=""))

sdat2<-sdat

# replace the "X" in the site name with "buf" (going to remove this later, but it's needed for merge)
site_assig$site<-paste("buf",substr(site_assig$site,2,nchar(site_assig$site)),sep="")
head(sdat2,2); head(site_assig, 2)

# should be zero:
length(unique(site_assig$site)[which(unique(site_assig$site) %in% as.character(unique(sdat2$site))!=T)])

# all assignment probs should appear in site data (these should all be TRUE)
table(site_assig$site %in% sdat2$site)

# site nine and site 4 were not sequenced, so there are four sites in the site data not in the assignment data (this probably doesn't matter if site_assig is the x data frame in the merge):
table(sdat2$site %in% site_assig$site)

head(site_assig, 2); dim(site_assig)

all_dat<-merge(site_assig,sdat2,by="site",all.x=T,all.y=F)

# order by longitude, with most westerly sites first (on left of plot); these are the smaller numbers in longitude, so order ascending. 

# However, when sorting by longitude, the burn or unburnt site is not consistently to the W, so it make the plots difficult to interpret. Since the block numbers are already sorted W to E, the way around this would be to sort by site, then re-arrange so that block buf11 is first
all_dat<-all_dat[order(all_dat$site,all_dat$indiv),]
all_dat<-tidy.df(all_dat)

all_dat<-all_dat[c(which(as.character(all_dat$block) %in% "buf11"), which(!as.character(all_dat$block) %in% "buf11")),]

all_dat<-tidy.df(all_dat)
head(all_dat,3); dim(all_dat)

# optionally change site names, removing 'buf':
all_dat$site<-substr(all_dat$site,4,nchar(all_dat$site))

# for plotting both on the same panel
# save K3 and re-run for K4
all_dat_K3<-all_dat
all_dat_K4<-all_dat

# Summarise groupings for downstream analysis:

# run all_dat for K3 first:
sum_dat<-all_dat
sum_dat$K3<-apply(all_dat[,grep("assig", colnames(all_dat))],1,function(x) which(x==max(x)))

# then re-run all_dat for K4:
sum_dat$K4<-apply(all_dat[,grep("assig", colnames(all_dat))],1,function(x) which(x==max(x)))

# remove K3 assig cols to avoid confusion:
sum_dat[,grep("assig", colnames(sum_dat))]<-NULL
head(sum_dat,3); dim(sum_dat)

# write.table(sum_dat, file="K_genetic_clusters.txt", sep="\t", quote=F, row.names=F)

#########################################
####  	       BAR PLOTS:    		 ####
#########################################

# The sequential palettes names are:
#  Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd

# The diverging palettes are 
# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

switch.col<-brewer.pal(K,"Accent")

### Order samples West - East:

### MAIN PLOT (K=3):
quartz(title="Fig",width=16,height=4,dpi=80,pointsize=10)
par(mfrow=c(1,1),oma=c(1,0,1,0))
str_plot_V10(K = 3,all_dat_K3,sdat2,las.opt=2,yaxs.loc=-3,cex.axis=0.7,col.pal="switch.col",site.lab="site")
mtext("West - East", side=1, line=2.8, cex=1.8)
par(xpd=NA)
arrows(x0=c(41,52),y0=c(-0.2,-0.2),x1=c(31,62),y1=c(-0.2,-0.2),code=2, length=0.2)
par(xpd=F)

### K=3 and K=4 PLOT, for SI:
quartz(title="Fig",width=16,height=8,dpi=80,pointsize=10)
par(mfrow=c(2,1),oma=c(2,0,1.5,0))
str_plot_V10(K = 3,all_dat_K3,sdat2,las.opt=2,yaxs.loc=-3,cex.axis=0.7,col.pal="switch.col",site.lab="site")
mtext("(a) K = 3", side=3, line=0.5, cex=1.8, at = 0, adj=0)
mtext("West - East", side=1, line=2.8, cex=1.8)
par(xpd=NA)
arrows(x0=c(41,52),y0=c(-0.2,-0.2),x1=c(31,62),y1=c(-0.2,-0.2),code=2, length=0.2)
par(xpd=F)

str_plot_V10(K = 4,all_dat_K4,sdat2,las.opt=2,yaxs.loc=-3,cex.axis=0.7,col.pal="switch.col",site.lab="site")
mtext("(b) K = 4", side=3, line=0.5, cex=1.8, at = 0, adj=0)
mtext("West - East", side=1, line=2.8, cex=1.8)
par(xpd=NA)
arrows(x0=c(41,52),y0=c(-0.2,-0.2),x1=c(31,62),y1=c(-0.2,-0.2),code=2, length=0.2)
par(xpd=F)


### Order samples along dendrogram:

ddir<-"00_Data/Filtered_DartSeq_format"
ddat<-read.table(paste(ddir, "dartseq_filt2.txt", sep="/"), header=T)
ghead(ddat)

# re-do distance matrix on raw data:
clust_dat<-ddat

# update individual names so they don't include the redundant X:
clust_dat$ind<-substr(clust_dat$ind,2,nchar(clust_dat$ind))
rownames(clust_dat)<-clust_dat$ind
clust_dat<-clust_dat[,3:length(clust_dat)]
ghead(clust_dat)

hc_dist<-dist(x = clust_dat, method="euclidean")
euc_clust<-hclust(hc_dist)
str(euc_clust)

# get name label order:
hc1_names<-data.frame(ind=hclust_name_order(euc_clust), hclust_order=1:length(hclust_name_order(euc_clust)))
head(hc1_names)

# for str_plot_V10 the ordering is done outside the function, so we can re-sort and use the same function:
K3_dend_order<-all_dat_K3

# remove X from indiv name:
K3_dend_order$indiv<-substr(K3_dend_order$indiv,2,nchar(K3_dend_order$indiv))
head(K3_dend_order); dim(K3_dend_order)
head(hc1_names,3); dim(hc1_names)

# check all names are present in both data sets:
table(hc1_names$ind %in% K3_dend_order$indiv)
table(K3_dend_order$indiv %in%  hc1_names$ind)

K3_dend_order<-merge(K3_dend_order, hc1_names, by.x = "indiv", by.y = "ind", all.x = T, all.y = F)
K3_dend_order<-K3_dend_order[order(K3_dend_order$hclust_order),]
K3_dend_order<-tidy.df(K3_dend_order)
head(K3_dend_order,3); dim(K3_dend_order)

# plot dendro & structure together:
quartz("",10,4,dpi=120)
par(mfrow=c(2,1),mar=c(0,4,1,0), mgp=c(2,0,-1),oma=c(1,0,1,0))

layout(matrix(c(1,1,1,1,2), 5, 1, byrow = TRUE))
layout.show(2)

plot(euc_clust, cex=0.8, xlab="", main="", cex.lab=1, las=1, ylab="Genetic distance (Euclidean)", sub="")
mtext("A", side=3, line=0.5, adj=0, at=-6.5)

str_plot_V11(K = 3,K3_dend_order,sdat2,las.opt=2,yaxs.loc=-2.5,cex.axis=0.7,col.pal="switch.col",site.lab="")
mtext("B", side=3, line=0.5, adj=0, at=-7.5)



