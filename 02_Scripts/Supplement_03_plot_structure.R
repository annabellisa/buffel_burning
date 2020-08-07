
# ------------------------------------ #
# ----------- SUPPLEMENT 03  --------- #
# ------------------------------------ #

# Plot structure results:
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

library("RColorBrewer")
library(rworldmap)
library("plotrix")

# --- *** DATA SET UP *** --- #

# SET WORKING DIRECTORIES, FILES AND TIDY DATA:

# The project dir is the location of the structure input files:
proj_dir<-"RESULTS/STRUCTURE/STRUCTURE_DIR/Cenchrus_filt1"
dir(proj_dir)

# The results dir contains the structure results, usually within proj_dir:
res_dir<-paste(proj_dir,"Cenchrus_filt1_results",sep="/")
dir(res_dir)

file_name<-"Cenchrus_filt1"

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
# don't need this for Cenchrus: site_assig$V1<-toupper(site_assig$V1)

head(site_assig)
head(sdat,2)

# this is the order for .fam files
colnames(site_assig)<-c("site","indiv",paste("assig",1:K,sep=""))

sdat2<-sdat

# should be zero:
length(unique(site_assig$site)[which(unique(site_assig$site) %in% as.character(unique(sdat2$site))!=T)])

# replace the "X" in the site name with "buf"
site_assig$site<-paste("buf",substr(site_assig$site,2,nchar(site_assig$site)),sep="")
head(sdat2,2); head(site_assig, 2)

# all assignment probs should appear in site data (these should all be TRUE)
table(site_assig$site %in% sdat2$site)

# site nine and site 4 were not sequenced, so there are four sites in the site data not in the assignment data (this probably doesn't matter if site_assig is the x data frame in the merge):
table(sdat2$site %in% site_assig$site)

head(site_assig, 2); dim(site_assig)

all_dat<-merge(site_assig,sdat2,by="site",all.x=T,all.y=F)
all_dat<-all_dat[order(all_dat$block,all_dat$site,all_dat$indiv),]
all_dat<-tidy.df(all_dat)
head(all_dat,3); dim(all_dat)

#########################################
####  	       BAR PLOTS:    		 ####
#########################################

# Single barplot (for main document):

# The sequential palettes names are:
#  Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd

# The diverging palettes are 
# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

switch.col<-brewer.pal(K,"Accent")

# Test colours:
quartz(title="Fig",width=4,height=4,dpi=80,pointsize=10)
barplot(matrix(data=6,nrow=6,ncol=1),col=brewer.pal(K,"Accent"),yaxt="n")
text(1,seq(1,35,length.out=6),c("non-native","cultivar","outgroup","nth central eur","atlantic","greece"))
text(0.5,seq(1,35,length.out=6),switch.col)

### MAIN PLOT:
quartz(title="Fig",width=16,height=4,dpi=80,pointsize=10)
par(mfrow=c(1,1),oma=c(3,0,2,0))
str_plot_V10(K,all_dat,sdat2,las.opt=2,yaxs.loc=-3,cex.axis=0.7,col.pal="switch.col",site.lab="indiv")































