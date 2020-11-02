
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### Diversity & distances
### Author: Annabel Smith & Binyin Di

# AS Load workspace:
load("03_workspaces/divdist_ALL.RData")

# BD Neutral Dataset: 20159 (20161)
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/NeturalWksp.RData")

# Non Netural Dataset: 3892 (3894)
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/NonNeturalWksp.RData")

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# Load libraries:
library(diveRsity);library(geosphere);library(hierfstat);library(adegenet); library(ecodist); library(AICcmodavg)

#  GENIND object & site data:    	# ----

# Load Genepop files
gp_dir<-"00_Data/Genepop_Files"
dir(gp_dir)

# NF
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/NF_Format"

# N-NF
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/N_NF_Format"
dir(gp_dir)

# Make genind objects:

# ~~ Neutral (~ 5 min for 20159 loci)
genind_neutral<-read.genepop(file=paste(gp_dir,"Genepop_Neutral_filt2.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_neutral

# ~~ Non-Neutral (~ 10 seconds for 3892 loci)
genind_nonneutral<-read.genepop(file=paste(gp_dir,"Genepop_NonNeutral.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_nonneutral

# Load site data:
sdat<-read.table(paste("00_data/Cenchrus_site_data.txt",sep=""),header=T)
sdat<-sdat[sdat$sequenced==1,]
sdat<-tidy.df(sdat)
head(sdat); dim(sdat)

# save.image("03_Workspaces/divdist_ALL.RData")

# close genind object ----

# RESULT:
# See parameter files in gp_dir for filters
genind_neutral # all filters + neutral markers only (20159 loci)
genind_nonneutral # all filters + non-neutral markers only  (3892 loci)
head(sdat,3); dim(sdat)

#  Genetic & enviro distances:    	# ----
  
### -- *** CALCULATE FST:

# Load Genepop files
gp_dir<-"00_Data/Genepop_Files"
dir(gp_dir)
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/Offline Winter Project/Shared/Genepop_Files/"
gp_dir<-"D:/Onedrive/OneDrive - The University of Queensland/Offline Winter Project/Shared/Genepop_Files"
dir(gp_dir)

# Get FST:
# ~ 2 min for Cenchrus 20159 loci
print(Sys.time())
fst<-diffCalc(paste(gp_dir,"Genepop_Neutral_filt2.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())
# save.image("03_Workspaces/divdist_ALL.RData")
# save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp.RData")

head(fst$pairwise$Fst)
head(fst$pairwise$gst)
head(fst$pairwise$Gst)
head(fst$pairwise$GGst)
head(fst$pairwise$D[,1:5])

# Convert square matrix to column matrix:
fstn<-rownames(fst$pairwise$Fst)
fstn<-substr(fstn,1,(nchar(fstn)-2))

ind_endnum<-which(unlist(gregexpr("[0-9]",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endnum<-fstn[ind_endnum]
fstn[ind_endnum]<-substr(with_endnum,1,nchar(with_endnum)-1)

ind_endund<-which(unlist(gregexpr("_",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endund<-fstn[ind_endund]
fstn[ind_endund]<-substr(with_endund,1,nchar(with_endund)-1)

fst_df<-data.frame(pop1=combn(fstn,2)[1,],pop2=combn(fstn,2)[2,],fst=fst$pairwise$Fst[lower.tri(fst$pairwise$Fst)],gst=fst$pairwise$gst[lower.tri(fst$pairwise$gst)],Gst=fst$pairwise$Gst[lower.tri(fst$pairwise$Gst)],GGst=fst$pairwise$GGst[lower.tri(fst$pairwise$GGst)],D=fst$pairwise$D[lower.tri(fst$pairwise$D)])
head(fst_df)

mean(fst_df$fst)
range(fst_df$fst,na.rm=T)

# write.table(fst_df,"fst_all_sites.txt",row.names=F,quote=F,sep="\t")

### -- *** ADD GEOGRAPHIC DISTANCE:

dat_dir<-"RESULTS/Diversity_and_Distance/FST"
dat_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/00_Data"
dir(dat_dir)

# load fst data:
pwpop<-read.table(paste(dat_dir,"fst_all_sites.txt",sep="/"),header=T)
head(pwpop)

# site data:

# add site code to match fst data:
sdat$pop<-sdat$site
sdat$pop<-paste("X",substr(x=sdat$pop,start = 4,stop = nchar(as.character(sdat$pop))), sep="")
head(sdat,3)

# Check all sites in data have site data:
unique(c(levels(pwpop$pop1),levels(pwpop$pop2))) %in% sdat$pop

sll<-sdat[,c("pop","lat","long")]
head(sll)
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

# write.table(m2,"m2.txt",row.names=F,quote=F,sep="\t")
# save.image("03_Workspaces/divdist_ALL.RData")

# Analyse FST:

dir(dat_dir)
# load fst data:
fst_all<-read.table(paste(dat_dir,"fst_and_distances_all_sites.txt",sep="/"),header=T)
head(fst_all,2)

# mantel test:
mant1<-mantel(formula = fst~geog_dist, data = fst_all)
mant1
mant2<-mantel(formula = fst~geog_dist+same_block, data = fst_all)
mant2

# plot FST:
quartz("",6,4,dpi=100)
par(mar=c(4,4,2,1), mgp=c(2.5,1,0))
plot(fst_all$geog_dist, fst_all$fst, pch=20, xlab="Geographic distance (m)", ylab="FST", las=1)
# pval2 = one-tailed p-value (null hypothesis: r >= 0).
mtext(paste("mean FST = ",round(mean(fst_all$fst),2),"; mantel r = ",round(mant1[1],2),"; p = ",round(mant1[3],2), sep=""), adj=0)

# close distances ----

#  Calculate genetic diversity:    	# ----

## -- ** POPULATION LEVEL GENETIC DIVERSITY:

genind_neutral # all filters + neutral markers only (20159 loci)
genind_nonneutral # all filters + non-neutral markers only  (3892 loci)
head(sdat,3); dim(sdat)

# Calculate genetic diversity per population in hierfstat:

# neutral:
gendiv_neutral <- basic.stats(genind_neutral, diploid = TRUE, digits = 2)
str(gendiv_neutral) # 1 min
head(gendiv_neutral$Ho)
tail(gendiv_neutral$Ho)
# save.image("03_workspaces/divdist_ALL.RData")

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

# save.image("03_workspaces/divdist_ALL.RData")

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

# save.image("03_workspaces/divdist_ALL.RData")

# NF
save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp_NF.RData")

# N-NF
save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp_N_NF.RData")

# Summarise

# Netural Theory of Molecular Evolution

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
# save.image("03_workspaces/divdist_ALL.RData")

# Binding
# dplyr from tidyverse to bind columns 

library(tidyverse)
sdatcoord<- sdat %>%
  select(burn_unburnt,lat,long)
colnames(sdatcoord)[1]<-"treatment"

# Neutral 
ar<-data.frame(ar= ar_res_rd$ar_default_rd)
joint_data_frame<-bind_cols(gd_filt1,sdatcoord, ar)
# write.table(joint_data_frame,"joint_netural_dataset.txt",sep="\t",row.names=F,quote=F)

# NonNetural
ar<-data.frame(ar = ar_res_adapt$ar_default_adapt)
joint_data_frame<-bind_cols(gd_filt1,sdatcoord,ar)
# write.table(joint_data_frame,"joint_nonnetural_dataset.txt",sep="\t",row.names=F,quote=F)

# close calc genetic diversity ----

#  ANALYSE Genetic diversity:    	# ----

gd_all<-read.table("RESULTS/Diversity_and_Distance/Genetic_Diversity/Genetic_Diversity_ALL.txt", header=T)
head(gd_all,3)

# simplify site data
sdt<-sdat[,c("block","site","pop","burn_unburnt","no_samples","lat","long")]
sdt$site<-paste("X",substr(sdt$site, 4, nchar(as.character(sdt$site))),sep="")

# merge with genetic diversity data:
sdt$site %in% gd_all$site
gd_all$site %in% sdt$site
gd_all<-merge(gd_all, sdt, by="site", all.x=T, all.y=F)

# add second treatment variable:
gd_all$treatment<-gd_all$burn_unburnt
gd_all$burn_unburnt<-NULL
gd_all$trt<-as.character(gd_all$treatment)
gd_all$trt[grep("b",gd_all$trt)]<-"b1"
gd_all$trt[grep("X11",gd_all$site)]<-"b2"
gd_all$trt<-factor(gd_all$trt, levels=c("u","b1","b2"))
head(gd_all,3); dim(gd_all)

# neutral models:
mod1.a<-lm(ar_neutral~1, data=gd_all)
mod2.a<-lm(ar_neutral~trt, data=gd_all)
mod3.a<-lm(ar_neutral~trt+long, data=gd_all)
mod4.a<-lm(ar_neutral~trt*long, data=gd_all)
AICc(mod1.a); AICc(mod2.a); AICc(mod3.a); AICc(mod4.a)
summary(mod2.a); anova(mod2.a)

# non-neutral models:
mod1.b<-lm(ar_nonneutral~1, data=gd_all)
mod2.b<-lm(ar_nonneutral~trt, data=gd_all)
mod3.b<-lm(ar_nonneutral~trt+long, data=gd_all)
mod4.b<-lm(ar_nonneutral~trt*long, data=gd_all)
AICc(mod1.b); AICc(mod2.b); AICc(mod3.b); AICc(mod4.b)
summary(mod2.b); anova(mod2.b)


# add genetic clusters from STRUCTURE analysis (this is only useful for the individual level analysis:
gclust<-read.table("00_Data/K_genetic_clusters.txt", header=T)
gclust$site %in% sdt$site
head(gclust,3); dim(gclust)


# BD:

joint_data_frame$trt<-as.character(joint_data_frame$treatment)

joint_data_frame$trt[grep("X11",joint_data_frame$site)]<-"b2"
joint_data_frame$trt<-as.factor(joint_data_frame$trt)
joint_data_frame$trt<-relevel(joint_data_frame$trt,"u")
levels(joint_data_frame$trt)

# Models
mod1<-lm(ar~treatment , data = joint_data_frame)
summary(mod1)
anova(mod1)

mod1.1<-lm(ar~trt, data = joint_data_frame)
summary(mod1.1)
anova(mod1.1)

plot(ar~trt, data = joint_data_frame)

mod2<-lm(He~treatment + long, data =joint_data_frame)
summary(mod2)
anova(mod2)


mod3<-lm(ar~treatment + long, data = joint_data_frame[-which(joint_data_frame$site%in%c("X11b1_0", "X11b2_0", "X11b3_0")),])
summary(mod3)
anova(mod3)



plot(as.factor(joint_data_frame$treatment), joint_data_frame$ar)



# joint_data_frame$treatment[which(joint_data_frame$site%in%c("X11b1_0", "X11b2_0", "X11b3_0"))]


new.dataframe<-data.frame(trt = factor(c("u", "b", "b2"),levels = c("u", "b", "b2")))



p1<-predict(mod1.1, new.dataframe, se.fit = TRUE)

p1

new.dataframe$se<- p1$se.fit
new.dataframe$fit<-p1$fit
print(new.dataframe)



# upper & lower confidence interval

# Question: why are the differences? 
p1$uci<-p1$fit+(1.96*p1$se.fit)
p1$lci<-p1$fit-(1.96*p1$se.fit)

# install.packages("Rmisc", dependencies = TRUE)
# library(Rmisc)
# p1$uci<-CI(p1$fit,ci=0.95)
# p1$lci<-CI(p1$fit,ci=0.95)


new.dataframe$uci<-p1$uci
new.dataframe$lci<-p1$lci

library("tidyverse")
ggplot(data = new.dataframe, mapping = aes(x = trt,
                                           y = fit)) +
  ylab("estimated allelic richness")+
  geom_point()+
  geom_errorbar(aes(ymin = lci, ymax = uci, position = "dodge"), width = 0.2)

# close analyse genetic diversity ----



