
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### Diversity & distances
### Author: Annabel Smith & Binyin Di

# Load workspace:
load("../04_workspaces/STEP04_divdist_wksp")
load("03_workspaces/divdist_wksp.RData")
dir("03_workspaces")


# Neutral Dataset: 20159 (20161)
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/NeturalWksp.RData")

# Non Netural Dataset: 3892 (3894)
load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/NonNeturalWksp.RData")

# load functions:
dir()
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))


# Load libraries:
install.packages(c("diveRsity","geosphere","hierfstat","adegenet"),repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com",dependencies=TRUE))

library(diveRsity);library(geosphere);library(hierfstat);library(adegenet)

# Update R
# install.packages("installr", dependencies = TRUE)
# library(installr)
# updateR()
# library(diveRsity) # require package "mnormt", not available in R 3.6, download under R 4.0. 



#########################################
##     GENIND object & site data:      ##
#########################################
{

gp_dir<-"../ANALYSIS_RESULTS/Genepop_DATA_FILES"

#NF
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/NF_Format"

#N-NF
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/N_NF_Format"


dir(gp_dir)
setwd(gp_dir)
getwd()

# Make genind objects:
# ~~
genind_filt1<-read.genepop(file=paste(gp_dir,"Genpop_Diversity_Original.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_filt1


save.image("../04_workspaces/STEP04_divdist_wksp")



} # close genind

# RESULT:
# See parameter files in gp_dir for filters
genind_filt1 # no OG or cultivars
head(sdat,3); dim(sdat)




#########################################
##     Genetic & enviro distances:     ##
#########################################
{

### -- *** CALCULATE FST:

gp_dir<-"../ANALYSIS_RESULTS/Genepop_DATA_FILES"
gp_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/Offline Winter Project/Shared/Genepop_Files/"
gp_dir<-"D:/Onedrive/OneDrive - The University of Queensland/Offline Winter Project/Shared/Genepop_Files"
dir(gp_dir)

# USE ALL POPULATIONS, including outgroups and cultivars (Genpop_Diversity_Original.gen); can subset this later, but we need all the FSTs:

# Get FST:
# <10min
print(Sys.time())
fst<-diffCalc(paste(gp_dir,"Genpop_Diversity_Original.gen",sep="/"),fst=T,pairwise=T) # 4.5 mins
print(Sys.time())
save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp.RData")

head(fst$pairwise$Fst)
head(fst$pairwise$gst)
head(fst$pairwise$Gst)
head(fst$pairwise$GGst)
head(fst$pairwise$D[,1:5])

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

og_lines<-c(grep("OG",fst_df$pop1),grep("OG",fst_df$pop2))
nonog_lines<-which(rownames(fst_df) %in% og_lines==F)

mean(fst_df$fst[og_lines])
mean(fst_df$fst[nonog_lines])

range(fst_df$fst[og_lines],na.rm=T)
range (fst_df$fst[nonog_lines],na.rm=T)

# write.table(fst_df,"fst.txt",row.names=F,quote=F,sep="\t")

### -- *** ADD GEOGRAPHIC DISTANCE:

dat_dir<-"../00_Data"
dat_dir<-"C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/00_Data"
dir(dat_dir)

# load fst data:
pwpop<-read.table(paste(dat_dir,"pw_pop_stats.txt",sep="/"),header=T)
head(pwpop)

# Update site names:
pwpop$pop1<-as.character(pwpop$pop1)
pwpop$pop2<-as.character(pwpop$pop2)
pwpop$pop1[which(pwpop$pop1=="VIR")]<-"VA"
pwpop$pop2[which(pwpop$pop2=="VIR")]<-"VA"
pwpop$pop1<-as.factor(pwpop$pop1)
pwpop$pop2<-as.factor(pwpop$pop2)

# site data:
head(sdat,3)

# Check all sites in data have site data:
unique(c(levels(pwpop$pop1),levels(pwpop$pop2))) %in% sdat$site_code

sll<-sdat[,c("site_code","latitude","longitude")]
head(sll)
m1<-merge(pwpop,sll,by.x="pop1",by.y="site_code",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("latitude","longitude")]<-c("lat1","lon1")
m1<-merge(m1,sll,by.x="pop2",by.y="site_code",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("latitude","longitude")]<-c("lat2","lon2")
m1<-m1[order(m1$pop1,m1$pop2),]
m1<-tidy.df(m1)
m1<-m1[,c(2,1,3:length(m1))]
m1$geog_dist<-distGeo(m1[,c("lon1","lat1")],m1[,c("lon2","lat2")])
m1$lat_dist<-abs(m1$lat1)-abs(m1$lat2)
head(m1)
head(sdat,4)

# Add native distance:
ndis<-sdat[,c("site_code","native")]
head(ndis)
m2<-merge(m1,ndis,by.x="pop1",by.y="site_code",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("native")]<-c("nat1")
m2<-merge(m2,ndis,by.x="pop2",by.y="site_code",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("native")]<-c("nat2")
m2<-m2[order(m2$pop1,m2$pop2),]
m2<-tidy.df(m2)
m2<-m2[,c(2,1,3:length(m2))]
head(m2)

m2$nat_dist<-m2$nat1==m2$nat2
m2$nat_dist<-ifelse(m2$nat_dist==T,0,1)
check.rows(m2)

# write.table(m2,"m2.txt",row.names=F,quote=F,sep="\t")

} # close distances

#########################################
##  	  Genetic diversity:  	       ##
#########################################
{

## -- ** POPULATION LEVEL GENETIC DIVERSITY:

# GenFilter1

# Calculate genetic diversity per population in hierfstat (conservative dataset):
gendiv_filt1 <- basic.stats(genind_filt1, diploid = TRUE, digits = 2)
str(gendiv_filt1) # a few mins
head(gendiv_filt1$Ho)
tail(gendiv_filt1$Ho)
# save.image("../04_workspaces/STEP04_divdist_wksp")

gd_filt1<-data.frame(site=names(apply(gendiv_filt1$Ho,2,mean,na.rm=T)),max_n=apply(gendiv_filt1$n.ind.samp,2,max,na.rm=T),Ho=apply(gendiv_filt1$Ho,2,mean,na.rm=T),He=apply(gendiv_filt1$Hs,2,mean,na.rm=T),Fis=apply(gendiv_filt1$Fis,2,mean,na.rm=T))
gd_filt1<-tidy.df(gd_filt1)
head(gd_filt1); dim(gd_filt1)

gd_filt1$site<-substr(gd_filt1$site,1,nchar(as.character(gd_filt1$site))-1)

gd_filt1$site[grep("_",substr(gd_filt1$site,nchar(gd_filt1$site),nchar(gd_filt1$site)))]<-substr(gd_filt1$site[grep("_",substr(gd_filt1$site,nchar(gd_filt1$site),nchar(gd_filt1$site)))],1,nchar(gd_filt1$site[grep("_",substr(gd_filt1$site,nchar(gd_filt1$site),nchar(gd_filt1$site)))])-1)

# write.table(gd_filt1,"gd_filt1.txt",sep="\t",row.names=F,quote=F)

## -- ** ALLELIC RICHNESS:

head(sdat,3); dim(sdat)
genind_filt1

library("hierfstat")

print(Sys.time())
ar_default_rd<- allelic.richness(genind_filt1, diploid = TRUE) # 2 mins
print(Sys.time())

str(ar_default_rd)
head(ar_default_rd$Ar,2)
ar_default_rd$Ar[,1]

save.image("../03_Workspaces/divdist_wksp.RData")

# NF
save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp_NF.RData")

# N-NF
save.image("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp_N_NF.RData")



load("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/03_Workspaces/divdist_wksp_NF.RData")


setwd("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/NF_Format")
setwd("C:/Users/s4467005/OneDrive - The University of Queensland/GitHub/Binyin_Winter/RESULTS/N_NF_Format")

# setwd("RESULTS/NF_Format")

# ar_dat<-read.table("C:/Users/s4467005/OneDrive - The University of Queensland/Offline Winter Project/Shared/Genepop_Files/joint_dataframe.txt", header = TRUE)
# head(ar_dat)


# Summarise

# Netural Theory of Molecular Evolution

# Neutral:
head(ar_default_rd$Ar,2)
ar_default_rd$min.all

ar_res_rd<-data.frame(
site=levels(genind_filt1@pop),
ar_default_rd=apply(ar_default_rd$Ar,2,mean,na.rm=T))
head(ar_res_rd); dim(ar_res_rd)

# combine allelic richness into table with other genetic diversity data:
head(gd_filt1); dim(gd_filt1)
head(ar_res_rd); dim(ar_res_rd)
# fix up the names
gd_filt1$site_code<-paste(gd_filt1$site, gd_filt1$max_n, sep = "")
library(tidyverse)
gd_filt1 %>% full_join(ar_res_rd, by = c("site_code" = "site"))
# Jump to Binding





# Non-neutral
# ar_default_adapt<-ar_default_rd


head(ar_default_adapt$Ar,2)
ar_default_adapt$min.all

ar_res_adapt<-data.frame(
site=levels(genind_filt1@pop),
ar_default_adapt=apply(ar_default_adapt$Ar,2,mean,na.rm=T))

test_df<-cbind(ar_res_adapt,ar_adapt=ar_res_adapt[,2])
head(test_df)
plot(test_df$ar_default_adapt, test_df$ar_adapt)

# write.table(ar_res_adapt,"ar_res_adapt.txt",sep="\t",row.names=F,quote=F)

save.image("../04_workspaces/STEP04_divdist_wksp")


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
  geom_jitter()+
  geom_errorbar(aes(ymin = lci, ymax = uci, position = "dodge"), width = 0.2)
  
  
  
  #stat_summary(geom = "bar", fun.y = fit, position = "dodge") +
  #stat_summary(geom = "errorbar", fun.data = se, position = "dodge")



# calculating std errors
# install.packages("plotrix") # Install plotrix R package
# library("plotrix")  
# std.error()



## -- ** Admixture Diversity Score:

### ---- level of admixture per population

# Working from equations in Harismendy et al. 2019, Journal of the American Medical Informatics Association, 26(5), 457–46

# Structure summary objects loaded from: Supplement_03_plot_structure.R

# 491 individuals including cultivars and outgroups
head(all_dat); dim(all_dat)
head(site_assig); dim(all_dat)

# Associated site data (61 pops including cultivars and outgroups):
head(sdat,2); dim(sdat)
head(sdat2); dim(sdat)

# The same data but with outgroups and cultivars removed (53 pops, 454 individuals):
head(sdpie,2); dim(sdpie)
head(cluster.data); dim(cluster.data)

# Aggregated data to give the assignment proportion for each population (53 pops):
head(r1); dim(r1)
head(agg_dat); dim(agg_dat)

# Use individual level data for the analysis (6 clusters with assignment probabilities for 454 individuals):
ind6<-cluster.data[,grep("assig",colnames(cluster.data))]
head(ind6); dim(ind6)

# The Diversity Score from Harismendy et al. We don't need eqn 1 (the cumulative admixture fraction), because ours already sum to 1. 

# The diversity score is the same as the shannon diversity index but with a scaling factor to account for the number of clusters (Hmax). That is, it makes the diversity score relative to complete evenness. 

DS<-function(x,K) {
	Hmax<-K*((1/K)*log(1/K))
	(-sum(x*log(x)))/-Hmax
	}

# SET K:
K<-6

# Calculate admixture diversity for all individuals:
ds_div<-apply(ind6, 1, DS, K=6)
range(ds_div)

ind6_dv<-data.frame(indiv=cluster.data$indiv, admix_div=ds_div)
ind_site<-cluster.data[,c("site_code","indiv","native","region","country","latitude","longitude")]
ind6_dv<-merge(ind6_dv, ind_site, by="indiv", all.x=T, all.y=F)
head(ind6_dv); dim(ind6_dv)

# Summarise site-level admixture diversity:
addiv_site<-aggregate(admix_div~site_code,FUN=mean,data=ind6_dv)
ex_dat<-ind6_dv[,c("site_code","native","region","country","latitude","longitude")]
ex_dat<-ex_dat[-which(duplicated(ex_dat$site_code)),]
ex_dat<-tidy.df(ex_dat)
head(ex_dat); dim(ex_dat)
addiv_site<-merge(addiv_site, ex_dat, by="site_code", all.x=T, all.y=F)
head(addiv_site); dim(addiv_site)

# Allelic richness:
ar_dat<-read.table(paste("../01_data/gen_div.txt",sep=""),header=T)

# Check all are true:
addiv_site$site_code %in% ar_dat$site

# Add admixture diversity
ar_dat<-merge(ar_dat, addiv_site, by.x="site", by.y="site_code", all.x=T, all.y=F)
head(ar_dat); dim(ar_dat)

# write.table(ar_dat, file="ardat.txt", quote=F, row.names=F, sep="\t")

# Is there a difference in admixture between the native and non-native ranges?
library(lmerTest)
library(lme4)
library(AICcmodavg)

head(ind6_dv); dim(ind6_dv)

dv_ind<-lmer(admix_div~native+(1|site_code), data=ind6_dv)
summary(dv_ind)
admix_df<-predictSE(dv_ind, newdata=data.frame(native=c("native","non_native")), se.fit=T)
admix_pred<-data.frame(native=c("native","non_native"),fit=admix_df$fit, se=admix_df$se.fit)
admix_pred$lci<-admix_pred$fit-(1.96*admix_pred$se)
admix_pred$uci<-admix_pred$fit+(1.96*admix_pred$se)

### PLOT for SI:
{

quartz("",4,4)
par(mfrow=c(2,2),mar=c(4,4,1,1), mgp=c(2.5,0.8,0))

plot(1:2,admix_pred$fit, xlim=c(0.5,2.5), ylim=c(min(admix_pred$lci),max(admix_pred$uci)), pch=20, xaxt="n", xlab="", ylab="admixture diversity", las=1,cex.axis=0.85,cex.lab=0.85)
arrows(1:2, admix_pred$lci, 1:2, admix_pred$uci, length=0.02, code=3, angle=90)
axis(side=1, at=c(1,2), labels=c("native","non-native"), cex.axis=0.85)
mtext(bquote("range "~italic("p = ")~.(round(summary(dv_ind)$coefficients[2,5],3))), side=3, line=0,cex=0.7, adj=0)
mtext("(a)",side=3, line=0, adj=0, at=-0.6, cex=0.8)

blankplot()

plot(ar_dat$ar[ar_dat$native=="native"], ar_dat$admix_div[ar_dat$native=="native"], pch=20, col="black", cex.axis=0.8, las=1,ylab="admixture diversity",xlab="", cex.lab=0.85)
points(ar_dat$ar[ar_dat$native=="non_native"], ar_dat$admix_div[ar_dat$native=="non_native"], pch=20, col="red", cex.axis=0.8, las=1,cex.lab=0.8)
title(xlab="allelic richness (neutral)", cex.lab=0.85, mgp=c(2,1,0))
mtext("(b)",side=3, line=0, adj=0, at=1.119, cex=0.8)

# Is there are r.ship between ar and admixture diversity?
ar_admix<-lm(admix_div~ar*native,data=ar_dat)
mtext(bquote("range x ar"~italic("p = ")~.(round(summary(ar_admix)$coefficients[4,4],3))), side=3, line=0,cex=0.7, adj=0)

plot(ar_dat$ar_adapt[ar_dat$native=="native"], ar_dat$admix_div[ar_dat$native=="native"], pch=20, col="black", cex.axis=0.8, las=1,ylab="admixture diversity",xlab="", cex.lab=0.8)
points(ar_dat$ar_adapt[ar_dat$native=="non_native"], ar_dat$admix_div[ar_dat$native=="non_native"], pch=20, col="red", cex.axis=0.8, las=1,cex.lab=0.8)
title(xlab="allelic richness (adaptive)", cex.lab=0.85, mgp=c(2,1,0))
mtext("(c)",side=3, line=0, adj=0, at=1.121, cex=0.8)

# Is there are r.ship between ar_adapt and admixture diversity?
adapt_admix<-lm(admix_div~ar_adapt*native,data=ar_dat)
mtext(bquote("range x ar"~italic("p = ")~.(round(summary(adapt_admix)$coefficients[4,4],3))), side=3, line=0,cex=0.7, adj=0)

} # close plot

} # close genetic diversity

























