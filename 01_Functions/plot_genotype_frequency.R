
# Author: Annabel Smith & Di Binyin

# Plot genotype frequencies for outlier loci by longitude and burnt/unburnt

# Modified by BD from AS's Plantago scripts for use on Cenchrus data: 

library(fields)

# There are two functions:

# plot geno frequencies by location (modified by BD, working on Windows only):

plot_freq_location<-function(loci,genotype_data,site_data,out.dir,number_to_plot){

# loci = a character vector of loci to plot
# gentoype_data = genetic data from which to calculate genotype frequencies
  # site_data = site data containing the variables for ordering (e.g. longitude)
  # out.dir = directory to put plots
# number_to_plot = the number of loci to plot, from 1:number_to_plot, in case you don't want to plot all of them

for(i in 1:number_to_plot){
  
loc.thisrun<-loc.toanalyse[i]
dat.thisrun<-data.frame(site=gt_data$site,locus=factor(gt_data[,loc.thisrun],levels=c(0,1,2)))
head(dat.thisrun)

df.thisrun<-as.data.frame.matrix(table(dat.thisrun))
head(df.thisrun)

prop.df<-as.data.frame.matrix(t(apply(df.thisrun,1,function(y) unlist(lapply(y,function(x) x/sum(y))))))
prop.df<-data.frame(site=rownames(prop.df),prop.df)
colnames(prop.df)<-c("site","0","1","2")
prop.df<-tidy.df(prop.df)
head(prop.df)

plr2<-site_data[,c("site_code","burn_unburnt","lat","long")] 
if(length(which(duplicated(plr2$site_code)))>0) plr2<-plr2[-which(duplicated(plr2$site_code)),]
prop.df<-merge(prop.df,plr2,by.x="site",by.y="site_code",all.x=T,all.y=F)

# Order by burnt/unburnt and longitude:
prop.df<-prop.df[order(prop.df$burn_unburnt, prop.df$long),] # rewriten by BD
prop.df<-tidy.df(prop.df)
head(prop.df)

jpeg(file=paste(out.dir,paste(loc.thisrun,".jpg",sep=""),sep="/"), width = 800, height = 200, units = "px",pointsize=10, res = NA)

# dev.new("",8,2,dpi=120,pointsize=10)
# quartz("",8,2,dpi=120,pointsize=10)
par(mfrow=c(1,1),mar=c(4,5,2,1.5),oma=c(0,0,0,0))

bk <- c(-100,seq(0,100,by=10))
mc<-c("grey50",colorRampPalette(colors = c("white","blue"))(length(bk)-2))
library("fields")
image.thisrun<-image.plot(as.matrix(prop.df[,which(colnames(prop.df) %in% c("0","1","2"))][,3:1]),xaxt="n",yaxt="n",col=mc)

seq.thisrun<-seq(par("usr")[1],par("usr")[2],length.out=nrow(prop.df)+1)
incr.thisrun<-(seq.thisrun[2]-seq.thisrun[1])/2

axis(side=1,at=seq.thisrun,labels=F,las=2,cex.axis=1)
axis(side=1,at=(seq.thisrun+incr.thisrun)[1:(length(seq.thisrun)-1)],labels=prop.df$site,las=2,cex.axis=0.9,tick=F)
axis(side=2,at=c(0,0.5,1),labels=c("het","hom_a","hom_b"),tick=F,las=2)
axis(side=2,at=c(-0.25,0.25,0.75,1.25),labels=F)
arrows(-1,c(0.25,0.75),2,c(0.25,0.75),code=3,length=0)
arrows(seq.thisrun,-1,seq.thisrun,2,col="grey40")

mtext(loc.thisrun,side=4,line=6.5)
dev.off()

} # close i

} # close plot_freq_location

# plot geno frequencies by longitude (modified by BD, working on Windows only):

plot_freq_long<-function(loci,genotype_data,site_data,out.dir,number_to_plot){
  
  # loci = a character vector of loci to plot
  # gentoype_data = genetic data from which to calculate genotype frequencies
  # site_data = site data containing the variables for ordering (e.g. longitude)
  # out.dir = directory to put plots
  # number_to_plot = the number of loci to plot, from 1:number_to_plot, in case you don't want to plot all of them
  
  for(i in 1:number_to_plot){
    
    loc.thisrun<-loc.toanalyse[i]
    dat.thisrun<-data.frame(site=gt_data$site,locus=factor(gt_data[,loc.thisrun],levels=c(0,1,2)))
    head(dat.thisrun)
    
    df.thisrun<-as.data.frame.matrix(table(dat.thisrun))
    head(df.thisrun)
    
    prop.df<-as.data.frame.matrix(t(apply(df.thisrun,1,function(y) unlist(lapply(y,function(x) x/sum(y))))))
    prop.df<-data.frame(site=rownames(prop.df),prop.df)
    colnames(prop.df)<-c("site","0","1","2")
    prop.df<-tidy.df(prop.df)
    head(prop.df)
    plr2<-site_data[,c("site_code","burn_unburnt","lat","long")] 
    if(length(which(duplicated(plr2$site_code)))>0) plr2<-plr2[-which(duplicated(plr2$site_code)),]
    prop.df<-merge(prop.df,plr2,by.x="site",by.y="site_code",all.x=T,all.y=F)

    # Order by longitude, then burnt/ unburnt:
    prop.df<-prop.df[order(prop.df$long, prop.df$burn_unburnt),] # rewriten by BD
    prop.df<-tidy.df(prop.df)
    head(prop.df)
    
    jpeg(file=paste(out.dir,paste(loc.thisrun,".jpg",sep=""),sep="/"), width = 800, height = 200, units = "px",pointsize=10, res = NA)
    
    # dev.new("",8,2,dpi=120,pointsize=10)
    # quartz("",8,2,dpi=120,pointsize=10)
    par(mfrow=c(1,1),mar=c(4,5,2,1.5),oma=c(0,0,0,0))
    
    bk <- c(-100,seq(0,100,by=10))
    mc<-c("grey50",colorRampPalette(colors = c("white","blue"))(length(bk)-2))
    library("fields")
    image.thisrun<-image.plot(as.matrix(prop.df[,which(colnames(prop.df) %in% c("0","1","2"))][,3:1]),xaxt="n",yaxt="n",col=mc)
    
    seq.thisrun<-seq(par("usr")[1],par("usr")[2],length.out=nrow(prop.df)+1)
    incr.thisrun<-(seq.thisrun[2]-seq.thisrun[1])/2
    
    axis(side=1,at=seq.thisrun,labels=F,las=2,cex.axis=1)
    axis(side=1,at=(seq.thisrun+incr.thisrun)[1:(length(seq.thisrun)-1)],labels=prop.df$site,las=2,cex.axis=0.9,tick=F)
    axis(side=2,at=c(0,0.5,1),labels=c("het","hom_a","hom_b"),tick=F,las=2)
    axis(side=2,at=c(-0.25,0.25,0.75,1.25),labels=F)
    arrows(-1,c(0.25,0.75),2,c(0.25,0.75),code=3,length=0)
    arrows(seq.thisrun,-1,seq.thisrun,2,col="grey40")
    
    mtext(loc.thisrun,side=4,line=6.5)
    dev.off()
    
  } # close i
  
} # close plot_freq_location


