
# Author: Annabel Smith

# Plot genotype frequencies for outlier loci by location (country / region and latitude), environment or principal components

library(fields)

# There are three functions:

# plot_freq_pc plots allele frequencies against principal components. It is fairly specific to LFMM as it plots only the sites used in the analysis for which PCs were obtained. 
# plot_freq_environ and plot_freq_location are more general and can be used to plot all of the sites which have environmental data and latitude, respectively, in the main site data

plot_freq_pc<-function(loci,genotype_data,site_data,out.dir,vector_to_plot,plot_by,plot_order){

# loci = a character vector of loci to plot
# gentoype_data = genetic data from which to calculate genotype frequencies
# site_data = site data containing the following variables: c("site_code","native","country","region","latitude")
# out.dir = directory to put PDF plots
# vector_to_plot = a numeric vector to select which loci to plot
# plot_by = what pc would you like to order sites by? Should be paste("Comp.",n) for n number of pcs
# plot_order = descending or ascending?

for(i in vector_to_plot){

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

# In this case we merge the allele frequencies with the PCS generated for lfmm and make all.x=F and all.y=T:
prop.df<-merge(prop.df,lfpcs,by.x="site",by.y="site",all.x=F,all.y=T)

# Then add the additional columns and re-merge. I've kept the other environ variables so this could be updated to make it more general. 
plr2<-site_data[,c("site_code","native","region","latitude","Altitude","mt","st","ap","sp","mm","sm")]
prop.df<-merge(prop.df,plr2,by.x="site",by.y="site_code",all.x=T,all.y=F)

order_col<-which(colnames(prop.df)==plot_by)

# Order by native / non_native and relevan:
if(plot_order=="descending") prop.df<-prop.df[order(prop.df[,which(colnames(prop.df)=="native")],prop.df[,which(colnames(prop.df)=="region")],-prop.df[,order_col]),] else prop.df<-prop.df[order(prop.df[,which(colnames(prop.df)=="native")],prop.df[,which(colnames(prop.df)=="region")],prop.df[,order_col]),]
prop.df<-tidy.df(prop.df)
head(prop.df)

# quartz("",8,2,dpi=120,pointsize=10,file=paste(out.dir,paste(plot_by,"_",plot_order,"_",loc.thisrun,".pdf",sep=""),sep="/"),type="pdf")
# quartz("",8,2,dpi=120,pointsize=10)
par(mfrow=c(1,1),mar=c(4,5,2,1.5),oma=c(0,0,0,0))

bk <- c(-100,seq(0,100,by=10))
mc<-c("grey50",colorRampPalette(colors = c("white","blue"))(length(bk)-2))
image.thisrun<-image.plot(as.matrix(prop.df[,which(colnames(prop.df) %in% c("0","1","2"))][,3:1]),xaxt="n",yaxt="n",col=mc)

seq.thisrun<-seq(par("usr")[1],par("usr")[2],length.out=nrow(prop.df)+1)
incr.thisrun<-(seq.thisrun[2]-seq.thisrun[1])/2

axis(side=1,at=seq.thisrun,labels=F,las=2,cex.axis=1)
axis(side=1,at=(seq.thisrun+incr.thisrun)[1:(length(seq.thisrun)-1)],labels=prop.df$site,las=2,cex.axis=0.9,tick=F)
axis(side=2,at=c(0,0.5,1),labels=c("het","hom_a","hom_b"),tick=F,las=2)
axis(side=2,at=c(-0.25,0.25,0.75,1.25),labels=F)
arrows(-1,c(0.25,0.75),2,c(0.25,0.75),code=3,length=0)
arrows(seq.thisrun,-1,seq.thisrun,2,col="grey40")

axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region))],labels=F)

axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region[-which(prop.df$region=="Africa")]))]+(incr.thisrun*3),labels=as.character(unique(prop.df$region))[-which(as.character(unique(prop.df$region))=="Africa")],line=-1,tick=F,hadj=0,cex.axis=1)

axis(side=3,at=seq.thisrun[which(prop.df$region=="Africa")]+incr.thisrun,labels=as.character(unique(prop.df$region))[which(as.character(unique(prop.df$region))=="Africa")],line=-0.7,tick=F,hadj=0,cex.axis=0.65,las=2)

arrows(seq.thisrun[which(!duplicated(prop.df$region))],-1,seq.thisrun[which(!duplicated(prop.df$region))],2,code=3,length=0)

mtext(loc.thisrun,side=4,line=6.5)
dev.off()

} # close i

} # close plot_freq_pc


plot_freq_environ<-function(loci,genotype_data,site_data,out.dir,vector_to_plot,plot_by,plot_order){

# loci = a character vector of loci to plot
# gentoype_data = genetic data from which to calculate genotype frequencies
# site_data = site data containing the following variables: c("site_code","native","country","region","latitude")
# out.dir = directory to put PDF plots
# vector_to_plot = a numeric vector to select which loci to plot
# plot_by = what environmental variable would you like to order sites by? Options are: "latitude", "Altitude","mt","st","ap","sp","mm","sm"
# plot_order = descending or ascending?

for(i in vector_to_plot){

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
plr2<-site_data[,c("site_code","native","country","region","latitude","Altitude","mt","st","ap","sp","mm","sm")]
if(length(which(duplicated(plr2$site_code)))>0) plr2<-plr2[-which(duplicated(plr2$site_code)),]
prop.df<-merge(prop.df,plr2,by.x="site",by.y="site_code",all.x=T,all.y=F)

order_col<-which(colnames(prop.df)==plot_by)

# Order by native / non_native and relevan:
if(plot_order=="descending") prop.df<-prop.df[order(prop.df[,which(colnames(prop.df)=="native")],prop.df[,which(colnames(prop.df)=="region")],-prop.df[,order_col]),] else prop.df<-prop.df[order(prop.df[,which(colnames(prop.df)=="native")],prop.df[,which(colnames(prop.df)=="region")],prop.df[,order_col]),]
prop.df<-tidy.df(prop.df)
head(prop.df)

# quartz("",8,2,dpi=120,pointsize=10,file=paste(out.dir,paste(plot_by,"_",plot_order,"_",loc.thisrun,".pdf",sep=""),sep="/"),type="pdf")
# quartz("",8,2,dpi=120,pointsize=10)
par(mfrow=c(1,1),mar=c(4,5,2,1.5),oma=c(0,0,0,0))

bk <- c(-100,seq(0,100,by=10))
mc<-c("grey50",colorRampPalette(colors = c("white","blue"))(length(bk)-2))
image.thisrun<-image.plot(as.matrix(prop.df[,which(colnames(prop.df) %in% c("0","1","2"))][,3:1]),xaxt="n",yaxt="n",col=mc)

seq.thisrun<-seq(par("usr")[1],par("usr")[2],length.out=nrow(prop.df)+1)
incr.thisrun<-(seq.thisrun[2]-seq.thisrun[1])/2

axis(side=1,at=seq.thisrun,labels=F,las=2,cex.axis=1)
axis(side=1,at=(seq.thisrun+incr.thisrun)[1:(length(seq.thisrun)-1)],labels=prop.df$site,las=2,cex.axis=0.9,tick=F)
axis(side=2,at=c(0,0.5,1),labels=c("het","hom_a","hom_b"),tick=F,las=2)
axis(side=2,at=c(-0.25,0.25,0.75,1.25),labels=F)
arrows(-1,c(0.25,0.75),2,c(0.25,0.75),code=3,length=0)
arrows(seq.thisrun,-1,seq.thisrun,2,col="grey40")

axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region))],labels=F)

axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region[-which(prop.df$region=="Africa")]))]+(incr.thisrun*3),labels=as.character(unique(prop.df$region))[-which(as.character(unique(prop.df$region))=="Africa")],line=-1,tick=F,hadj=0,cex.axis=1)

axis(side=3,at=seq.thisrun[which(prop.df$region=="Africa")]+incr.thisrun,labels=as.character(unique(prop.df$region))[which(as.character(unique(prop.df$region))=="Africa")],line=-0.7,tick=F,hadj=0,cex.axis=0.65,las=2)

arrows(seq.thisrun[which(!duplicated(prop.df$region))],-1,seq.thisrun[which(!duplicated(prop.df$region))],2,code=3,length=0)

mtext(loc.thisrun,side=4,line=6.5)
dev.off()

} # close i

} # close plot_freq_environ


plot_freq_location<-function(loci,genotype_data,site_data,out.dir,number_to_plot){

# loci = a character vector of loci to plot
# gentoype_data = genetic data from which to calculate genotype frequencies
# site_data = site data containing the following variables: c("site_code","native","country","region","latitude")
# out.dir = directory to put PDF plots
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
# plr2<-site_data[,c("site_code","native","country","region","latitude")]
plr2<-site_data[,c("site_code","burn_unburnt","lat","long")] 
if(length(which(duplicated(plr2$site_code)))>0) plr2<-plr2[-which(duplicated(plr2$site_code)),]
prop.df<-merge(prop.df,plr2,by.x="site",by.y="site_code",all.x=T,all.y=F)
# prop.df<-prop.df[order(prop.df$native,prop.df$region,prop.df$country),]

# Order by native / non_native and latitude:
# prop.df<-prop.df[order(prop.df$native,prop.df$region,-prop.df$latitude),]
prop.df<-prop.df[order(prop.df$burn_unburnt, prop.df$long),]# rewriten by BD
prop.df<-tidy.df(prop.df)
head(prop.df)

dev.new("",8,2,dpi=120,pointsize=10,file=paste(out.dir,paste(loc.thisrun,".png",sep=""),sep="/"),type="pdf")
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

# axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region))],labels=F)

# axis(side=3,at=seq.thisrun[which(!duplicated(prop.df$region[-which(prop.df$region=="Africa")]))]+(incr.thisrun*3),labels=as.character(unique(prop.df$region))[-which(as.character(unique(prop.df$region))=="Africa")],line=-1,tick=F,hadj=0,cex.axis=1)

# axis(side=3,at=seq.thisrun[which(prop.df$region=="Africa")]+incr.thisrun,labels=as.character(unique(prop.df$region))[which(as.character(unique(prop.df$region))=="Africa")],line=-0.7,tick=F,hadj=0,cex.axis=0.65,las=2)

# arrows(seq.thisrun[which(!duplicated(prop.df$region))],-1,seq.thisrun[which(!duplicated(prop.df$region))],2,code=3,length=0)

mtext(loc.thisrun,side=4,line=6.5)
dev.off()

} # close i

} # close plot_freq_location

