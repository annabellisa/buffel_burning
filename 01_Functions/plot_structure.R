
# Author: Annabel Smith

# V10 was updated from the PLANTPOPNET V8 

str_plot_V10<-function(K,cluster.data,site.data,las.opt,yaxs.loc,col.pal,site.lab,...){

head(cluster.data,3)

# Organise columns:
cluster.data$max_p<-apply(cluster.data[,3:(2+K)],1,function(x) which(x==max(x)))

cluster.data<-cbind(cluster.data[,1:2],cluster.data[,3:(K+2)][,unique(cluster.data$max_p)],cluster.data[,(K+3):length(cluster.data)])

# make assignment data matrix:
dat.thisrun<-apply(cluster.data[grep("assig",colnames(cluster.data))],1,rbind)
head(dat.thisrun)

# make site data:
sdat.thisrun<-cluster.data[,c("site",colnames(cluster.data)[which(colnames(cluster.data)=="block"):length(cluster.data)])]

# Set colour scheme:
if(col.pal %in% rownames(brewer.pal.info)) cols.thisrun<-brewer.pal(K,col.pal) else cols.thisrun<-get(col.pal)[1:K]

# Plot:
par(mar=c(5,5,1,0),mgp=c(1,yaxs.loc+1,yaxs.loc),xpd=F)
p1<-barplot(dat.thisrun,cex.axis=1.5,col=cols.thisrun,space=0,border=cols.thisrun,xaxt="n",las=las.opt,ylab="Probability \nof assignment",cex.lab=1.5)

# Add lines for individuals and sites:
arrows(which(!duplicated(sdat.thisrun$site))-1,-0.08,which(!duplicated(sdat.thisrun$site))-1,0.995,code=0,lwd=0.8)
arrows(nrow(sdat.thisrun),-0.08,nrow(sdat.thisrun),0.995,code=0,lwd=0.8)
arrows(1:nrow(sdat.thisrun),0,1:nrow(sdat.thisrun),0.995,code=0,lwd=0.05)

# Draw a box:
arrows(0,1,nrow(sdat.thisrun),1,code=0,lwd=0.8)
arrows(0,0,nrow(sdat.thisrun),0,code=0,lwd=0.8)

head(sdat.thisrun,3)
head(cluster.data,3)
dim(dat.thisrun)

# Add labels for sites
axis(1,which(!duplicated(sdat.thisrun$site))+1,labels=unique(cluster.data$site),tick=F,line=2.5,las=1,cex.axis=1.5)

# Add longer tick marks between sites:
axis(side=1, at=which(!duplicated(sdat.thisrun$block))-1, labels=F, tick = T, line=0, tck=-0.1, lwd=0, lwd.ticks = 1)

} # close structure plot function V10






