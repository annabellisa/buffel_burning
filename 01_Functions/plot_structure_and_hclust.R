
# Author: Annabel Smith

# STRUCTURE PLOT: V10 was updated from the PLANTPOPNET V8 

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

# Name order for hclust output:

hclust_name_order<-function(hclust_output){
  
  # Author = Annabel Smith
  
  # Arguments:
  # hclust_output = the named output from an hclust fit
  
  # Details:
  
  # This function extracts the sample names in the order they appear on the plot of the hclust output, from left to right. It was created to access the index for ordering a structure plot, in the same order as the hclust output. 
  
  # There are three options with merge:
  
  # 1. Double individual (both negative)
  # 2. Single individual + 1 group (one negative, one positive)
  # 3. Double group
  
  # For option # 1, where both clusters are single individuals, they are ordered by their position in the original data (from hclust? 'merges involving two observations place them in order by their observation sequence number')
  
  # For option # 2, the negative (single individual) is ALWAYS on the left, i.e. the first column in row i of merge. So when a single observation merges with an already merged group, the single is ALWAYS on the left.
  
  # For option # 3, they are ordered in terms of the order they were merged. Check the data in merge and you see that the lower number is always on the left (i.e. the first column). This means they can be ordered simply by their group number. 
  
  # These options dictate three different ways of extracting the order data.
  
  # The process proceeds iteratively, over each step in merge (i.e. n-1, where n is the number of samples, see hclust help file)
  
  ## testing space:
  ## hclust_ouput<-euc_clust
  ## end test space
  
  # create data frame with the name labels and their order (original and permuted):
  hc_order<-data.frame(sample=hclust_ouput$labels, orig_pos=1:length(hclust_ouput$order), order=hclust_ouput$order)
  
  # create data frame with the merge data:
  hc_merge<-data.frame(hclust_ouput$merge)
  hc_merge$height<-hclust_ouput$height
  
  # create a list to store the results of each step:
  res_store<-list()
  
  for (i in 1:nrow(hc_merge)){
    
    # determine the type of merge, from the three options above:
    step_thisrun<-hc_merge[i,]
    clust1<-step_thisrun[,1]
    clust2<-step_thisrun[,2]
    
    # OPTION # 1:
    if (clust1<0 & clust2<0){
      # if double individual, assign both clusters to vector, in their named order:
      res_store[[i]]<-c(as.character(hc_order[hc_order$orig_pos==as.numeric(abs(clust1)),]$sample),as.character(hc_order[hc_order$orig_pos==as.numeric(abs(clust2)),]$sample))
    } # close if option # 1
    
    # OPTION # 2:
    if (clust1<0 & clust2>0){
      # if single individual + group, assign individual to vector and pull in group from previous step:
      res_store[[i]]<-c(as.character(hc_order[hc_order$orig_pos==as.numeric(abs(clust1)),]$sample))
      
      # the group is always in the second column of the current row:
      group.thisrun<-abs(clust2)
      
      # "Single observations are the tightest clusters possible"
      # A single observation merging with an already merged group ALWAYS places the single observation on the left, so we can put the single observation automatically on the left:
      res_store[[i]]<-c(res_store[[i]], res_store[[group.thisrun]])
      
    } # close if option # 2
    
    # OPTION # 3:
    if (clust1>0 & clust2>0){
      
      # if double group, pull in both groups, and order them by their merge order:
      group1.thisrun<-abs(clust1)
      group2.thisrun<-abs(clust2)
      res_store[[i]]<-c(res_store[[group1.thisrun]], res_store[[group2.thisrun]])
    } # close if option # 3
    
    if (i == nrow(hc_merge)) return(res_store[[nrow(hc_merge)]])
    
  } # close for i row of merge
  
} # close hclust_name_order function





