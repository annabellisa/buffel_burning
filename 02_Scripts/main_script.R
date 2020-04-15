
# Binyin winter project

# by Di Binyin & Annabel Smith
dir()
dir("01_Functions")

# load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))


# load workspace with relevant data:
load("binyin_winter.RData")

# see this page for how to change the origin:
# https://stackoverflow.com/questions/39435240/rstudio-changing-origin-for-git-version-control-of-project
# in terminal:
# > git remote -v
# if an origin is already specified:
# > git remote rm origin
# > git remote addorigin https://github.com/annabellisa/Binyin_Winter.git


#git problem 14 April
# git remote -v is showing the new origin (tamielleB), while Tools > Version Control > Project settings is showing the original Binyin origin - so far so good
# updating after reset

# reset git, added EVERYTHING to gitignore in root dir
# final check after fixing tam's repo


# in terminal

cd /Users/annabelsmith/OneDrive\ -\ The\ University\ of\ Queensland/TEACHING/UQ_Masters_and_special_topic/Binyin_Winter/Analysis/Binyin_analysis

git init

git commit -m "First commit"



#----------Script 04/04/2020---------#
library(vcfR)
install.packages("vcfR")
library(adegenet)
install.packages("adegenet")
library(adegraphics)
install.packages("adegraphics")
library(pegas)
install.packages("pegas")
library(StAMPP)
install.packages("STAMPP")
library(lattice)
library(gplots)
install.packages("gplots")
library(ape)
library(ggmap) 
install.packages("ggmap")


options(max.print = 100000)
ghead(filtered_data); dim(filtered_data)
ls(); View(filtered_data)
str(filtered_data)

#----------Script 14/04/2020---------#
#cr<-0.5 @#42

#quality rate: 0.05 from 
#"C:/Users/s4467005/OneDrive - The University of Queensland/Smith Lab/Schilling_etal_2014_PlosOne_HWE.pdf"
#Call rate 0.8-0.9:0.8/0.85/0.9
#---0.9-0.8:
#"There were no loci with more than 90 % missing data; none removed"
#0.5:"6164 loci with more than 50 % missing data removed"


# Remarks derived from filter_SNPs:
# Cannot use DartSeq CallRate filter because it was calculated on the whole dataset. Individuals from the RCH and FS populations will cause higher proportion of missing data in the rest of the data:
# This metric is highly correlated with the DartSeq metric (99%), but it filters loci biased by the presence of other populations, e.g. removing the firescape samples results in approx. 500 extra loci filtered.
# cr<-linf[,c("locus","CallRate")]
# miss<-merge(miss,cr,by="locus",all.x=T,all.y=F)
# plot(miss$missing_data,miss$CallRate)
# cor.test(miss$missing_data,miss$CallRate)
# callrate50<-miss$locus[which(miss$missing_data>0.5)]
# callrate50_dart<-miss$locus[which(linf$CallRate<0.5)]


#ra<-0.98 @#50
#0.95 instead? 


#malim<-0.01 @#65
#0.01*1 v 0.05*2--> 0.05 instead?


#> LD_dir<-"../../ANALYSIS_RESULTS/LINKAGE_DISEQUILIBRIUM/LD_parameters"
#> dir(LD_dir)
#character(0)????
#Need a file: LD_r75_over5pop_LOCI_FOR_REMOVAL.txt



#Not Sure:
# Filter loci with were in LD in > 5 populations with a correlation of 0.75 (see Supplement_03_LD_tests.R for details):
#@73 Error
#> LD_dir<-"../../ANALYSIS_RESULTS/LINKAGE_DISEQUILIBRIUM/LD_parameters"
#> dir(LD_dir)
#character(0)
#> ld_loc<-read.table(paste(LD_dir, "LD_r75_over5pop_LOCI_FOR_REMOVAL.txt",sep="/"),header=T)
#Error in file(file, "rt") : cannot open the connection
#In addition: Warning message:
  #In file(file, "rt") :
  #cannot open file '../../ANALYSIS_RESULTS/LINKAGE_DISEQUILIBRIUM/LD_parameters/LD_r75_over5pop_LOCI_FOR_REMOVAL.txt': No such file or directory


