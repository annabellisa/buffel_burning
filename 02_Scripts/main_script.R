
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
#cr<-0.5 @#42
#ra<-0.98 @#50
#malim<-0.01 @#65

#End @#54
#site     ind L6 L7 L8 L11 L12 L26 L41 L44
#1 X01b X01b_01  2  1  2   0   0   0   0   1
#2 X01b X01b_02  2  1  2   0   0   0  NA   1
#3 X01b X01b_03  2  1  2   0   0   0   0   1
#4 X01b X01b_04  2  1  2   0   0   0   0   1
#5 X01u X01u_01  1  2  0   0   0   0   1   2
#6 X01u X01u_02  2  1  2   0   0   0  NA   1
#[1]    93 30119

#End @#67
#site     ind L6 L7 L8 L11 L12 L26 L41 L44
#1 X01b X01b_01  2  1  2   0   0   0   0   1
#2 X01b X01b_02  2  1  2   0   0   0  NA   1
#3 X01b X01b_03  2  1  2   0   0   0   0   1
#4 X01b X01b_04  2  1  2   0   0   0   0   1
#5 X01u X01u_01  1  2  0   0   0   0   1   2
#6 X01u X01u_02  2  1  2   0   0   0  NA   1
#[1]    93 29009


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


