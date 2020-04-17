
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

#----------Script 15/04/2020---------#
#> ghead(filtered_data); dim(filtered_data)
#site     ind L6 L7 L8 L9 L12 L26 L37 L47
#1 X01b X01b_01  2  1  2  0   0   0   2   1
#2 X01b X01b_02  2  1  2  0   0   0   2   1
#3 X01b X01b_03  2  1  2  0   0   0   2   1
#4 X01b X01b_04  2  1  2  0   0   0   2   1
#5 X01u X01u_01  1  2  0  0   0   0   0   1
#6 X01u X01u_02  2  1  2  0   0   0   2   1
#[1]    93 19120

#> ghead(dat_test); dim(dat_test)
#1 2 3 4 5 6 7 8 9 10
#L6  2 2 2 2 1 2 2 2 2  2
#L7  1 1 1 1 2 1 1 1 1  1
#L8  2 2 2 2 0 2 2 2 2  2
#L9  0 0 0 0 0 0 0 0 0  0
#L12 0 0 0 0 0 0 0 0 0  0
#L26 0 0 0 0 0 0 0 0 0  0
#[1] 19118    93

#> head(df_test); dim(df_test)
#loc1 loc2           r2
#1   L6   L7 0.0001216915
#2   L6   L8 0.4972949160
#3   L6   L9 0.2656598485
#4   L6  L12 0.6284829721
#5   L6  L26 0.6263439360
#6   L6  L37 0.6318159116
#[1] 182739403         3

#18/04/2020
#ld_loc<-read.table(paste(LD_dir, "LD_r75_over5pop_LOCI_FOR_REMOVAL.txt",sep="/"),header=T)
#head(ld_loc)

