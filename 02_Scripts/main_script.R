
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


#20/04/2020
#> ghead(filtered_data); dim(filtered_data)
#site     ind L6 L7 L8 L9 L12 L26 L37 L47
#1 X01b X01b_01  2  1  2  0   0   0   2   1
#2 X01b X01b_02  2  1  2  0   0   0   2   1
#3 X01b X01b_03  2  1  2  0   0   0   2   1
#4 X01b X01b_04  2  1  2  0   0   0   2   1
#5 X01u X01u_01  1  2  0  0   0   0   0   1
#6 X01u X01u_02  2  1  2  0   0   0   2   1
#[1]    93 19120

#"no loci before ld filt = 19120"
#"no loci after ld filt = 11332"
#> head(hwe.res); dim(hwe.res)
#locus            p        p.adj
#1    L6 1.000000e+00 1.000000e+00
#2    L7 2.314964e-01 1.000000e+00
#3    L8 6.324177e-02 1.000000e+00
#4    L9 1.423518e-23 1.546795e-19
#5   L12 6.206069e-02 1.000000e+00
#6   L26 4.986915e-02 1.000000e+00
#[1] 11330     3

#"no loci before ld filt = 11332"
#"no loci after ld filt = 9371"
#> ghead(filtered_data); dim(filtered_data)
#site     ind L9 L37 L47 L54 L59 L164 L169 L176
#1 X01b X01b_01  0   2   1   2   0    2    1    1
#2 X01b X01b_02  0   2   1   2   0    2    1    1
#3 X01b X01b_03  0   2   1   2   0    2    1    1
#4 X01b X01b_04  0   2   1   2   0    2    1    1
#5 X01u X01u_01  0   0   1   0   1    0    0   NA
#6 X01u X01u_02  0   2   1   0   0    2    1    1
#[1]   93 9371

#Not Sure what's wrong with 
#> filtered_data<-get(data_name)
#Error in get(data_name) : object 'snp_onerow' not found
