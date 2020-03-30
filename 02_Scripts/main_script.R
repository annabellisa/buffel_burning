
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


# in terminal

cd /Users/annabelsmith/OneDrive\ -\ The\ University\ of\ Queensland/TEACHING/UQ_Masters_and_special_topic/Binyin_Winter/Analysis/Binyin_analysis

git init

git commit -m "First commit"




#Hello World
#Welcome
#Changed by BD


