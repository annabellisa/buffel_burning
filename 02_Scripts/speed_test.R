# Conclusions of speed test:

# 1. No difference between OneDrive and Desktop
# 2. No difference between R Studio and R app
# 3. FRESH sessions (i.e. cleared of data stored in current workspace) are always faster
# 4. Big mac's speed is coming from the faster processor, not the multiple cores (16 vs 4 cores). 
# 5. Speed difference between laptop and big mac are not huge
# Only way to speed this up is to use parallel functions

# Results (mins):

# Big mac, OneDrive, R Studio (fresh)
5000 = 0.116
15000 = 1.85
35000 = 14.116

# Big mac, OneDrive, R Studio (not fresh)
5000 = NA
15000 = 2.75
35000 = 16

# Big mac, OneDrive, R app (fresh)
5000 = 0.116
15000 = 1.86

# Big mac, Desktop, R app (fresh)
5000 = 0.783
15000 = 1.83

# Laptop, Desktop, R app (fresh)
5000 = 0.15
15000 = 2.73
35000 = 16.96


library(parallel)
numCores<-detectCores()

# the mac is slow, find out why:

test_number<-35000
test_names<-sample(letters[],size=test_number,replace=T)
start_time <- gsub("[: -]", "" , Sys.time(), perl=TRUE)
test_combn<-combn(test_names,2)
end_time <- gsub("[: -]", "" , Sys.time(), perl=TRUE)

print(paste("names test START: ", start_time, sep=""))
print(paste("There were this many names: ", test_number, sep=""))
print(paste("There were this many pair wise combinations: ", length(test_combn)))
print(paste("names test END: ", end_time, sep=""))
paste("it took this many minutes: ",(as.numeric(end_time)-as.numeric(start_time))/60,sep="")

