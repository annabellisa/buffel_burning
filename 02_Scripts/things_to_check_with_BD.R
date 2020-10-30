# Things to check with Binyin
# 30th Oct 2020

# Hi Binyin, a couple of things to follow up on. 

# PCAdapt questions:
# we have identified 7040 outliers from PCAdapt, based on q values and a threshold of 0.05
# the values in the outliers file for PCAdapt (i.e. before it was combined with the other methods) range from 0:3. What are the zeros? These sum to 7040, so basically all of these were identified as outliers, but I assumed this was the PC on which the locus was associated. We set K to 3, so how can there be four here?

# LFMM questions:
# Why did we choose K=3 for PCAdapt but K=5 for LFMM?

dir("00_Data/")
xx<-read.table("00_Data/outliers_all.txt", header=T)
head(xx); dim(xx)
table(xx$pca_outl)

dir("RESULTS/PCAdapt/PCAdapt_results/")
xx<-read.table("RESULTS/PCAdapt/PCAdapt_results/outliers_PCAdapt_from_snp_pc.txt", header=T)
head(xx); dim(xx)
table(xx$PC)



