# Things to check with Binyin
# 30th Oct 2020

# Hi Binyin, a couple of things to follow up on. 

# General questions

# "RESULTS/NF_Format" Neutral format == NF? 

# NF = neutral format
# N-NF = non-neutral format

# is this ^ correct?

# we need to remove ALL workspaces in the 03_Workspaces folder that are not being used
# we need a shareable key that describes the filtering that was used for each workspace and what it is used for

# We should no longer be using the files in Offline_Results/Genepop_Files (e.g. Genpop_Diversity_Original.gen) because they contain an old subset of data (24051 loci). I've moved those to Offline_Results/01_old_results/Genepop_Files

# I can see in Binyin_Winter/RESULTS there are two folders:
# NF_Format
# and
# N_NF_Format

# both of these contain genepop files called "Genpop_Diversity_Original.gen", but which were created using different filters. Here lies the necessity of creating the param files when doing filtering. 

# first of all, it's good practice to name files that contain different data sets with different file names

# HOWEVER, you've put two param files in each of these folders - each param file describes the same number of loci, but with DIFFERENT filters

# it's very important that we get this right - we need to know exactly what filtering was done to create each data set. With different multiple files, for the same data set, but with different filters - how can we know which filters were applied?




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



