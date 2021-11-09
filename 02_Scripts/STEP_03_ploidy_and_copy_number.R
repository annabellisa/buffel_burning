
# ------------------------------------ #
# ------------- STEP 03  ------------- #
# ------------------------------------ #

### Estimate ploidy of individuals (Gompert and Mock 2017) and copy number of loci (McKinney et al. 2017)

# Both of these methods require data on read counts for each individual locus across all individuals. We only have the average counts; have emailed Andrzej for more info. Will decide how to proceed after hearing from him. 

### Updated following peer-review

### Code author: Annabel Smith (except where indicated)

# Load workspace including data and functions:
load("03_Workspaces/STEP03_ploidy_and_CN.RData")

# The data (filtered_data) include all 40,711 loci after only the most basic fileters were applied:
# Duplicate sequences removed
# Monomorphic loci removed
# Call rate > 50% 

ghead(filtered_data); dim(filtered_data)

library(gbs2ploidy)
?gbs2ploidy

# The simulated data is for 200 individuals at 10,000 loci
# The first two components are N (number of individuals) by P (number of SNPs) matrixes with allele counts for the first and second allele at each locus, respectively. The third component is a numeric vector that gives the true ploidy for each individual (2 = diploid, 4 = tetraploid).

# From the Gompert & Mock paper: "Let yij denote the number of sequence reads containing the nonreference allele for individual j at heterozygous SNP i (i.e. we only consider SNPs where individual j is heterozygous), and nij denote the total number of reads covering that SNP for individual j"

# So these are the allele counts, which we don't have from the dartseq data. We only have the average allele counts for each locus. 

head(linf,2)
ghead(filtered_data)
data(dat)
str(dat)
head(dat[[1]][,1:10],3); dim(dat[[1]])
head(dat[[2]][,1:10],3); dim(dat[[2]])
head(dat[[3]][,1:10],3)



# HDplot to identify paralogs: 

# This function, written by McKinney et al. (2017), was downloaded on 4/11/21 from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cm08m

# McKinney, G.J., Waples, R.K., Seeb, L.W. & Seeb, J.E. 2017. Paralogs are revealed by proportion of heterozygotes and deviations in read ratios in genotyping-by-sequencing data from natural populations. Mol Ecol Resour 17: 656-669.

#HDplot
HDplotData<-read.delim("HDplot_genericInput.txt",header=TRUE)
names(HDplotData)
#calculate percentage of heterozygotes (H)
HDplotData$hetPerc<-HDplotData$num_hets/HDplotData$num_samples
#calculate allele ratio
HDplotData$totalCounts<-HDplotData$depth_a+HDplotData$depth_b
HDplotData$ratio<-HDplotData$depth_a/HDplotData$totalCounts
#calculate read-ratio deviation (D)
HDplotData$std<-sqrt(HDplotData$totalCounts*0.5*0.5)
#calculate z-score for each locus
HDplotData$z<- -(HDplotData$totalCounts/2-HDplotData$depth_a)/HDplotData$std

#plot HD plot results
plot(z~hetPerc,data=HDplotData,xlab='H',ylab='D')
#plot H vs allele ratios
plot(ratio~hetPerc,data=HDplotData,xlab='H',ylab='Allele Ratio')





