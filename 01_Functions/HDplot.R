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

