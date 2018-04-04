library(cgdsr)
library(reshape)
library(RTCGA)


'''
Definitions: Simply put, a z-score is the number 
of standard deviations from the mean a data point is. But more technically its
a measure of how many standard deviations below or above the population mean 
a raw score is. A z-score is also known as a standard score and it can be placed 
on a normal distribution curve.


'''


# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")

# Get list of cancer studies at server
cancerstudies<-getCancerStudies(mycgds)
write.csv(cancerstudies,"cancerstudies.csv")

#Get Complete Gene List
prad_study <- getCancerStudies(mycgds)[117,1]
pca_case_list <- getCaseLists(mycgds,prad_study)[1,1]
Full.Gene.List <- read.csv("GeneNameList.csv")
#rows.gene.list <-Full.Gene.List$V1
length(Full.Gene.List)
i <- 1
while(i <= length(Full.Gene.List)){
  if(i == 1){
    z_pca <- matrix()
  }
  k <- i+499
  temp <- getProfileData(mycgds, Full.Gene.List[i:k,1], "prad_tcga_rna_seq_v2_mrna_median_Zscores",pca_case_list)
  message(paste("done for genes -->",i," to ",k))
  if(dim(temp)[1] != 0){
    z_pca <- cbind(z_pca,  temp)
  }
  i <- k+1
}

z_pca <- t(z_pca)
rem_pca <- as.numeric(which(apply(z_pca,1,function(x) sumis.na(x))) == ncol(z_pca))
z_pca <- z_pca[-rem_pca,]
dim(z_pca)

colnames(z_pca) <- substr(colnames(z_pca),1,12)

#GBM
GBM.study = getCancerStudies(mycgds)[51,1]
find.GBM.caselists = getCaseLists(mycgds,GBM.study)
View(findRNAexp)


get.GBM.caselist = getCaseLists(mycgds,GBM.study)
allprofiles <- getGeneticProfiles(mycgds,GBM.study)
#Get mRNA expression Z-scores (see Z-score def)

#Get Mutation Data
#getMutationData(x, caseList, geneticProfile, genes, ...)

getMutationData(mycgds,GBM.study,get.GBM.caselist[7,1])


#View(allprofiles)
GBM.mRNA.exp.profile = getGeneticProfiles(mycgds,GBM.study)[1,1]
GBM.mutation.profile = getGeneticProfiles(mycgds,GBM.study)[7,1]
GBM.genes.rnaseq <- getProfileData(mycgds,Full.Gene.List,GBM.mRNA.exp.profile,GBM.caselist)
scaled.GBM.genes <- scale(GBM.genes.rnaseq)
df.GBM<-data.frame(scaled.GBM.genes)
#df.GBM<-data.frame(GBM.genes.rnaseq)

df.GBM["CancerType"] <- "GBM" 

cancerdata<-rbind(df.BRCA, df.GBM, df.LAML, df.LGG, df.PAAD, df.PRAD)

write.csv(cancerdata, file="/Volumes/Macintosh HD/Users/louisecabansay/Dropbox (Personal)/CancerHackers/CancerGeneData.csv")