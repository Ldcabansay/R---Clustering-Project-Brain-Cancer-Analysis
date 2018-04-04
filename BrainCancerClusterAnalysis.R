#load cluster libraries
library(cluster)
library(NbClust)
library(flexclust)
library(rattle)
library(ggplot)


library(readxl)
#Loading data files from excel and subsetting glioblastoma and glioma cancer gene sets
precogdata <- read_excel("recog_MetaZ.xlsx")
gbm.oncolnc <- read_excel("GBM_mrna.xlsx")
gli.oncolnc <- read_excel("LowerGradeGlioma.xlsx")
gbm.precog <-  data.frame(precogdata$Gene, precogdata$Name, precogdata$`Unweighted_meta-Z_of_all_cancers`, precogdata$Brain_cancer_Glioblastoma)
gli.precog <-  data.frame(precogdata$Gene, precogdata$Name, precogdata$`Unweighted_meta-Z_of_all_cancers`, precogdata$Brain_cancer_Glioma)

#examine data
names(gbm.oncolnc)[1] <- "GeneName"
names(gli.oncolnc)[1] <- "GeneName"
names(gbm.precog)[1] <-"GeneName"
names(gli.precog)[1] <-"GeneName"
str(gbm.precog)
str(gbm.oncolnc)
str(gli.oncolnc)
describe(gbm.precog)
describe(gbm.oncolnc)
describe(gli.oncolnc)

#combine dataframes to create one Glioblastoma master dataframe
gbm.df <- merge(gbm.precog, gbm.oncolnc, by = "GeneName")
gli.df <- merge(gli.precog, gli.oncolnc, by = "GeneName")
head(gbm.df)
head(gli.df)
names(gbm.df)[2:4] <- c("GeneFunction", "Unweighted_meta_Z_allcancers", "precog_Z")
names(gli.df)[2:4] <- c("GeneFunction", "Unweighted_meta_Z_allcancers", "precog_Z")
sig.gbm<-gbm.df[gbm.df$`Raw p-value` <0.05,]
sig.gli<-gli.df[gli.df$`Raw_p-value`<0.05,]

gli.genes<-sig.gli[sig.gli$GeneName %in% sig.gbm$GeneName,]
gbm.genes<-sig.gbm[sig.gbm$GeneName %in% gli.genes$GeneName,]

#standardize num variables to a mean of 50 with a standard deviation of 10 for ease of comparison
gbm.scaled <- scale(gbm.genes[3:9])
gli.scaled <- scale(gli.genes[3:9])
glioblastoma<-cbind(gbm.genes[1:2], gbm.scaled)
glioma<-cbind(gli.genes[1:2], gli.scaled)
names(glioblastoma)[6:7] <-c("pvalue", "bh_pvalue")
names(glioma)[6:7] <-c("pvalue", "bh_pvalue")
describe(glioblastoma)
describe(glioma)

#Average-Linkage Clustering
#best w/small datasets (less than 150) so will select random groups

#w/50 randomly sampled genes
r50.glioblastoma <- glioblastoma[sample(1:nrow(glioblastoma), 50, replace=FALSE),]
r50.glioma <- glioma[glioma$GeneName %in% r50.glioblastoma$GeneName,]
row.names(r50.glioblastoma) <- r50.glioblastoma$GeneName
row.names(r50.glioma) <- r50.glioma$GeneName   

#save set to CSV to replicate later
write.csv(r50.glioblastoma, "r50glioblastoma.csv")
write.csv(r50.glioma, "r50glioma.csv")



opar <- par(no.readonly=TRUE) #save current state

par(mfrow=c(2,1)) #layout dendogram hierarchy plots 2 rows, 1 column

#make glioblastoma cluster dendogram
nc.r50.glioblastoma <- r50.glioblastoma[3:9]                       
d.r50.glioblastoma <- dist(r50.glioblastoma)                                          
fit.gbm.average50 <- hclust(d.r50.glioblastoma, method="average")                          
plot(fit.gbm.average50, hang=-1, cex=.8, main="50 gene Average Linkage Clustering for Glioblastoma")

#make glioma cluster dendogram
nc.r50.glioma <- r50.glioma[3:9]                       
d.r50.glioma <- dist(r50.glioma)                                          
fit.gli.average50 <- hclust(d.r50.glioma, method="average")                          
plot(fit.gli.average50, hang=-1, cex=.8, main="50 gene Average Linkage Clustering for Glioma")

par(opar)

#Selecting best number of clusters for both glioma50 and glioblastoma50

nc.gli50 <- NbClust(nc.r50.glioma, distance="euclidean", 
                  min.nc=2, max.nc=10, method="average")
table(nc.gli50$Best.n[1,])

nc.gbm50 <- NbClust(nc.r50.glioblastoma, distance="euclidean", 
                    min.nc=2, max.nc=10, method="average")
table(nc.gbm50$Best.n[1,])

par(opar)
par(mfrow=c(1,2)) #layout cluster bar plots next to each other plots 1rows, 2column
barplot(table(nc.gli50$Best.n[1,]), 
        xlab="Numer of Clusters", ylab="Number of Criteria",
        main="Number of GLI Clusters Chosen by 7 Criteria")

barplot(table(nc.gbm50$Best.n[1,]), 
        xlab="Numer of Clusters", ylab="Number of Criteria",
        main="Number of GBM Clusters Chosen by 7 Criteria") 

par(opar)
#Create 5 Cluster dendogram for 50 glioma genes
gli.50.clusters <- cutree(fit.gli.average50, k=5) 
table(gli.50.clusters) #makes table showing how many genes in each cluster
aggregate(nc.r50.glioma, by=list(cluster=gli.50.clusters), median) 
aggregate(as.data.frame(nc.r50.glioma), by=list(cluster=gli.50.clusters),
          median)

#plot clusters w/ red boxes indicating clusters
plot(fit.gli.average50, hang=-1, cex=.8,  
     main="GLI Average Linkage Clustering\n5 Cluster Solution") 
#overlay red rectangle on cluster
rect.hclust(fit.gli.average50, k=5)

par(opar)
#Create 5 Cluster dendogram for 50 glioblastoma genes
gbm.50.clusters <- cutree(fit.gbm.average50, k=5) 
table(gbm.50.clusters)
aggregate(nc.r50.glioblastoma, by=list(cluster=gbm.50.clusters), median) 
aggregate(as.data.frame(nc.r50.glioblastoma), by=list(cluster=gbm.50.clusters),
          median)
plot(fit.gbm.average50, hang=-1, cex=.8,  
     main="GBM Average Linkage Clustering\n5 Cluster Solution")
rect.hclust(fit.gbm.average50, k=5)


par(opar)

gli.50.clist <- lapply(sort(unique(gli.50.clusters)), function(x) r50.glioma[which(gli.50.clusters==x),])

gli.cluster.50.1 <- cbind(gli.50.clist[[1]], Cluster = "Gli.Cluster1")
gli.cluster.50.2 <- cbind(gli.50.clist[[2]], Cluster = "Gli.Cluster2")
gli.cluster.50.3 <- cbind(gli.50.clist[[3]], Cluster = "Gli.Cluster3")
gli.cluster.50.4 <- cbind(gli.50.clist[[4]], Cluster = "Gli.Cluster4")
gli.cluster.50.5 <- cbind(gli.50.clist[[5]], Cluster = "Gli.Cluster5")
#gli.cluster.50.6 <- cbind(gli.50.clist[[6]], Cluster = "Gli.Cluster6")
#gli.cluster.50.7 <- cbind(gli.50.clist[[7]], Cluster = "Gli.Cluster7")
#gli.cluster.50.8 <- cbind(gli.50.clist[[8]], Cluster = "Gli.Cluster8")
Glioma_50_Genes_C <- rbind(gli.cluster.50.1, gli.cluster.50.2, gli.cluster.50.3, gli.cluster.50.4,
                                 gli.cluster.50.5)# gli.cluster.50.6,gli.cluster.50.7, gli.cluster.50.8)

Glioma_50_Genes_C[,-c(3:9)]


#Create table for clustering results
gbm.50.clist <- lapply(sort(unique(gbm.50.clusters)), function(x) r50.glioblastoma[which(gbm.50.clusters==x),])

gbm.cluster.50.1 <- cbind(gbm.50.clist[[1]], Cluster = "Gbm.Cluster1")
gbm.cluster.50.2 <- cbind(gbm.50.clist[[2]], Cluster = "Gbm.Cluster2")
gbm.cluster.50.3 <- cbind(gbm.50.clist[[3]], Cluster = "Gbm.Cluster3")
gbm.cluster.50.4 <- cbind(gbm.50.clist[[4]], Cluster = "Gbm.Cluster4")
gbm.cluster.50.5 <- cbind(gbm.50.clist[[5]], Cluster = "Gbm.Cluster5")
#gbm.cluster.50.6 <- cbind(gbm.50.clist[[6]], Cluster = "Gbm.Cluster6")
#gbm.cluster.50.7 <- cbind(gbm.50.clist[[7]], Cluster = "Gbm.Cluster7")
#gbm.cluster.50.8 <- cbind(gbm.50.clist[[8]], Cluster = "Gbm.Cluster8")
glioblastoma_50_Genes_C <- rbind(gbm.cluster.50.1, gbm.cluster.50.2, gbm.cluster.50.3, gbm.cluster.50.4,
                               gbm.cluster.50.5) #gbm.cluster.50.6,gbm.cluster.50.7, gbm.cluster.50.8)

#save cluster tables results into csvs(optional)
write.csv(glioblastoma_50_Genes_C, "Table_Glioblastoma_50genes_5Clust.csv")



#Repeat Ward Clustering w/ 100 genes
#generate 100 random genes that are significant in both gli and gbm
#grabs from glioblastoma & glioma dfs (contains same genes)

#r100.glioblastoma <- glioblastoma[sample(1:nrow(glioblastoma), 100, replace=FALSE),]
#r100.glioma <- glioma[glioma$GeneName %in% r100.glioblastoma$GeneName,]

#save randomly generated gene df into csv for reproducibility 

#note: for markdown file, loaded previous randomly generated df for report analysis
gbm.oncolnc <- read.csv("GBM_mrna.xlsx")
gli.oncolnc <- read.csv("LowerGradeGlioma.xlsx")

write.csv(r100.glioblastoma, "r100glioblastoma.csv")
write.csv(r100.glioma, "r100glioma.csv")

row.names(r100.glioblastoma) <- r100.glioblastoma$GeneName
row.names(r100.glioma) <- r100.glioma$GeneName


par(opar)
nc.r100.glioblastoma <- r100.glioblastoma[3:9]
d.r100.glioblastoma <- dist(nc.r100.glioblastoma)                                          
fit.gbm.ward100 <- hclust(d.r100.glioblastoma, method="ward")                          
plot(fit.gbm.ward100, hang=-1, cex=.8, main="100 Gene Ward Method Clustering: Glioblastoma")

nc.r100.glioma <- r100.glioma[3:9]
d.r100.glioma <- dist(nc.r100.glioma)                                          
fit.gli.ward100 <- hclust(d.r100.glioma, method="ward")                          
plot(fit.gli.ward100, hang=-1, cex=.8, main="100 Gene Ward Method Clustering: Glioma")
par(opar)


#find best number of clusters for 100 genes hierarchial clustering

#find for 100 gli genes
nc.gli.100 <- NbClust(nc.r100.glioma, distance="euclidean", 
              min.nc=2, max.nc=20, method="ward.D")
table(nc.gli.100$Best.n[1,])
barplot(table(nc.gli.100$Best.n[1,]), 
        xlab="Number of Clusters", ylab="Number of Criteria",
        main="Number of Ward Clusters Chosen by 7 Criteria: 100 GLI") 
par(opar)


#find for 100 gbm genes
nc.gbm.100 <- NbClust(nc.r100.glioblastoma, distance="euclidean", 
                  min.nc=2, max.nc=20, method="ward.D")
table(nc.gbm.100$Best.n[1,])
barplot(table(nc.gbm.100$Best.n[1,]), 
        xlab="Numer of Clusters", ylab="Number of Criteria",
        main="Number of Ward Clusters Chosen by 7 Criteria: 100 GBM") 
par(opar)


#plot 100 gene glioma dendogram w/red boxes around 10 clusters
gli.100.clusters <- cutree(fit.gli.ward100, k=10) 
table(gli.100.clusters)
aggregate(nc.r100.glioma, by=list(cluster=gli.100.clusters), median) 
aggregate(as.data.frame(nc.r100.glioma), by=list(cluster=gli.100.clusters),
          median)
plot(fit.gli.ward100, hang=-1, cex=.8,  
     main="Glioma 100 Gene Ward Method Clustering\n10 Cluster Solution")
#outline k clusters in red boxes on dendogram
rect.hclust(fit.gli.ward100, k=10)

par(opar)
#plot 100 gene glioblastoma dendogram w/redboxes around 10 clusters
gbm.100.clusters <- cutree(fit.gbm.ward100, k=10) 
table(gbm.100.clusters)
aggregate(nc.r100.glioblastoma, by=list(cluster=gbm.100.clusters), median) 
aggregate(as.data.frame(nc.r100.glioblastoma), by=list(cluster=gbm.100.clusters),
          median)
plot(fit.gbm.ward100, hang=-1, cex=.8,  
     main="Glioblastoma 100 Gene Ward Method Clustering\n10 Cluster Solution")
#outline k clusters in red boxes on dendogram
rect.hclust(fit.gbm.ward100, k=10)


#Generate tables for ward clustering results

#100 gene glioma 10 cluster table
gli.100.clist <- lapply(sort(unique(gli.100.clusters)), function(x) r100.glioma[which(gli.100.clusters==x),])
gli.cluster.100.1 <- cbind(gli.100.clist[[1]], Cluster = "gli.Cluster1")
gli.cluster.100.2 <- cbind(gli.100.clist[[2]], Cluster = "gli.Cluster2")
gli.cluster.100.3 <- cbind(gli.100.clist[[3]], Cluster = "gli.Cluster3")
gli.cluster.100.4 <- cbind(gli.100.clist[[4]], Cluster = "gli.Cluster4")
gli.cluster.100.5 <- cbind(gli.100.clist[[5]], Cluster = "gli.Cluster5")
gli.cluster.100.6 <- cbind(gli.100.clist[[6]], Cluster = "gli.Cluster6")
gli.cluster.100.7 <- cbind(gli.100.clist[[7]], Cluster = "gli.Cluster7")
gli.cluster.100.8 <- cbind(gli.100.clist[[8]], Cluster = "gli.Cluster8")
gli.cluster.100.9 <- cbind(gli.100.clist[[9]], Cluster = "gli.Cluster9")
gli.cluster.100.10 <- cbind(gli.100.clist[[10]], Cluster = "gli.Cluster10")
glioma_100_Genes_C <- rbind(gli.cluster.100.1, gli.cluster.100.2, gli.cluster.100.3, gli.cluster.100.4,
                           gli.cluster.100.5,gli.cluster.100.6,gli.cluster.100.7, gli.cluster.100.8, 
                           gli.cluster.100.9, gli.cluster.100.10)

glioma_100_Genes_C
#100 gene glioblastoma 10 cluster table
gbm.100.clist <- lapply(sort(unique(gbm.100.clusters)), function(x) r100.glioblastoma[which(gbm.100.clusters==x),])
gbm.cluster.100.1 <- cbind(gbm.100.clist[[1]], Cluster = "gbm.Cluster1")
gbm.cluster.100.2 <- cbind(gbm.100.clist[[2]], Cluster = "gbm.Cluster2")
gbm.cluster.100.3 <- cbind(gbm.100.clist[[3]], Cluster = "gbm.Cluster3")
gbm.cluster.100.4 <- cbind(gbm.100.clist[[4]], Cluster = "gbm.Cluster4")
gbm.cluster.100.5 <- cbind(gbm.100.clist[[5]], Cluster = "gbm.Cluster5")
gbm.cluster.100.6 <- cbind(gbm.100.clist[[6]], Cluster = "gbm.Cluster6")
gbm.cluster.100.7 <- cbind(gbm.100.clist[[7]], Cluster = "gbm.Cluster7")
gbm.cluster.100.8 <- cbind(gbm.100.clist[[8]], Cluster = "gbm.Cluster8")
gbm.cluster.100.9 <- cbind(gbm.100.clist[[9]], Cluster = "gbm.Cluster9")
gbm.cluster.100.10 <- cbind(gbm.100.clist[[10]], Cluster = "gbm.Cluster10")
glioblastoma_100_Genes_C <- rbind(gbm.cluster.100.1, gbm.cluster.100.2, gbm.cluster.100.3, gbm.cluster.100.4,
                            gbm.cluster.100.5,gbm.cluster.100.6,gbm.cluster.100.7, gbm.cluster.100.8, 
                            gbm.cluster.100.9, gbm.cluster.100.10)

View(Glioblastoma_100_Genes_C)
#save cluster table results to csvs (optional)
write.csv(glioblastoma_100_Genes_C, "Table_Glioblastoma_100genes_10WClust.csv")
write.csv(glioma_100_Genes_C, "Table_Glioma_100genes_10WClust.csv")

