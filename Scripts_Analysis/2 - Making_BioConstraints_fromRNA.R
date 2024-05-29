###### getting the host expression data in a Useable form 
###### 
#Library necessary packages
library(data.table)
#
#import Modelgenelist

Model_genelist<-read.csv("Data/Model_genelist.csv",header = F)

#any necessary tair translations
Tair_translation<-read.csv("Data/Tair_translation.csv")



#create a genelist we will use the trim the RNA dataset

Enzyme_geneNames<-as.data.frame(unlist(c(Model_genelist,Tair_translation[,2])))
#hostRNA 
HOST_RNA<-read.csv("Data/Host_RNAseq_Counts.csv")



#delete non matching columns host camalexin was previously screened for this step
HOST_RNA<-HOST_RNA[,c("gene",intersect(colnames(HOST_RNA), HOST_Camalexin$SampleName))]


#segment Gene names to query rows. 

#segment Gene names to query rows. 

TranscriptNames<-HOST_RNA[,1]
#remove the .1 variations

TranscriptNames<-as.data.frame(sub("\\..*", "", HOST_RNA[,1]))
colnames(TranscriptNames)<-"Gene-id"
#find gene names in model not in host_rna
#missing model genes
MissingTranslations<-MissingModelGenese<-setdiff(Tair_translation$OriginalName,Model_genelist$V1)
#this should be 0 all the genes have a tair translation sometimes the RNA data is funky we need to change the host RNAdata to the new translation
##find gene names in model not in host_rna
MissingModelGenese<-setdiff(Enzyme_geneNames[,1],TranscriptNames[,1])
#if missing genes are present you need to change the RNA data names

MissingModelGenes_translation<- Tair_translation[which(Tair_translation$OriginalName %in% MissingModelGenese),]

#replace HostRNA names & transcript names
library(plyr)
TranscriptNames[,1] <- mapvalues(TranscriptNames[,1], from=MissingModelGenes_translation$TairName, to=MissingModelGenes_translation$OriginalName)

#change HostRNA names
HOST_RNA$gene<-TranscriptNames$`Gene-id`

#subset host RNA to only biosynthetic enzyme in your gene list
BiosyntheticRNA<-HOST_RNA[which(TranscriptNames[,1] %in% Enzyme_geneNames[,1]),]



#consolidate all isoforms into a single gene
BiosyntheticRNA_ready<-as.data.table(BiosyntheticRNA)[, lapply(.SD, sum), by = gene]

#change NA to 0 
#
BiosyntheticRNA_ready[is.na(BiosyntheticRNA_ready)] <- 0




#export data 
write.csv(BiosyntheticRNA_ready,"Data/BiosyntheticConstraints.csv",  row.names = F)

####Creating Camalexin Constraints

ListofSamples<-as.data.frame(colnames(BiosyntheticRNA)[-1])



#hostCamalexin
HOST_Camalexin<-read.csv("Data/Camalexin_Constraints.csv")

View(HOST_Camalexin)

#Trim RNA based on camlexin?

#should be 0
setdiff( unique(HOST_Camalexin$SampleName),colnames(BiosyntheticRNA_ready))
#should be Gene
setdiff(colnames(BiosyntheticRNA_ready),unique(HOST_Camalexin$SampleName))
#if these are fine you are good to go