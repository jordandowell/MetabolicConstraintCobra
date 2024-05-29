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

#hostCamalexin
HOST_Camalexin<-read.csv("Data/Camalexin_Constraints.csv")

#delete non matching columns host camalexin was previously screened for this step
HOST_RNA<-HOST_RNA[,c("gene",intersect(colnames(HOST_RNA), HOST_Camalexin$SampleName))]


#segment Gene names to query rows. 



#segment Gene names to query rows. 

TranscriptNames<-HOST_RNA[,1]
#remove the .1 variations

TranscriptNames<-as.data.frame(sub("\\..*", "", HOST_RNA[,1]))
colnames(TranscriptNames)<-"Gene-id"


#subset host RNA to only biosynthetic enzyme in your gene list
BiosyntheticRNA<-HOST_RNA[which(TranscriptNames[,1] %in% Enzyme_geneNames[,1]),]

#remake the gene-id label so that we remove isoform
#
BiosyntheticRNA$gene<-as.data.frame(sub("\\..*", "", BiosyntheticRNA[,1]))

colnames(BiosyntheticRNA)[1]<-"gene"

#consolidate all isoforms into a single gene
BiosyntheticRNA_ready<-as.data.table(BiosyntheticRNA)[, lapply(.SD, sum), by = gene]

#change NA to 0 
#
BiosyntheticRNA_ready[is.na(BiosyntheticRNA_ready)] <- 0

#convert data to model approved names
#
#
#Replace column value with another based on Tair condition
#
i<-1
for (i in 1:length(Tair_translation)) {
  BiosyntheticRNA_ready$gene[BiosyntheticRNA_ready$gene == Tair_translation[i,2]] <- Tair_translation[i,1]
  
}

#export data 
write.csv(BiosyntheticRNA_ready,"Data/BiosyntheticConstraints.csv",  row.names = F)

