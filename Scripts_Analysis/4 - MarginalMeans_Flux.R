
library(lme4)
library(emmeans)


#import experimental Design
#
ExperimentalDesign<-read.csv('Data/Experiment_Metadata.csv',row.names = 1)
#971,972,981,982 are all controls give them the same name
listofIsolates<-as.character(ExperimentalDesign$Isolate)


controls<-c("971", "972","981", "982")
ChangeTo <- c("controls","controls","controls","controls")
ids <- match(listofIsolates, controls)
ExperimentalDesign$Isolate<-replace(listofIsolates, which(!is.na(ids)), ChangeTo[ids[!is.na(ids)]])


#get a list of necessary priors to specifcy based on data

cols<-colnames(ExperimentalDesign)

#ensure metadata are factors
ExperimentalDesign[cols] <- lapply(ExperimentalDesign[cols], factor)  ## as.factor() could also be used



#relevel genotypes to make WT col.0 number 1 & controls 

ExperimentalDesign$HostGenotype <- relevel(ExperimentalDesign$HostGenotype, "col.0")  ## as.factor() could also be used


ExperimentalDesign$Isolate <- relevel(ExperimentalDesign$Isolate, "controls")  ## as.factor() could also be used



#import model fluxes
#
MODEL_FLUX<-as.data.frame(t(read.csv('Data/Model_fluxes.csv',header = T, row.names = 1)))
#Erase columns wiht 0 flux 
MODEL_FLUX<-MODEL_FLUX[, colSums(MODEL_FLUX) != 0]

#check distribution of biomass
hist(MODEL_FLUX$bio1_biomass)


#merge model flux with experimental Design 
#
Model_Experiment_dataframe <- merge(ExperimentalDesign, MODEL_FLUX, by = 0)


Model_Experiment_dataframe2<- Model_Experiment_dataframe[,-c(1:5,7,9:11)]

All_Flux_Average<-aggregate(. ~ HostGenotype + Isolate, data=Model_Experiment_dataframe2, FUN=mean)

dim(All_Flux_Average)

View(All_Flux_Average)
####################################

All_Flux_Means<-data.frame()

i<-12


FLUX_lm1 <- lm(rxn00001_c ~ Isolate + HostGenotype +0, data = Model_Experiment_dataframe)

#get marginal Means
Flux_Means.emm.n <- emmeans(FLUX_lm1, c("Isolate", "HostGenotype"))


Flux_Means <- as.data.frame(print(Flux_Means.emm.n))


Flux_Means_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),3])
Flux_Means_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),3])
Flux_Means_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),3])

paste0(unique(Flux_Means$HostGenotype)[1])
paste0(unique(Flux_Means$HostGenotype)[2])
paste0(unique(Flux_Means$HostGenotype)[3])

c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))

All_Flux_Average<-aggregate(levels_together ~ group_by_var + kernel_type, data=df, FUN=mean)
#create dataframe of isolates we will add genotype name to trait 
All_Flux_Means<-Flux_Means[,1]
All_Flux_SE<-Flux_Means[,1]
All_Flux_df<-Flux_Means[,1]
All_Flux_lower.CL<-Flux_Means[,1]
All_Flux_upper.CL<-Flux_Means[,1]

i<-12

for (i in 12:ncol(Model_Experiment_dataframe)) {
  
  #FLUX_lm1 <- lmer(Model_Experiment_dataframe[,i] ~  (1|Experiment/GrowingFlat/AgarFlat) + Isolate + HostGenotype +0, data = Model_Experiment_dataframe)
FLUX_lm1 <- lmer(Model_Experiment_dataframe$DM_cpd11416_c ~  (1|Experiment/GrowingFlat/AgarFlat) + Isolate + HostGenotype +0, data = Model_Experiment_dataframe)
  
  Model_Experiment_dataframe$bio1_biomass
  #get marginal Means
  Flux_Means.emm.n <- emmeans(FLUX_lm1, c("Isolate", "HostGenotype"))
  
  
  Flux_Means <- as.data.frame(print(Flux_Means.emm.n))
  
  Needdata<-Flux_Means%>%group_by(HostGenotype)%>%summarise(Average=mean(emmean))
  hist(Needdata$Average)
  
  #mean dataset #separate by host genotype 
  #
  Flux_Means_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),3])
  #rename the column
  colnames(Flux_Means_col.0) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[1])))
  All_Flux_Means <- cbind(All_Flux_Means, Flux_Means_col.0)
    
    
  Flux_Means_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),3])
  colnames(Flux_Means_coil.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[2])))
  All_Flux_Means <- cbind(All_Flux_Means, Flux_Means_coil.1)
  
  Flux_Means_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),3])
  colnames(Flux_Means_npr.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))
  All_Flux_Means <- cbind(All_Flux_Means, Flux_Means_npr.1)
  
  
  #
  #Standard error dataset
  
  #
  Flux_SE_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),4])
  #rename the column
  colnames(Flux_SE_col.0) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[1])))
  All_Flux_SE <- cbind(All_Flux_SE, Flux_SE_col.0)
  
  
  Flux_SE_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),4])
  colnames(Flux_SE_coil.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[2])))
  All_Flux_SE <- cbind(All_Flux_SE, Flux_SE_coil.1)
  
  Flux_SE_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),4])
  colnames(Flux_SE_npr.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))
  All_Flux_SE <- cbind(All_Flux_SE, Flux_SE_npr.1)
  
  
  
  
  
  
  #
  #
  #df dataset
  #

  #
  Flux_df_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),5])
  #rename the column
  colnames(Flux_df_col.0) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[1])))
  All_Flux_df <- cbind(All_Flux_df, Flux_df_col.0)
  
  
  Flux_df_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),5])
  colnames(Flux_df_coil.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[2])))
  All_Flux_df <- cbind(All_Flux_df, Flux_df_coil.1)
  
  Flux_df_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),5])
  colnames(Flux_df_npr.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))
  All_Flux_df <- cbind(All_Flux_df, Flux_df_npr.1)
  
  
  
  
  #
  #
  #lower.CL dataset
  #
  Flux_lower.CL_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),6])
  #rename the column
  colnames(Flux_lower.CL_col.0) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[1])))
  All_Flux_lower.CL <- cbind(All_Flux_lower.CL, Flux_lower.CL_col.0)
  
  
  Flux_lower.CL_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),6])
  colnames(Flux_lower.CL_coil.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[2])))
  All_Flux_lower.CL <- cbind(All_Flux_lower.CL, Flux_lower.CL_coil.1)
  
  Flux_lower.CL_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),6])
  colnames(Flux_lower.CL_npr.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))
  All_Flux_lower.CL <- cbind(All_Flux_lower.CL, Flux_lower.CL_npr.1)
  
  #
  #
  #upper.CL dataset
  
  Flux_upper.CL_col.0<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("col.0"),6])
  #rename the column
  colnames(Flux_upper.CL_col.0) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[1])))
  All_Flux_upper.CL <- cbind(All_Flux_upper.CL, Flux_upper.CL_col.0)
  
  
  Flux_upper.CL_coil.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("coil.1"),6])
  colnames(Flux_upper.CL_coil.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[2])))
  All_Flux_upper.CL <- cbind(All_Flux_upper.CL, Flux_upper.CL_coil.1)
  
  Flux_upper.CL_npr.1<- data.frame(Flux_Means[Flux_Means$HostGenotype == c("npr.1"),7])
  colnames(Flux_upper.CL_npr.1) <-c(paste0(gsub('_','',colnames(Model_Experiment_dataframe)[i]),"_",paste0(unique(Flux_Means$HostGenotype)[3])))
  
  #
  #
  #
  #
  
  
}

system("say Just finished!")






#bind the isolate names
#
#
IsolateKey<-read.csv("GenomicData/IsolateKey.csv")
#ensure things are characters
IsolateKey$Sample<- as.character(IsolateKey$Sample)

All_Flux_Means$All_Flux_Means<- as.character(All_Flux_Means$All_Flux_Means)
colnames(All_Flux_Means)[1]<-'Isolate'
#swapnames
All_Flux_Means$Isolate<-IsolateKey$Isolate[match(All_Flux_Means$Isolate,IsolateKey$Sample)]


write.csv(All_Flux_Means,"Data/arabidopsis_flux_means.csv")




All_Flux_SE$All_Flux_SE<- as.character(All_Flux_SE$All_Flux_SE)
colnames(All_Flux_SE)[1]<-'Isolate'
#swapnames
All_Flux_SE$Isolate<-IsolateKey$Isolate[match(All_Flux_SE$Isolate,IsolateKey$Sample)]


write.csv(All_Flux_SE,"Data/arabidopsis_All_Flux_SE.csv")


############################
All_Flux_df$All_Flux_df<- as.character(All_Flux_df$All_Flux_df)
colnames(All_Flux_df)[1]<-'Isolate'
#swapnames
All_Flux_df$Isolate<-IsolateKey$Isolate[match(All_Flux_df$Isolate,IsolateKey$Sample)]


write.csv(All_Flux_df,"Data/arabidopsis_All_Flux_df.csv")
###########################
###########################
###########################
All_Flux_lower.CL$All_Flux_lower.CL<- as.character(All_Flux_lower.CL$All_Flux_lower.CL)
colnames(All_Flux_lower.CL)[1]<-'Isolate'
#swapnames
All_Flux_lower.CL$Isolate<-IsolateKey$Isolate[match(All_Flux_lower.CL$Isolate,IsolateKey$Sample)]


write.csv(All_Flux_lower.CL,"Data/arabidopsis_All_Flux_lower.CL.csv")
###########################
###########################
###########################

All_Flux_upper.CL$All_Flux_upper.CL<- as.character(All_Flux_upper.CL$All_Flux_upper.CL)
colnames(All_Flux_upper.CL)[1]<-'Isolate'
#swapnames
All_Flux_upper.CL$Isolate<-IsolateKey$Isolate[match(All_Flux_upper.CL$Isolate,IsolateKey$Sample)]


write.csv(All_Flux_upper.CL,"Data/arabidopsis_All_Flux_upper.CL.csv")



write.csv(All_Flux_Average_22,"Data/arabidopsis_All_Flux_Average.csv")
All_Flux_Average_22<-reshape(All_Flux_Average, idvar = "Isolate", timevar = "HostGenotype", direction = "wide" ,sep = "...")


All_Flux_Average_22$Isolate<-IsolateKey$Isolate[match(All_Flux_Average_22$Isolate,IsolateKey$Sample)]
colnames(All_Flux_Average_22)<-sub("_", "", colnames(All_Flux_Average_22), fixed = TRUE)
colnames(All_Flux_Average_22)<-sub("...", "_", colnames(All_Flux_Average_22), fixed = TRUE)

View(data.frame(All_Flux_Average_22[1:5,1:5]))
length(All_Flux_Average_22$HostGenotype)
