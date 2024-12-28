
#connecting the working directory

setwd("G:/Data")


#reading phytoplankton species abundance data: File_S1.tsv - tab-delimited text file

data<-read.table("phyto_species_abundance_class_names.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]

Abundance<-data #Abundance - data frame with the phytoplankton species abundance


#reading phytoplankton species biomass data: File_S2.tsv - tab-delimited text file

data<-read.table("phyto_species_biomass_class_names.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]

Biomass<-data #Biomass - data frame with the phytoplankton species biomass


#reading phytoplankton ASV read counts data: File_S3.tsv - tab-delimited text file

data<-read.table("ASV_phyto_class_names.tsv",header=TRUE, sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]


#converting phytoplankton ASV read counts to relative abundance of phytoplankton ASV

S<-colSums(data)
for(i in 1:ncol(data)) data[,i]<-data[,i]/S[i]
colSums(data)

ASV<-data #ASV - data frame with the relative abundance of phytoplankton ASV


#/////////////////////////////////////////////////////////////////////
#converting relative abundance of ASV into absolute abundance of ASV

ASV_Abundance<-ASV


#Суммирование общей численности видов фитопланктона в образцах
S<-colSums(Abundance)

for(i in 1:nrow(ASV_Abundance)) ASV_Abundance[,i]<-S[i]*ASV_Abundance[,i]

colSums(Abundance)
colSums(ASV_Abundance)

#ASV_Abundance - data frame with absolute abundance of ASV

#saving data of absolute abundance of ASV

write.table(cbind(id=rownames(ASV_Abundance), ASV_Abundance), "ASV_Abundance.tsv", sep="\t", row.names=F)


#/////////////////////////////////////////////////////////////////////
#converting relative abundance of ASV into absolute biomass of ASV

ASV_Biomass<-ASV

#Суммирование общей биомассы видов фитопланктона в образцах
S<-colSums(Biomass)

for(i in 1:nrow(ASV_Biomass)) ASV_Biomass[,i]<-S[i]*ASV_Biomass[,i]

#ASV_Biomass - data frame with absolute biomass of ASV

colSums(Biomass)
colSums(ASV_Biomass)

#saving data of absolute biomass of ASV

write.table(cbind(id=rownames(ASV_Biomass), ASV_Biomass), "ASV_Biomass.tsv", sep="\t", row.names=F)
