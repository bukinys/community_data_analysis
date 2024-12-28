
#connecting the working directory
setwd("G:/Data")


#reading phytoplankton species abundance data: phyto_species_abundance_class_names.tsv - tab-delimited text file

data<-read.table("phyto_species_abundance_class_names.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]

Abundance<-as.data.frame(data[,-ncol(data)])
Abundance_class<-data
Abundance_class<-aggregate(.~Class, Abundance_class, sum)

#Abundance_class - data frame with abundance of phytoplankton classes


#reading phytoplankton species biomass data: phyto_species_biomass_class_names.tsv - tab-delimited text file

data<-read.table("phyto_species_biomass_class_names.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]

Biomass<-as.data.frame(data[,-ncol(data)])
Biomass_class<-data
Biomass_class<-aggregate(.~Class, Biomass_class, sum)
#Biomass_classs - data frame with biomass of phytoplankton classes


#reading phytoplankton ASV read counts data: ASV_phyto_class_names.tsv - tab-delimited text file

data<-read.table("ASV_phyto_class_names.tsv",header=TRUE, sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]


#converting phytoplankton ASV read counts to relative abundance of phytoplankton ASV per class

S<-colSums(data[,-ncol(data)])
for(i in 1:ncol(data[,-ncol(data)])) data[,i]<-data[,i]/S[i]
colSums(data[,-ncol(data)])

ASV<-data
ASV<-ASV[order(ASV$Class, decreasing=F),]

Class<-unique(ASV$Class)

for(i in 1:length(Class))
 {
  n<-which(ASV$Class %in% Class[i])
  S<-colSums(ASV[n, -ncol(ASV)])
  for(i in 1:ncol(ASV[n, -ncol(ASV)])) ASV[n, i]<-ASV[n, i]/S[i]
  colSums(ASV[n,-ncol(ASV)])
 }

ASV[is.na(ASV)]<-0

#ASV - data frame with the relative abundance of phytoplankton ASV per class


#/////////////////////////////////////////////////////////////////////
#converting relative abundance of ASV into class-specific absolute abundance of ASV

ASV_Abundance<-ASV

length(unique(ASV$Class))
length(unique(Abundance_class$Class))

n<-which(Abundance_class$Class %in% ASV_Abundance$Class)
length(n)

for(i in 1:nrow(Abundance_class))
 {
  n<-which(ASV_Abundance$Class %in% Abundance_class$Class[i])
  S<-Abundance_class[i, -1]
  for(j in 1:length(n)) ASV_Abundance[n[j], -ncol(ASV_Abundance)]<-S*ASV_Abundance[n[j], -ncol(ASV_Abundance)]
 }

colSums(Abundance_class[,-1])
colSums(ASV_Abundance[, -ncol(ASV_Abundance)])

#ASV_Abundance - data frame with class-specific absolute abundance of ASV

#saving data of class-specific absolute abundance of ASV
write.table(cbind(id=rownames(ASV_Abundance), ASV_Abundance), "ASV_Abundance_Class.tsv", sep="\t", row.names=F)


#/////////////////////////////////////////////////////////////////////
#converting relative abundance of ASV into class-specific absolute biomass of ASV

ASV_Biomass<-ASV

length(unique(ASV$Class))
length(unique(Biomass_class$Class))

n<-which(Biomass_class$Class %in% ASV_Biomass$Class)
length(n)

for(i in 1:nrow(Biomass_class))
 {
  n<-which(ASV_Biomass$Class %in% Biomass_class$Class[i])
  S<-Biomass_class[i, -1]
  for(j in 1:length(n)) ASV_Biomass[n[j], -ncol(ASV_Biomass)]<-S*ASV_Biomass[n[j], -ncol(ASV_Biomass)]
 }

colSums(Biomass_class[,-1])
colSums(ASV_Biomass[, -ncol(ASV_Biomass)])

#ASV_Biomass - data frame with class-specific absolute biomass of ASV

#saving data of class-specific absolute biomass of ASV
write.table(cbind(id=rownames(ASV_Biomass), ASV_Biomass), "ASV_Biomass_Class.tsv", sep="\t", row.names=F)
