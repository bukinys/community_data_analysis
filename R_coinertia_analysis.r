
library(vegan)
library(ade4)
library(compositions)

#Setting the working directory and path

setwd("D:\\ivan\\eu_3south\\coinertia")

#///////////////////////////////////////

#reading file with read counts of ASV phytoplankton

data<-read.table("ASV_phyto.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]

#transformation of read counts of ASV phytoplankton to relative abundance of ASV

ASV_relative<-data
colSums(ASV_relative)

S<-colSums(ASV_relative)
for(i in 1:ncol(ASV_relative)) ASV_relative[,i]<-ASV_relative[,i]/S[i]
colSums(ASV_relative)
ASV_relative<-as.matrix(ASV_relative)

#clr-transformation of the relative abundance of ASV

ASV_clr<-ASV_relative

for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])

#///////////////////////////////////////

#reading a file with the absolute abundance of phytoplankton species

data<-read.table("phyto_south_species_abundance.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]

Abundance<-as.data.frame(data)
Abundance<-as.matrix(Abundance)
colSums(Abundance)

#transformation of absolute abundance of phytoplankton species into their relative abundance

Abundance_relative<-Abundance
S<-colSums(Abundance_relative)
Sm<-mean(S)
Sm<-mean(Abundance_relative)
for(i in 1:ncol(Abundance_relative)) Abundance_relative[,i]<-Abundance_relative[,i]/S[i]
colSums(Abundance_relative)

#clr-transformation of relative abundance of phytoplankton species

Abundance_clr<-Abundance_relative

for(i in 1:ncol(Abundance_clr)) Abundance_clr[,i]<-clr(Abundance_clr[,i])

#///////////////////////////////////////

#reading the phytoplankton species biomass file

data<-read.table("phyto_south_species_biomass.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-c(1, ncol(data))]

Biomass<-as.data.frame(data)
Biomass<-as.matrix(Biomass)
colSums(Biomass)

#transformation of absolute biomass of phytoplankton species in their relative biomass

Biomass_relative<-Biomass
S<-colSums(Biomass_relative)
Sm<-mean(S)
Sm<-mean(Biomass_relative)
for(i in 1:ncol(Biomass_relative)) Biomass_relative[,i]<-Biomass_relative[,i]/S[i]
colSums(Biomass_relative)


#clr-transformation of relative biomass of phytoplankton species

Biomass_clr<-Biomass_relative

for(i in 1:ncol(Biomass_clr)) Biomass_clr[,i]<-clr(Biomass_clr[,i])

#transformation of relative abundance of ASV to absolute abundance of ASV based on total phytoplankton abundance from microscopy

ASV_Abundance<-ASV_relative
S<-colSums(Abundance)
S
for(i in 1:ncol(ASV_Abundance)) ASV_Abundance[,i]<-S[i]*ASV_Abundance[,i]
colSums(ASV_Abundance)

#transformation of relative abundance of ASV to absolute biomass of ASV based on total phytoplankton biomass by microscopy

ASV_Biomass<-ASV_relative
S<-colSums(Biomass)
S
for(i in 1:ncol(ASV_Biomass)) ASV_Biomass[,i]<-S[i]*ASV_Biomass[,i]
colSums(ASV_Biomass)


#co-inertia analysis

#co-inertial analysis for relative values

dudi_ASV_rel <- dudi.pca(t(ASV_relative), scale = FALSE, scan = FALSE, nf = 3)
dudi_Abundance_rel<- dudi.pca(t(Abundance_relative), scale = FALSE, scan = FALSE, nf = 2)
dudi_Biomass_rel<- dudi.pca(t(Biomass_relative), scale = FALSE, scan = FALSE, nf = 2)

coin<- coinertia(dudi_ASV_rel, dudi_Abundance_rel, scan = FALSE, nf = 2)
coin$RV

coin<- coinertia(dudi_ASV_rel, dudi_Biomass_rel, scan = FALSE, nf = 2)
coin$RV

#co-inertia analysis for clr-transformed relative values

dudi_ASV_clr <- dudi.pca(t(ASV_clr), scale = FALSE, scan = FALSE, nf = 3)
dudi_Abundance_clr<- dudi.pca(t(Abundance_clr), scale = FALSE, scan = FALSE, nf = 2)
dudi_Biomass_clr<- dudi.pca(t(Biomass_clr), scale = FALSE, scan = FALSE, nf = 2)

coin<- coinertia(dudi_ASV_clr, dudi_Abundance_clr, scan = FALSE, nf = 2)
coin$RV

coin<- coinertia(dudi_ASV_clr, dudi_Biomass_clr, scan = FALSE, nf = 2)
coin$RV

#coineration analysis for relative values transformed by ranking from 0 to 1

dudi_ASV_rel <- dudi.pca(decostand(t(ASV_relative), method="range", MARGIN=2), scale = FALSE, scan = FALSE, nf = 3)
dudi_Abundance_rel<- dudi.pca(decostand(t(Abundance_relative), method="range", MARGIN=2), scale = FALSE, scan = FALSE, nf = 2)
dudi_Biomass_rel<- dudi.pca(decostand(t(Biomass_relative), method="range", MARGIN=2), scale = FALSE, scan = FALSE, nf = 2)

coin<- coinertia(dudi_ASV_rel, dudi_Abundance_rel, scan = FALSE, nf = 2)
coin$RV

coin<- coinertia(dudi_ASV_rel, dudi_Biomass_rel, scan = FALSE, nf = 2)
coin$RV

#Coinertia analysis for absolute values of ASV and species

dudi_ASV_Abundance <- dudi.pca(t(ASV_Abundance), scale = FALSE, scan = FALSE, nf = 3)
dudi_ASV_Biomass <- dudi.pca(t(ASV_Biomass), scale = FALSE, scan = FALSE, nf = 3)
dudi_Abundance<- dudi.pca(t(Abundance), scale = FALSE, scan = FALSE, nf = 2)
dudi_Biomass<- dudi.pca(t(Biomass), scale = FALSE, scan = FALSE, nf = 2)


coin<- coinertia(dudi_ASV_Abundance, dudi_Abundance, scan = FALSE, nf = 2)
coin$RV

coin<- coinertia(dudi_ASV_Biomass, dudi_Biomass, scan = FALSE, nf = 2)
coin$RV
