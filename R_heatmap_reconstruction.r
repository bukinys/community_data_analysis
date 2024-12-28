
library(vegan)
library(gplots)
library(compositions)


setwd("D:\\ivan\\eu_3south\\heatmap")


#///////////////////////////////////////

#reading a file with the absolute abundance of phytoplankton classes (or species)

data<-read.table("phyto_Class_abundance.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
ncol(data)

Abundance<-as.data.frame(data)
Abundance<-as.matrix(Abundance)

#transformation of absolute abundance of phytoplankton classes into their relative abundance

Abundance_relative<-Abundance
S<-colSums(Abundance_relative)
Sm<-mean(S)
Sm<-mean(Abundance_relative)
for(i in 1:ncol(Abundance_relative)) Abundance_relative[,i]<-Abundance_relative[,i]/S[i]
colSums(Abundance_relative)

#clr-transformation of relative abundance of phytoplankton classes

Abundance_clr<-Abundance_relative
for(i in 1:ncol(Abundance_clr)) Abundance_clr[,i]<-clr(Abundance_clr[,i])

#///////////////////////////////////////

#reading the phytoplankton species biomass file

data<-read.table("phyto_Class_biomass.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
ncol(data)

Biomass<-as.data.frame(data)
Biomass<-as.matrix(Biomass)

#transformation of absolute biomass of phytoplankton classes in their relative biomass

Biomass_relative<-Biomass
S<-colSums(Biomass_relative)
Sm<-mean(S)
Sm<-mean(Biomass_relative)
for(i in 1:ncol(Biomass_relative)) Biomass_relative[,i]<-Biomass_relative[,i]/S[i]
colSums(Biomass_relative)

#clr-transformation of relative biomass of phytoplankton classes

Biomass_clr<-Biomass_relative
for(i in 1:ncol(Biomass_clr)) Biomass_clr[,i]<-clr(Biomass_clr[,i])

#///////////////////////////////////////

#reading file with read counts of ASV phytoplankton

data<-read.table("ASV_class_phyto.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
ncol(data)

S<-colSums(data)
Sm<-mean(S)
for(i in 1:ncol(data)) data[,i]<-data[,i]/S[i]
colSums(data)

#transformation of read counts of ASV phytoplankton to relative abundance of ASV

ASV_relative<-as.data.frame(data)
ASV_relative<-as.matrix(ASV_relative)

#clr-transformation of relative abundance of ASV

ASV_clr<-ASV_clr
for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])


#//////////////////////////////////construction of heat maps based on correlation coefficients

#ASV_relative - Abundance_relative; Selecting a pair of tables for calculating pairwise correlation coefficients and transferring these tables to intermediate variables

d_row<-ASV_relative       
d_col<-Abundance_relative 

#creating an empty matrix to store correlation coefficients
DD=matrix(nrow=nrow(d_row), ncol=nrow(d_col))

rownames(DD)<-rownames(d_row)
colnames(DD)<-rownames(d_col)

#creating an empty matrix to store p_value of correlation coefficients 
DP<-DD

#Calculation of correlation coefficients with replacement of unreliable correlation coefficients by 0
for(i in 1:nrow(d_row))
 {
   for(j in 1:nrow(d_col))
    {
     R=cor.test(d_row[i,], d_col[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-R$p.value
	 DP[i,j]<-p.adjust(DP[i,j], method="BH", n=nrow(d_row))
	 if(DP[i,j]>0.05) DD[i,j]=0
    }
 }

#heat map construction with row and column clustering(Euclidean distance metric and the complete-link clustering method)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=25)

dfun<-function(x)
 {
  return(vegdist(x, method="euclidean"))
 }
heatmap.2(DD, Rowv=T, Colv=T, distfun=dfun, trace="none", col=my_palette, margins = c(10, 10), cexCol=0.8, cexRow=0.8)

#saving correlation coefficients to a text table
write.table(cbind(sp=rownames(DD),as.data.frame(DD)), "cor_class_ASV_relative_Abundance_relative.tsv", sep="\t", row.names=F)