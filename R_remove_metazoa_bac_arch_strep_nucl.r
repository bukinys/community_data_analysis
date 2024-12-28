#Script to remove Metazoa, Bacteria, Archaea, Streptophyta, Eukaryota:nucl, NA


library(tibble)

#Setting the working directory and path

setwd("D:\\ivan\\eu_3south\\data_filtr")


#reading, preparing and transforming eu_3south_seqtab.nochim.txt. Before reading, manually add the "sample" identifier to the first line of eu_3south_seqtab.nochim.txt

data<-read.csv("eu_3south_seqtab.nochim.txt",header=TRUE,sep="\t")

rownames(data)<-data[,1]
data<-data[,-1]

col_names<-paste("ASV_", c(1:ncol(data)), sep="")

colnames(data)<-col_names
data<-t(data)
data<-as.data.frame(data)
total<-rowSums(data)
data<-cbind(data, total=total)
data<-add_column(data, ASV_id=col_names, .before = 1)


#reading and preparing taxa_eu_3south_pr2.txt. Before reading, manually add the identifier "representative_sequence" to the first line of taxa_eu_3south_pr2.txt

tax_pr2<-read.csv("taxa_eu_3south_pr2.txt", header=TRUE,sep="\t") 
tax_pr2<-add_column(tax_pr2, ASV_id=col_names, .before = 1)
colnames(tax_pr2)[1]<-"ASV_id"


#reading and preparing taxa_eu_3south_silva.txt. Before reading, manually add the identifier "representative_sequence" to the first line of taxa_eu_3south_silva.txt

tax_silva<-read.csv("taxa_eu_3south_silva.txt", header=TRUE,sep="\t")
tax_silva<-add_column(tax_silva, ASV_id=col_names, .before = 1)
colnames(tax_silva)[1]<-"ASV_id"


#filtering the dataset by the number of reads. ASVs with read count less than or equal to 2 are deleted

data<-data[total>=3,]
tax_pr2<-tax_pr2[total>=3,]
tax_silva<-tax_silva[total>=3,]


#selection of ASV Metazoa and saving data on them in read tables, taxonomy tables and sequences in fasta format and storing ASV Metazoa indices
 

tax_line<-paste(tax_pr2$Domain, tax_pr2$Supergroup, tax_pr2$Division, tax_pr2$Subdivision, tax_pr2$Class, tax_pr2$Order, tax_pr2$Family, tax_pr2$Genus, tax_pr2$Species, sep=";")

n_Metazoa<-which(grepl("Metazoa", tax_line))
if(length(n_Metazoa)==0) print("no metazoa data")

if(length(n_Metazoa)!=0)
 {
  data_Metazoa<-data[n_Metazoa,]
  tax_pr2_Metazoa<-tax_pr2[n_Metazoa,]
  tax_silva_Metazoa<-tax_silva[n_Metazoa,]

  wcon<-file("rep_ASV_Metazoa.fasta", open ="w")
  for(i in 1:nrow(data_Metazoa))
   {
    writeLines(paste(">", data_Metazoa$ASV_id[i]), wcon)
    writeLines(tax_pr2_Metazoa$representative_sequence[i], wcon)
   }
  close(wcon)

  write.table(data_Metazoa, "ASV_Metazoa_table.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
  write.table(tax_pr2_Metazoa, "ASV_Metazoa_tax_pr2.txt", quote=F, row.names=F, col.names=T, sep = "\t")
  write.table(tax_silva_Metazoa, "ASV_Metazoa_tax_silva.txt", quote=F, row.names=F, col.names=T, sep = "\t")
 } 
 

#selection of ASV Bacteria and saving data on them in read tables, taxonomy tables and sequences in fasta format and storing ASV Bacteria indices
 

tax_line<-paste(tax_pr2$Domain, tax_pr2$Supergroup, tax_pr2$Division, tax_pr2$Subdivision, tax_pr2$Class, tax_pr2$Order, tax_pr2$Family, tax_pr2$Genus, tax_pr2$Species, sep=";")

n_bacteria<-which(grepl("Bacteria", tax_line))
if(length(n_bacteria)==0) print("no bacteria data")

if(length(n_bacteria)!=0)
 {
  data_bacteria<-data[n_bacteria,]
  tax_pr2_bacteria<-tax_pr2[n_bacteria,]
  tax_silva_bacteria<-tax_silva[n_bacteria,]

   wcon<-file("rep_ASV_bacteria.fasta", open ="w")
  for(i in 1:nrow(data_bacteria))
   {
    writeLines(paste(">", data_bacteria$ASV_id[i]), wcon)
    writeLines(tax_pr2_bacteria$representative_sequence[i], wcon)
   }
  close(wcon)

  write.table(data_bacteria, "ASV_bacteria_table.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
  write.table(tax_pr2_bacteria, "ASV_bacteria_tax_pr2.txt", quote=F, row.names=F, col.names=T, sep = "\t")
  write.table(tax_silva_bacteria, "ASV_bacteria_tax_silva.txt", quote=F, row.names=F, col.names=T, sep = "\t")
 } 
 


#selection of ASV Archaea and saving data on them in read tables, taxonomy tables and sequences in fasta format and storing ASV Archaea indices


tax_line<-paste(tax_pr2$Domain, tax_pr2$Supergroup, tax_pr2$Division, tax_pr2$Subdivision, tax_pr2$Class, tax_pr2$Order, tax_pr2$Family, tax_pr2$Genus, tax_pr2$Species, sep=";")

n_archaea<-which(grepl("Archaea", tax_line))
if(length(n_archaea)==0) print("no archaea data")

if(length(n_archaea)!=0)
 {
  data_archaea<-data[n_archaea,]
  tax_pr2_archaea<-tax_pr2[n_archaea,]
  tax_silva_archaea<-tax_silva[n_archaea,]

   wcon<-file("rep_ASV_archaea.fasta", open ="w")
  for(i in 1:nrow(data_archaea))
   {
    writeLines(paste(">", data_archaea$ASV_id[i]), wcon)
    writeLines(tax_pr2_archaea$representative_sequence[i], wcon)
   }
  close(wcon)

  write.table(data_archaea, "ASV_archaea_table.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
  write.table(tax_pr2_archaea, "ASV_archaea_tax_pr2.txt", quote=F, row.names=F, col.names=T, sep = "\t")
  write.table(tax_silva_archaea, "ASV_archaea_tax_silva.txt", quote=F, row.names=F, col.names=T, sep = "\t")
 } 


#selection of ASV Streptophyta and saving data on them in read tables, taxonomy tables and sequences in fasta format and storing ASV Streptophyta indices
 

tax_line<-paste(tax_pr2$Domain, tax_pr2$Supergroup, tax_pr2$Division, tax_pr2$Subdivision, tax_pr2$Class, tax_pr2$Order, tax_pr2$Family, tax_pr2$Genus, tax_pr2$Species, sep=";")

n_streptophyta<-which(grepl("Streptophyta", tax_line))
if(length(n_streptophyta)==0) print("no streptophyta data")

if(length(n_streptophyta)!=0)
 {
  data_streptophyta<-data[n_streptophyta,]
  tax_pr2_streptophyta<-tax_pr2[n_streptophyta,]
  tax_silva_streptophyta<-tax_silva[n_streptophyta,]

  wcon<-file("rep_ASV_streptophyta.fasta", open ="w")
  for(i in 1:nrow(data_streptophyta))
   {
    writeLines(paste(">", data_streptophyta$ASV_id[i]), wcon)
    writeLines(tax_pr2_streptophyta$representative_sequence[i], wcon)
   }
  close(wcon)

  write.table(data_streptophyta, "ASV_streptophyta_table.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
  write.table(tax_pr2_streptophyta, "ASV_streptophyta_tax_pr2.txt", quote=F, row.names=F, col.names=T, sep = "\t")
  write.table(tax_silva_streptophyta, "ASV_streptophyta_tax_silva.txt", quote=F, row.names=F, col.names=T, sep = "\t")

  }   


#selection of ASV Eukaryota:nucl (Nucleomorph) and saving data on them in read tables, taxonomy tables and sequences in fasta format and storing ASV Eukaryota:nucl indices


tax_line<-paste(tax_pr2$Domain, tax_pr2$Supergroup, tax_pr2$Division, tax_pr2$Subdivision, tax_pr2$Class, tax_pr2$Order, tax_pr2$Family, tax_pr2$Genus, tax_pr2$Species, sep=";")

n_eukaryota_nucl<-which(grepl("Eukaryota:nucl", tax_line))
if(length(n_eukaryota_nucl)==0) print("no eukaryota:nucl data")

if(length(n_eukaryota_nucl)!=0)
 {
  data_eukaryota_nucl<-data[n_eukaryota_nucl,]
  tax_pr2_eukaryota_nucl<-tax_pr2[n_eukaryota_nucl,]
  tax_silva_eukaryota_nucl<-tax_silva[n_eukaryota_nucl,]

  wcon<-file("rep_ASV_eukaryota_nucl.fasta", open ="w")
  for(i in 1:nrow(data_eukaryota_nucl))
   {
    writeLines(paste(">", data_eukaryota_nucl$ASV_id[i]), wcon)
    writeLines(tax_pr2_eukaryota_nucl$representative_sequence[i], wcon)
   }
  close(wcon)

  write.table(data_eukaryota_nucl, "ASV_eukaryota_nucl_table.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
  write.table(tax_pr2_eukaryota_nucl, "ASV_eukaryota_nucl_tax_pr2.txt", quote=F, row.names=F, col.names=T, sep = "\t")
  write.table(tax_silva_eukaryota_nucl, "ASV_eukaryota_nucl_tax_silva.txt", quote=F, row.names=F, col.names=T, sep = "\t")

  }      



#NA removal

n_na<-which(is.na(tax_pr2$Domain))


#microeukaryotic data extraction and saving data on them in read tables, taxonomy tables and sequences in fasta format without Metazoa, Bacteria, Archaea, Streptophyta, Eukaryota:nucl, NA

n<-c(n_Metazoa, n_bacteria, n_archaea, n_streptophyta, n_eukaryota_nucl, n_na)

data<-data[-n,]
tax_pr2<-tax_pr2[-n,]
tax_silva<-tax_silva[-n,]

write.table(data, "ASV_table_filtr.txt", quote=F, row.names=F, col.names=T, sep = "\t") 
write.table(tax_pr2, "ASV_tax_pr2_filtr.txt", quote=F, row.names=F, col.names=T, sep = "\t")
write.table(tax_silva, "ASV_tax_silva_filtr.txt", quote=F, row.names=F, col.names=T, sep = "\t")

wcon<-file("rep_ASV_filtr.fasta", open ="w")
for(i in 1:nrow(data))
 {
  writeLines(paste(">", data$ASV_id[i]), wcon)
  writeLines(tax_pr2$representative_sequence[i], wcon)
 }
close(wcon)