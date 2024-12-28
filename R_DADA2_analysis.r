#Sequence (read) analysis using DADA2


library(dada2)


#Setting the working directory and path

setwd("D:/ivan/eu_3south/fastq_data")
path <- "D:/ivan/eu_3south/fastq_data"

#Read the names of the fastq files and get matched lists of the forward and reverse fastq files
#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Place filtered files in filtered/ subdirectory 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


#Filter and trim the reads 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(21, 20), truncLen=c(250,165), maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)

#trimLeft=c(21, 20) - the number of nucleotides to remove from the start of each read (remove primer regions). Forward primer - 21, reverse primer - 20.
#truncLen=c(250, 165) - truncate reads after truncLen bases (remove low quality tails). Reads shorter than this are discarded. Forward reads is 250 nucleotides and reverse reads is 165 nucleotides
#maxN=0 - Default 0. After truncation, sequences with more than maxN Ns will be discarded (remove read pairs with unidentified nucleotides)
#maxEE=c(2, 3) - after truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#truncQ=2 - Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.
#rm.phix=TRUE - Default TRUE. If TRUE, discard reads that match against the phiX (bacteriophage) genome, as determined by isPhiX.


#Learn the error rates 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


#Visualizе of the estimated error rates

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#Correction of reads taking into account error rates

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

                                                                                                                                                                                                                                                                                                                                                                                          sequence
#Construct sequence table

seqtab <- makeSequenceTable(mergers)


#Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


#Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab.nochim)))


#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


#make a table from a data frame. sep="\t" - tab separator

write.table(seqtab.nochim, file='D:\\ivan\\eu_3south\\eu_3south_seqtab.nochim.txt', sep="\t")
write.table(track, file='D:\\ivan\\eu_3south\\track_eu_3south.txt', sep="\t")


#Assign taxonomy by SILVA

taxa_silva <- assignTaxonomy(seqtab.nochim, "D:\\ivan\\eu_3south\\silva_132.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)


#make a table from a data frame

write.table(taxa_silva, file='D:\\ivan\\eu_3south\\taxa_eu_3south_silva.txt', sep="\t")


#Assign taxonomy by PR2

taxa_pr2 <- assignTaxonomy(seqtab.nochim, taxLevels = c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"), "D:\\ivan\\eu_3south\\pr2_version_5.0.0_SSU_dada2.fasta.gz", multithread=TRUE)


#make a table from a data frame

write.table(taxa_pr2, file='D:\\ivan\\eu_3south\\taxa_eu_3south_pr2.txt', sep="\t")
