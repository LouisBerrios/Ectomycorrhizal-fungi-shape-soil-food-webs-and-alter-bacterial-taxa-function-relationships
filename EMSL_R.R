### EMSL 16S Dataset ###

# unload file folder folders to deposit FASTQ files

find ~/Documents/EMSL_Project/FASTQ_Files/ -type f -name '*L001*' -exec mv {} ~/Documents/EMSL_R/ \;

# DADA2 workflow for processing 16S raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(ShortRead)
library(Biostrings)


#############################

path <- "~/Documents/EMSL_R"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_L001_R1_001.fastq.gz","",basename(fnFs))
basefilenames_Rs <- sub("_L001_R2_001.fastq.gz","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_L001_R1_001.fastq.gz"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_L001_R2_001.fastq.gz"))]
}


# Identify primers - used 515F & 806R from Hiro's spreadsheet
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer seq
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...


# Get all orientations of primers, just to be safe
# Note - changed this for the dimensions project due to different sequencing primers
# Used

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 1, multithread = TRUE)

# (From tutorial) We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations. 
# Identifying and counting the primers on one set of paired end FASTQ files is
# sufficient, assuming all the files were created using the same library preparation,
# so we’ll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# If you see the reverse-complement of the forward primer in the reverse reads (cells [2,4] and [3,4]),
# it's because the ITS region is short and it is reading part of the forward primer.

# Remove primers using cutadapt

cutadapt <- "/Users/louisberrios/Documents/Cutadapt/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  # for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files; fnFs.filtN replaced by fnFs.filtN, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Count primers in first post-cutadapt sample (should all be 0):
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[20]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[20]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[20]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[20]]))

# Since they are zero, skip step to remove other orientations of primers

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), split="_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles of forward reads #1-2
plotQualityProfile(cutFs[1:12])

# Inspect read quality profiles of reverse reads #1-2
plotQualityProfile(cutRs[1:12])

# Filter and trim

# Assigning the filenames for the output of the filtered reads 
# to be stored as fastq.gz files.
filtFs <- file.path(path, "filtered", basename(fnFs.filtN))
filtRs <- file.path(path, "filtered", basename(fnRs.filtN))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

filtFs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R1_001.fastq.gz", full.names=TRUE)
filtRs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R2_001.fastq.gz", full.names=TRUE)

# Learn the error rates
errF <- learnErrors(filtFs.out, multithread = TRUE)
errR <- learnErrors(filtRs.out, multithread = TRUE)

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs.out, verbose = TRUE)
derepRs <- derepFastq(filtRs.out, verbose = TRUE)
# Name the derep-class objects by the sample names
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(filtFs.out, get.sample.name))

# DADA2's core sample inference algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,trimOverhang = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

saveRDS(seqtab,"EMSL_16S_seqtab.bac.rds")
saveRDS(seqtab.nochim,"EMSL_16S_seqtab.bac.nochim.rds")

# Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))

#format out to accommodate dropped samples
raw.sample.names <- unname(sapply(row.names(out), get.sample.name))

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

track2<-cbind(out,track[match(row.names(out),row.names(track)),])

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                      "nonchim")
write.csv(track2,"EMSL_16S_summary.csv")

rownames(track) <- sample.names
head(track2)


# Assign taxonomy using the UNITE database
silva.ref<-"silva_nr99_v138.1_train_set.fa.gz"
silva.species<-"silva_species_assignment_v138.1.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = TRUE, tryRC = TRUE)

taxa <- addSpecies(taxa, silva.species)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "EMSL_16S_seqtab_nochim_taxa.rds")
saveRDS(taxa, "EMSL_16S_taxa.rds")

########Construct phyloseq object with sample data variables#########

# make sample table
library(phyloseq)
sample <- as.data.frame(rownames(seqtab.nochim))
library(dplyr)
library(tibble) 
sample2 <- sample %>%
  mutate(Condition = ifelse(grepl("SP", Sample), "EcM", "No-EcM"))
sample3 <- sample2 %>% mutate(Library = ifelse(grepl("DNA", Sample), "DNA", "RNA"))
sample4 <- column_to_rownames(sample3, var = "Sample")
SAM <- sample_data(sample4, errorIfNULL = T)

#read in environmental sample table (this is used when building sample data outside of R)

sample <- readRDS("EMSL_16S_seqtab_nochim_taxa.rds")
write.csv(SAM, "EMSL_SampleData.csv") # removed controls
SAM <- read.csv("EMSL_SampleData.csv", row.names = 1)
SAM <- sample_data(SAM, errorIfNULL = T)
taxa <- readRDS("EMSL_16S_taxa.rds")

#create a phyloseq object
bac.ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(SAM), tax_table(taxa))

#filter out unwanted taxa (e.g., mitochondira and chloroplast sequences)
bac.ps.filt<-subset_taxa(bac.ps,Family!="Mitochondria")
bac.ps.filt<-subset_taxa(bac.ps.filt,Genus!="Chloroplast")

# Removing sequence rownames for display only
taxa.print <- tax_table(bac.ps.filt)
rownames(taxa.print) <- NULL
head(taxa.print)

#save the filtered dataset 
saveRDS(bac.ps.filt,"EMSL.16S.filtered.NoCon.rds")

#filter out low abundant sequences
bac.ps.filt2 = prune_taxa(taxa_sums(bac.ps.filt) > 10, bac.ps.filt) 
bac.ps.filt2 = prune_samples(sample_sums(bac.ps.filt2)>1000, bac.ps.filt2)

#save the filtered+pruned dataset
saveRDS(bac.ps.filt2, "EMSL.16S.filt-prune.NoCon.rds")

#rarefy the dataset
bac.ps.rare <- rarefy_even_depth(bac.ps.filt2, sample.size = min(sample_sums(bac.ps.filt2)),
                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
# remove controls
rare2 <- subset_samples(bac.ps.rare, Condition != "NA") # use this file for the DESeq2 analysis
raredf <- psmelt(rare2)
#save rarefied dataset
saveRDS(bac.ps.rare, "Rarefied-EMSL-Bacteria.NoCon.rds")
rare <- readRDS("Rarefied-EMSL-Bacteria.NoCon.rds")
#create ASV-TAX table
taxa_names(bac.ps.rare) <- paste0("Seq", seq(ntaxa(bac.ps.rare)))
ASV_BAC <- data.frame(otu_table(bac.ps.rare))
TAX_ASV <- data.frame(tax_table(bac.ps.rare))
ASV_BAC_T <- t(ASV_BAC)
MERGER2 <- merge(ASV_BAC_T, TAX_ASV, by = "row.names")
write.csv(MERGER2, "ASV_Bacteria_EMSL_Rare_Final.csv")

#create the combined OTU table plus taxonomy string metacoder needs

#read it in using phyloseq object
tax_table(bac.ps.rare)<-tax_table(bac.ps.filt2)[,1:7]
x1<-parse_phyloseq(bac.ps.rare)

#transform data [compositional]
pseq <- microbiome::transform(bac.ps.rare, "compositional")

#save as RDS
saveRDS(pseq, "EMSL-Bacteria-PS_filter+rare+tran.NoCon.rds")

# alpha diversity
# Remove NA first
rare2 <- subset_samples(rare, Condition != "NA")
bacaDIV <- estimate_richness(rare2)
write.csv(bacaDIV, "EMSL_Richness.AlphaDiv.csv")
# plot using phyloseq
plot_richness(rare2, x = "Condition", color = "Condition", measures = c("Observed")) + theme_bw() + xlab("") + labs(title = "Bacterial Community") + geom_boxplot() 
# plot using ggplot
library(tibble)
bacaDIV2 <- rownames_to_column(bacaDIV)
colnames(bacaDIV2)[1] <- "Sample"
library(dplyr)
bacaDIV3 <- bacaDIV2 %>% 
mutate(Condition = ifelse(grepl("SP", Sample), "EcM", "No-EcM"))
# add facet variable
bacaDIV3$title <- "Bacterial Richness"
# stats
anova_alpha2 <- aov(Observed ~ Condition, bacaDIV3)
summary(anova_alpha2)
# plot w/ stats
ggplot(bacaDIV3, aes(x=Condition, y=Observed, color=Condition)) + 
  geom_boxplot() + theme_bw() + scale_color_manual(values=c("black", "red")) + 
  xlab("") + theme(axis.text.x = element_text(size = 12, face = "bold")) + theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + facet_grid(~title) + theme(axis.title.x = element_text(size = 16, face="bold")) +
  theme(strip.text = element_text(size = 14, face = "bold", color="white")) + theme(axis.title.y = element_text(size = 16, vjust = 1.5, face="bold")) + 
  annotate("text", x = 1.5, y = 250, label = "Condition:F = 5.577; p = 0.026", size = 4, fontface="bold") +
  theme(strip.background = element_rect(fill = "burlywood4")) + scale_x_discrete(labels=c("EcMF", "No EcMF"))

# Alternative plot with adjusted colors --> colors moving forward with be from the wesanderson package (Chevalier1)
ggplot(bacaDIV3, aes(x=Condition, y=Observed, color=Condition)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + scale_color_manual(values=c("#446455","#C7B19C")) + 
  xlab("") + theme(axis.text.x = element_text(size = 12, face = "bold")) + theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + facet_grid(~title) + theme(axis.title.x = element_text(size = 16, face="bold")) +
  theme(strip.text = element_text(size = 14, face = "bold", color="white")) + theme(axis.title.y = element_text(size = 16, vjust = 1.5, face="bold")) + 
  annotate("text", x = 1.5, y = 250, label = "Condition:F = 5.577; p = 0.026", size = 4, fontface="bold") +
  theme(strip.background = element_rect(fill = "black")) + scale_x_discrete(labels=c("EcMF", "No EcMF")) + 
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), shape=1, size = 3, aes(color=factor(Condition)), show.legend = F)

# save plot 1
ggsave(
  filename = "EMSL_16S_Richness.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5.5,
  units = c("in"),
  dpi = 300)

# save plot 2
ggsave(
  filename = "EMSL_16S_Richness_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5.5,
  units = c("in"),
  dpi = 300)

#create dataframe 
pseqDF <- psmelt(pseq2)
# subset out the controls
df_subset <- pseqDF[grepl("NM|SP", pseqDF$Sample), ]
library(RColorBrewer)
library(ggplot2)
#get appropriate number of colors for data
getPalette <- colorRampPalette(brewer.pal(13, "Set1"))
colourCount = length(unique(df_subset$Phylum))

#plot the relative abundance by Phylum
ggplot(df_subset, aes(fill=Phylum, y=Abundance, x=Condition)) + geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = c("red3", "black", "goldenrod3", "darkslategray4", "burlywood4", "cornsilk4", "lightgoldenrod3", "lightpink3", "royalblue3", "coral3", "antiquewhite3", "darkolivegreen", "cyan4")) + ylab("Relative Abundance") + scale_y_continuous(labels = scales::percent) + xlab("") + theme(axis.text.x = element_text(size=12, face="bold")) + theme(axis.text.y = element_text(face = "bold", size = 12)) + theme(axis.title.y = element_text(size = 16, face="bold"))
# save
ggsave(
  filename = "EMSL_16S_by_ConditionBarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5.5,
  units = c("in"),
  dpi = 300)
# create dataframe by Phylum
Agg.phy <- aggregate(Abundance ~ Phylum + Condition, data=df_subset, FUN=sum)
# create aggregated by genus dataframe
AGG <- aggregate(Abundance ~ Genus + Condition, data = raredf, FUN=sum)
# subset top 20 most abundant genera for each condition (i.e., EcM + No EcM)
EcM <- subset(AGG, Condition == "EcM")
noEcM <- subset(AGG, Condition == "No-EcM")
top20EcM <- subset(EcM, Abundance > 450) # top 20 
topnoECM <- subset(noEcM, Abundance > 480) # top 20
# merge dataframes 
top20 <- merge(top20EcM, topnoECM, all=TRUE)
# save dataframe 
write.csv(top20, "Merged_Top20_Genus-Condition.csv")
# plot barplot
ggplot(top20, aes(x=Genus, y=Abundance)) + xlab("") + ylab("No. of Sequences") + geom_bar(stat="identity", fill = c("black")) + coord_flip() + theme_bw() + theme(strip.text = element_text(size = 8, face = "bold")) + theme(axis.text.x = element_text(size = 6, face = "bold")) + theme(axis.text.y = element_text(size = 6, face = "bold.italic")) + theme(axis.title.x = element_text(size = 7, face = "bold", vjust = 0.5)) + facet_grid(~Condition)
ggsave(
  filename = "EMSL_16S_Top20_by_ConditionBarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5.5,
  units = c("in"),
  dpi = 300)
####### NMDS ordination ###########
# subset out controls
pseq2 <- subset_samples(pseq, Condition != "NA")
saveRDS(pseq2, "EMSL_Bac_Phylo_NoCon.rds")
pseq2 <- readRDS("EMSL_Bac_Phylo_NoCon.rds")
#perform an NMDS ordination
ord.bac <- ordinate(pseq2, "NMDS", "bray", k=2)

#plot NMDS ordination [axes 1:2] w/ PERMANOVA results from below.
plot_ordination(pseq2, ord.bac, type="sample", color="Condition", axes=1:2) + theme_bw() + scale_color_manual(values = c("#446455","#C7B19C")) + geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) + annotate("text", x = -1.91, y = 1.7, label = "DNA/RNA:F = 0.889; r2 = 0.028; p = 0.530", size = 4, fontface="bold") + annotate("text", x = -2, y = 1.8, label = "Condition:F = 5.56; r2 = 0.18; p = 0.001", size = 4, fontface="bold")
# add new column for faceting 
sample_data(pseq2)$Community <- rep("Bacterial Community", nrow(sample_data(pseq2)))
plot_ordination(pseq2, ord.bac, type="sample", color="Condition", axes=1:2) + 
  theme_bw() + scale_color_manual(values = c("#446455","#C7B19C"), labels = c("EcMF", "No EcMF")) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) + theme(axis.title = element_text(size =16, face = "bold")) +
  annotate("text", x = -0.9, y = 1.7, label = "DNA/RNA:F = 0.889; r2 = 0.028; p = 0.530", size = 3, fontface="bold") + 
  annotate("text", x = -1, y = 1.8, label = "Condition:F = 5.56; r2 = 0.18; p = 0.001", size = 3, fontface="bold") + 
  facet_grid(~Community) + theme(strip.text = element_text(face = "bold", size=12, color="white")) + theme(strip.background = element_rect(fill = "black"))
########PERMANOVA###########

#calculate bray-curtis distance matrix
library(vegan)
ord.bac.bray <- phyloseq::distance(pseq2, method = "bray", k=2)


#make a dataframe from the sample_data
ord.bac.bray.sampleDF <- data.frame(sample_data(pseq2))

#adonis test
PERMANOVA <- adonis2(ord.bac.bray ~ Condition, data = ord.bac.bray.sampleDF)
PERMANOVA2 <- adonis2(ord.bac.bray ~ Condition + Library, data = ord.bac.bray.sampleDF)
#write PERMANOVA results to table CSV
list <- PERMANOVA[["aov.tab"]]
list2 <- PERMANOVA2[["aov.tab"]]
write.csv(list, "EMSL_Bac_Permanova.csv")
write.csv(list2, "EMSL_Bac_Permanova2.csv")

# save plot 1
ggsave(
  filename = "EMSL_16S_BrayNMDS_PERMANOVA.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

# save plot2 
ggsave(
  filename = "EMSL_16S_BrayNMDS_PERMANOVA_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

# betadisper
library(vegan)
ord.bac.bray <- phyloseq::distance(pseq2, method = "bray", k=3)
beta_adonis <- betadisper(ord.bac.bray, ord.bac.bray.sampleDF$Library, bias.adjust = TRUE)
beta_adonis2 <- betadisper(ord.bac.bray, ord.bac.bray.sampleDF$Condition, bias.adjust = TRUE)
#stats
stat_disp_anova <- anova(beta_adonis)
stat_disp_anova2 <- anova(beta_adonis2)
stat_disp_tukey <- TukeyHSD(stat_disp_anova)
stat_disp_tukey2 <- TukeyHSD(stat_disp_anova2)

plot(beta_adonis, hull=FALSE, ellipse=TRUE)
plot(beta_adonis2, hull=FALSE, ellipse=TRUE)

#pairwise permutation test for homogeneity of multivariate dispersions
perm.test <- permutest(beta_adonis, pairwise = TRUE, permutations = 1000)
perm.test2 <- permutest(beta_adonis2, pairwise = TRUE, permutations = 1000)


######### DESeq2 Analysis ###########

library(DESeq2)

#convert phyloseq to DESeq2 format
dds <- phyloseq_to_deseq2(rare2, ~ Condition)

#calculate size factors using edgeR
library(edgeR)
sizeFactors(dds) <- calcNormFactors(counts(dds))

#run DESeq function
dds = DESeq(dds, test="Wald", fitType="parametric")

res <- results(dds, cooksCutoff = FALSE)
alpha = 0.1
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(rare2)[rownames(sigtab), ], "matrix"))

#prepare data to plot DESeq2 output
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
colourCount = length(unique(sigtab$Phylum))

library(ggplot2)
#plot DESeq2 output w/ ggplot2 [EcM vs. No EcM]
sigtab$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium", sigtab$Genus)
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4, alpha = 0.5) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme_bw() + coord_flip() + scale_color_manual(values = getPalette(colourCount)) + xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 7)) + labs(title = "Differentially Abundant Bacterial Genera [EcMF vs. No EcMF]") + theme(plot.title = element_text(size = 7)) + geom_hline(yintercept=0, col="red", linetype = 2)
# modified 1
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4, alpha = 0.7) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme_bw() + coord_flip() + theme(axis.title = element_text(size=12, face="bold")) +
  scale_color_manual(values = c("black", "darkslategray4", "red3", "darkgoldenrod", "burlywood4", "orchid", "dodgerblue3", "ivory4", "lightblue")) + 
  xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 10)) + theme(plot.title = element_text(size = 7)) + 
  geom_hline(yintercept=0, col="red", linetype = 2) + theme(axis.text.y = element_text(face="bold.italic")) + theme(axis.text.x = element_text(face="bold"))

#modified (don't use)
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5, alpha = 0.5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme_bw() + coord_flip() + 
  scale_color_manual(values = getPalette(colourCount)) + xlab("") + ylab("[Log2] Fold-Change") + 
  theme(axis.text.y = element_text(size = 7, face="bold")) + labs(title = "Differentially Abundant Bacterial Genera [Bulk soil vs. Endosphere]") + 
  theme(plot.title = element_text(size = 10))  + theme(axis.text.y = element_text(size = 7, face="bold")) + 
  theme(axis.text.x = element_text(size = 8, face = "bold")) + scale_x_discrete(labels = c("Aquaspirillum", "Acidibacter", "Acidicapsa", "Acidipilia", "Acidothermus", "Actinocatenispora", "Actinospica", "[Unclassified] Verrucomicrobia", "Afipia", "Aliidongia", "Rhizobium", "Aminobacter", "Amycolatopsis", "Anaeromyxobacter", "Aquisphaera", "Asanoa", "Bradyrhizobium", "Bryobacter", "Burkholderia", "Candidatus Koribacter", "Candidatus Solibacter", "Candidatus Udaeobacter", "Catenulispora", "Conexibacter", "Dongia", "Edaphobacter", "Edaphobaculum", "Frankia", "Gaiella", "Gemmata", "Gemmatimonas", "Haliangium", "Inquilinus", "[Unclassified] IS-44", "Jatrophihabitans", "Labrys", "Massila", "Mesorhizobium", "Pseudomonas sp. MND1", "Mucilaginibacter", "Mycobacterium", "Novosphingobium", "Occallatibacter", "Pseudonocardia", "Puia", "Reyranella", "Rhodanobacter", "Rhodoplanes", "Roseiarcus", "Rubrobacter", "Sphaerotilus", "Sphingomonas", "Sporocytophaga", "Streptomyces", "Terriglobus", "Tundrisphaera")) + geom_hline(yintercept=c(-2,2), linetype="dashed", color=c("blue", "red"))

# save plot 
ggsave(
  filename = "EMSL_16S_DESEq2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6.5,
  height = 5,
  units = c("in"),
  dpi = 300)

# modified plot
ggsave(
  filename = "EMSL_16S_DESEq2_Final.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 8.5,
  units = c("in"),
  dpi = 300)

### Metabolite Data ###

# manually formatted in Excel
met <- read.csv("~/Documents/EMSL_R/Metabolomics/EMSL_Met_DataFrame.csv")
library(dplyr)
met_known <- met %>%
  filter(!grepl("Unknown", Metabolite)) # filter out all 'Unknown' metabolites
library(wesanderson)
my_colors <- wes_palette("Darjeeling1", 3)
my_colors2 <- c("red", "black")
my_colors3 <- c("darkslategray2", "black")
my_colors4 <- c("white", "dodgerblue4")
# add Z-score data
met_known$zscore <- ave(met_known$Abundance, met_known$Metabolite, met_known$Condition, FUN=scale)
#plot
library(ggplot2)
ggplot(met_known, aes(x=Condition, y=Metabolite, fill=zscore)) + geom_tile(color = "black") + theme_bw() + xlab("") + ylab("") + scale_fill_gradientn(colors = my_colors2, na.value = "white")
# save plot 
ggsave(
  filename = "EMSL_Met_HeatMap.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6.5,
  units = c("in"),
  dpi = 300)

# add classification (loose), subset, and plot
met_known$Class <- ifelse(grepl("ose|lse", met_known$Metabolite), "Sugar",
                   ifelse(grepl("^L-", met_known$Metabolite), "Amino Acid", "Other"))

amino <- subset(met_known, Class == "Amino Acid")
sugar <- subset(met_known, Class == "Sugar")
other <- subset(met_known, Class == "Other")
# aminos
ggplot(amino, aes(x=Sample, y=Metabolite, fill=zscore)) + geom_tile(color = "white") + theme_bw() + xlab("") + ylab("") + scale_fill_gradientn(colors = my_colors3, na.value = "white") + ggtitle("Amino Acids") + theme(axis.text = element_text(size=10, face="bold")) + theme(axis.text.x = element_text(angle=45, hjust=1))
                                                                                                                                                                                                                                                                                                      
ggsave(
  filename = "EMSL_Aminos_HeatMap2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 4,
  units = c("in"),
  dpi = 300)
# sugars
ggplot(sugar, aes(x=Sample, y=Metabolite, fill=zscore)) + geom_tile(color = "white") + theme_bw() + xlab("") + ylab("") + scale_fill_gradientn(colors = my_colors4, na.value = "white") + ggtitle("Sugars") + theme(axis.text = element_text(size=12, face="bold")) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
ggsave(
  filename = "EMSL_Sugars_HeatMap2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 4.5,
  units = c("in"),
  dpi = 300)
# others
ggplot(other, aes(x=Condition, y=Metabolite, fill=zscore)) + geom_tile(color = "black") + theme_bw() + xlab("") + ylab("") + scale_fill_gradientn(colors = my_colors2, na.value = "white") + ggtitle("Other Metabolites") + theme(axis.text = element_text(size=12, face="bold"))
ggsave(
  filename = "EMSL_Others_HeatMap2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10.5,
  height = 8,
  units = c("in"),
  dpi = 300)


##### NMDS Ordinations #######

met.data <- read.csv("Met_Data.csv", row.names = 1, header = TRUE)

met.data2<-t(met.data)
sp.metab.norm<-decostand(met.data2,"log",margin=2) #1=rows, 2=cols
sp.metab.norm<-wisconsin(sp.metab.norm); rowSums(sp.metab.norm); colSums(sp.metab.norm)
sp.metab.norm[1:5,1:5]; colSums(sp.metab.norm)

metab.mds<-metaMDS(sp.metab.norm,k=2,trymax=50,autotransform = T)
plot(metab.mds, display="sites")

treatments<-rep("NM",length(row.names(sp.metab)))
treatments[grep("Sp.",row.names(sp.metab))]<-"SP"

mds.scores<-data.frame(metab.mds$points)
mds.scores$Sample <- row.names(mds.scores)
mds.scores$Treatment <- ifelse(grepl("Sp", mds.scores$Sample), "EcMF","No EcMF")

library(viridis)


ggplot(mds.scores,aes(x=MDS1,y=MDS2,color=Treatment)) + 
  geom_point(size=3) +
  theme_bw() + scale_color_manual(values=c("black", "red3"))
library(vegan)
adonis2(met.data~Treatment, data=mds.scores)

# Change Axis names and annotate plot

colnames(mds.scores)[1] <- "NMDS1"
colnames(mds.scores)[2] <- "NMDS2"
mds.scores$Title <- "Non-Lipid Metabolites"

ggplot(mds.scores,aes(x=NMDS1,y=NMDS2,color=Treatment)) + 
  geom_point(size=3) + theme_bw() + scale_color_manual(values=c("#446455","#C7B19C")) + 
  annotate("text", x = 0.08, y = 0.15, label = "F = 1.63; r2 = 0.03581; p = 0.132", size = 3, fontface="bold") +
  stat_ellipse(type = "norm", linetype=2) + facet_grid(~Title) +
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(size=12, face="bold", color="white"), axis.text = element_text(size=10, face="bold"), axis.title = element_text(size=12, face="bold"))

# Save plot 1
ggsave(
  filename = "EMSL_Metabolites_NMDS.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

# save plot 2
ggsave(
  filename = "EMSL_Metabolites_NMDS_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

################################

#### DESeq2 of Metabolites ######

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(met.data2), # round the count data since DESeq2 needs integers
                              colData = mds.scores,
                              design = ~ Treatment)

dds = DESeq(dds, test="Wald", fitType="parametric")


res <- results(dds, cooksCutoff = FALSE)
alpha = 0.1
sigtab.met <- res[which(res$padj < alpha), ]
sigtab.met2 <- cbind(as(sigtab.met, "data.frame")) # All negative fold changes are EcMF enriched. All positives are Rhizosphere-enriched.

write.csv(sigtab.met2, "Metabolite_DESeq2.table.csv") # non-labelled table 

sigtab.met2$Metabolite <- row.names(sigtab.met2)

library(dplyr)
sigtab.met3 <- sigtab.met2 %>% mutate(Treatment = case_when(log2FoldChange <= 0 ~ "EcMF", 
                                                            log2FoldChange >= 0 ~ "No EcMF"))

sigtab.met3$Title <- "Differentially Abundant Non-Lipid Metabolites"
write.csv(sigtab.met3, "Metabolite_DESeq2.table.Treatments.csv")
sigtab.met3 <- read.csv("Metabolite_DESeq2.table.Treatments.csv", row.names = 1)
# Plot DESeq2 Output 

library(ggplot2)
ggplot(sigtab.met3, aes(x=factor(Metabolite, level = c("Unknown 137", "Unknown 125", "Sucrose", "Sorbose", "Scyllo-Inositol", "Malonic Acid", "Maleic Acid", "L-Threonine", "L-Pyroglutamic Acid", "L-Phenylalanine", "L-Ornithine", "L-Glutamic Acid", "Galactitol", "D-Xylitol", "D-Ribose", "D-Mannitol", "D-Arabitol", "Acetol", "4-Aminobutyric Acid (GABA)")), y= log2FoldChange, color = Treatment)) +
                                              geom_point(size=6) + theme_bw() + coord_flip() + scale_color_manual(values = c("#446455","#C7B19C")) + xlab("") + facet_grid(~Title) + 
                                              theme(strip.text = element_text(color = "white", size = 9, face = "bold")) +
                                              theme(strip.background = element_rect(fill = "black")) + labs(colour = "Treatment") + 
                                              theme(axis.text = element_text(face = "bold")) + 
                                              ylab("[Log2] Fold-Change") + theme(axis.title = element_text(face = "bold")) + 
                                              theme(legend.title = element_text(face = "bold")) + theme(axis.text.y = element_text(size = 12, face = "bold")) + 
                                              geom_segment(aes(x = Metabolite, xend = Metabolite, y = 0, yend = log2FoldChange))
# Save plot 1 
ggsave(
  filename = "EMSL_Metabolites-DESeq2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5,
  units = c("in"),
  dpi = 300)

# Save plot 2
ggsave(
  filename = "EMSL_Metabolites-DESeq2_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5,
  units = c("in"),
  dpi = 300)

### Prep for PiCrust and MIMOSA ################################################################

asv_seqs <- colnames(rare)
asv_headers <- vector(dim(rare)[2], mode="character")
for (i in 1:dim(rare)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
asv_otu <- t(rare)
row.names(asv_otu) <- sub(">", "", asv_headers)
asv_tax <- tax_asv
row.names(asv_tax) <- sub(">", "", asv_headers)
OTU_TAX_table <- merge(asv_otu, asv_tax, by=0)
write(asv_fasta, "asv_fasta.fa")
write.table(asv_otu, "asv_otu.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax.csv", sep=",", quote=F, col.names=NA)
write.table(OTU_TAX_table, "OTU_TAX_table.csv", sep=",", quote=F, col.names=NA)

################################################################################################

otu_seqs <- getSequences(ps)
otu_ids <- rownames(otu_table(ps))
otu_refs <- refseq(ps)

# Create mapping between sequence names/ids/ASV ids and OTU IDs
id_map <- data.frame(id = rownames(otu_seqs), otu_id = otu_ids, refseq_id = otu_refs)
id_map <- merge(id_map, tax_table(ps), by = "otu_id")

# Format OTU sequences as a FASTA file with matching IDs
fasta_file <- "my_otu_seqs.fasta"
write.fasta(otu_seqs, names = id_map$refseq_id, file.out = fasta_file)


########### Picrust2 #########################################################################

library(metagMisc); library(seqinr)

Canteen_clean <- readRDS('Rarefied-EMSL-Bacteria.NoCon.rds')

sample_names(Canteen_clean) <- gsub(pattern = '-', replacement = '__', x = sample_names(Canteen_clean))
cc_table <- phyloseq_to_df(Canteen_clean, addtax = T, addtot = F, addmaxrank = F, sorting = 'abundance')
colnames(cc_table) <- gsub(pattern = '__', replacement = '-', x = colnames(cc_table))

counts_table <- cc_table[,-c(1:8)]

asv_headers <- vector(dim(cc_table)[1], mode = 'character')
for (i in 1:dim(cc_table)[1]) {
  asv_headers[i] <- paste('ASV', i, sep = '_')
}

row.names(counts_table) <- sub('>', '', asv_headers)

write.table(counts_table, 'ASVs_counts.tsv', sep = '\t', quote = F, col.names = NA)

asv_seqs <- cc_table[,1]

write.fasta(as.list(asv_seqs), asv_headers, 'ASVs.fa', open = 'w', nbchar = 60, as.string = FALSE)

############################################################################################################


### Bacterial Characteristics Analysis ###

#Gram stain + motility characteristics obtained from Google searches

# load data
bac.char <- read.csv("DESeq2_BacterialCharacteristics.csv", header = T)
# only EcMF
bac.char.ecmf <- subset(bac.char, Enrichment.Type == "EcMF")
sum(bac.char.ecmf$Ratio)
# only No EcMF
bac.char.no.ecmf <- subset(bac.char, Enrichment.Type == "No EcMF")
sum(bac.char.no.ecmf$Ratio)
# Gram stain
ggplot(bac.char, aes(x=Enrichment.Type,y=Ratio,fill=Enrichment.Type)) + geom_col() + 
  theme_bw() + facet_grid(~Cell.Wall) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold")) + scale_y_continuous(labels = scales::percent)

library(ggpubr) 
my_comps <- list(c("EcMF", "No EcMF"))
#plot motility
ggplot(bac.char, aes(x=Enrichment.Type,y=Ratio,fill=Enrichment.Type)) + geom_col() + 
  theme_bw() + facet_grid(~Motility) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + theme(strip.background = element_rect(fill = "black")) +
  theme(strip.text = element_text(size=12, face="bold", color = "white")) + scale_y_continuous(labels = scales::percent) + stat_compare_means(comparisons = my_comps, label.y =0.82)

# save plot
ggsave(
  filename = "EMSL_Bacterial.Characteristics.Motility.DESeq2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

# plot cell wall type
ggplot(bac.char, aes(x=Enrichment.Type,y=Ratio,fill=Enrichment.Type)) + geom_col() + 
  theme_bw() + facet_grid(~Cell.Wall) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + theme(strip.background = element_rect(fill = "black")) +
  theme(strip.text = element_text(size=12, face="bold", color = "white")) + scale_y_continuous(labels = scales::percent) + stat_compare_means(comparisons = my_comps, label.y =0.9)

# save
ggsave(
  filename = "EMSL_Bacterial.Characteristics.CellWall.DESeq2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

