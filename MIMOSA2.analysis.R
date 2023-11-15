##### MIMOSA2 Analysis #####

# The files used to run this analysis are:
#       MIMOSA2_DNA-RNA_table.tsv
#       MIMOSA2_metabolite_table.tsv

# The pre/post processing of these data in R (produced by Andressa Venturini):

### MIMOSA2 ###
library(sqldf); library(tidyverse)

getwd(); setwd('/Users/andressa/Documents/PD_Stanford/Parcerias/Louis/EMSL/PICRUSt2_strat/picrust2_out_pipeline/KO_metagenome_out/')
data_pic <- read.delim('pred_metagenome_contrib.tsv', header = T)

unique(data_pic$sample)

head(data_pic)
names(data_pic)[names(data_pic) == 'sample'] <- 'Sample'
names(data_pic)[names(data_pic) == 'function.'] <- 'Gene'
names(data_pic)[names(data_pic) == 'taxon'] <- 'OTU'
names(data_pic)[names(data_pic) == 'genome_function_count'] <- 'GeneCountPerGenome'
names(data_pic)[names(data_pic) == 'taxon_abun'] <- 'OTUAbundanceInSample'
names(data_pic)[names(data_pic) == 'taxon_function_abun'] <- 'CountContributedByOTU'

data_pic <- data_pic[,-c(5,8,9)]
unique(data_pic$Sample)

data_pic_noM <- data_pic %>% filter(!grepl('Bac-Mock_S37_L001_R1_001.fastq.gz', Sample))
data_pic_noM <- data_pic_noM %>% filter(!grepl('NM21-RNA_S23_L001_R1_001.fastq.gz', Sample))
data_pic_noM <- data_pic_noM %>% filter(!grepl('SP35-RNA_S13_L001_R1_001.fastq.gz', Sample))
unique(data_pic_noM$Sample)

data_pic_noM <- data.frame(lapply(data_pic_noM, function(x) {gsub('-.*', '', x)}))
data_pic_noM <- data.frame(lapply(data_pic_noM, function(x) {gsub('NM', 'NM.', x)}))
data_pic_noM <- data.frame(lapply(data_pic_noM, function(x) {gsub('SP', 'Sp.', x)}))
unique(data_pic_noM$Sample)

getwd(); setwd('/Users/andressa/Documents/PD_Stanford/Parcerias/Louis/EMSL/Mimosa/')
data_met1 <- read.csv('EMSL_Met_DataFrame.csv', header = TRUE, sep = ',')
data_met2 <- read.csv('Met_Data.csv', header = TRUE, sep = ',')
data_met3 <- read.csv('name_map.csv', header = TRUE, sep = ',')

names(data_met3)[names(data_met3) == 'Query'] <- 'Metabolite'
data_met <- merge(data_met2, data_met3, by = 'Metabolite')
data_met <- data_met[,c(51,2:47)]

unique(data_pic_noM$Sample)
colnames(data_met)
data_met <- select(data_met, KEGG, NM.20, NM.21, NM.26, NM.35, NM.38, NM.39, NM.42, NM.44, Sp.01, Sp.07, Sp.08, Sp.11, Sp.12,
                   Sp.16, Sp.20, Sp.28, Sp.35, Sp.37, Sp.38, Sp.43, Sp.44)

write.table(data_pic_noM, 'DNA+RNA/MIMOSA2_DNA-RNA_table.tsv', sep = '\t', quote = F, col.names = NA)
write.table(data_met, 'DNA+RNA/MIMOSA2_metabolite_table.tsv', sep = '\t', quote = F, col.names = NA)

### Note that there are 21 samples total; there are 28 bacterial amplicon samples;
# And two of the samples NM21-RNA and SP35-RNA were removed since there was a DNA
# equivalent sample AND there were no significant differences between RNA and DNA communities

# load the communityMetabolicPotentialScores.csv file (from the allResults_log folder)
# renamed to CMP.txt <- convert to .csv file

CMP <- read.csv("CMP.csv", header = T)

#create ASV-TAX table
taxa_names(rare2) <- paste0("ASV_", seq(ntaxa(rare2)))
ASV_BAC <- data.frame(otu_table(rare2))
TAX_ASV <- data.frame(tax_table(rare2))
ASV_BAC_T <- t(ASV_BAC)
MERGER2 <- merge(ASV_BAC_T, TAX_ASV, by = "row.names")
write.csv(MERGER2, "ASV_EMSL_MIMOSA2.csv")

# load the ASV table to match taxa to ASV numbers and CMPs
MERGER2 <- read.csv("ASV_EMSL_MIMOSA2.csv", header = T)
# remove seq abundances from the table
asv_tax_only <- MERGER2[,c(1,30:36)] 
# save data frame
write.csv(asv_tax_only, "EMSL_ASV_TAX_ONLY.csv")
# match the CMP data frame to the asv_tax_only data frame
# first change the column name of the asv_tax_only data frame from Row.names to Species
colnames(asv_tax_only)[1] <- "Species"
library(dplyr)

# Add new column 'Condition' to the CMP data frame
library(dplyr)

# Add new column based on sample names (i.e., SP or NM)
new_CMP <- CMP %>%
  mutate(Condition = case_when(
    grepl("^Sp", Sample) ~ "EcMF",  # For each Sp sample, label as "EcMF"
    grepl("^NM", Sample) ~ "No EcMF",  # For each Nm sample, label as "No EcMF"
    TRUE ~ NA_character_  # For other cases, you can set a default value (e.g., NA) or remove this line
  ))
CMP.agg <- aggregate(CMP ~ Condition + Species + compound, FUN = mean, data = new_CMP)

# match the compound names to their KEGG codes

### read in the compound name + codes ###

comp.names <- read.csv("name_map.csv", header = T)
# subset to only include needed information
comps.subs <- comp.names[,c(1,2,5)]
unique(comps.subs$KEGG) # there are 61 compounds
unique(CMP.agg$compound) # there are 16 compounds

# Add code - compound matches 

compound_mapping <- c(
  "C00025" = "L-Glutamic acid",
  "C00037" = "Glycine",
  "C00042" = "Succinic acid",
  "C00122" = "Fumaric acid",
  "C00183" = "L-Valine",
  "C00184" = "Dihydroxyacetone",
  "C00186" = "L-Lactic acid",
  "C00188" = "L-Threonine",
  "C00249" = "Palmitic acid",
  "C00258" = "Glyceric acid",
  "C01384" = "Maleic acid", 
  "C01530" = "Stearic acid",
  "C01606" = "Phthalic acid",
  "C01879" = "Pyroglutamic acid",
  "C06104" = "Adipic acid",
  "C06423" = "Caprylic acid"
)
# load dplyr package
library(dplyr)
# This is for the aggregated dataset 
CMP.agg.comps <- CMP.agg %>%
  mutate(Compound_name = ifelse(compound %in% names(compound_mapping), compound_mapping[compound], "Compound code not found"))
colnames(CMP.agg.comps)[2] <- "ASV"
CMP.final <- merge(CMP.agg.comps, TAX_ASV, by="ASV")

# save this data frame 
write.csv(CMP.agg.comps, "EMSL_MIMOSA_AggByCondition_withCompoundNames.csv")
write.csv(CMP.final, "EMSL_MIMOSA_SuppTable_MetsASV_TAXA.csv")
# plot the results
library(ggplot2)

ggplot(CMP.agg.comps, aes(x=Compound_name, y=CMP)) + facet_grid(~Condition) +
  geom_col() + theme_bw()

# add the compound names to the MIMOSA CMP data that has not been aggregated 
CMP.comps <- new_CMP %>%
  mutate(Compound_name = ifelse(compound %in% names(compound_mapping), compound_mapping[compound], "Compound code not found"))

ggplot(CMP.comps, aes(x=Sample, y=CMP, color=Compound_name)) + facet_grid(~Condition) +
  geom_col() + theme_bw()

# try analyzing each metabolite (n=16) separately (T test significant comparisons will be listed)
glut.acid <- subset(CMP.comps, Compound_name == "L-Glutamic acid") # L-Glutamic acid 
ggplot(glut.acid, aes(x=Condition, y=log(CMP))) + facet_grid(~Compound_name) +
  geom_boxplot() + theme_bw()
t.test(log(CMP) ~ Condition, data = glut.acid) # not sig. results 
# add stats the the plot with ggpubr
library(ggpubr)
my_comps = c("EcMF", "No EcMF")
ggplot(glut.acid, aes(x=Condition, y=log(CMP))) + facet_grid(~Compound_name) +
  geom_boxplot() + theme_bw() + stat_compare_means(comparisons = my_comps)
# L-glutmaic acid does not seem to differ substantially between conditions
# try glycine
glycine <- subset(CMP.comps, Compound_name == "Glycine")
ggplot(glycine, aes(x=Condition, y=log(CMP))) + facet_grid(~Compound_name) +
  geom_boxplot() + theme_bw() + stat_compare_means(comparisons = my_comps)
t.test(log(CMP) ~ Condition, data = glycine) # not sig. resutls 
# try succinic acid
succ.acid <- subset(CMP.comps, Compound_name == "Succinic acid")
# plot w/ wesanderson colors #446455" "#FDD262" "#D3DDDC" "#C7B19C"

P <- ggplot(succ.acid, aes(x=Condition, y=log(CMP), fill=Condition)) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("#446455", "#C7B19C"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + theme(axis.text = element_text(size = 12, face="bold"),
                                                                                            axis.title.x = element_blank(),
                                                                                            strip.text = element_text(face = "bold", size=16, color="white"),
                                                                                            strip.background = element_rect(fill = "black"),
                                                                                            axis.title.y = element_text(size = 12, face="bold"))
P + annotate("text", x=1.5, y=5, label="p=0.0189", size=6)

t.test(log(CMP) ~ Condition, data = succ.acid) # p = 0.01894

ggsave(
  filename = "EMSL_MIMOSA_Succinic-Acid.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

####Welch Two Sample t-test

#### data:  log(CMP) by Condition
#### t = -2.3493, df = 1510.1, p-value = 0.01894
#### alternative hypothesis: true difference in means between group EcMF and group No EcMF is not equal to 0
#### 95 percent confidence interval:
####  -0.28470267 -0.02561088
#### sample estimates:
####  mean in group EcMF mean in group No EcMF 
#### 2.162748              2.317905 

# Aggregate Succinic acid CMPs by ASVs and match ASV to taxa classes
succ.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = succ.acid)
colnames(succ.agg)[1] <- "ASV"
TAX_ASV$ASV <- row.names(TAX_ASV)
succ.tax <- merge(succ.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition
succ.tax.em <- subset(succ.tax, Condition == "EcMF")
succ.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = succ.tax.em)
succ.tax.em.agg.top20 <- subset(succ.tax.em.agg, CMP > 620)
succ.tax.em.agg.top20$Condition <- "EcMF"
succ.tax.noem <- subset(succ.tax, Condition == "No EcMF")
succ.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = succ.tax.noem)
succ.tax.noem.agg.top20 <- subset(succ.tax.noem.agg, CMP > 300)
succ.tax.noem.agg.top20$Condition <- "No EcMF"

# aggregate these data frames 
succ.agg.top20 <- merge(succ.tax.em.agg.top20, succ.tax.noem.agg.top20, all=T)
succ.agg.top20$Metabolite <- "Succinic acid"
succ.agg.top20$Genus[succ.agg.top20$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
succ.agg.top20$Genus[succ.agg.top20$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"

write.csv(succ.agg.top20, "SuccinicAcidTOP20.csv")

# plot

ggplot(succ.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("black", "red")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# plot 2
ggplot(succ.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))


#save
ggsave(
  filename = "EMSL_MIMOSA_Succinic-Acid.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

#save plot 2
ggsave(
  filename = "EMSL_MIMOSA_Succinic-Acid.TOP20_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

##############################################################################

# try fumaric acid 
fum.acid <- subset(CMP.comps, Compound_name == "Fumaric acid") 
fum.acid.pos <- subset(fum.acid, CMP > 0)
t.test(log(CMP) ~ Condition, data = fum.acid) # the non-transformed do not have sig. results
t.test(CMP ~ Condition, data = fum.acid)
# plot
PP <- ggplot(fum.acid, aes(x=Condition, y=log(CMP), fill=Condition)) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("#446455","#C7B19C"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + theme(axis.text = element_text(size = 12, face="bold"),
                                      axis.title.x = element_blank(),
                                      strip.text = element_text(face = "bold", size=16, color="white"),
                                      strip.background = element_rect(fill = "black"),
                                      axis.title.y = element_text(size = 12, face="bold"))
PP + annotate("text", x=1.5, y=5, label="p=0.043", size=6)
# save plot 
ggsave(
  filename = "EMSL_MIMOSA_Fumaric-Acid.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

#### Welch Two Sample t-test ####

#### data:  log(CMP) by Condition
#### t = -2.0254, df = 1495.8, p-value = 0.04301
#### alternative hypothesis: true difference in means between group EcMF and group No EcMF is not equal to 0
#### 95 percent confidence interval:
####   -0.252765441 -0.004044957
#### sample estimates:
####   mean in group EcMF mean in group No EcMF 
#### 1.901426              2.029831 

# Aggregate Fumaric acid CMPs by ASVs and match ASV to taxa classes
fum.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = fum.acid)
colnames(fum.agg)[1] <- "ASV"
fum.tax <- merge(fum.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition
fum.tax.em <- subset(fum.tax, Condition == "EcMF")
fum.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = fum.tax.em)
fum.tax.em.agg.top20 <- subset(fum.tax.em.agg, CMP > 400)
fum.tax.em.agg.top20$Condition <- "EcMF"
fum.tax.noem <- subset(fum.tax, Condition == "No EcMF")
fum.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = fum.tax.noem)
fum.tax.noem.agg.top20 <- subset(fum.tax.noem.agg, CMP > 270)
fum.tax.noem.agg.top20$Condition <- "No EcMF"

# aggregate these data frames 
fum.agg.top20 <- merge(fum.tax.em.agg.top20, fum.tax.noem.agg.top20, all=T)
fum.agg.top20$Metabolite <- "Fumaric Acid"
fum.agg.top20$Genus[fum.agg.top20$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
fum.agg.top20$Genus[fum.agg.top20$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"
write.csv(fum.agg.top20, "FumaricAcidTOP20.csv")

# plot

ggplot(fum.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_Fumaric-Acid.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)
#save plot 2
ggsave(
  filename = "EMSL_MIMOSA_Fumaric-Acid.TOP20_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# L-Valine 
l.val <- subset(CMP.comps, Compound_name == "L-Valine")
t.test(log(abs(CMP)) ~ Condition, data = l.val) # for some reason, the log transformed comparison gives an error? Maybe because these data are not normally dist.
# plot NOTE THAT these data are the absolute values because plotting would not work otherwise
PP.1 <- ggplot(l.val, aes(x=Condition, y=-log(abs(CMP))), fill=Condition) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black", "red"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + scale_y_continuous(labels = abs) + theme(axis.text = element_text(size = 12, face="bold"),
                                      axis.title.x = element_blank(),
                                      strip.text = element_text(face = "bold", size=16, color="white"),
                                      strip.background = element_rect(fill = "black"),
                                      axis.title.y = element_text(size = 12, face="bold")) # color? # Also, the log trans does not plot???

PP.1 + annotate("text", x=1.5,y=3, label="p=0.0290", size=6)
# save plot 
ggsave(
  filename = "EMSL_MIMOSA_L-Valine.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

#### Welch Two Sample t-test

#### data:  CMP by Condition
#### t = 2.1858, df = 1284.6, p-value = 0.02901
#### alternative hypothesis: true difference in means between group EcMF and group No EcMF is not equal to 0
#### 95 percent confidence interval:
####   0.324426 6.007759
#### sample estimates:
####   mean in group EcMF mean in group No EcMF 
#### -16.41328             -19.57937 

# Aggregate l-valine CMPs by ASVs and match ASV to taxa classes
lval.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = l.val)
colnames(lval.agg)[1] <- "ASV"
lval.tax <- merge(lval.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition
lval.tax.em <- subset(lval.tax, Condition == "EcMF")
lval.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = lval.tax.em)
lval.tax.em.agg.top20 <- subset(lval.tax.em.agg, CMP < -400)
lval.tax.em.agg.top20$Condition <- "EcMF"
lval.tax.noem <- subset(lval.tax, Condition == "No EcMF")
lval.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = lval.tax.noem)
lval.tax.noem.agg.top20 <- subset(lval.tax.noem.agg, CMP < -230)
lval.tax.noem.agg.top20$Condition <- "No EcMF"
# aggregate these data frames 
lval.agg.top20 <- merge(lval.tax.em.agg.top20, lval.tax.noem.agg.top20, all=T)
lval.agg.top20$Metabolite <- "L-Valine"
lval.agg.top20$Genus[lval.agg.top20$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
lval.agg.top20$Genus[lval.agg.top20$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"
write.csv(lval.agg.top20, "L-ValineTOP20.csv")

# plot

ggplot(lval.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("black", "red")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_L-Valine.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Dihydroxyacetone
dihydro <- subset(CMP.comps, Compound_name == "Dihydroxyacetone")
t.test(log(abs(CMP)) ~ Condition, data = dihydro) # p = 0.02601
# plot
PP.2 <- ggplot(dihydro, aes(x=Condition, y=-log(abs(CMP))), fill=Condition) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black", "red"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + 
  scale_y_continuous(labels = abs) + theme(axis.text = element_text(size = 12, face="bold"),
                                                       axis.title.x = element_blank(),
                                                       strip.text = element_text(face = "bold", size=16, color="white"),
                                                       strip.background = element_rect(fill = "black"),
                                                       axis.title.y = element_text(size = 12, face="bold")) 

PP.2 + annotate("text", x=1.5,y=3, label="p=0.0260", size=6)
# save plot 
ggsave(
  filename = "EMSL_MIMOSA_Dihydroxyacetone.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Aggregate Dihydroxyacetone CMPs by ASVs and match ASV to taxa classes
dihydro.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = dihydro)
colnames(dihydro.agg)[1] <- "ASV"
dihydro.tax <- merge(dihydro.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition
dihydro.tax.em <- subset(dihydro.tax, Condition == "EcMF")
dihydro.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = dihydro.tax.em)
dihydro.tax.em.agg.top20 <- subset(dihydro.tax.em.agg, CMP < -60)
dihydro.tax.em.agg.top20$Condition <- "EcMF"
dihydro.tax.noem <- subset(dihydro.tax, Condition == "No EcMF")
dihydro.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = dihydro.tax.noem)
dihydro.tax.noem.agg.top20 <- subset(dihydro.tax.noem.agg, CMP < -60)
dihydro.tax.noem.agg.top20$Condition <- "No EcMF"
# aggregate these data frames 
dihydro.agg.top20 <- merge(dihydro.tax.em.agg.top20, dihydro.tax.noem.agg.top20, all=T)
dihydro.agg.top20$Metabolite <- "Dihydroxyacetone"
dihydro.agg.top20$Genus[dihydro.agg.top20$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
dihydro.agg.top20$Genus[dihydro.agg.top20$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"
write.csv(dihydro.agg.top20, "DihydroxyacetoneTOP20.csv")

# plot

ggplot(dihydro.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("black", "red")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_Dihdroxyacetone.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

#L-Lactic acid
lac.acid <- subset(CMP.comps, Compound_name == "L-Lactic acid")
t.test(log(abs(CMP)) ~ Condition, data = lac.acid) # neither the log-transformed nor the standard data give significant results 

# L-Threonine
l.threo <- subset(CMP.comps, Compound_name == "L-Threonine") 
l.threo.pos <- subset(l.threo, CMP > 0)
t.test(CMP ~ Condition, data = l.threo) # only the non-log transformed data produce significant results; p = 0.007216 
t.test(CMP ~ Condition, data = l.threo.pos)
# plot
PP2 <- ggplot(l.threo, aes(x=Condition, y=log(CMP))) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black", "red"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + theme(axis.text = element_text(size = 12, face="bold"),
                                     axis.title.x = element_blank(),
                                     strip.text = element_text(face = "bold", size=16, color="white"),
                                     strip.background = element_rect(fill = "black"),
                                     axis.title.y = element_text(size = 12, face="bold")) # should be tan to light brown color theme
PP2 + annotate("text", x=1.5, y=5, label="p=0.0072", size=6)
PP2.1 <- ggplot(l.threo.pos, aes(x=Condition, y=log(CMP))) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black", "red"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + theme(axis.text = element_text(size = 12, face="bold"),
                                                                                            axis.title.x = element_blank(),
                                                                                            strip.text = element_text(face = "bold", size=16, color="white"),
                                                                                            strip.background = element_rect(fill = "black"),
                                                                                            axis.title.y = element_text(size = 12, face="bold"))
PP2.1 + annotate("text", x=1.5, y=5, label="p=0.0289", size=6)

# save plot
ggsave(
  filename = "EMSL_MIMOSA_L-Threonine.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Aggregate L-Threonine CMPs by ASVs and match ASV to taxa classes
l.threo.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = l.threo)
colnames(l.threo.agg)[1] <- "ASV"
l.threo.tax <- merge(l.threo.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition
l.threo.tax.em <- subset(l.threo.tax, Condition == "EcMF")
l.threo.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = l.threo.tax.em)
l.threo.tax.em.agg.top20 <- subset(l.threo.tax.em.agg, CMP > 45)
l.threo.tax.em.agg.top20$Condition <- "EcMF"
l.threo.tax.noem <- subset(l.threo.tax, Condition == "No EcMF")
l.threo.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = l.threo.tax.noem)
l.threo.tax.noem.agg.top20 <- subset(l.threo.tax.noem.agg, CMP > 30)
l.threo.tax.noem.agg.top20$Condition <- "No EcMF"

# aggregate these data frames 
l.threo.agg.top20 <- merge(l.threo.tax.em.agg.top20, l.threo.tax.noem.agg.top20, all=T)
l.threo.agg.top20$Metabolite <- "L-Threonine"
l.threo.agg.top20$Genus[l.threo.agg.top20$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
l.threo.agg.top20$Genus[l.threo.agg.top20$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"
write.csv(l.threo.agg.top20, "L-ThreonineTOP20.csv")

# plot

ggplot(l.threo.agg.top20, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("#446455","#C7B19C")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10,angle=45,hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_L-Threonine.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Palmitic acid
palm.acid <- subset(CMP.comps, Compound_name == "Palmitic acid")
t.test(log(CMP) ~ Condition, data = palm.acid) # p = 0.07227 /0.1055

# Glyceric acid
gly.acid <- subset(CMP.comps, Compound_name == "Glyceric acid")
t.test(CMP ~ Condition, data = gly.acid) # p = 0.368

# Maleic acid
mal.acid <- subset(CMP.comps, Compound_name == "Maleic acid")
t.test(log(CMP) ~ Condition, data = mal.acid) # p = 0.5193 / 0.625

# Stearic acid 
ster.acid <- subset(CMP.comps, Compound_name == "Stearic acid")
t.test(log(CMP) ~ Condition, data = ster.acid) # p = 0.863 / p = 0.6265

# Phthalic acid
phth.acid <- subset(CMP.comps, Compound_name == "Phthalic acid") 
t.test(CMP ~ Condition, data = phth.acid) # p = 0.05084 # if data are normally distributed, then do not log transform
# plot
PP3 <- ggplot(phth.acid, aes(x=Condition, y=log(CMP))) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black", "red"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") + theme(axis.text = element_text(size = 12, face="bold"),
                                      axis.title.x = element_blank(),
                                      strip.text = element_text(face = "bold", size=16, color="white"),
                                      strip.background = element_rect(fill = "black"),
                                      axis.title.y = element_text(size = 12, face="bold")) # should be white color theme
PP3 + annotate("text", x=1.5,y=3.5, label="p=0.051", size=6)

#plot 
ggsave(
  filename = "EMSL_MIMOSA_Phthalic-Acid.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Aggregate Phthalic acid CMPs by ASVs and match ASV to taxa classes
phth.acid.agg <- aggregate(CMP ~ Species + Condition, FUN = sum, data = phth.acid)
colnames(phth.acid.agg)[1] <- "ASV"
phth.acid.tax <- merge(phth.acid.agg, TAX_ASV, by ="ASV")

# subset based on the top 20 bacterial taxa for each condition # not 20 
phth.acid.tax.em <- subset(phth.acid.tax, Condition == "EcMF")
phth.acid.tax.em.agg <- aggregate(CMP ~ Genus, FUN = sum, data = phth.acid.tax.em)
# phth.acid.tax.em.agg.top20 <- subset(phth.acid.tax.em.agg, CMP > 45) # ignore
phth.acid.tax.em.agg$Condition <- "EcMF"
phth.acid.tax.noem <- subset(phth.acid.tax, Condition == "No EcMF")
phth.acid.tax.noem.agg <- aggregate(CMP ~ Genus, FUN = sum, data = phth.acid.tax.noem)
#phth.acid.tax.noem.agg.top20 <- subset(phth.acid.tax.noem.agg, CMP > 30) # ignore
phth.acid.tax.noem.agg$Condition <- "No EcMF"

# aggregate these data frames 
phth.acid.agg <- merge(phth.acid.tax.em.agg, phth.acid.tax.noem.agg, all=T)
phth.acid.agg$Metabolite <- "Phthalic acid"
phth.acid.agg$Genus[phth.acid.agg$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
phth.acid.agg$Genus[phth.acid.agg$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"
write.csv(phth.acid.agg, "Phthalic.Acid.ALL.csv")

# plot

ggplot(phth.acid.agg, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("black", "red")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_Phthalic-Acid.TOP20.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)

# Pyroglutamic acid
pyro.acid <- subset(CMP.comps, Compound_name == "Pyroglutamic acid")
t.test(log(CMP) ~ Condition, data = pyro.acid) # p = 0.2087

# Adipic acid
adip.acid <- subset(CMP.comps, Compound_name == "Adipic acid")
t.test(log(CMP) ~ Condition, data = adip.acid) # p = 0.7644 / 0.8289

# Caprylic acid
capry.acid <- subset(CMP.comps, Compound_name == "Caprylic acid")
PP4 <- ggplot(capry.acid, aes(x=Condition, y=log(CMP))) + facet_wrap(~Compound_name) +
  geom_boxplot(fill=c("black"), alpha=0.6) + theme_bw() + labs(y="Log[CMP]") +
  theme(axis.text = element_text(size = 12, face="bold"),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size=16, color="white"),
        strip.background = element_rect(fill = "black"),
        axis.title.y = element_text(size = 12, face="bold")) 
# plot
ggsave(
  filename = "EMSL_MIMOSA_Caprylic-Acid.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 3,
  height = 6,
  units = c("in"),
  dpi = 300)
t.test(CMP ~ Condition, data = capry.acid) # this is only present in the EcMF condition, so no comparison is needed / available

colnames(capry.acid)[2] <- "ASV"
capry.acid.agg <- aggregate(CMP ~ ASV, FUN=sum, data=capry.acid)
capry.acid.agg.tax <- merge(capry.acid.agg, TAX_ASV, by="ASV")
capry.acid.agg.tax$Condition <- "EcMF"
capry.acid.agg.tax$Metabolite <- "Caprylic acid"
capry.acid.agg.tax2 <- capry.acid.agg.tax[,c(2,8,10,11)]

# plot

ggplot(capry.acid.agg.tax2, aes(x=Genus,y=CMP,fill=Condition)) + geom_col() + 
  theme_bw() + coord_flip() + facet_grid(~Condition) + scale_fill_manual(values=c("black", "red")) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size=10, face="bold.italic"), axis.text.x = element_text(size=10)) +
  theme(axis.title.y = element_blank()) +
  theme(strip.text = element_text(size=12, face="bold"))

# save
ggsave(
  filename = "EMSL_MIMOSA_Caprylic-Acid.TOP20.ONLY3.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 6,
  units = c("in"),
  dpi = 300)
# merge all top metabolite data together 

# succinic acid, fumaric acid, caprylic acid, dihydroxyacetone, l-threonine, l-valine, phthalic acid

master.met <- merge(succ.agg.top20,l.threo.agg.top20, all=T)
master.met2 <- merge(master.met, phth.acid.agg, all=T)
master.met3 <- merge(master.met2, dihydro.agg.top20, all=T)
master.met4 <- merge(master.met3, fum.agg.top20, all=T)
master.met5 <- merge(master.met4, lval.agg.top20, all=T)
master.met6 <- merge(master.met5, capry.acid.agg.tax2, all=T)
master.met6$Genus[master.met6$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
master.met6$Genus[master.met6$Genus == "Methylobacterium-Methylorubrum"] <- "Methylobacterium"

# save
write.csv(master.met6, "Master.Met.Top.MIMOSA.csv")

# tidy up all code prior to GitHub push
# plot

ggplot(master.met6, aes(x=Genus, y=CMP)) + facet_grid(~Metabolite) + theme_bw() + geom_col(fill=c("black", "red")) + coord_flip()
