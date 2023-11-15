##### Processing the discriminant ion EMSL data (from roots) #####

# Data obtained: 6/27/2023

# Format: These data are in an excel spreadsheet. 
# Purpose: These data show the intensity (or amount) of the metabolites
#          in the roots of Bishop pine plants either not colonized by EcM
#          or colonized by EcM (Suillus pungens). 


# Import data frame

met.data <- read.csv("metabolite.root.data.csv")

# This dataframe needs to be reformatted. 
# Removed the unncessary columns and blank rows in the spreadsheet 
#               renamed to 'metabolite.root.data.adj.csv'

# Import adjusted data frame
met.data.adj <- read.csv("metabolite.root.data.adj.csv")

# subset only the significant comparisons (AUC > 0.75 --> NM enriched; 
#                                          AUC < 0.25 --> SP enriched)


NM.met.adj.sub <- subset(met.data.adj, AUC >= 0.75)
SP.met.adj.sub <- subset(met.data.adj, AUC <= 0.35) # AUC , 0.35 also enriched in SP

# none of the AUC values are equal to or less than 0.25. 

# Save data frames

write.csv(NM.met.adj.sub, "NM_Root_Enriched.csv")
write.csv(SP.met.adj.sub, "SP_Root_Enriched.csv")

# Create a new data frame that only includes ions that have a single corresponding
# metabolite (for both SP and NM conditions)

# NM metabolites
single.metabolite <- subset(NM.met.adj.sub, molecule == "N-Monomethyl-2-aminoethylphosphonate" | 
                              molecule == "Xanthopterin-B2" | molecule == "Rhein" | molecule == "Riccionidin A" | 
                              molecule == "Indoleglycerol phosphate" | molecule == "N-Feruloylglycine" | molecule == "Dyphylline" | 
                              molecule == "5'-Oxoinosine" | molecule == "8-Oxodeoxycoformycin" | molecule == "Inosine" | molecule == "Abacavir" | 
                              molecule == "Cyclopiazonic acid" | molecule == "Pyridoxamine phosphate")

# Add condition column (No EcMF)
single.metabolite$Condition <- "No EcMF"

# Add condition column for EcMF
SP.met.adj.sub$Condition <- "EcMF"

# Merge data frames
single.met.SP.NM <- merge(single.metabolite, SP.met.adj.sub, all=TRUE)

# Add new column for facetting
single.met.SP.NM$Title <- "Root Metabolites"

# Plot these data
library(ggplot2)

ggplot(single.met.SP.NM, aes(x=reorder(molecule, AUC), y=AUC, color=Condition)) + geom_point(size=6) + theme_bw() + coord_flip() + scale_color_manual(values = c("#446455","#C7B19C")) + xlab("") + facet_grid(~Title) + 
  theme(strip.text = element_text(color = "white", size = 12, face = "bold")) +
  theme(strip.background = element_rect(fill = "black")) + labs(colour = "Condition") + 
  theme(axis.text = element_text(face = "bold")) + 
  ylab("AUC") + theme(axis.title = element_text(face = "bold")) + 
  theme(legend.title = element_text(face = "bold")) + theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  geom_segment(aes(x = molecule, xend = molecule, y = 0, yend = AUC))

# Save plot

ggsave(
  filename = "EMSL_Root_Metabolites.SingleMetabolites_Revised.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 6,
  units = c("in"),
  dpi = 300)




