################################################################################
############## This script contains classification of mapped reads #############
################################################################################
pacman::p_load(readxl, dplyr, ggplot2, tidyr, wrappedtools, RColorBrewer, forcats,
               ggrepel, patchwork, pals)

# Leaf
leaf_summary <- read_excel("PATH_HERE/MappedRegion_summary.xlsx")

leaf_summary <-  leaf_summary |> 
  separate(col = "Input_1",
           into = c("Leaf 1", "Leaf 1 (%)"),
           sep = "\\(") |> 
  separate(col = "Input_2",
           into = c("Leaf 2", "Leaf 2 (%)"),
           sep = "\\(") |> 
  separate(col = "Input_3",
           into = c("Leaf 3", "Leaf 3 (%)"),
           sep = "\\(")

for (i in seq_along(leaf_summary$Samples)) {
  leaf_summary$`Leaf 1 (%)`[i] <- gsub("%\\)", "", leaf_summary$`Leaf 1 (%)`[i])
  leaf_summary$`Leaf 2 (%)`[i] <- gsub("%\\)", "", leaf_summary$`Leaf 2 (%)`[i])
  leaf_summary$`Leaf 3 (%)`[i] <- gsub("%\\)", "", leaf_summary$`Leaf 3 (%)`[i])
}

columns <- cn(leaf_summary)
columns_to_convert <- columns[-1]
leaf_summary <- leaf_summary |> 
  mutate(across(all_of(columns_to_convert), as.double))

leaf_mean <- leaf_summary |> 
  select(c("Leaf 1", "Leaf 2", "Leaf 3")) |> 
  rowMeans()
leaf_100 <- leaf_summary |> 
  select(c("Leaf 1 (%)", "Leaf 2 (%)", "Leaf 3 (%)")) |> 
  rowMeans()

leaf_summary <- leaf_summary |> 
  mutate("Samples"= as_factor(Samples),
         "Counts_Mean"= leaf_mean,
         "Percentage_Mean"= leaf_100)

leaf_summary <- leaf_summary[
   leaf_summary$Samples %in% c("Others",
                               "protein_coding",
                               "sense_intronic",
                               "snoRNA",
                               "SRP_RNA"), ]
 

################################################################################
# Phloem
phloem_summary <- read_excel("PATH_HERE/MappedRegion_summary.xlsx")

phloem_summary <-  phloem_summary |> 
  separate(col = "Phloem1",
           into = c("Phloem 1", "Phloem 1 (%)"),
           sep = "\\(") |> 
  separate(col = "Phloem2",
            into = c("Phloem 2", "Phloem 2 (%)"),
            sep = "\\(") |> 
  separate(col = "Phloem3",
           into = c("Phloem 3", "Phloem 3 (%)"),
           sep = "\\(")

for (i in seq_along(phloem_summary$Samples)) {
   phloem_summary$`Phloem 1 (%)`[i] <- gsub("%\\)", "", phloem_summary$`Phloem 1 (%)`[i])
   phloem_summary$`Phloem 2 (%)`[i] <- gsub("%\\)", "", phloem_summary$`Phloem 2 (%)`[i])
   phloem_summary$`Phloem 3 (%)`[i] <- gsub("%\\)", "", phloem_summary$`Phloem 3 (%)`[i])
}

columns <- cn(phloem_summary)
columns_to_convert <- columns[-1]
phloem_summary <- phloem_summary |> 
  mutate(across(all_of(columns_to_convert), as.double))

phloem_mean <- phloem_summary |> 
  select(c("Phloem 1", "Phloem 2", "Phloem 3")) |> 
  rowMeans()
phloem_100 <- phloem_summary |> 
  select(c("Phloem 1 (%)", "Phloem 2 (%)", "Phloem 3 (%)")) |> 
  rowMeans()
  

phloem_summary <- phloem_summary |> 
  mutate("Samples"= as_factor(Samples),
         "Counts_Mean"= phloem_mean,
         "Percentage_Mean"= phloem_100)

phloem_summary <- phloem_summary[
  phloem_summary$Samples %in% c("Others",
                              "protein_coding",
                              "sense_intronic",
                              "snoRNA",
                              "SRP_RNA"), ]


################################################################################
################## ##### Plotting: Pie chart ###################################
################################################################################
## Leaf
leaf_plot <- ggplot(leaf_summary,
                    aes(x="", y=Percentage_Mean,
                        fill=(paste0(Samples, " (", round(Percentage_Mean, digits=2), "%)")))) +
  geom_col() +
  coord_polar("y")+
  #scale_fill_brewer(palette="Blues")+
  scale_fill_manual(values=c("#A8D0E6", "#23305E", "#39424E", "#F99797", "#F76C6C"))+
  guides(fill=guide_legend(title=""))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 20, b = 20),
        text=element_text(size=35),
        legend.key.size = unit(1, 'cm'))

ggsave(filename = "mapping/Leaf_Mapped_Classes.png",
       plot = leaf_plot,
       width = 10, height = 8,
       dpi = 600)

## Phloem
phloem_plot <- ggplot(phloem_summary,
                      aes(x="", y=Percentage_Mean,
                          fill=paste0(Samples, " (", round(Percentage_Mean, digits=2), "%)"))) +
  geom_col() +
  coord_polar("y")+
  #scale_fill_brewer(palette="Accent")+
  scale_fill_manual(values=c("#A8D0E6", "#23305E", "#39424E", "#F99797", "#F76C6C"))+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title=""))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 20, b = 20),
        text=element_text(size=35),
        legend.key.size = unit(1, 'cm'))

ggsave(filename = "mapping/Phloem_Mapped_Classes.png",
       plot = phloem_plot,
       width = 10, height = 8,
       dpi = 600)

## Combining plots with patchwork
# combined_plots <-  (leaf_plot|phloem_plot) +
#   plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 24, face='bold'))
# 
# ## Saving into pdf
# pdf("Mapping_AllReadsClassification.pdf")
# print(combined_plots)
# dev.off()
