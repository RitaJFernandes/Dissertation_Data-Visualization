################################################################################
## This script contains the FPKM distribution of leaf and phloem samples for  ##:
##   -All transcripts                                                         ##
##   -mRNAs                                                                   ##
##   -lncRNAs                                                                 ##
################################################################################
pacman::p_load(readxl, dplyr, ggplot2, tidyr, stringr, RColorBrewer)

###################### Across biological replicates ###########################
# All transcripts
leaf_transcripts_FPKM <- read_excel("PATH_HERE/transcripts.FPKM.xlsx")
phloem_transcripts_FPKM <- read_excel("PATH_HERE/transcripts.FPKM.xlsx")

leaf_log_transcripts <- leaf_transcripts_FPKM |> 
  mutate(Input_1=log10(leaf_transcripts_FPKM$Input_1+1),
         Input_2=log10(leaf_transcripts_FPKM$Input_2+1),
         Input_3=log10(leaf_transcripts_FPKM$Input_3+1)) |>
  rename(`Leaf 1`=Input_1,
         `Leaf 2`=Input_2,
         `Leaf 3`=Input_3)

phloem_log_transcripts <- phloem_transcripts_FPKM |> 
  mutate(Phloem1=log10(phloem_transcripts_FPKM$Phloem1+1),
         Phloem2=log10(phloem_transcripts_FPKM$Phloem2+1),
         Phloem3=log10(phloem_transcripts_FPKM$Phloem3+1)) |>
  rename(`Phloem 1`=Phloem1,
         `Phloem 2`=Phloem2,
         `Phloem 3`=Phloem3)


leaf_log_transcripts_long <- pivot_longer(data=leaf_log_transcripts,
                                                 cols=-transcript_id,
                                                 names_to = "Sample",
                                                 values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))

phloem_log_transcripts_long <- pivot_longer(data=phloem_log_transcripts,
                                            cols=-transcript_id,
                                            names_to = "Sample",
                                            values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))


transcripts_FPKM_plot <- ggplot()+
  geom_boxplot()+
  geom_boxplot(data = leaf_log_transcripts_long,
               aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample), alpha=.7)+
  geom_boxplot(data = phloem_log_transcripts_long,
               aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample), alpha=.7)+
  theme(axis.text = element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=35),
        legend.position = "none",
        aspect.ratio=1/1)+
  scale_fill_manual(values=c("#A8D0E6","#474386", "#23305E", "#F0D0D5","#F99797", "#F76C6C"))

ggsave(filename = "quantification/transcripts_FPKM.png",
       plot = transcripts_FPKM_plot,
       width = 10, height = 8,
       dpi = 600)


# lncRNAs
leaf_lnc_transcripts_FPKM <- read_excel("PATH_HERE/lncRNA.transcripts.FPKM.xlsx")
phloem_lnc_transcripts_FPKM <- read_excel("PATH_HERE/lncRNA.transcripts.FPKM.xlsx")

leaf_log_lnc_transcripts <- leaf_lnc_transcripts_FPKM |> 
  mutate(Input_1=log10(leaf_lnc_transcripts_FPKM$Input_1+1),
         Input_2=log10(leaf_lnc_transcripts_FPKM$Input_2+1),
         Input_3=log10(leaf_lnc_transcripts_FPKM$Input_3+1)) |>
  rename(`Leaf 1`=Input_1,
         `Leaf 2`=Input_2,
         `Leaf 3`=Input_3)

phloem_log_lnc_transcripts <- phloem_lnc_transcripts_FPKM |> 
  mutate(Phloem1=log10(phloem_lnc_transcripts_FPKM$Phloem1+1),
         Phloem2=log10(phloem_lnc_transcripts_FPKM$Phloem2+1),
         Phloem3=log10(phloem_lnc_transcripts_FPKM$Phloem3+1)) |>
  rename(`Phloem 1`=Phloem1,
         `Phloem 2`=Phloem2,
         `Phloem 3`=Phloem3)


leaf_log_lnc_transcripts_long <- pivot_longer(data=leaf_log_lnc_transcripts,
                                          cols=-transcript_id,
                                          names_to = "Sample",
                                          values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))

phloem_log_lnc_transcripts_long <- pivot_longer(data=phloem_log_lnc_transcripts,
                                            cols=-transcript_id,
                                            names_to = "Sample",
                                            values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))


lncRNA_transcripts_FPKM_plot <- ggplot()+
  geom_boxplot()+
  geom_boxplot(data = leaf_log_lnc_transcripts_long,
               aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample), alpha=.7)+
  geom_boxplot(data = phloem_log_lnc_transcripts_long,
               aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample), alpha=.7)+
  theme(axis.text = element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=35),
        legend.position = "none",
        aspect.ratio=1/1)+
  scale_fill_manual(values=c("#A8D0E6","#474386", "#23305E", "#F0D0D5","#F99797", "#F76C6C"))

ggsave(filename = "quantification/lncRNA.transcripts_FPKM.png",
       plot = lncRNA_transcripts_FPKM_plot,
       width = 10, height = 8,
       dpi = 600)

# leaf_transcripts_FPKM_plot <- ggplot(data = leaf_log_transcripts_long,
#                                        aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
#   geom_boxplot(alpha=0.5)+
#   theme(axis.text = element_text(size=15),
#         axis.title = element_text(size=18),
#         legend.position = "none")
# 
# phloem_transcripts_FPKM_plot <- ggplot(data = phloem_log_transcripts_long,
#                                   aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
#   geom_boxplot(alpha=0.5)+
#   theme(axis.text = element_text(size=15),
#         axis.title = element_text(size=18),
#         legend.position = "none")

# pdf("transcripts_FPKM.pdf")
# print(leaf_transcripts_FPKM_plot)
# print(phloem_transcripts_FPKM_plot)
# dev.off()
################################################################################
# mRNAS
leaf_mRNA_FPKM <- read_excel("PATH_HERE/mRNA.transcripts.FPKM.xlsx")
phloem_mRNA_FPKM <- read_excel("PATH_HERE/mRNA.transcripts.FPKM.xlsx")

leaf_log_mRNA <- leaf_mRNA_FPKM |> 
  mutate(Input_1=log10(leaf_mRNA_FPKM$Input_1+1),
         Input_2=log10(leaf_mRNA_FPKM$Input_2+1),
         Input_3=log10(leaf_mRNA_FPKM$Input_3+1)) |>
  rename(`Leaf 1`=Input_1,
         `Leaf 2`=Input_2,
         `Leaf 3`=Input_3)

phloem_log_mRNA <- phloem_mRNA_FPKM |> 
  mutate(Phloem1=log10(phloem_mRNA_FPKM$Phloem1+1),
         Phloem2=log10(phloem_mRNA_FPKM$Phloem2+1),
         Phloem3=log10(phloem_mRNA_FPKM$Phloem3+1)) |>
  rename(`Phloem 1`=Phloem1,
         `Phloem 2`=Phloem2,
         `Phloem 3`=Phloem3)


leaf_log_mRNA_long <- pivot_longer(data=leaf_log_mRNA,
                                          cols=-transcript_id,
                                          names_to = "Sample",
                                          values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))

phloem_log_mRNA_long <- pivot_longer(data=phloem_log_mRNA,
                                            cols=-transcript_id,
                                            names_to = "Sample",
                                            values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))


leaf_mRNA_FPKM_plot <- ggplot(data = leaf_log_mRNA_long,
                         aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
  geom_boxplot(alpha=0.5)+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.position = "none")

phloem_mRNA_FPKM_plot <- ggplot(data = phloem_log_mRNA_long,
                           aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
  geom_boxplot(alpha=0.5)+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.position = "none")

pdf("mRNA_FPKM.pdf")
print(leaf_mRNA_FPKM_plot)
print(phloem_mRNA_FPKM_plot)
dev.off()
################################################################################
# lncRNAs
leaf_lncRNA_FPKM <- read_excel("PATH_HERE/lncRNA.transcripts.FPKM.xlsx")
phloem_lncRNA_FPKM <- read_excel("PATH_HERE/lncRNA.transcripts.FPKM.xlsx")

leaf_log_lncRNA <- leaf_lncRNA_FPKM |> 
  mutate(Input_1=log10(leaf_lncRNA_FPKM$Input_1+1),
         Input_2=log10(leaf_lncRNA_FPKM$Input_2+1),
         Input_3=log10(leaf_lncRNA_FPKM$Input_3+1)) |>
  rename(`Leaf 1`=Input_1,
         `Leaf 2`=Input_2,
         `Leaf 3`=Input_3)

phloem_log_lncRNA <- phloem_lncRNA_FPKM |> 
  mutate(Phloem1=log10(phloem_lncRNA_FPKM$Phloem1+1),
         Phloem2=log10(phloem_lncRNA_FPKM$Phloem2+1),
         Phloem3=log10(phloem_lncRNA_FPKM$Phloem3+1)) |>
  rename(`Phloem 1`=Phloem1,
         `Phloem 2`=Phloem2,
         `Phloem 3`=Phloem3)


leaf_log_lncRNA_long <- pivot_longer(data=leaf_log_lncRNA,
                                     cols=-transcript_id,
                                     names_to = "Sample",
                                     values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))

phloem_log_lncRNA_long <- pivot_longer(data=phloem_log_lncRNA,
                                            cols=-transcript_id,
                                            names_to = "Sample",
                                            values_to = "log10(FPKM+1)") |>
  mutate(Sample=as.factor(Sample))


leaf_lncRNA_FPKM_plot <- ggplot(data = leaf_log_lncRNA_long,
                           aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
  geom_boxplot(alpha=0.5)+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.position = "none")

phloem_lncRNA_FPKM_plot <- ggplot(data = phloem_log_lncRNA_long,
                             aes(x=Sample, y=`log10(FPKM+1)`, fill=Sample))+
  geom_boxplot(alpha=0.5)+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        legend.position = "none")

pdf("lncRNA_FPKM.pdf")
print(leaf_lncRNA_FPKM_plot)
print(phloem_lncRNA_FPKM_plot)
dev.off()

################################################################################
######################## Between leaf and phloem ###############################
################################################################################
leaf_transcripts_mean <- as_tibble(rowMeans(leaf_transcripts_FPKM[-1])) #119716
leaf_mRNA_mean <- as_tibble(rowMeans(leaf_mRNA_FPKM[-1])) #111423, missing 8293
leaf_lncRNA_mean <- as_tibble(rowMeans(leaf_lncRNA_FPKM[-1])) #1213, missing 118503

empty_leafmRNA <- tibble(value = NA_real_, .rows = 8293)
empty_leaflncRNA <- tibble(value = NA_real_, .rows = 118503) # difference of rows between tables

temp_leafmRNA <- leaf_mRNA_mean |> rbind(empty_leafmRNA) |> rename(mRNA = value)
temp_leaflncRNA <- leaf_lncRNA_mean |> rbind(empty_leaflncRNA) |> rename(lncRNA = value)

leaf_combined <- leaf_transcripts_mean |>  rename("All transcripts" = value) |>
  cbind(temp_leafmRNA, temp_leaflncRNA) |> 
  pivot_longer(cols= c("All transcripts", mRNA, lncRNA),
               values_to = "FPKM" ,
               names_to = "Leaf") |> 
  mutate(Leaf = factor(Leaf, levels= c("All transcripts", "mRNA", "lncRNA"))) |> 
  filter(!is.na(FPKM)) |> 
  mutate(`log10(FPKM+1)` = log10(FPKM+1)) #232352 rows.


ggplot(leaf_combined, aes(x=Leaf, y=`log10(FPKM+1)`, fill=Leaf))+
  geom_boxplot()


phloem_transcripts_mean <- as_tibble(rowMeans(phloem_transcripts_FPKM[-1])) #130836
phloem_mRNA_mean <- as_tibble(rowMeans(phloem_mRNA_FPKM[-1])) #118902, missing 11934
phloem_lncRNA_mean <- as_tibble(rowMeans(phloem_lncRNA_FPKM[-1])) #2397, missing 128439

empty_phloemmRNA <- tibble(value = NA_real_, .rows = 11934)
empty_phloemlncRNA <- tibble(value = NA_real_, .rows = 128439) # difference of rows between tables

temp_phloemmRNA <- phloem_mRNA_mean |> rbind(empty_phloemmRNA) |> rename(mRNA = value)
temp_phloemlncRNA <- phloem_lncRNA_mean |> rbind(empty_phloemlncRNA) |> rename(lncRNA = value)

phloem_combined <- phloem_transcripts_mean |>  rename("All transcripts" = value) |>
  cbind(temp_phloemmRNA, temp_phloemlncRNA) |> 
  pivot_longer(cols= c("All transcripts", mRNA, lncRNA),
               values_to = "FPKM" ,
               names_to = "Phloem") |> 
  mutate(Phloem = factor(Phloem, levels= c("All transcripts", "mRNA", "lncRNA"))) |> 
  filter(!is.na(FPKM)) |> 
  mutate(`log10(FPKM+1)` = log10(FPKM+1)) #252135 rows

ggplot(phloem_combined, aes(x=Phloem, y=`log10(FPKM+1)`, fill=Phloem))+
  geom_boxplot()

# Combining leaf and phloem into one graph
empty_combined <- tibble(Leaf = NA_real_,
                         `Leaf FPKM` = NA_real_,
                         value = NA_real_, .rows = 19783) |> #232352+19783=252135 --> phloem_combined row#
  rename(`Leaf log10(FPKM+1)`=value)

temp_leaf_combined <- leaf_combined |> 
  rename(`Leaf FPKM` = FPKM,
         `Leaf log10(FPKM+1)`= `log10(FPKM+1)`) |>
  rbind(empty_combined)

both_combined <- phloem_combined |> 
  rename(`Phloem FPKM` = FPKM,
         `Phloem log10(FPKM+1)`= `log10(FPKM+1)`) |> 
  cbind(temp_leaf_combined) |> 
  pivot_longer(cols= c(`Phloem log10(FPKM+1)`, `Leaf log10(FPKM+1)`),
               values_to = "log10(FPKM+1)" ,
               names_to = "Tissue") |> 
  select(-`Phloem FPKM`, -`Leaf FPKM`) |> 
  filter(!is.na(`log10(FPKM+1)`)) |> 
  mutate(Tissue= gsub(" log10\\(FPKM\\+1\\)",
                      "", 
                      as.character(Tissue)))

All_FPKM_Dist <- ggplot(both_combined, aes(x=Tissue, y=`log10(FPKM+1)`, fill=Phloem))+
  geom_boxplot(outlier.alpha = 0.1)+
  #theme_bw()+
  labs(x=NULL)+
  theme(text = element_text(size=40),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #axis.line = element_line(colour = "black", linewidth = 2),
        #axis.ticks.length = unit(0.3, "cm"),
        #axis.ticks = element_line(linewidth = 2),
        legend.key.size = unit(2, 'cm'),
        aspect.ratio=1/1,
        legend.position="bottom",
        legend.title = element_text(size = 29))+
  guides(fill=guide_legend(title=""))+
  scale_fill_brewer(palette="Set1")

ggsave(filename = "quantification/All_FPKM_Distribution.png",
       plot = All_FPKM_Dist,
       width = 10, height = 8,
       dpi = 600)
#saveRDS(All_FPKM_Dist, file='quantification/All_FPKM_Distribution.RDS')
