pacman::p_load(readxl, dplyr, tidyr, writexl, rtracklayer, ggplot2)

################################## Leaf ########################################
leaf_locations <- readGFF("PATH_HERE/novel_lncRNA.gtf")

leaf_exons <- leaf_locations |> 
  group_by(transcript_id) |> 
  summarize(leaf_exons = n()) |> #1142
  mutate(leaf_exons = as.factor(leaf_exons))

leaf_exons_counts <- as_tibble(table(leaf_exons$leaf_exons), .name_repair = "universal") |> 
  rename("...1" = "Class",
         "n" = "Leaf") |> 
  rbind(tibble("Class" = 12,
               "Leaf" = 0)) |> 
  mutate(Class = factor(Class, levels=(c(2:12))))

################################# Phloem #######################################
phloem_locations <- readGFF("PATH_HERE/novel_lncRNA.gtf")

phloem_exons <- phloem_locations |> 
  group_by(transcript_id) |> 
  summarize(phloem_exons = n()) |> #2326
  mutate(phloem_exons = as.factor(phloem_exons))

phloem_exons_counts <- as_tibble(table(phloem_exons$phloem_exons), .name_repair = "universal") |> 
  rename("...1" = "Class",
         "n" = "Phloem") |> 
  mutate(Class = factor(Class, levels=(c(2:12))))

################################################################################
graph_table <- as_tibble(merge(phloem_exons_counts, leaf_exons_counts)) |> 
  arrange(Class)

graph_table_normalized <- graph_table |> 
  mutate_at(vars(Leaf, Phloem), function(x) x / sum(x) * 100) |> 
  pivot_longer(cols = c(Leaf, Phloem), names_to = "Tissue", values_to = "Percentage") |> 
  mutate(Tissue = factor(Tissue, levels= c("Phloem", "Leaf")))

strand_bar <- ggplot(graph_table_normalized, aes(x = Class, y = Percentage, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Exon number",
       y = "Percentage",
       fill= NULL) +
  scale_fill_manual(values=c("#23305E",  "#F76C6C"))+
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(size=40),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(linewidth = 2),
        legend.key.size = unit(2, 'cm'),
        aspect.ratio=1/1)+
  scale_y_continuous(expand = c(0, 0))

ggsave(filename = "lncRNAs/lncRNAs_ExonDistribution.png",
       plot = strand_bar,
       width = 10, height = 8,
       dpi = 600)
saveRDS(strand_bar, file='lncRNAs/lncRNAs_ExonDistribution.RDS')

