pacman::p_load(dplyr, tidyr, writexl, rtracklayer, ggplot2)

leaf_strands <- readGFF("PATH_HERE/novel_lncRNA.gtf") |> 
  select(strand) |> 
  rename("strand"="Leaf")

phloem_strands <- readGFF("PATH_HERE/novel_lncRNA.gtf") |> 
  select(strand) |> 
  rename("strand"="Phloem")

leaf_length_sense <- length(which(leaf_strands$Leaf=='+'))
leaf_length_antisense <- length(which(leaf_strands$Leaf=='-'))
phloem_length_sense <- length(which(phloem_strands$Phloem=='+'))
phloem_length_antisense <- length(which(phloem_strands$Phloem=='-'))   
                            
                            
graph_table <- tibble(
  "Leaf" = c(leaf_length_sense, leaf_length_antisense),
  "Phloem" = c(phloem_length_sense, phloem_length_antisense),
  "Strand" = factor(c("+", "-")))

# Normalize percentages within each tissue
graph_table_normalized <- graph_table %>%
  mutate_at(vars(Leaf, Phloem), function(x) x / sum(x) * 100) %>%
  pivot_longer(cols = c(Leaf, Phloem), names_to = "Tissue", values_to = "Percentage") |> 
  mutate(Tissue = factor(Tissue, levels= c("Phloem", "Leaf")))

# Create the bar graph
strand_bar <- ggplot(graph_table_normalized, aes(x = Strand, y = Percentage, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tissue",
       y = "Percentage",
       fill = NULL) +
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

ggsave(filename = "lncRNAs/lncRNAs_StrandDistribution.png",
       plot = strand_bar,
       width = 10, height = 8,
       dpi = 600)
saveRDS(strand_bar, file='lncRNAs/lncRNAs_StrandDistribution.RDS')
