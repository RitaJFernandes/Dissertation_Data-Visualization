################################################################################
################ This script contains the sizes of lncRNAs #####################
################################################################################
pacman::p_load(readxl, dplyr, tidyr, ggplot2, rtracklayer, scales)

################################### Leaf #######################################
# Getting sizes
leaf <- readGFF("PATH_HERE/novel_lncRNA.gtf")
leaf<- leaf |> 
  select(seqid, start, end, transcript_id, gene_id)

## By transcript (I used this one in the thesis)
leaf_size_by_transcript <- leaf |> 
  group_by(transcript_id) |>
  summarise(leaf_total_exon_length = sum(end - start +1),
            leaf_num_exons = n()) #1142 transcripts

temp1 <- leaf_size_by_transcript |> 
  rename("leaf_total_exon_length"="Size")

temp1$grouped_values <- ifelse(temp1$Size <= 2000, temp1$Size, 2000 + 1)

leaf_size_plot <- ggplot(temp1, aes(x = grouped_values)) +
  geom_histogram(binwidth = 100, fill = "#F76C6C", color = "black") +
  scale_x_continuous(breaks = c(seq(200, 1800, by = 200), 2001), 
                     labels = c(seq(200, 1800, by = 200), ">2000"),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "Size", y = "Count")+
  theme_bw()+
  theme(axis.text.x = element_text(size=35, angle = 45, hjust = 1),
        axis.text.y = element_text(size=35),
        axis.title = element_text(size= 40),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        aspect.ratio=1/1)

ggsave(filename = "lncRNAs/Leaf_lncRNAs_SizeDistribution.png",
       plot = leaf_size_plot,
       width = 10, height = 8,
       dpi = 600)

################################## Phloem #######################################
# Getting sizes
phloem <- readGFF("PATH_HERE/novel_lncRNA.gtf")
phloem <- phloem |> 
  select(seqid, start, end, transcript_id, gene_id)

phloem_size_by_transcript <- phloem |> 
  group_by(transcript_id) |>
  summarise(phloem_total_exon_length = sum(end - start +1),
            phloem_num_exons = n()) #2326 transcripts

temp2 <- phloem_size_by_transcript |> 
  rename("phloem_total_exon_length"="Size")

temp2$grouped_values <- ifelse(temp2$Size <= 2000, temp2$Size, 2000 + 1)

phloem_size_plot <- ggplot(temp2, aes(x = grouped_values)) +
   geom_histogram(binwidth = 100, fill = "#23305E", color = "black") +
   scale_x_continuous(breaks = c(seq(200, 1800, by = 200), 2001), 
                      labels = c(seq(200, 1800, by = 200), ">2000"),
                      expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "Size", y = "Count")+
  theme_bw()+
  theme(axis.text.x = element_text(size=35, angle = 45, hjust = 1),
        axis.text.y = element_text(size=35),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        axis.title = element_text(size= 40),
        aspect.ratio=1/1)

ggsave(filename = "lncRNAs/Phloem_lncRNAs_SizeDistribution.png",
       plot = phloem_size_plot,
       width = 10, height = 8,
       dpi = 600)


###################### Combining results in one table ##########################
# empty_sizes <- tibble(leaf_total_exon_length = NA_real_,
#                       .rows = 1184) # difference of rows between phloem and leaf size tables
# transcripts_length <- leaf_size_by_transcript |> 
#   select(leaf_total_exon_length) |> 
#   rbind(empty_sizes) |> 
#   rename(leaf_total_exon_length = "Leaf") |> 
#   cbind(phloem_size_by_transcript[2]) |>
#   rename(phloem_total_exon_length = "Phloem") |> 
#   pivot_longer(cols= c(Leaf, Phloem),
#                values_to = "Size" ,
#                names_to = "Tissue") |> 
#   mutate(Tissue = factor(Tissue, levels= c("Phloem", "Leaf"))) |> 
#   filter(!is.na(Size))
# ggsave(filename = "lncRNAs/lncRNAs_SizeDistribution.png",
#        plot = density_plot,
#        width = 10, height = 8,
#        dpi = 600)


############################### Combined Plot ##################################
# size_combined <- ggplot(data = transcripts_length, aes(x = Size, fill = Tissue)) +
#   geom_histogram(position = "dodge")+
#   theme_bw()+
#   theme(text = element_text(size=40),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "black", linewidth = 2),
#         axis.ticks.length = unit(0.3, "cm"),
#         axis.ticks = element_line(linewidth = 2),
#         legend.key.size = unit(2, 'cm'),
#         aspect.ratio=1/1)+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_x_continuous(limits= c(200, 5000), expand = c(0, 0))+
#   scale_fill_manual(values=c("#23305E",  "#F76C6C"))
# 
# ggsave(filename = "lncRNAs/lncRNAs_SizeDistribution.png",
#        plot = density_plot,
#        width = 10, height = 8,
#        dpi = 600)
# saveRDS(density_plot, file='lncRNAs/lncRNAs_SizeDistribution.RDS')

# Bar plot
# ggplot(data = phloem_sizes, aes(x = Size)) + 
#   geom_histogram() + 
#   xlim(0, 5000) +
#   #scale_x_continuous(breaks = seq(0, 2000, by = 100)) +
#   stat_bin(aes(y=..count.., label=..count..), geom="text", angle = 90, vjust = -.5)+
#   scale_fill_brewer(palette="Set1")


############################# Sanity check #####################################
# ggplot()+
#   geom_density()+
#   geom_density(data=leaf_sizes, aes(x=Size), color='green')+
#   geom_density(data=phloem_sizes, aes(x=Size), color='red')+
#   xlim(0, 2000)
