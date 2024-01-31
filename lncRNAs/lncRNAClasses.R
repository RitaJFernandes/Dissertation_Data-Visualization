################################################################################
############## This script contains classification of lncRNAs ##################
################################################################################

pacman::p_load(readxl, dplyr, ggplot2, tidyr, wrappedtools, RColorBrewer, forcats,
               ggrepel, patchwork, pals, forcats)

# Loading files
leaf_lncRNA_class <- read.table("PATH_HERE/lncRNA_classification.txt",
                                 header = TRUE, 
                                 sep = "\t")
phloem_lncRNA_class <- read.table("PATH_HERE/lncRNA_classification.txt",
                                header = TRUE, 
                                sep = "\t")

# Getting percentages and casting classes into factor
leaf_lncRNA_class <- leaf_lncRNA_class |> 
  mutate(Percentage = round((count / sum(count)) * 100, digits=2),
         lncRNA_class = as_factor(lncRNA_class))

phloem_lncRNA_class <- phloem_lncRNA_class |> 
  mutate(Percentage = round((count / sum(count)) * 100, digits=2),
         lncRNA_class = as_factor(lncRNA_class))

# Pie charts
leaf_plot <- ggplot(leaf_lncRNA_class,
                    aes(x="", y=Percentage,
                        fill=(paste0(lncRNA_class, " (", Percentage, "%)")))) +
  geom_col() +
  coord_polar("y")+
  scale_fill_manual(values=as.vector(cols25(11)))+
  guides(fill=guide_legend(title="lncRNA classes"))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 20, b = 20),
        text=element_text(size=24))

ggsave(filename = "lncRNAs/Leaf_lncRNAs_Classes.png",
       plot = leaf_plot,
       width = 10, height = 8,
       dpi = 600)
saveRDS(leaf_plot, file='lncRNAs/Leaf_lncRNAs_Classes.RDS')




phloem_plot <- ggplot(phloem_lncRNA_class,
                    aes(x="", y=Percentage,
                        fill=(paste0(lncRNA_class, " (", Percentage, "%)")))) +
  geom_col() +
  coord_polar("y")+
  scale_fill_manual(values=as.vector(cols25(11)))+
  guides(fill=guide_legend(title="lncRNA classes"))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 20, b = 20),
        text=element_text(size=24))

ggsave(filename = "lncRNAs/Phloem_lncRNAs_Classes.png",
       plot = phloem_plot,
       width = 10, height = 8,
       dpi = 600)
saveRDS(phloem_plot, file='lncRNAs/Phloem_lncRNAs_Classes.RDS')
