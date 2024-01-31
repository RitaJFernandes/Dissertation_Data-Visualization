################################################################################
########## This script contains the coding potential filter results ############
################################################################################
pacman::p_load(readxl, dplyr, tidyr, ggvenn, patchwork)

# Coding Potential Files
leafCP <- read_excel("PATH_HERE/Coding_potential.xlsx")
phloemCP <- read_excel("PATH_HERE/Coding_potential.xlsx")

# Removing duplicates
leafCP <- distinct(leafCP)
phloemCP <- distinct(phloemCP)

# Getting non-coding RNAS
## Leaf
### To confirm numbers
# leaf_lncRNAs <- leafCP |> filter(CNCI=='noncoding' & CPC=='noncoding' & PFAM=='noncoding')
# leaf_lncRNAs_TCONS <- c(distinct(leaf_lncRNAs, transcript_id))
# CNCI_coding <- rep("CNCI", times=sum(leafCP$CNCI=="coding"))
# CPC_coding <- rep("CPC", times=sum(leafCP$CPC=="coding"))
# PFAM_coding <- rep("PFAM", times=sum(leafCP$PFAM=="coding"))
# CNCI_noncoding <- rep("CNCI", times=sum(leafCP$CNCI=="noncoding"))
# CPC_noncoding <- rep("CPC", times=sum(leafCP$CPC=="noncoding"))
# PFAM_noncoding <- rep("PFAM", times=sum(leafCP$PFAM=="noncoding"))

# Checking if all lncRNAs IDs are present only in inputs
# leaf_inputsID <- read_excel("PATH_HERE/transcripts.FPKM.xlsx") |> 
#   select(transcript_id)
# check <- merge(leafCP, leaf_inputsID)

sets_leaf <- leafCP[-2:-5] |>
  mutate(CNCI= case_when(CNCI %in% c("noncoding") ~ TRUE,
                         CNCI %in% c("coding") ~ FALSE),
         CPC= case_when(CPC %in% c("noncoding") ~ TRUE,
                        CPC %in% c("coding") ~ FALSE),
         PFAM= case_when(PFAM %in% c("noncoding") ~ TRUE,
                         PFAM %in% c("coding") ~ FALSE))

LeafVennDiagramCP <- ggvenn(sets_leaf,
                            c("CNCI", "CPC", "PFAM"),
                            show_percentage = FALSE,
                            text_size = 8,
                            set_name_size = 8)

ggsave(filename = "lncRNAs/LeafCP_VennDiagram.png",
       plot = LeafVennDiagramCP,
       width = 10, height = 8,
       dpi = 600)
saveRDS(LeafVennDiagramCP, file = "lncRNAs/LeafCP_VennDiagram.RDS")
        
## Phloem
# phloem_lncRNAs <- phloemCP |> filter(CNCI=='noncoding' & CPC=='noncoding' & PFAM=='noncoding')
# phloem_lncRNAs_TCONS <- c(distinct(phloem_lncRNAs, transcript_id))

sets_phloem <- phloemCP[-2:-5] |> 
  mutate(CNCI= case_when(CNCI %in% c("noncoding") ~ TRUE,
                         CNCI %in% c("coding") ~ FALSE),
         CPC= case_when(CPC %in% c("noncoding") ~ TRUE,
                        CPC %in% c("coding") ~ FALSE),
         PFAM= case_when(PFAM %in% c("noncoding") ~ TRUE,
                         PFAM %in% c("coding") ~ FALSE))

PhloemVennDiagramCP <- ggvenn(sets_phloem,
                              c("CNCI", "CPC", "PFAM"),
                              show_percentage = FALSE,
                              text_size = 8,
                              set_name_size = 8)

ggsave(filename = "lncRNAs/PhloemCP_VennDiagram.png",
       plot = PhloemVennDiagramCP,
       width = 10, height = 8,
       dpi = 600)
saveRDS(PhloemVennDiagramCP, file = "lncRNAs/PhloemCP_VennDiagram.RDS")
