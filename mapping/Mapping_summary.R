pacman::p_load(readxl, writexl, dplyr, tidyr, forcats, flextable)

leaf_summary <- read_excel("PATH_HERE/mapping_summary.xlsx") |> 
  rename(Leaf1=Input_1,
         Leaf2=Input_2,
         Leaf3=Input_3)
phloem_summary <- read_excel("PATH_HERE/mapping_summary.xlsx")

mapping_sum_combined <- merge(phloem_summary, leaf_summary, sort=F) 
write_xlsx(mapping_sum_combined, "mapping/mapping_summary.xlsx")

  flextable(leaf_summary) |>
  fontsize(size = 12, part = "header") |>
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") |>
  fontsize(size = 11, part = "body") |> 
  align_text_col(align = "center", header = T, footer = T) |>
  align_nottext_col(align = "center") |>
  autofit() |>
  save_as_image(path = "mapping/Leaf_mapping_summary.png")
  
  flextable(phloem_summary) |>
    fontsize(size = 12, part = "header") |>
    bold(i = NULL, j = NULL, bold = TRUE, part = "header") |>
    fontsize(size = 11, part = "body") |> 
    align_text_col(align = "center", header = T, footer = T) |>
    align_nottext_col(align = "center") |>
    autofit() |>
    save_as_image(path = "mapping/Phloem_mapping_summary.png")

  