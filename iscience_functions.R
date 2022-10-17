################################################################################
# Functions used for analysis in the research paper of Kim et al., (2022)
# Relative infectivity of the SARS-CoV-2 Omicron variant in human alveolar cells Biorxiv 
################################################################################


############################################################################################################################################
## Figure 1E: average coverage of 10X and SMART-seq3
############################################################################################################################################
gene_draw <- function(site_short_vocs, voc, index){
  test <- site_short_vocs %>% dplyr::filter(!is.na(.data[[voc]])) %>% dplyr::mutate(POS2 = POS + nchar(ALT) - nchar(REF)) 
  color_values <- colorset_final
  print(voc)
  print(index)
  ggplot() + geom_gene_arrow(gtf %>% dplyr::filter(type == "transcript") %>% dplyr::mutate(seqnames = voc),
                             arrowhead_height = unit(1, "cm"), arrowhead_width = unit(0.13, "cm"), arrow_body_height = unit(1,"cm"),
                             mapping = aes(xmin = start, xmax = end, y = seqnames))  +
    geom_segment(data = test, aes(x = POS, xend = POS2, y = 0.5, yend = 1.5, color = "blue")) + scale_color_manual(values = color_values[index]) +
    theme_void() + theme(legend.position = "none") 
}
############################################################################################################################################
## Figure 1F-G: Cosine similarity between two methods 
############################################################################################################################################
extract_value <- function(line, name){ line %>% dplyr::pull(name)}
cosine_sim_func <- function(data_frame){
  nmf_colnames <- c("GR_ratio", "Alpha_ratio", "Delta_ratio", "Omicron_ratio")
  avg_vaf_colnames <- c("gr_assume", "alpha_assume", "delta_assume", "omicron_assume")
  nmf <- lapply(nmf_colnames, function(x,dat = data_frame){ extract_value(dat, x)}) %>% unlist()
  avg_vaf <- lapply(avg_vaf_colnames, function(x,dat = data_frame){ extract_value(dat, x)}) %>% unlist()
  na <- cosine(nmf, avg_vaf)
  avg_vaf[is.na(avg_vaf)] <- 0
  nazero <- cosine(nmf, avg_vaf)
  return(data.frame(cell = data_frame[["cell"]], na2zero= nazero, na_original = na))
}
