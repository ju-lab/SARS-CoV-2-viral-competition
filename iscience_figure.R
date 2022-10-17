################################################################################
# this script is used for analysis and drawing figures in Kim et al., (2022)
# Relative infectivity of the SARS-CoV-2 Omicron variant in human alveolar cells Biorxiv 
################################################################################

################################################################################
# Required source files 
################################################################################
source("iscience_functions.R") 
library(data.table)
library(readr)
library(RColorBrewer)
library(VennDiagram)
library(ggplot2)
library(gggenes)
library(ggpubr)
library(lsa)
library(gridExtra)
library(cowplot)
library(EMT)
library(lessR)
library(scales)

################################################################################
# color code 
################################################################################
colorset_final <- brewer.pal(9,"Set1")[c(2,4,3,1,7)] # blue, purple, green, red, brown
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# unique mutations  
################################################################################
vocs <- c("gr","alpha","delta","omicron")

unique_mutations <- "unique_mutations.csv"
site_short_vocs <- read_csv(unique_mutations)

################################################################################
# Figure S1A : clonal and subclonal mutations  
################################################################################
# read annotated files for 
heatmapfile <- "heatmap_vaf.tsv"

dt_heat2_tmp2 <- read_csv(heatmapfile)

dt_heat2_tmp2 <- dt_heat2_tmp2 %>% dplyr::distinct(.,CHROM, POS, REF, ALT,.keep_all = T) %>% arrange(POS)  %>% tidyr::unite("mut",CHROM:ALT) 
rownames(dt_heat2_tmp2) <- dt_heat2_tmp2[["mut"]] # name of mutation's position and REF ALT
dt_heat2<- data.matrix(dt_heat2_tmp2, rownames.force = T)
dt_heat2 <- dt_heat2[,-1]

ha = HeatmapAnnotation(
  viral_variant = anno_block(gp = gpar(fill = colorset_final[c(3,1,2,4)]), labels = c("gr","alpha","delta","omicron"), labels_gp = gpar(col = "white", fontsize = 1)), 
  type = rep(c(rep("batches",12), "gisaid", "unique"), 4), 
  col = list(type = c("batches" = colorBlindBlack8[3], "gisaid" = colorBlindBlack8[5], "unique" = colorBlindBlack8[7])),
  border = T, show_legend = c("type" = FALSE))

split = rep(1:4, each =14 )

col_fun = colorRamp2(c(0, 1), c("white", "blue"))

Heatmap(dt_heat2,
        na_col = "grey",  col = col_fun, 
        cluster_rows = F, cluster_columns = F, 
        show_row_names = T, show_column_names = F, show_heatmap_legend = F,
        column_split = split, 
        rect_gp = gpar(col = "grey", lwd = 2),
        top_annotation = ha)


############################################################################################################################################
# Figure 1B : unique mutations by venn diagram (it is not possible to draw with phylogenetic tree)
############################################################################################################################################
grid.newpage()
draw.quad.venn( # shown as 1 3 4 2
  area1 = site_short_vocs %>% dplyr::filter(!is.na(gr)) %>% nrow(.), 
  area3 = site_short_vocs %>% dplyr::filter(!is.na(alpha)) %>% nrow(.),
  area4 = site_short_vocs %>% dplyr::filter(!is.na(delta)) %>% nrow(.),
  area2 = site_short_vocs %>% dplyr::filter(!is.na(omicron)) %>% nrow(.), 
  n13 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(alpha) & !is.na(gr)) %>% nrow(.),
  n14 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(delta) & !is.na(gr)) %>% nrow(.),
  n12 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(gr) & !is.na(omicron)) %>% nrow(.),
  n34 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(delta) & !is.na(alpha)) %>% nrow(.),
  n23 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(alpha) & !is.na(omicron)) %>% nrow(.),
  n24 = site_short_vocs %>% dplyr::filter(vocs >= 2 & !is.na(delta) & !is.na(omicron)) %>% nrow(.),
  n134 = site_short_vocs %>% dplyr::filter(vocs >= 3 & !is.na(alpha) & !is.na(delta) & !is.na(gr)) %>% nrow(.),
  n123 = site_short_vocs %>% dplyr::filter(vocs >= 3 & !is.na(alpha) & !is.na(gr) & !is.na(omicron)) %>% nrow(.), 
  n124 = site_short_vocs %>% dplyr::filter(vocs >= 3 & !is.na(delta) & !is.na(gr) & !is.na(omicron)) %>% nrow(.),
  n234 = site_short_vocs %>% dplyr::filter(vocs >= 3  & !is.na(delta) & !is.na(alpha) & !is.na(omicron)) %>% nrow(.),
  n1234 = site_short_vocs %>% dplyr::filter(vocs == 4) %>% nrow(),
  category = c('Alpha', 'Omicron','Delta', 'GR'),
  fill = colorset_final[c(3,4,1,2)],
  lty = 'solid',
  lwd = rep(1,4),
  alpha = rep(0.6,4),
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica",
  cex = rep(1,15),
  cat.pos = c(-7,7,0,0)
)



############################################################################################################################################
# Figure 1C : expression levels of human alveolar type 2 cells 
############################################################################################################################################
exp_csv <- "20221004_fig_dat_per10k.tsv"

umi_dat <- read_tsv(exp_csv)  %>% 
  gather('SFTPB','NKX2-1','ACE2','TMPRSS2', key = 'gene', value = 'expression') %>%
  dplyr::mutate(gene = factor(gene, levels = c('SFTPB','NKX2-1','ACE2','TMPRSS2')))

fig1c_exp <- ggplot(umi_dat) +  
  geom_violin(aes(x= gene, y = expression + 1), scale= "width") +
  geom_jitter(aes(x= gene, y = expression + 1), width = 0.25, size = 1.5) +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  theme_pubr() + 
  theme( axis.text.x = element_text(family = "Helvetica",size = 15), axis.title.x = element_blank(), 
         axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20)) 

############################################################################################################################################
# Metadata : Infectivity information for all cells 
############################################################################################################################################
metadata <- "infected_cells_final_revision.csv"
norm_cell_voc <- read_csv(metadata)

### fix sample order
norm_cell_voc_whole <- gather(data = norm_cell_voc, 1:5, key = "vocs", value = "fraction")
norm_cell_voc_whole$vocs <- factor(norm_cell_voc_whole$vocs, levels =c("gr_assume","alpha_assume","delta_assume","omicron_assume","unknown_assume"))
norm_cell_voc_whole$sample_name <- factor(norm_cell_voc_whole$sample_name, levels = sample_name_order)

infected_norm_cell_voc <- norm_cell_voc_whole %>% dplyr::filter( !QC_fail & Infected & MOI == 2.5  ) %>% 
  dplyr::mutate(fraction = 100 * fraction, 
                virus_batch = factor(virus_batch, levels= c("O","A","B","C","D","E","F","G","H","I","J","K")),
                time_incubation_min = factor(time_incubation_min, levels = c(5,60)),
                vocs = ifelse(vocs == "gr_assume", "GR", 
                              ifelse(vocs == "alpha_assume", "Alpha", 
                                     ifelse(vocs == "delta_assume", "Delta", 
                                            ifelse(vocs == "omicron_assume", "Omicron",
                                                   ifelse(vocs == "unknown_assume", "Unknown", vocs)))))) 

# the number of cells passed infected criteria
ncells <- infected_norm_cell_voc %>% nrow(.) /5 # 191 cells  -> 244 cells (additional 53 cells)

# value of viral reads per cell
infected_norm_cell_voc %>% dplyr::pull(viral_percent_bwa) %>% min(., na.rm = T) # 0.00085
infected_norm_cell_voc %>% dplyr::pull(viral_percent_bwa) %>% max(., na.rm = T) # 60.70

############################################################################################################################################
## Figure 1D: viral and cellular reads information
############################################################################################################################################
## adding subplot over the barplots using gradient barplots 
viral_reads_dat <- gather(infected_norm_cell_voc %>% 
                            dplyr::mutate(percent_viral_bwa = total_viral_reads_bwa / total_reads_bwa * 100 ) %>% 
                            distinct(sample_name, .keep_all = T), "total_reads_bwa","total_viral_reads_bwa", key = "reads_type",value = "reads") %>%  
  arrange(sample_name) %>% dplyr::mutate(cell_number = rep(1:ncells, each = 2)) 

fig1_viral_reads_updated <- ggplot(viral_reads_dat, aes(x = cell_number, y = reads, fill = reads_type)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), 
                limits = c(1, 100000000)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values =colorBlindBlack8[c(6,2)]) + 
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20),
        legend.text =  element_text(family = "Helvetica",size = 15) , legend.title = element_text(family = "Helvetica",size =15), legend.justification = "center",
        panel.background = element_blank(), legend.position = "top") + 
  facet_grid( ~ boundaries, scales = "free_x") +
  theme(panel.spacing = unit(0.05, "cm"))


fig1_viral_heatmap <- ggplot(viral_reads_dat, aes(x = cell_number, y = 1, fill = percent_viral_bwa )) + geom_bar(stat ="identity") +
  scale_fill_gradient(name = "viral_read_percent(%)", trans = "log",
                      breaks = trans_breaks("log10", function(x) 10^x), 
                      labels = trans_format("log10", math_format(10^.x)),
                      limits = c(0.0001,100),low = "white",high ="brown") +
  theme(legend.position = "top", axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_blank(),
        legend.text =  element_text(family = "Helvetica",size = 15) , legend.title = element_text(family = "Helvetica",size =15),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(2, 'cm')) + 
  facet_grid( ~ boundaries, scales = "free_x") +
  theme(panel.spacing = unit(0.05, "cm"))

legend <- get_legend(fig1_viral_heatmap)   # get the legend of the first one plot

prow <- plot_grid( fig1_viral_reads_updated + theme(legend.position = "none", axis.title.y = element_blank()),
                   fig1_viral_heatmap + theme(legend.position="none"),
                   align = 'hv',
                   nrow = 2, rel_heights = c( 11, 1.5))
p <- plot_grid( prow, legend, rel_widths = c(3, .3))


############################################################################################################################################
## Figure S1C :  virus fraction by plaque assay and RNA sequencing  
############################################################################################################################################
virus_assay_dt <- "live_virus_assay.xlsx"
live_virus <- read_xlsx(path = virus_assay_dt) %>%
  dplyr::mutate(virus_batch = factor(virus_batch, levels = c("K","J","I","H")),
                variant = factor(variant, levels = c("GR","Alpha","Delta","Omicron")))

virus_fraction_by_method <- ggplot(live_virus , aes(x = virus_batch, y = fraction, fill = variant)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colorset_final[c(3,1,2,4)]) +
  facet_wrap( ~ method) + coord_flip() + 
  theme( axis.text.x = element_blank(),axis.title.x = element_blank(),
         axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20),
         legend.text =  element_text(family = "Helvetica",size = 15) , legend.title = element_text(family = "Helvetica",size =15), legend.justification = "center",
         panel.background = element_blank(), legend.position = "none",
         strip.text.x = element_text(size = 15), strip.text=element_text(vjust= 0), strip.switch.pad.grid = unit('1', "cm"), strip.placement = "outside")

############################################################################################################################################
## Figure 1E: average coverage of 10X and SMART-seq3
############################################################################################################################################
# draw genes in the ggplot
gtffile <-"ncbiGenes_non_nestedexon.gtf"
gtf <- as.data.frame(rtracklayer::import(gtffile)) %>% dplyr::mutate(molecule = seqnames, 
                                                                     gene_id = factor(gene_id, levels = c("ORF1a","ORF1ab","S","ORF3a",
                                                                                                          "E","M","ORF6",
                                                                                                          "ORF7a","ORF7b","ORF8",
                                                                                                          "N","ORF10"))) 

alpha_gene <- gene_draw(site_short_vocs, "alpha",1)
delta_gene <- gene_draw(site_short_vocs, "delta",2)
gr_gene <- gene_draw(site_short_vocs, "gr",3)
omicron_gene <- gene_draw(site_short_vocs, "omicron",4)
fig1e <- ggarrange(gr_gene, alpha_gene,delta_gene,omicron_gene, main_viral_depth,  ncol = 1, align = "v", heights = c(1,1,1,1,10))

infected_cells_viral_read <- infected_norm_cell_voc %>% distinct(sample_name, .keep_all = T) #%>%   dplyr::filter(total_viral_reads_bwa < 9893 & total_viral_reads_bwa > 33 ) 
samtools_depth_smart <- norm_avg_cell_depth(data_frame = infected_cells_viral_read)

# draw coverage of 10X and SMART-seq3 
depth_file <- "samtools_depth_final.csv"
samtools_depth_final <- read_csv(depth_file, col_names = c('CHR','POS','depth', 'method'))

main_viral_depth <- ggplot(samtools_depth_final, aes(x = POS , y = depth)) + 
  geom_line(aes(color = method)) + 
  theme_pubclean() + 
  theme( axis.text.x = element_text(family = "Helvetica",size = 15), axis.title.x = element_blank(), 
         axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20)) +
  labs(fill = "Viruses",  y = "average coverage") +   
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-7,1e-2))



############################################################################################################################################
## Figure 1F-G: Cosine similarity between two methods 
############################################################################################################################################
# since infected_norm_cell_voc gathers each value of fraction from viruses . So I use norm_cell_voc which already filters QC and viral infection
cos_info <- lapply(1:nrow(norm_cell_voc), function(x){cosine_sim_func(norm_cell_voc %>% dplyr::slice(x))}) %>% rbindlist() #%>% do.call("rbind")
cossim <- full_join(norm_cell_voc, cos_info,by = "cell")

# Figure 1F
cossim_plot <- ggpar(plot, xlim = c(0,1), yticks.by = 50, xticks.by = 0.1, ylim = c(0,250))

# Figure 1G
cossim_plot_sub <- ggplot(cossim) + 
  geom_point(aes(x = nmf_unique_sites, y = na2zero , color = as.factor(voc_stat)), size = 4) + xlim(0,120) + ylim(0,1) +
  scale_color_brewer(palette = "Greys") +
  theme_pubr() + 
  theme( axis.text.x = element_text(family = "Helvetica",size = 15), axis.title.x = element_blank(), 
         axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20)) 



############################################################################################################################################
## Figure 2A: All infected cells' decomposed results 
############################################################################################################################################
sub_a <-  ggplot(infected_norm_cell_voc %>% arrange(sample_name) %>%  dplyr::mutate(cell_number = rep(1:ncells, each= 5)) , aes(x = cell_number, y = fraction, fill = vocs)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = colorset_final) +
  theme(legend.text =  element_text(family = "Helvetica",size = 15) ,
        panel.background = element_blank(), legend.position = "right", plot.margin = unit(c(-0.35,1,1,1), units = "cm")) +
  labs(fill = "Viruses", x = "Cell", y = "Virus Proportion (%)") + 
  facet_grid( ~ boundaries, scales = "free_x") +
  theme(panel.spacing = unit(0.05, "cm"))

sub_b <- ggplot(infected_norm_cell_voc %>% arrange(sample_name) %>%  dplyr::mutate(cell_number = rep(1:ncells, each= 5)), aes(x=cell_number, y= 1, fill = time_post_infection_min)) + 
  geom_bar(stat = "identity") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_brewer(palette = "Oranges") +
  facet_grid( ~ boundaries, scales = "free_x") +
  theme(panel.spacing = unit(0.05, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

sub_c <- ggplot(infected_norm_cell_voc %>% arrange(sample_name) %>%  dplyr::mutate(cell_number = rep(1:ncells, each= 5)), aes(x=cell_number, y= 1, fill = time_incubation_min)) + geom_tile(stat = "identity") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_fill_brewer(palette= "Greys")+
  facet_grid( ~ boundaries, scales = "free_x") +
  theme(panel.spacing = unit(0.05, "cm"))

legend <- get_legend(sub_a)   # get the legend of the first one plot

prow <- plot_grid( sub_a+ theme(legend.position="none", axis.title.y = element_blank()),
                   sub_b + theme(legend.position="none"),
                   sub_c + theme(legend.position = "none"),
                   align = 'v',
                   nrow = 3, rel_heights = c( 18, 1, 1))

p <- plot_grid( prow, legend, rel_widths = c(3, .3))

############################################################################################################################################
## Figure 2B : infectivity across experimental batches 
############################################################################################################################################
fig2b_dat <- infected_norm_cell_voc %>%
  dplyr::filter(time_incubation_min == 5 ) %>% #& virus_batch != "A") %>% 
  arrange(sample_name) %>% dplyr::mutate(cell_number = rep(1: (nrow(.)/5), each = 5))

fig2b <- ggplot(fig2b_dat, aes(x = sample_name, y = fraction, fill = vocs)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorset_final) +
  facet_grid( ~ virus_batch, scales = 'free', space = 'free') +
  theme( axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
         axis.text.y = element_text(family = "Helvetica",size = 15), axis.title.y = element_text(family = "Helvetica",size = 20),
         legend.text =  element_text(family = "Helvetica",size = 15) , legend.title = element_text(family = "Helvetica",size =15), legend.justification = "center",
         panel.background = element_blank(), legend.position = "none",
         strip.text.x = element_text(size = 15), strip.text=element_text(vjust= 0), strip.switch.pad.grid = unit('1', "cm"), strip.placement = "outside") 


fig2b_sub_b <- ggplot(fig2b_dat, aes(x = sample_name, y= 1, fill = time_post_infection_min)) + 
  geom_bar(stat = "identity") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_fill_brewer(palette = "Oranges") +
  facet_grid( ~ virus_batch, scales = "free", space = 'free') 

fig2b_sub_c <- ggplot(fig2b_dat, aes(x =sample_name, y= 1, fill = time_incubation_min)) + geom_tile(stat = "identity") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_fill_brewer(palette= "Greys")+
  facet_grid( ~ virus_batch, scales = "free", space = 'free') 

legend <- get_legend(fig2b)   # get the legend of the first one plot

prow <- plot_grid( fig2b + theme(legend.position="none", axis.title.y = element_blank()),
                   fig2b_sub_b + theme(legend.position="none"),
                   fig2b_sub_c + theme(legend.position = "none"),
                   align = 'v',
                   nrow = 3, rel_heights = c( 18, 1, 1))

############################################################################################################################################
## Figure 2C: infectivity over time course
############################################################################################################################################
phist <- gghistogram(
  infected_norm_cell_voc %>% dplyr::filter(time_post_infection_min  == 5  ) %>% dplyr::filter(vocs =="Omicron"),
  x = "fraction", fill = "time_incubation_min", 
  add = "mean", rug = FALSE, 
  palette = colorBlindBlack8[c(5,6)], 
  bins = 40 )

pdensity <- ggdensity(
  infected_norm_cell_voc %>% dplyr::filter(time_post_infection_min  == 5  ) %>% dplyr::filter(vocs =="Omicron"),
  x = "fraction", color = "time_incubation_min" ,
  palette = colorBlindBlack8[c(5,6)],
  alpha = 0
) + scale_y_continuous( position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")

# 3. Align the two plots and then overlay them.
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])


############################################################################################################################################
## Analysis : infectivity calculation
############################################################################################################################################
# the number of infected cells by single, double, and multiple infections
norm_cell_voc %>% group_by(voc_stat, unknown) %>% dplyr::summarize(n=n())

# value for how many cells were infected by single virus only
infected_norm_cell_voc %>% distinct(sample_name, .keep_all = T) %>% dplyr::filter(voc_stat == 1 & unknown == 0) %>% 
  group_by(gr_stat,alpha_stat,delta_stat, omicron_stat) %>% dplyr::summarize(n=n()) #%>% View(.)

# value for how many cells were infected by total 
infect<- infected_norm_cell_voc %>% distinct(sample_name, .keep_all = T) %>% dplyr::filter( unknown == 0) %>% 
  group_by(gr_stat,alpha_stat,delta_stat, omicron_stat) %>% dplyr::summarize(n=n()) 
infect %>% dplyr::filter(omicron_stat) %>% pull(n) %>% sum() 
infect %>% dplyr::filter(alpha_stat) %>% pull(n) %>% sum() 
infect %>% dplyr::filter(delta_stat) %>% pull(n) %>% sum() 
infect %>% dplyr::filter(gr_stat) %>% pull(n) %>% sum()

# calculate the area under VOCs by integral
integral_voc <- function(dat, voc){return(dat %>% dplyr::filter(vocs == voc) %>% dplyr::pull(fraction) %>% sum(., na.rm = T))}
integrals <- c("Alpha","Delta","GR","Omicron")
integral_list <- lapply(1:length(integrals), function(x){
  integral_voc(infected_norm_cell_voc %>%  dplyr::filter( unknown == 0) , integrals[x])})
names(integral_list) <- integrals

### calculate the ratio of omicron's dominance
# single infection cases only 
tot_single <- 89
omi_single <- 56
single_omi_test_single <- prop.test(omi_single, tot_single, p = 0.25)
single_omi_test_ratio <- single_omi_test_single$estimate/0.25
single_omi_test_conf <- single_omi_test_single$conf.int / 0.25
single_omi_test_pvalue <- single_omi_test_single$p.value

omi_single / (tot_single - omi_single) / (1/4)

# test all variant's proportion are not same 
library(lessR)
variants <-c(199,114,75,71)# c(161,106,70,48)
names(variants) <- c("omicron","alpha","delta","gr")
tot_variants <- rep(244,4) 
multinomial_multi <- lessR::Prop_test(n_succ = variants, n_tot = tot_variants)
multinomial_multi$p.value
#https://rpubs.com/raysunau/compare_proportions

# multi infection cases by considering proportion in each cell
tot_multi <- 244
omi_multi <- 142.8
tot_omi_test_single <- prop.test(omi_multi, tot_multi, p = 0.25)
tot_omi_test_ratio <- tot_omi_test_single$estimate/0.25
tot_omi_test_conf <- tot_omi_test_single$conf.int / 0.25
tot_omi_test_pvalue <- tot_omi_test_single$p.value


## odds ratio 
(omi_multi/(244-omi_multi)) / (1/4) # for average infection
(omi_multi/(244-omi_multi)) / (55.9 / (244 - 55.9)) # 4.8 for alpha
(omi_multi/(244-omi_multi)) / (11.7 / (244 -11.7)) # 28
(omi_multi/(244-omi_multi)) / (30.6 / (244 - 30.6)) # 10 for GR
