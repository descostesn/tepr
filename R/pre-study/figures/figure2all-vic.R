# Figure 2 final
# 2B Universe_metagene
# 2C AUC_HS_vs_control_pvalue
# 2D MARCHF6_plot
# 2E EGFR_plot
# 2F_left AUC_HS_vs_control_groups
# 2F_right_top Attenuated_metagene
# 2F_right_bottom Outgroup_metagene
# 2I histo_per
# 2J histo_kb


####################
## LOADING TABLES ##
####################

## Loading
tst_df <- read.delim("/g/romebioinfo/Projects/tepr/downloads/Cugusi2022_AttenuationScores_10_200.tsv", header=T, sep="\t")
working_directory <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/RAW_TTseq/deduplicated/MAPQ255/stranded/bedgraph255"
extension <- "*.bg"

PATH <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/RAW_TTseq/"
rounding <- 10
window_n <- 200
PaperYear <- "Cugusi2022"

mean_value_control_full <- "MeanValueFull_ctrl"
mean_value_stress <- "MeanValueFull_HS"
AUC_ctrl <- "AUC_ctrl"
AUC_stress <- "AUC_HS"
p_value_KStest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
p_value_theoritical<- "adjFDR_p_AUC_ctrl" 

# loading the summary table and applying threshold for attenuated vs non attenuated
tst_df <- tst_df %>%
  mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 & !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
                             !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%  relocate(Universe, .before = 1) 
tst_df <- tst_df %>%
  mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest))>2, "Attenuated", NA),
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group)  
  ) %>%  relocate(Group, .before = 2)

## Saving list of attenuated and outgroup transcripts for downstream analysis
write.table(tst_df, 
            file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/LOLA_analysis/Cugusi_analysis/20240818/All_transcripts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
write.table(tst_df %>% filter(Group=="Attenuated"), 
            file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/LOLA_analysis/Cugusi_analysis/20240818/Attenuation_HS_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
write.table(tst_df %>% filter(Group=="Outgroup"), 
            file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/LOLA_analysis/Cugusi_analysis/20240818/Outgroup_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
write.table(tst_df %>% filter(Universe==T), 
            file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/LOLA_analysis/Cugusi_analysis/20240818/Universe_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)

write.table(tst_df %>% filter(Universe==T) %>%  mutate(score_bed=".") %>% select(c(chr, coor1, coor2, transcript, score_bed, strand )), 
            file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/LOLA_analysis/Cugusi_analysis/20240818/Universe_transcripts.nosorted.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

# For metagenes
#concat_dfFX_res created from results_ECDF with calculates_meanFx function
concat_dfFX_res <- read.delim("/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/RAW_TTseq/Big_intermediary_file/concat_dfFX_res_Cugusi2rep.tsv", header=T, sep="\t")



##################
## SCATTER PLOT ##
##################

## Fig 2C. AUC in HS vs AUC in control with pValue
AUC_HS_vs_control_pvalue <- ggplot(tst_df %>% arrange(-log10(!!sym(p_value_KStest))), aes( AUC_ctrl, AUC_HS, color= -log10(!!sym(p_value_KStest)))) +
  geom_point(size=0.5) +
  geom_density_2d()+
  geom_label_repel(data = subset(tst_df, gene %in% c("EGFR","DAP","FLI1","MARCHF6", "LINC01619")), aes(label = gene),
   box.padding   = 0.55,
   point.padding = 0,
   segment.color = 'black', max.overlaps = 50, color="red") +
  scale_color_gradient2(midpoint=0,  low="white", mid="grey", high = "darkgreen") +
  xlim(-10,100) + ylim(-10,100)+
  labs(x="AUC in Control", y="AUC in HS", legend="-log10 p-value", color="-log10 p-value") +
  coord_fixed(ratio = 1) +   # Set aspect ratio to 1:1
  theme_classic() +
  theme(legend.position = "bottom" )

## Fig 2F Groups of genes
# palette: "#264653"  "#2a9d8f" "#e9c46a" "#f4a261" "#e76f51"
AUC_HS_vs_control_groups <- ggplot(tst_df %>% filter(Universe==F), aes( AUC_ctrl, AUC_HS) ) +
  geom_point(size=0.5, color="grey") +
  geom_point(data=tst_df %>% filter(Group=="Attenuated"), aes( x=AUC_ctrl , y=AUC_HS), colour="black" , size=1.3) +
  geom_point(data=tst_df %>% filter(Group=="Attenuated"), aes( x=AUC_ctrl , y=AUC_HS), color="#e76f51" , size=1) +
  geom_point(data=tst_df %>% filter(Group=="Outgroup"), aes( x=AUC_ctrl , y=AUC_HS), colour="black" , size=1.3) +
  geom_point(data=tst_df %>% filter(Group=="Outgroup"), aes( x=AUC_ctrl , y=AUC_HS), color="#e9c46a" , size=1) +
  xlim(-10,100) + ylim(-10,100)+
  labs(x="AUC in Control", y="AUC in HS", legend="-log10 p-value", title = "AUC Control vs HS", subtitle = "Genes selected for Unibind") +
  coord_fixed(ratio = 1) +   # Set aspect ratio to 1:1
  theme_classic() +
  theme(legend.position = "none" )

###############
## METAGENES ##
###############

#This work also for single genes
Attenuation_list<- tst_df %>% filter(Group=="Attenuated") %>% pull(transcript)
Outgroup_list<- tst_df %>% filter(Group=="Outgroup") %>% pull(transcript)
Universe_list<- tst_df %>% filter(Universe==T) %>% pull(transcript)
All_list<- tst_df  %>% pull(transcript)

AUC_allcondi <- tst_df %>% select(transcript, gene, strand, contains("Full"), dAUC_Diff_meanFx_HS_ctrl, AUC_ctrl, AUC_HS, -contains(c("UP","DOWN")), window_size)

normalize_and_summarize <- function(list) {
  result <- concat_dfFX_res %>%
    filter(transcript %in% list) %>%
    left_join(., AUC_allcondi, by = c("transcript", "gene")) %>%
    # mutate(
    #  across(contains("mean_value"), ~ . / mean_value_control_full, .names = "Normal_{.col}")
    # ) %>%
    select(transcript, gene, coord, contains("mean_value"), -contains("Full"))  %>%
    dplyr::group_by(coord) %>%
    #   summarise(across(contains("Normal"), ~ mean(., na.rm = TRUE))) %>%
    
    summarise(across(contains("mean_value"), ~ mean(., na.rm = TRUE))) # %>%
  
  return(result)
}

# Usage
Attenuation_concat_df <- normalize_and_summarize(Attenuation_list)
Outgroup_concat_df <- normalize_and_summarize(Outgroup_list)
Universe_concat_df <- normalize_and_summarize(Universe_list)
ALL_concat_df <- normalize_and_summarize(All_list)


## plotting:
y="Transcription density"
Attenuated_metagene <-  ggplot() +
  geom_line(data=Attenuation_concat_df, aes(x=coord/2, y= mean_value_ctrl), color="#00AFBB", size=1.5)+
  geom_line(data=Attenuation_concat_df, aes(x=coord/2, y= mean_value_HS), color="#FC4E07", size=1.5)+
  theme_bw() +
  ylim(0,7)+
  #scale_x_continuous(limits = c(0, 200), breaks = seq(0, 100, by = 25)) +
  labs(x="TSS to TTS")+
  labs(title = "Attenuated genes", subtitle = length(Attenuation_list), y=y) +
  theme(legend.position = "none", legend.box = "vertical")

Outgroup_metagene <-   ggplot() +
  geom_line(data=Outgroup_concat_df, aes(x=coord/2, y= mean_value_ctrl), color="#00AFBB", size=1.5)+
  geom_line(data=Outgroup_concat_df, aes(x=coord/2, y= mean_value_HS), color="#FC4E07", size=1.5)+
  theme_bw() +
  ylim(0,7)+
  labs(x="TSS to TTS")+
  labs(title = "Outgroup genes", subtitle = length(Outgroup_list), y=y) +
  theme(legend.position = "none", legend.box = "vertical")

Universe_metagene <-   ggplot() +
  geom_line(data=Universe_concat_df, aes(x=coord/2, y= mean_value_ctrl), color="#00AFBB", size=1.5)+
  geom_line(data=Universe_concat_df, aes(x=coord/2, y= mean_value_HS), color="#FC4E07", size=1.5)+
  theme_bw() +
  ylim(0,7)+
  labs(x="TSS to TTS")+
  labs(title = "Universe genes", subtitle = length(Universe_list), y=y) +
  theme(legend.position = "none", legend.box = "vertical")

All_metagene <-   ggplot() +
  geom_line(data=ALL_concat_df, aes(x=coord/2, y= mean_value_ctrl), color="#00AFBB", size=1.5)+
  geom_line(data=ALL_concat_df, aes(x=coord/2, y= mean_value_HS), color="#FC4E07", size=1.5)+
  theme_bw() +
  ylim(0,7)+
  labs(x="TSS to TTS")+
  labs(title = "Universe genes", subtitle = length(All_list), y=y) +
  theme(legend.position = "none", legend.box = "vertical")

###############
## ECDF PLOT ##
###############

## see Fig2_D_plotECDF_function.R for the function 
#pattern_colors  "#10AFBB" "#FC4E07"
MARCHF6_plot <- plot_gene_ECDF("MARCHF6", concat_dfFX_res, tst_df, extension, working_directory, pattern_colors, 2)
EGFR_plot <- plot_gene_ECDF("EGFR", concat_dfFX_res, tst_df, extension, working_directory, pattern_colors, 2)

################
## HISTOGRAMS ##
################

histo_per <- ggplot(tst_df %>% filter(Universe==TRUE & Group=="Attenuated"), aes(x=knee_AUC_HS/2)) +
  geom_histogram(binwidth=5, fill = "grey", color = "black", boundary=0)+
  xlim(0,100) + xlab("Distance TSS to knee (% of the gene)") +
  theme_classic()

histo_kb <- ggplot(tst_df %>% filter(Universe==TRUE & Group=="Attenuated"), aes(x=(knee_AUC_HS*window_size)/1000)) +
  geom_histogram(binwidth=10, fill = "grey", color = "black", boundary=0)+
  xlim(0, 350) + 
  #scale_x_continuous(breaks = c(seq(0, 150, by = 50), 200), 
  #                   labels = c(seq(0, 150, by = 50), ">200")) 
  labs(x = "Distance TSS to knee (kb)", y = "Count", title = "TSS to knee position in kb for attenuated genes") +
  theme_classic()   

##################
## SAVING FILES ##
##################

PATH_ggsave <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/R_script/plots_paper_function/Plot_V2_20240807/Figure_2/" 

# 2B Universe_metagene
# 2C AUC_HS_vs_control_pvalue
# 2D MARCHF6_plot
# 2E EGFR_plot
# 2F_left AUC_HS_vs_control_groups
# 2F_right_top Attenuated_metagene
# 2F_right_bottom Outgroup_metagene
# 2I histo_per
# 2J histo_kb

ggsave(paste0(PATH_ggsave,"Fig2_","Universe_metagene",".pdf"), plot = Universe_metagene, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","AUC_HS_vs_control_pvalue",".pdf"), plot = AUC_HS_vs_control_pvalue, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","MARCHF6_plot",".pdf"), plot = MARCHF6_plot, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","EGFR_plot",".pdf"), plot = EGFR_plot, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","AUC_HS_vs_control_groups",".pdf"), plot = AUC_HS_vs_control_groups, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","Attenuated_metagene",".pdf"), plot = Attenuated_metagene, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","Outgroup_metagene",".pdf"), plot = Outgroup_metagene, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","histo_per",".pdf"), plot = histo_per, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"Fig2_","histo_kb",".pdf"), plot = histo_kb, width = 8, height = 6, dpi=300)


