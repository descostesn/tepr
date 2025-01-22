library(ggpubr)

###############
## Plot list ##
###############
knee_AUC_DRB_plot
attenuted_gene_knee_ctrl
attenuted_gene_knee_HS
metagene_attenuated_genes
LRP5_DRB
RSF1_DRB
LRP5
RSF1

#############
## LOADING ##
#############

AUC_knee_DRB <- read.delim( file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/DRB_TTseq/Big_intermediary_file/AUC_knee_DRB.tsv",
                            sep = "\t",  header = T)

concat_dfFX_res_CugusiDRB <- read.delim( file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/DRB_TTseq/Big_intermediary_file/concat_dfFX_res_CugusiDRB.tsv",
                            sep = "\t", header = T)
# ## 2rep
concat_dfFX_res_Cugusi2rep <- read.delim( file = "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/RAW_TTseq/Big_intermediary_file/concat_dfFX_res_Cugusi2rep.tsv",
                            sep = "\t",  header = T)

############################
## ELONGATION RATE Fig 5A ##
############################


# Subset the data
data_universe <- AUC_knee_DRB %>% filter(Universe==T, 
                                         Group=="Outgroup",
                                         adjFDR_p_AUC_DRB_HS_30<0.1 , 
                                         adjFDR_p_AUC_DRB_HS_40<0.1 , 
                                         adjFDR_p_AUC_DRB_ctrl_30<0.1 , 
                                         adjFDR_p_AUC_DRB_ctrl_40 <0.1
) 
# Reshape the data into a long format
long_data <- data_universe %>% 
  select(gene, transcript, window_size, matches(".*10"), matches(".*20"), matches(".*30"), matches(".*40")) %>%
  #select(gene, transcript, window_size, contains("knee_AUC_DRB"), -matches("knee_AUC_DRB_.*_0")) %>%
  pivot_longer(
    cols = starts_with("knee_AUC_DRB_"),
    names_to = "condition",
    values_to = "knee_AUC_DRB"
  ) %>%
  mutate(knee_AUC_DRB = knee_AUC_DRB * window_size / 1000) %>%
  mutate(condition = factor(condition, levels = c("knee_AUC_DRB_ctrl_10", "knee_AUC_DRB_HS_10", 
                                                  "knee_AUC_DRB_ctrl_20", "knee_AUC_DRB_HS_20",
                                                  "knee_AUC_DRB_ctrl_30", "knee_AUC_DRB_HS_30",
                                                  "knee_AUC_DRB_ctrl_40", "knee_AUC_DRB_HS_40")))

knee_AUC_DRB_plot <- ggplot(long_data, aes(x = knee_AUC_DRB, y = condition, fill=condition, color=condition)) +
  geom_density_ridges(alpha = 0.8,  linewidth=1, scale=2) +
  scale_color_manual(values = c(
    "knee_AUC_DRB_HS_10" = "#FC4E07",
    "knee_AUC_DRB_HS_20" = "#FC4E07",
    "knee_AUC_DRB_HS_30" = "#FC4E07",
    "knee_AUC_DRB_HS_40" = "#FC4E07",
    "knee_AUC_DRB_ctrl_10" = "#00AFBB",
    "knee_AUC_DRB_ctrl_20" = "#00AFBB",
    "knee_AUC_DRB_ctrl_30" = "#00AFBB",
    "knee_AUC_DRB_ctrl_40" = "#00AFBB")) +
  scale_fill_manual(values = c(
    "knee_AUC_DRB_HS_10" = "#FC4E07",
    "knee_AUC_DRB_HS_20" = "#FC4E07",
    "knee_AUC_DRB_HS_30" = "#FC4E07",
    "knee_AUC_DRB_HS_40" = "#FC4E07",
    "knee_AUC_DRB_ctrl_10" = "#00AFBB",
    "knee_AUC_DRB_ctrl_20" = "#00AFBB",
    "knee_AUC_DRB_ctrl_30" = "#00AFBB",
    "knee_AUC_DRB_ctrl_40" = "#00AFBB")) +
  # scale_fill_viridis_d(direction = -1, guide = "none")+
  theme_classic() +
  xlim(0, 200) +
  xlab("Knee position (kb) in HS with DRB") +
  theme(legend.position = "none")

############################################
## MULTIMODE CALCULATION AND SPEED FIGURE ##
############################################
library("multimode")

HS_10 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_HS_10") %>% pull(knee_AUC_DRB) ),display=F)[[1]]
HS_20 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_HS_20") %>% pull(knee_AUC_DRB) ),display=F)[[1]] 
HS_30 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_HS_30") %>% pull(knee_AUC_DRB) ),display=F)[[1]]
HS_40 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_HS_40") %>% pull(knee_AUC_DRB) ),display=F)[[1]]
        
ctrl_10 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_ctrl_10") %>% pull(knee_AUC_DRB) ),display=F)[[1]]
ctrl_20 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_ctrl_20") %>% pull(knee_AUC_DRB) ),display=F)[[1]] 
ctrl_30 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_ctrl_30") %>% pull(knee_AUC_DRB) ),display=F)[[1]]
ctrl_40 <- locmodes(data= (long_data %>% filter(condition=="knee_AUC_DRB_ctrl_40") %>% pull(knee_AUC_DRB) ),display=F)[[1]]

# Sample data
time <- c(0, 10, 20, 30, 40)
HS <- c(0,HS_10,HS_20,HS_30,HS_40) # HS wave positions
ctrl <- c(0,ctrl_10,ctrl_20,ctrl_30,ctrl_40)# ctrl wave positions
 

# Create a data frame for ggplot
df <- data.frame(
  Time = rep(time, 2),
  Wave_position = c(HS, ctrl),
  Condition = rep(c("HS", "ctrl"), each = length(time))
)
rm(speed_knee_Outgroup)

# Plotting WITH linear regression
speed_knee_Outgroup_linear <- ggplot(df, aes(x = Time, y = Wave_position, color = Condition)) +
  geom_point(size = 3) +  # Points
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = F, size = 1) +  # Linear regression lines through (0,0)  
  scale_color_manual(values = c("#00AFBB","#FC4E07")) +  # Manual colors
  labs(x = "Time (min)", y = "Wave position (kb)", color = "") + xlim(0,42) +
 theme_classic()

# Linear regression models without intercept (forcing line through (0,0))
lm_HS <- lm(Wave_position ~ 0 + Time, data = df[df$Condition == "HS", ])
lm_ctrl <- lm(Wave_position ~ 0 + Time, data = df[df$Condition == "ctrl", ])

# Extract R-squared values
r2_HS <- summary(lm_HS)$r.squared
r2_ctrl <- summary(lm_ctrl)$r.squared

# Extract coefficients for the regression equations
coef_HS <- summary(lm_HS)$coefficients
coef_ctrl <- summary(lm_ctrl)$coefficients

# Create the regression equation as text (note: no intercept)
eq_HS <- paste0("HS: y = ", round(coef_HS[1, 1], 2), "x, R² = ", round(r2_HS, 4))
eq_ctrl <- paste0("Ctrl: y = ", round(coef_ctrl[1, 1], 2), "x, R² = ", round(r2_ctrl, 4))

# Plotting without linear regression
# speed_knee_Outgroup <- ggplot(df, aes(x = Time, y = Wave_position, color = Condition)) +
#   geom_line(size = 1) +  # Line plot
#   geom_point(size = 3) +  # Points
#   scale_color_manual(values = c("#00AFBB","#FC4E07")) +  # Manual colors
#   labs(x = "Time (min)", y = "Wave position (kb)", color = "") + xlim(0,42) +
#   theme_classic()



# Adding text annotations for the slopes
speed_knee_Outgroup_linear <- speed_knee_Outgroup_linear +
  annotate("text", x = 12, y = 12, label = paste0(round(ctrl_10/10,digits = 1), " kb/min"), color = "#00AFBB", size = 5) +
  annotate("text", x = 22, y = 40, label = paste0(round(ctrl_20/20,digits = 1), " kb/min"), color = "#00AFBB", size = 5) +
  annotate("text", x = 32, y = 70, label = paste0(round(ctrl_30/30,digits = 1), " kb/min"), color = "#00AFBB", size = 5) +
  annotate("text", x = 40, y = 100, label = paste0(round(ctrl_40/40,digits = 1), " kb/min"), color = "#00AFBB", size = 5) +
  annotate("text", x = 8, y = 38, label = paste0(round(HS_10/10,digits = 1), " kb/min"), color = "#FC4E07", size = 5) +
  annotate("text", x = 18, y = 78, label = paste0(round(HS_20/20,digits = 1), " kb/min"), color = "#FC4E07", size = 5) +
  annotate("text", x = 28, y = 125, label = paste0(round(HS_30/30,digits = 1), " kb/min"), color = "#FC4E07", size = 5) +
  annotate("text", x = 36, y = 160, label = paste0(round(HS_40/40,digits = 1), " kb/min"), color = "#FC4E07", size = 5) +
  annotate("text", x = 5, y = 150, label = eq_HS, color = "#FC4E07", size = 5, hjust = 0) + ## for linear regression
  annotate("text", x = 5, y = 100, label = eq_ctrl, color = "#00AFBB", size = 5, hjust = 0) + ## for linear regression
  theme(legend.position = "top")

########################################
## SCATTER PLOT ATTENUTED CTRL Fig 5B ##
########################################

## Adjust colors manually for higher saturation
library(scales)
intense_orange <- "orange"  # Keep as is
intense_lightgreen <- scales::muted("lightgreen", l = 70, c = 100)  # Increase chroma (saturation)
intense_pink <- scales::muted("pink", l = 70, c = 100)  # Increase chroma
intense_lightblue <- scales::muted("lightblue", l = 70, c = 100)  # Increase chroma

## subseting data
data_sub <- AUC_knee_DRB %>% filter(Universe==T, Group=="Attenuated", #knee_AUC_HS*window_size < 100000, size > 70000, 
                                    adjFDR_p_AUC_DRB_HS_30<0.1 ,
                                    adjFDR_p_AUC_DRB_HS_40 <0.1 ,
                                    adjFDR_p_AUC_DRB_ctrl_30<0.1 , 
                                    adjFDR_p_AUC_DRB_ctrl_40 <0.1,
                                    window_size>200) 

condi <- "ctrl"
# Create the scatter plot
Main <- ggplot(data_sub) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_10"))*window_size/1000), color=intense_lightblue) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_20"))*window_size/1000), color=intense_pink) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_30"))*window_size/1000), color=intense_lightgreen) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="orange") + 
  
  geom_abline(color="red")+  xlab("Knee position (kb) in HS ") +ylab("Knee position (kb) in HS with DRB") + labs(title=paste0("Attenuated genes in ", condi)) +
  ylim(0,170)+ xlim(0,170) +
  theme_classic()

# Create the density plot aligned to the right
yplot <- ggplot(subset(data_sub )) + #, knee_AUC_HS*window_size/1000 > 25 & knee_AUC_HS*window_size/1000 < 45)) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_10"))*window_size/1000), fill = intense_lightblue, alpha = 0.5 ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_20"))*window_size/1000), fill = intense_pink, alpha = 0.5 ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_30"))*window_size/1000), fill = intense_lightgreen, alpha = 0.5, ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), fill = "orange", alpha = 0.9, ) +
  coord_flip() + theme_classic() +
  xlim(0,170) + xlab("") +
  theme(legend.position = "right")


attenuted_gene_knee_ctrl <- ggarrange(Main, yplot, 
          ncol = 2, nrow = 1,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)


################################
## SCATTER PLOT ATTENUATED HS ##
################################

condi <- "HS"

Main <- ggplot(data_sub) +
   # geom_point(data=subset(data_sub, knee_AUC_DRB_ctrl_40<knee_AUC_DRB_ctrl_30),  ## genesize
  #            aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="black", size=3) +
  # geom_point(data=subset(data_sub, window_size*200/1000<40*4), # gene that are short enough to go through one round of transcription in 40min
  #            aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="green", size=3.5) +
# geom_point(data=subset(data_sub, knee_AUC_DRB_HS_40>knee_AUC_HS), ## genes in forbidden zone
#            aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="grey", size=2.5) + 

geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_10"))*window_size/1000), color=intense_lightblue) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_20"))*window_size/1000), color=intense_pink) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_30"))*window_size/1000), color=intense_lightgreen) +
  
  geom_point(data=subset(data_sub, knee_AUC_DRB_HS_40<knee_AUC_DRB_HS_20 | knee_AUC_DRB_HS_40<knee_AUC_DRB_HS_10), ## genes attenuated again
             aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="black", size=3) +
  geom_point(aes(knee_AUC_HS*window_size/1000, !!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), color="orange") + 
  
  geom_abline(color="red")+  xlab("Knee position (kb) in HS ") +ylab("Knee position (kb) in HS with DRB") + labs(title=paste0("Attenuated genes in ", condi)) +
  ylim(0,170)+ xlim(0,170) +
  theme_classic()


# Create the density plot aligned to the right
yplot <- ggplot(subset(data_sub )) + #, knee_AUC_HS*window_size/1000 > 25 & knee_AUC_HS*window_size/1000 < 45)) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_10"))*window_size/1000), fill = intense_lightblue, alpha = 0.5 ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_20"))*window_size/1000), fill = intense_pink, alpha = 0.5 ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_30"))*window_size/1000), fill = intense_lightgreen, alpha = 0.5, ) +
  geom_density(aes(!!sym(paste0("knee_AUC_DRB_",condi,"_40"))*window_size/1000), fill = "orange", alpha = 0.9, ) +
  coord_flip() + theme_classic() +
  xlim(0,170) + xlab("") +
  theme(legend.position = "right")


attenuted_gene_knee_HS <-ggarrange(Main, yplot, 
          ncol = 2, nrow = 1,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)
  
###################
## METAGENE PLOT ##
###################

### metagene analysis
#concat_dfFX_res created from results_ECDF with calculates_meanFx function
#This work also for single genes
Attenuation_list<- AUC_knee_DRB %>% filter(Group=="Attenuated", knee_AUC_HS*window_size/1000>25 &  knee_AUC_HS*window_size/1000<45 ) %>% pull(transcript)
Outgroup_list<- AUC_knee_DRB %>% filter(Group=="Outgroup", window_size>500 & window_size <700) %>% pull(transcript)

AUC_allcondi <- AUC_knee_DRB %>% select(transcript, gene, window_size)

normalize_and_summarize <- function(list,concat_dfFX_res ) {
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
#Attenuation_concat_df <- normalize_and_summarize(Attenuation_list, concat_dfFX_res_CugusiDRB)
Attenuation_concat_df <- normalize_and_summarize(Attenuation_list, concat_dfFX_res_CugusiDRB)
Attenuation2rep_concat_df <- normalize_and_summarize(Attenuation_list, concat_dfFX_res_Cugusi2rep)
Outgroup_concat_df <- normalize_and_summarize(Outgroup_list, concat_dfFX_res_CugusiDRB)
Outgroup2rep_concat_df <- normalize_and_summarize(Outgroup_list, concat_dfFX_res_Cugusi2rep)


metagene_attenuated_genes <- ggplot() +
  geom_line(data=Attenuation_concat_df, aes(x=coord, y= mean_value_DRB_HS_40), color="#FC4E07",  size=1.4)+
  geom_line(data=Attenuation_concat_df, aes(x=coord, y= mean_value_DRB_ctrl_40), color="#00AFBB",   size=1.4)+
  geom_line(data=Attenuation2rep_concat_df, aes(x=coord, y= mean_value_HS), color="#FC4E07",  size=0.7)+
  geom_line(data=Attenuation2rep_concat_df, aes(x=coord, y= mean_value_ctrl), color="#00AFBB",  size=0.7)+
  theme_bw() +   ylim(0,NA)+
  labs(x="")+
  labs(title = "Attenuated genes with knee between 25 and 45kb", subtitle = length(Attenuation_list), y=y) +
  theme(legend.position = "bottom", legend.box = "vertical")


##########
## ECDF ##
##########
#pattern_colors <- ask_for_colors(  "*.bg", "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/DRB_TTseq/FASTQ_FILES/folder_clean/BAM/deduplicated/MAPQ255/stranded/bedgraph255/" )

pattern_colors <- c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#10AFBB","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FC4E07")
#pattern_colors <- c("#FFFFFF",intense_lightblue,intense_pink,intense_lightgreen,"orange","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")

gene_name <- "LRP5"  # Example gene

LRP5_DRB <- 
  plot_gene_ECDF("LRP5", concat_dfFX_res_CugusiDRB, AUC_knee_DRB, "*.bg", 
                          "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/DRB_TTseq/FASTQ_FILES/folder_clean/BAM/deduplicated/MAPQ255/stranded/bedgraph255/", 
                 pattern_colors, 2)

RSF1_DRB <- 
  plot_gene_ECDF("RSF1", concat_dfFX_res_CugusiDRB, AUC_knee_DRB, "*.bg", 
                 "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Cugusi2022/DRB_TTseq/FASTQ_FILES/folder_clean/BAM/deduplicated/MAPQ255/stranded/bedgraph255/", 
                 pattern_colors, 2)

LRP5 <- 
  plot_gene_ECDF("LRP5", concat_dfFX_res, tst_df, extension, working_directory, c("#10AFBB","#FC4E07"), 2)
RSF1 <- 
  plot_gene_ECDF("RSF1", concat_dfFX_res, tst_df, extension, working_directory, c("#10AFBB","#FC4E07"), 2)

############
## SAVING ##
############
knee_AUC_DRB_plot
attenuted_gene_knee_ctrl
attenuted_gene_knee_HS
metagene_attenuated_genes
LRP5_DRB
RSF1_DRB
LRP5
RSF1
speed_knee_Outgroup
speed_knee_Outgroup_linear

PATH_ggsave <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/R_script/plots_paper_function/Plot_V2_20240807/Figure_5/" 

ggsave(paste0(PATH_ggsave,"knee_AUC_DRB_plot",".pdf"), plot = knee_AUC_DRB_plot, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"attenuted_gene_knee_ctrl",".pdf"), plot = attenuted_gene_knee_ctrl, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"attenuted_gene_knee_HS",".pdf"), plot = attenuted_gene_knee_HS, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"metagene_attenuated_genes",".pdf"), plot = metagene_attenuated_genes, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"LRP5_DRB",".pdf"), plot = LRP5_DRB, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"RSF1_DRB",".pdf"), plot = RSF1_DRB, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"LRP5",".pdf"), plot = LRP5, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"RSF1",".pdf"), plot = RSF1, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"speed_knee_Outgroup",".pdf"), plot = speed_knee_Outgroup, width = 8, height = 6, dpi=300)
ggsave(paste0(PATH_ggsave,"speed_knee_Outgroup_linear",".pdf"), plot = speed_knee_Outgroup_linear, width = 8, height = 6, dpi=300)
