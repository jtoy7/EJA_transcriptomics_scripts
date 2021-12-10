#Enriched Pathways Overlap Analysis (Venn Diagrams)
#Jason A. Toy
#2021-12-10

library(tidyverse)


#read in fgsea results

up_2015 <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/2015_experiment/gsea/brain/fgsea/fgsea_results_brain_all_dbs_pos_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)

down_2015 <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/2015_experiment/gsea/brain/fgsea/fgsea_results_brain_all_dbs_neg_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)



up_2017 <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/lowstat_v_modstat/fgsea_results_lowstat_v_modstat_all_dbs_pos_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)

down_2017 <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/lowstat_v_modstat/fgsea_results_lowstat_v_modstat_all_dbs_neg_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)


#all
library(ggVennDiagram)
ggVennDiagram(list("Exp 1 Down" = down_2015,
                   "Exp 1 Up" = up_2015,
                   "Exp 2 Up" = up_2017,
                   "Exp 2 Down" = down_2017),
              color = "black",
              label_alpha = 0,
              label_color = "black",
              label = "count") +
  scale_fill_gradient(low="white",high = "red") +
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000')) +
  labs(fill="Count")


#save figure
#ggsave(
  "overlap_venn_2015_vs_2017_.jpg",
  plot = last_plot(),
  device = "jpg",
  path = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/",
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)


#2015 vs 2017 up
ggVennDiagram(list("2015 Up" = up_2015,
                   "2017 Up" = up_2017))

#2015 vs 2017 down
ggVennDiagram(list("2015 Down" = down_2015,
                   "2017 Down" = down_2017))



#Mod Stat/Var up vs Low Stat/Var up

#load data
up_mod <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/modstat_v_modvar/fgsea_results_modstat_v_modvar_all_dbs_pos_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)

down_mod <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/modstat_v_modvar/fgsea_results_modstat_v_modvar_all_dbs_neg_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)


up_low <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/lowstat_v_lowvar/fgsea_results_lowstat_v_lowvar_all_dbs_pos_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)

down_low <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/lowstat_v_lowvar/fgsea_results_lowstat_v_lowvar_all_dbs_neg_sig.tsv", col_names = TRUE) %>% 
  pull(pathway)


#diagram
ggVennDiagram(list("7.85 Up" = up_mod,
                   "7.70 Up" = up_low))


#Mod Stat/Var down vs. Low Stat/Var down
ggVennDiagram(list("7.85 Down" = down_mod,
                   "7.70 Down" = down_low))


#all
ggVennDiagram(x = list("7.85-Var Down" = down_mod,
                       "7.85-Var Up" = up_mod,
                       "7.70-Var Up" = up_low,
                       "7.70-Var Down" = down_low),
              color = "black",
              label_alpha = 0,
              label_color = "black",
              label = "count",
              set_size = 3.5) +
  scale_fill_gradient(low="white",high = "red") +
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000')) +
  labs(fill="Count")


#save figure
#ggsave(
  "overlap_venn_2017_stat_vs_var.jpg_",
  plot = last_plot(),
  device = "jpg",
  path = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/",
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)





### Extract overlapping pathways


#2015 vs 2017 comparisons

a <- list("Exp 1 Down" = down_2015,
          "Exp 1 Up" = up_2015,
          "Exp 2 Up" = up_2017,
          "Exp 2 Down" = down_2017)


venn_a <- Venn(a)


overlap_a <- process_region_data(venn_a)


exp_1_2_up <- overlap_a$item[8] %>% as.data.frame()
exp_1_2_down <- overlap_a$item[7] %>% as.data.frame()


#write to tsv
#write_tsv(exp_1_2_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_exp_1_2_up_.tsv", col_names = FALSE)
#write_tsv(exp_1_2_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_exp_1_2_down_.tsv", col_names = FALSE)



#Stat vs Var comparisons
b <- list("7.85 Down" = down_mod,
          "7.85 Up" = up_mod,
          "7.70 Up" = up_low,
          "7.70 Down" = down_low)  


venn_b <- Venn(b)


overlap_b <- process_region_data(venn_b)


mod_up_low_up <- overlap_b$item[8] %>% as.data.frame()
mod_down_low_down <- overlap_b$item[7] %>% as.data.frame()
mod_up_low_down <- overlap_b$item[9] %>% as.data.frame()
mod_down_low_up <- overlap_b$item[6] %>% as.data.frame()


#write to tsv
#write_tsv(mod_up_low_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_mod_up_low_up_.tsv", col_names = FALSE)
#write_tsv(mod_down_low_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_mod_down_low_down_.tsv", col_names = FALSE)
#write_tsv(mod_up_low_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_mod_up_low_down_.tsv", col_names = FALSE)
#write_tsv(mod_down_low_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/overlap_pathways_mod_down_low_up_.tsv", col_names = FALSE)



