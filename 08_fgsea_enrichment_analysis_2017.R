#fgsea - fast preranked gene set enrichment analysis (GSEA)
#2021-12-10
#Jason Toy

library(fgsea)
library(tidyverse)
library(edgeR)




#LOW STAT vs MOD STAT (Mod Stat = baseline)

#Alternatively, import edgeR DE results file from Trinity
result_table <- read.delim("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/edgeR_out_dir/RSEM.gene.counts.matrix.Low_Stat_vs_Mod_Stat.edgeR.DE_results", header = TRUE, sep = "\t") %>% 
  rownames_to_column(var = "gene_id")




#Annotate genes in result_table using Trinotate sprot annotation
#import annotation database
annot <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/trinotate_annotation/annotated_genes_sprot_blastx_clean_final.tsv", col_names = TRUE) %>% 
  dplyr::select(c(1:4))


#Add gene symbols and lengths to result_table by merging the data frames
gsea_table <- merge(result_table, annot, by="gene_id")

#Sort result_table by gene symbol, then by FDR, then gene length, then filter out unique values.
#The distinct() function keeps the first row of values (which is now the row with the lowest FDR/longest gene length) for each duplicate gene symbol.
gsea_table_uniq <- 
  gsea_table %>% 
  arrange(gene_symbol, FDR, desc(length)) %>% 
  distinct(gene_symbol, .keep_all = TRUE)

#create GSEA rnk file (list of ranked genes)
rank <- gsea_table_uniq %>% 
  mutate(score = sign(logFC) * -log10(PValue))

rnk_dat <- rank %>% 
  select(c(gene_symbol, score)) %>% 
  arrange(desc(score))

#write to tsv for use in cytoscape
#write_tsv(rnk_dat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/2017_lowstat_modstat.rnk")

#create scores vector
rnk <- pull(rnk_dat, score)
names(rnk) <- rnk_dat$gene_symbol


#plot ranks
barplot(rnk)

#load pathways (downloaded from gsea site on 2021-06-06)
kegg <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c2.cp.kegg.v7.4.symbols.gmt")
go <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c5.go.v7.4.symbols.gmt")
hallmark <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/h.all.v7.4.symbols.gmt")

pathways <- c(go, kegg, hallmark)

#writeGmtPathways(pathways, "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/combined.gmt")

#run fGSEA
fgsea_results <- fgsea(pathways, rnk, minSize=15, maxSize = 500)

#tidy results
fgsea_tidy <- fgsea_results %>% 
  as_tibble() %>%
  arrange(desc(NES))

#plot results
library(ggplot2)
ggplot(fgsea_tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO, KEGG, and Hallmark NES from GSEA")

fgsea_tidy_flat <- fgsea_tidy %>% 
  mutate(leadingEdge = as.character(leadingEdge)) #%>% 
#  rename("NOM p-val" = pval) %>% 
#  rename("FDR q-val" = padj) %>% 
#  rename(SIZE = size) %>% 
#  rename(NAME = pathway) %>% 
#  select(NAME, SIZE, ES, NES, "NOM p-val", "FDR q-val")

# write_tsv(fgsea_tidy_flat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_modstat_all_dbs_.tsv", col_names = TRUE)
# write_tsv(filter(fgsea_tidy_flat, NES > 0), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_modstat_all_dbs_pos_.tsv", col_names = TRUE)
# write_tsv(filter(fgsea_tidy_flat, NES < 0) %>% arrange(NES), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_modstat_all_dbs_neg_.tsv", col_names = TRUE)


#filter pathways with padj < 0.05
sig_up <- fgsea_tidy_flat %>% filter(NES > 0) %>% filter(padj < 0.05)
sig_down <- fgsea_tidy_flat %>% filter(NES < 0) %>% filter(padj < 0.05)

# write_tsv(sig_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_modstat_all_dbs_pos_sig_.tsv", col_names = TRUE)
# write_tsv(sig_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_modstat_all_dbs_neg_sig_.tsv", col_names = TRUE)





### MOD STAT vs. MOD VAR ### (Mod Stat = baseline)

# import edgeR DE results file from Trinity
result_table <- read.delim("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/edgeR_out_dir/RSEM.gene.counts.matrix.Mod_Stat_vs_Mod_Var.edgeR.DE_results", header = TRUE, sep = "\t") %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(logFC = logFC * -1)


#Annotate genes in result_table using Trinotate sprot annotation
#import annotation database
annot <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/trinotate_annotation/annotated_genes_sprot_blastx_clean.tsv", col_names = TRUE) %>% 
  dplyr::select(c(1:4))


#Add gene symbols and lengths to result_table by merging the data frames
gsea_table <- merge(result_table, annot, by="gene_id")

#Sort result_table by gene symbol, then by FDR, then gene length, then filter out unique values.
#The distinct() function keeps the first row of values (which is now the row with the lowest FDR/longest gene length) for each duplicate gene symbol.
gsea_table_uniq <- 
  gsea_table %>% 
  arrange(gene_symbol, FDR, desc(length)) %>% 
  distinct(gene_symbol, .keep_all = TRUE)

#create GSEA rnk file (list of ranked genes)
rank <- gsea_table_uniq %>% 
  mutate(score = sign(logFC) * -log10(PValue))

rnk_dat <- rank %>% 
  select(c(gene_symbol, score)) %>% 
  arrange(desc(score))

#write to tsv for use in cytoscape
#write_tsv(rnk_dat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/2017_modstat_modvar.rnk")

rnk <- pull(rnk_dat, score)
names(rnk) <- rnk_dat$gene_symbol


#plot ranks
barplot(rnk)

#load pathways (downloaded from gsea site on 2021-06-06)
kegg <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c2.cp.kegg.v7.4.symbols.gmt")
go <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c5.go.v7.4.symbols.gmt")
hallmark <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/h.all.v7.4.symbols.gmt")

pathways <- c(go, kegg, hallmark)

#writeGmtPathways(pathways, "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/combined.gmt")

#run fGSEA
fgsea_results <- fgsea(pathways, rnk, minSize=15, maxSize = 500)

#tidy results
fgsea_tidy <- fgsea_results %>% 
  as_tibble() %>%
  arrange(desc(NES))

#plot results
library(ggplot2)
ggplot(fgsea_tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO, KEGG, and Hallmark NES from GSEA")

fgsea_tidy_flat <- fgsea_tidy %>% 
  mutate(leadingEdge = as.character(leadingEdge)) 


#write_tsv(fgsea_tidy_flat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_modstat_v_modvar_all_dbs_.tsv", col_names = TRUE)
#write_tsv(filter(fgsea_tidy_flat, NES > 0), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_modstat_v_modvar_all_dbs_pos_.tsv", col_names = TRUE)
#write_tsv(filter(fgsea_tidy_flat, NES < 0) %>% arrange(NES), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_modstat_v_modvar_all_dbs_neg_.tsv", col_names = TRUE)


#filter pathways with padj < 0.05
sig_up <- fgsea_tidy_flat %>% filter(NES > 0) %>% filter(padj < 0.05)
sig_down <- fgsea_tidy_flat %>% filter(NES < 0) %>% filter(padj < 0.05)

#write_tsv(sig_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_modstat_v_modvar_all_dbs_pos_sig_.tsv", col_names = TRUE)
#write_tsv(sig_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_modstat_v_modvar_all_dbs_neg_sig_.tsv", col_names = TRUE)




### LOW STAT VS LOW VAR ### (Low Stat = baseline)

#Alternatively, import edgeR DE results file from Trinity
result_table <- read.delim("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/edgeR_out_dir/RSEM.gene.counts.matrix.Low_Stat_vs_Low_Var.edgeR.DE_results", header = TRUE, sep = "\t") %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(logFC = logFC * -1)


#Annotate genes in result_table using Trinotate sprot annotation
#import annotation database
annot <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/trinotate_annotation/annotated_genes_sprot_blastx_clean.tsv", col_names = TRUE) %>% 
  dplyr::select(c(1:4))


#Add gene symbols and lengths to result_table by merging the data frames
gsea_table <- merge(result_table, annot, by="gene_id")

#Sort result_table by gene symbol, then by FDR, then gene length, then filter out unique values.
#The distinct() function keeps the first row of values (which is now the row with the lowest FDR/longest gene length) for each duplicate gene symbol.
gsea_table_uniq <- 
  gsea_table %>% 
  arrange(gene_symbol, FDR, desc(length)) %>% 
  distinct(gene_symbol, .keep_all = TRUE)

#create GSEA rnk file (list of ranked genes)
rank <- gsea_table_uniq %>% 
  mutate(score = sign(logFC) * -log10(PValue))

rnk_dat <- rank %>% 
  select(c(gene_symbol, score)) %>% 
  arrange(desc(score))

#write to tsv for use in cytoscape
#write_tsv(rnk_dat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/2017_lowstat_lowvar.rnk")

rnk <- pull(rnk_dat, score)
names(rnk) <- rnk_dat$gene_symbol


#plot ranks
barplot(rnk)

#load pathways (downloaded from gsea site on 2021-06-06)
kegg <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c2.cp.kegg.v7.4.symbols.gmt")
go <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/c5.go.v7.4.symbols.gmt")
hallmark <- gmtPathways("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/h.all.v7.4.symbols.gmt")

pathways <- c(go, kegg, hallmark)

#writeGmtPathways(pathways, "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/gsea_databases/combined.gmt")

#run fGSEA
fgsea_results <- fgsea(pathways, rnk, minSize=15, maxSize = 500)

#tidy results
fgsea_tidy <- fgsea_results %>% 
  as_tibble() %>%
  arrange(desc(NES))

#plot results
library(ggplot2)
ggplot(fgsea_tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO, KEGG, and Hallmark NES from GSEA")

fgsea_tidy_flat <- fgsea_tidy %>% 
  mutate(leadingEdge = as.character(leadingEdge)) 


#write_tsv(fgsea_tidy_flat, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_lowvar_all_dbs_.tsv", col_names = TRUE)
#write_tsv(filter(fgsea_tidy_flat, NES > 0), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_lowvar_all_dbs_pos_.tsv", col_names = TRUE)
#write_tsv(filter(fgsea_tidy_flat, NES < 0) %>% arrange(NES), file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_lowvar_all_dbs_neg_.tsv", col_names = TRUE)


#filter pathways with padj < 0.05
sig_up <- fgsea_tidy_flat %>% filter(NES > 0) %>% filter(padj < 0.05)
sig_down <- fgsea_tidy_flat %>% filter(NES < 0) %>% filter(padj < 0.05)

#write_tsv(sig_up, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_lowvar_all_dbs_pos_sig_.tsv", col_names = TRUE)
#write_tsv(sig_down, file = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/enrichment_analysis/fgsea/fgsea_results_lowstat_v_lowvar_all_dbs_neg_sig_.tsv", col_names = TRUE)

