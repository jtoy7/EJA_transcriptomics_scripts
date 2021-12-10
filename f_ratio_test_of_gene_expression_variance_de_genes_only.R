#f-ratio variance test of gene expression (DE genes only)
#fratio = higher variance over lower variance
#calculate f ratios, plot distribution of ratios, would expect it to be centered around 1. Color code
#Jason Toy
#2021-12-10

library(tidyverse)
library(vegan)
library(grid)
library(edgeR)

#set working directory
setwd("~/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2")

#import raw read counts matrix
raw <- read.table("RSEM.gene.counts.matrix", header = TRUE, row.names = 1)
#raw <- raw[,-(1:4)] #uncomment to remove ambient replicates


#create DGEList object (object type used by edgeR)

#assign treatment groups (with or without Amb)
trt <- c("Amb", "Amb", "Amb", "Amb", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Var", "Low_Var", "Low_Var", "Low_Var", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Var", "Mod_Var", "Mod_Var", "Mod_Var")
#trt <- c("Low_Stat", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Var", "Low_Var", "Low_Var", "Low_Var", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Var", "Mod_Var", "Mod_Var", "Mod_Var")

#form DGEList object using raw read counts matrix and trt list
y <- DGEList(raw, group = trt)

#display library sizes for each sample
summary(y$samples$lib.size)

#filter out genes with very low expression levels (filters list of 39,258 genes down to 33,597 genes)
keep <- filterByExpr(y)         #logical (T/F) list of genes to keep
yf <- y[keep, , keep.lib.sizes = FALSE]         #create new DGEList with only genes with sufficiently high counts; keep.lib.sizes = FALSE recalculates new library sizes



#display new library sizes for each sample
summary(yf$samples$lib.size)

#TPM normalize raw read counts using edgeR's calcNormFactors()
tpm <- calcNormFactors(yf)

#log2 transform using edgeR's cpm (prior.count = 2 is default value)
log2cpm <- cpm(tpm, log= TRUE, prior.count = 2)

#write TPM-normalized and centered data to file
#write_tsv(as.data.frame(log2cpm) %>% rownames_to_column(var = "gene_id"), "/home/jason/Documents/ucsc/projects/black_perch_exp_2017/de_analysis_round_2/nmds_and_permanova/log2cpm_tpm_expr_matrix", col_names = TRUE)

#import list of DE genes
de_genes <- read.csv("edgeR_out_dir/diffExpr.P0.001_C1.5.list_of_de_genes.csv", header = FALSE) %>% 
  mutate(V1 = as.character(V1))

#filter annotation database to keep only DE genes
log2cpm_de <- rownames_to_column(as.data.frame(log2cpm), var = "gene_id")

log2cpm_de <- filter(log2cpm_de, log2cpm_de$gene_id %in% de_genes$V1)

log2cpm_de <- column_to_rownames(log2cpm_de, var = "gene_id")

log2cpm_de <- as.matrix(log2cpm_de)


#convert from wide to long format
tidy = log2cpm_de %>%
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather("Amb_1", "Amb_2", "Amb_3", "Amb_4", "Low_Stat_1", "Low_Stat_2", "Low_Stat_3", "Low_Stat_4", "Low_Var_1", "Low_Var_2", "Low_Var_3", "Low_Var_4", "Mod_Stat_1", "Mod_Stat_2", "Mod_Stat_3", "Mod_Stat_4", "Mod_Var_1", "Mod_Var_2", "Mod_Var_3", "Mod_Var_4", key = rep, value = log2cpm)


#create treatment column
tidy_trt <- tidy %>% 
  mutate(trt = word(rep, 1, sep = "_\\d") %>% as.factor)

#calculate variance for each gene and treatment
variance <- tidy_trt %>% group_by(gene_id, trt) %>% summarise(var = var(log2cpm))

sum_var <- variance %>% group_by(trt) %>% summarise(mean_var = mean(var), sd = sd(var))

#boxplot of variances by treatment

boxplot(var ~ trt, data = variance,
        notch  = TRUE,
        ylim = c(0,20),
        ylab = expression("Variance of log "[2]*"CPM"),
        xlab = NULL
        )
  points(sum_var$mean_var, col = "red", pch = 19)
  
#with ggplot

pdata <- variance

#change treatment names
pdata$order <- recode_factor(pdata$trt, Amb = "1",
                           Low_Stat = "4",
                           Low_Var = "5",
                           Mod_Stat = "2",
                           Mod_Var = "3"
                           )
pdata <- pdata %>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = as.numeric(order))

pdata$trt <- recode_factor(pdata$trt, Amb = "Ambient",
                           Low_Stat = "7.70 static",
                           Low_Var = "7.70 variable",
                           Mod_Stat = "7.85 static",
                           Mod_Var = "7.85 variable"
)


ggplot(pdata, aes(x = reorder(trt, order), y = var)) +
  geom_boxplot(notch = TRUE, aes(color = trt)) +
  scale_color_manual(values = c("#FB9A99", "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=rel(4))
        ) +
  ylab(expression("Variance of log "[2]*"CPM"))

#without ambient
ggplot(pdata %>% filter(trt != "Ambient"), aes(x = reorder(trt, order), y = var)) +
  geom_boxplot(notch = TRUE, aes(color = trt), outlier.shape = NA) +
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=rel(4))
  ) +
  ylab(expression("Variance of log"["2"]*"CPM")) +
  ylim(c(0,10))

#horizontal
ggplot(pdata, aes(x = reorder(trt, -order), y = var)) +
  geom_boxplot(notch = TRUE, aes(color = trt)) +
  scale_color_manual(values = c("#FB9A99", "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  coord_flip() +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size=rel(4)),
        axis.title.y = element_blank()
  ) +
  ylab(expression("Variance of log "[2]*"CPM"))
  
  


#save plot
#  ggsave(
    "expression_variance_de_genes_only_outliers_removed.jpg",
    plot = last_plot(),
    device = "jpg",
    path = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/",
    width = 8,
    height = 6,
    units = "in",
    dpi = 300,
    limitsize = TRUE
  )
  


#separate data by treatment
amb <- variance %>% ungroup() %>%  filter(trt == "Amb")
lowstat <- variance %>% ungroup() %>%  filter(trt == "Low_Stat")
lowvar <- variance %>% ungroup() %>%  filter(trt == "Low_Var")
modstat <- variance %>% ungroup() %>%  filter(trt == "Mod_Stat")
modvar <- variance %>% ungroup() %>%  filter(trt == "Mod_Var")

#calculate f-ratios, then log transform them
fratio_low <- data.frame(gene_id = lowstat$gene_id, fratio = lowvar$var/lowstat$var) %>% mutate(logf = log(fratio))
fratio_mod <- data.frame(gene_id = modstat$gene_id, fratio = modvar$var/modstat$var) %>% mutate(logf = log(fratio))


#boxplots to get sense of data distribution
bplot_low <- boxplot(fratio_low$logf)
bplot_low$stats


bplot_mod <- boxplot(fratio_mod$logf)
bplot_mod$stats


#plot histograms
ggplot(fratio_low, aes(x = logf)) +
  geom_histogram(bins = 160) +
  ylim(c(0,70))

ggplot(fratio_mod, aes(x = logf)) +
  geom_histogram(bins = 160) +
  ylim(c(0,70))

#Run t-test
t.test(x = fratio_low$logf, mu = 0, alternative = "greater")
#t = 1.8234, df = 199, p-value = 0.03487

t.test(x = fratio_mod$logf, mu = 0, alternative = "greater")
#t = 3.7205, df = 199, p-value = 0.0001293

