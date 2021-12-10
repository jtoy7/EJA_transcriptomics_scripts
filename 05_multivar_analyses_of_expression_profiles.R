#nMDS Plot and PERMANOVA of Global Expression Profiles using vegan/adonis()
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


#create DGEList object (object type used by edgeR)

#assign treatment groups
trt <- c("Amb", "Amb", "Amb", "Amb", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Var", "Low_Var", "Low_Var", "Low_Var", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Var", "Mod_Var", "Mod_Var", "Mod_Var")

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

#flip matrix for use in downstream functions
transpose <- as.data.frame(t(log2cpm))



#create treatment label vector
trt <- c("Amb", "Amb", "Amb", "Amb", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Var", "Low_Var", "Low_Var", "Low_Var", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Var", "Mod_Var", "Mod_Var", "Mod_Var")
ph <- c("Amb", "Amb", "Amb", "Amb", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod")
var <- c("Amb", "Amb", "Amb", "Amb", "Stat", "Stat", "Stat", "Stat", "Var", "Var", "Var", "Var", "Stat", "Stat", "Stat", "Stat", "Var", "Var", "Var", "Var")

trt_no_amb <- c("Low_Stat", "Low_Stat", "Low_Stat", "Low_Stat", "Low_Var", "Low_Var", "Low_Var", "Low_Var", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Stat", "Mod_Var", "Mod_Var", "Mod_Var", "Mod_Var")
ph_no_amb <- c("Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod", "Mod")
var_no_amb <- c("Stat", "Stat", "Stat", "Stat", "Var", "Var", "Var", "Var", "Stat", "Stat", "Stat", "Stat", "Var", "Var", "Var", "Var")

#run metaMDS
mds <- metaMDS(transpose, distance = "manhattan")

#extract MDS scores
data.scores <- as.data.frame(scores(mds))

#add treatments columns to scores dataframe
data.scores$rep <- rownames(data.scores)
data.scores$treatment <- trt
data.scores$ph <- ph
data.scores$var <- var


#plot scores
library(ggpubr)
library(RColorBrewer)

#create plot data
pdata <- data.scores %>% 
  mutate(Treatment = ifelse(treatment == "Low_Stat", "7.70 static", treatment)) %>% 
  mutate(Treatment = ifelse(treatment == "Low_Var", "7.70 variable", Treatment)) %>% 
  mutate(Treatment = ifelse(treatment == "Mod_Stat", "7.85 static", Treatment)) %>% 
  mutate(Treatment = ifelse(treatment == "Mod_Var", "7.85 variable", Treatment)) %>% 
  mutate(Treatment = ifelse(treatment == "Amb", "ambient", Treatment))

ggplot(data = pdata) + 
  geom_point(aes(x = NMDS1, y = NMDS2, shape = Treatment, colour = Treatment), size = 4) +
  scale_shape_manual(values = c(17,17,15,15,16)) +
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A", "#FB9A99")) +
  #stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = treatment), level = 0.95) +      #data (prediction interval) ellipses
  stat_conf_ellipse(aes(x = NMDS1, y = NMDS2, color = Treatment), level = 0.95) + #confidence interval ellipses
  theme_bw()

#save plot
#ggsave(
  "nmds_global_expression_manhattan_conf_ellipse_95_log2cpm.jpg",
  plot = last_plot(),
  device = "jpg",
  path = "C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/manuscript/figs_final/",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)


# get Manhattan distances between samples
dist <- vegdist(transpose, method = "manhattan")

transpose_no_amb <- tail(transpose, -4)
dist_no_amb <- vegdist(transpose_no_amb, method = "manhattan")

# Cluster samples using hierarchical agglomerative clustering with complete linkage (default)
clust <- hclust(dist)

# Plot dendrogram based on clustering
plot(clust)

#----
#run PERMANOVA with trt as only factor
adonis1.results <- adonis(transpose ~ trt, method="manhattan", perm=1000000)
print(adonis1.results)

#run PERMANOVA with pH and variability level as factors
adonis2.results <- adonis(transpose_no_amb ~ ph_no_amb * var_no_amb, method="manhattan", perm=1000000)
print(adonis2.results)

#run pairwise.adonis for pairwise comparisons
library(pairwiseAdonis)

#factor = trt
pairwise_adonis_results_trt <- pairwise.adonis(transpose, factors = trt, sim.function = "vegdist", sim.method = "manhattan", perm = 1000000)

#factor = ph
pairwise_adonis_results_ph <- pairwise.adonis(transpose_no_amb, factors = ph_no_amb, sim.function = "vegdist", sim.method = "manhattan", perm = 1000000)


