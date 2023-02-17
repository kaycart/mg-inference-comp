### Carter et al.
### Metagenome inference comparison script
### Written by Kayla A. Carter
### February 17, 2023






##############################
### LOAD PACKAGES AND DATA ###
##############################



library(plyr)
library(gmodels)
library(tidyverse)
library(epiDisplay)
library(philentropy)
library(clValid)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
library(multgee)
library(geepack)
library(MASS)
library(nnet)
library(phyloseq)
library(DataCombine)
library(data.table)
library(ggpattern)
library(ggdendro)
library(dendsort)
library(rgr)
library(ggnewscale)
library(patchwork)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(scales)

# set working directory or include path in file names
data <-readRDS("Carter et al observed dataframe.rds")
# data should have 72 obs of 4425 vars
pi<-readRDS("Carter et al PICRUSt2 dataframe.rds")
# pi should have 72 obs of 5883 vars
tax<-readRDS("Carter et al Tax4Fun2 dataframe.rds")
# tax should have 72 obs of 7050 vars
ko <- readRDS("Carter et al KO category dataframe.rds")
# ko should have 2548 obs of 4 vars

# create function that is opposite of %in%
'%!in%' <- function(x,y)!('%in%'(x,y))






####################################
### GENERATE 16S PHYLOSEQ OBJECT ###
####################################



### generate phyloseq otu table

# generate transposed dataframe of 16S counts
rownames(data) <- data$sample
counts <- data.frame(t(data[, 2:116]))
# counts should have 115 obs of 72 vars

# generate otu table for phyloseq object
phylo_otu <- otu_table(counts, taxa_are_rows = T)

# assign sample names in phylo_sample
sample_names(phylo_otu) <- data$sample



### generate phyloseq sample data

# convert case, race, nugent score category to factor variables
data$C_CASE2 <- factor(data$C_CASE2)
data$C_MRACE <- factor(data$C_MRACE)
data$C_BVCAT <- factor(data$C_BVCAT)

# generate sample data for phyloseq object
phylo_sample <- phyloseq::sample_data(data[, c(1, 4420:4425)])
# phylo_sample should have 72 obs of 7 vars

# assign sample names in phylo_sample
sample_names(phylo_sample) <- data$sample



### generate phyloseq taxon table

# generate taxon names that have ";" as spacer between taxonomy levels (not ".")
taxa_names <- rownames(counts)
taxa_names <- str_replace(taxa_names, ".p_", ";p_")
taxa_names <- str_replace(taxa_names, ".c_", ";c_")
taxa_names <- str_replace(taxa_names, ".o_", ";o_")
taxa_names <- str_replace(taxa_names, ".f_", ";f_")
taxa_names <- str_replace(taxa_names, ".g_", ";g_")
taxa_names <- str_replace(taxa_names, ".s_", ";s_")

# split taxon names to generate level-specific taxonomy assignments, 
# replace empty cells with "Unassigned"
taxa_levels <- as.data.frame(str_split_fixed(taxa_names, ";", 7))
taxa_levels[taxa_levels == ""] <- "Unassigned"

# drop leading uninformative characters
taxa_levels$V1 <- str_replace(taxa_levels$V1, "d__", "")
taxa_levels$V2 <- str_replace(taxa_levels$V2, "p__", "")
taxa_levels$V3 <- str_replace(taxa_levels$V3, "c__", "")
taxa_levels$V4 <- str_replace(taxa_levels$V4, "o__", "")
taxa_levels$V5 <- str_replace(taxa_levels$V5, "f__", "")
taxa_levels$V6 <- str_replace(taxa_levels$V6, "g__", "")
taxa_levels$V7 <- str_replace(taxa_levels$V7, "s__", "")

# drop trailing uninformative characters
for(each in 1:5){
  taxa_levels$V2 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V2 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V2 <- str_replace(taxa_levels$V2, "\\.$", "")
  
  taxa_levels$V3 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V3 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V3 <- str_replace(taxa_levels$V2, "\\.$", "")
  
  taxa_levels$V4 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V4 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V4 <- str_replace(taxa_levels$V2, "\\.$", "")
  
  taxa_levels$V5 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V5 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V5 <- str_replace(taxa_levels$V2, "\\.$", "")
  
  taxa_levels$V6 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V6 <- str_replace(taxa_levels$V2, "_$", "")
  taxa_levels$V6 <- str_replace(taxa_levels$V2, "\\.$", "")
}

# set full taxon names as row names
rownames(taxa_levels) <- rownames(counts)

# generate tax_table for phyloseq object
phylo_tax <- phyloseq::tax_table(as.matrix(taxa_levels))



### generate phyloseq object
phylo <- phyloseq(phylo_otu, phylo_sample, phylo_tax)






###################################
### 16S HIERARCHICAL CLUSTERING ###
###################################



# run clValid to identify optimal number of clusters
valid <- clValid(counts, 
                 nClust = 2:10, 
                 clMethods = "hierarchical", 
                 validation = "internal", 
                 maxitems = 115, 
                 metric = "euclidean", 
                 method = "ward",
                 neighbSize = 10)
summary(valid)
optimalScores(valid)
# using 3 clusters to have L. iners dominated and mixed as separate clusters

# generate JSD matrix
jsd <- phyloseq::distance(phylo_otu, method = "jsd")
# jsd should have 2556 elements

# hierarchical clustering to generate clusters and dendrogram sorted by average distance
hclust <- dendsort(hclust(jsd,
                          method = "ward.D",
                          members = NULL),
                   type = "average")

# plot and save dendrogram, add rectangles around clusters, save order of leaf labels
dendro <- as.dendrogram(hclust)
plot(hclust)
rect.hclust(hclust, 
            k = 3,
            border = "red")
leaf <- dendro_data(dendro)$labels[3]

# generate cluster variable, add to data, generate factor version
cluster <- cutree(hclust,
                  k = 3)
cluster <- cluster[order(names(cluster))] 
data <- data[order(data$sample),] 
data <- cbind(data, cluster)
data$cluster_f <- factor(data$cluster,
                         levels = c(3, 2, 1),
                         labels = c("crisp", "iners", "mixed"))
# data should have 72 obs of 4427 vars






####################################
### ESTIMATE 16S ALPHA DIVERSITY ###
####################################



# estimate alpha diversity
alpha <- estimate_richness(phylo,
                           split = T,
                           measures = c("Observed",
                                        "Shannon",
                                        "Simpson",
                                        "InvSimpson"))
# alpha should have 72 obs of 4 vars

# add alpha to data
rownames(alpha) <- data$sample
alpha <- alpha[order(rownames(alpha)),] 
data <- cbind(data, alpha)
# data should have 72 obs of 4431 vars






##################################################
### UPDATE PHYLOSEQ OBJECT WITH CLUSTER, ALPHA ###
##################################################



# generate updated sample data for phyloseq object
phylo_sample <- phyloseq::sample_data(data[, c(1, 4420:ncol(data))])
# phylo_sample should have 72 obs of 13 vars

# assign sample names in phylo_sample
sample_names(phylo_sample) <- data$sample

# generate updated phyloseq object
phylo <- phyloseq(phylo_otu, phylo_sample, phylo_tax)






##############################
### DESCRIPTIVE STATISTICS ###
##############################



### summarize distributions of following variables overall and by case status:
### case status, race, age, gestational age at sample collection,
### gestational age at delivery, Nugent score category, cluster

tab1(data$C_CASE2)

tab1(data$C_MRACE)
tab1(data$C_MRACE[data$C_CASE2 == 1])
tab1(data$C_MRACE[data$C_CASE2 == 0])

summary(data$C_MAGE24)
summary(data$C_MAGE24[data$C_CASE2 == 1])
summary(data$C_MAGE24[data$C_CASE2 == 0])

summary(data$C_SAMGS2)
summary(data$C_SAMGS2[data$C_CASE2 == 1])
summary(data$C_SAMGS2[data$C_CASE2 == 0])

summary(data$C_DELGS2)
summary(data$C_DELGS2[data$C_CASE2 == 1])
summary(data$C_DELGS2[data$C_CASE2 == 0])

tab1(data$C_BVCAT)
tab1(data$C_BVCAT[data$C_CASE2 == 1])
tab1(data$C_BVCAT[data$C_CASE2 == 0])

tab1(data$cluster_f)
tab1(data$cluster_f[data$C_CASE == 1])
tab1(data$cluster_f[data$C_CASE == 0])



### summarize distributions of total 16S reads, total observed metagenome reads

summary(data$total_read_16s)
summary(data$total_read_meta_ko)






##############################################
### GENERATE FIGURE 2 ALPHA DIVERSITY PLOT ###
##############################################



### make long dataset for plotting
plot <- data[, colnames(data) %in% c("sample",
                                     "C_CASE2",
                                     "cluster_f",
                                     "Observed", 
                                     "Shannon", 
                                     "Simpson", 
                                     "InvSimpson")]
plot <- reshape2::melt(plot,
                       id.vars = c("sample", 
                                   "C_CASE2",
                                   "cluster_f"),
                       measure.vars = c("Observed",
                                        "Shannon", 
                                        "Simpson", 
                                        "InvSimpson"))
colnames(plot) <- c("sample", "C_CASE2", "cluster_f", "metric", "alpha")
# plot should have 288 obs of 5 vars

### generate plot
ggplot(data = plot, aes(x = cluster_f, y = alpha, color = cluster_f)) + 
  geom_violin(lwd = 0.75) + 
  scale_color_manual(values = c("#ca0020","#ffd92f","#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus")," dominated")),
                                expression(paste(italic("L. iners")," dominated")),
                                "Mixed"),
                     "Cluster") +
  geom_jitter(width = 0.2, height = 0.2, shape = 1) +
  scale_shape(solid = F) + 
  facet_wrap(plot$metric, scales = "free") + 
  xlab("") + 
  ylab("Alpha diversity") + 
  theme_light(base_size = 11) +
  theme(legend.text.align = 0,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(margin = margin(t = 0,
                                                    r = 15,
                                                    b = 0,
                                                    l = 0)))






####################################################################
### GENERATE FIGURE 1 CASE STATUS, 16S, OBSERVED METAGENOME PLOT ###
####################################################################



### generate dataframe for 16S stacked bar plot

# subset to sample, 16S relative abundance variables, cluster
plot <- data[, c(1, 118:232, 4426)]
# plot should have 72 obs of 117 vars

# subset by cluster
plot_1 <- plot[plot$cluster == 1,]
# plot_1 should have 24 obs of 117 vars
plot_2 <- plot[plot$cluster == 2,]
# plot_2 should have 31 obs of 117 vars
plot_3 <- plot[plot$cluster == 3,]
# plot_3 should have 17 obs of 117 vars



# generate dataframe of mean relative abundance of each taxon by cluster
otu_relabund <- c()
cluster <- c()
otu <- c()
for(x in 2:116){
  cluster1_x <- plot_1[, x]
  cluster2_x <- plot_2[, x]
  cluster3_x <- plot_3[, x]
  otu_x_means <- c(mean(cluster1_x), 
                   mean(cluster2_x), 
                   mean(cluster3_x))
  otu_x <- c(x, x, x)
  cluster_x <- c(1, 2, 3)
  otu_relabund <- append(otu_relabund, otu_x_means)
  otu <- append(otu, otu_x)
  cluster <- append(cluster, cluster_x)
}
prelim_plotdata <- data.frame(cluster, otu, otu_relabund)
prelim_plotdata$cluster <- factor(prelim_plotdata$cluster,
                                  levels = c(3, 2, 1))
# prelim_plotdata should have 345 obs of 3 vars



# retrieve top 10 taxa in each cluster
prelim_plotdata <- prelim_plotdata[order(prelim_plotdata$otu_relabund, 
                                         decreasing=T),]

prelim_plotdata_1 <- prelim_plotdata[prelim_plotdata$cluster == 1,]
# prelim_plotdata_1 should have 115 obs of 3 vars
top10_1 <- prelim_plotdata_1$otu[1:10]
(top10_1_names <- colnames(plot_1)[top10_1])

prelim_plotdata_2 <- prelim_plotdata[prelim_plotdata$cluster == 2,]
# prelim_plotdata_2 should have 115 obs of 3 vars
top10_2 <- prelim_plotdata_2$otu[1:10]
(top10_2_names <- colnames(plot_2)[top10_2])

prelim_plotdata_3 <- prelim_plotdata[prelim_plotdata$cluster == 3,]
# prelim_plotdata_3 should have 115 obs of 3 vars
top10_3 <- prelim_plotdata_3$otu[1:10]
(top10_3_names <- colnames(plot_3)[top10_3])

(taxa_plot <- unique(c(top10_1_names, 
                       top10_2_names, 
                       top10_3_names)))
# use for stacked bar plot:
# L. antri (otu 47 in prelim_plot_data), L. crispatus (49), L. iners (50), L. jensenii (51), L. paragasseri (52), 
# Lactobacillus (46), Lactobacillus metagenome (54)
# A. vaginae (17), Bifidobacteriaceae (8), Gardnerella (10), Gardnerella uncultured (11),  
# Megasphaera unidentified (94), S. sanguinegens (100), Streptococcus agalactiae (56),
# Ureaplasma (60)



# subset plot_1, plot_2, plot_3 to taxa of interest; add variable for "other" relative abundance, sample;
# rename columns
taxa_names <- colnames(plot_1)[c(47, 49, 50, 51, 52, 46, 54, 17, 8, 10, 11, 94, 100, 56, 60)]

plot_1 <- plot_1[, colnames(plot_1) %in% taxa_names]
plot_1$other <- 1 - rowSums(plot_1)
plot_1$sample <- rownames(plot_1)
colnames(plot_1) <- c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", "lpara",
                      "lmet", "strepa", "urea", "mega", "ssang", "other", "sample")
# plot_1 should have 24 obs of 17 vars

plot_2 <- plot_2[, colnames(plot_2) %in% taxa_names]
plot_2$other <- 1 - rowSums(plot_2)
plot_2$sample <- rownames(plot_2)
colnames(plot_2) <- c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", "lpara",
                      "lmet", "strepa", "urea", "mega", "ssang", "other", "sample")
#plot_2 should have 31 obs of 17 vars

plot_3 <- plot_3[, colnames(plot_3) %in% taxa_names]
plot_3$other <- 1 - rowSums(plot_3)
plot_3$sample <- rownames(plot_3)
colnames(plot_3) <- c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", "lpara",
                      "lmet", "strepa", "urea", "mega", "ssang", "other", "sample")
#plot_3 should have 17 obs of 17 vars



# make long dataset for plotting
plot_1 <- melt(data = plot_1, 
               id.vars = c("sample"),
               measure.vars = c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", 
                                "lpara", "lmet", "strepa", "urea", "mega", "ssang", "other"))
plot_1$cluster <- 1
# plot_1 should have 384 obs of 4 vars

plot_2 <- melt(data = plot_2, 
               id.vars = c("sample"),
               measure.vars = c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", 
                                "lpara", "lmet", "strepa", "urea", "mega", "ssang", "other"))
plot_2$cluster <- 2
# plot_2 should have 496 obs of 4 vars

plot_3 <- melt(data = plot_3, 
               id.vars = c("sample"),
               measure.vars = c("bif", "gard", "gard_un", "avag", "lacto", "lantri", "lcrisp", "liners", "ljens", 
                                "lpara", "lmet", "strepa", "urea", "mega", "ssang", "other"))
plot_3$cluster <- 3
# plot_3 should have 272 obs of 4 vars

plotdata_16s <- rbind(plot_1, plot_2, plot_3)
plotdata_16s <- plotdata_16s[order(plotdata_16s$cluster,
                           decreasing = T),]
# plotdata_16s should have 1152 obs of 4 vars



# make taxon name variable, factor taxon order variable, sample order variable
plotdata_16s <- plotdata_16s %>% mutate(name =
                                          case_when(
                                            variable == "lantri" ~ "Lactobacillus antri",
                                            variable == "lcrisp" ~ "Lactobacillus crispatus",
                                            variable == "liners" ~ "Lactobacillus iners",
                                            variable == "ljens" ~ "Lactobacillus jensenii",
                                            variable == "lpara" ~ "Lactobacillus paragasseri",
                                            variable == "lacto" ~ "Lactobacillus",
                                            variable == "lmet" ~ "Lactobacillus metagenome",
                                            variable == "avag" ~ "Fannyhessea vaginae",
                                            variable == "bif" ~ "Bifidobacteriaceae",
                                            variable == "gard" ~ "Gardnerella",
                                            variable == "gard_un" ~ "Gardnerella uncultured",
                                            variable == "mega" ~ "Megasphaera unidentified",
                                            variable == "urea" ~ "Ureaplasma",
                                            variable == "ssang" ~ "Sneathia sanguinegens",
                                            variable == "strepa" ~ "Streptococcus agalactiae",
                                            variable == "other" ~ "Other"
                                  ))

plotdata_16s$name <- factor(plotdata_16s$name, 
                            levels = c("Lactobacillus antri",
                                       "Lactobacillus crispatus",
                                       "Lactobacillus iners",
                                       "Lactobacillus jensenii",
                                       "Lactobacillus paragasseri",
                                       "Lactobacillus",
                                       "Lactobacillus metagenome",
                                       "Bifidobacteriaceae",
                                       "Fannyhessea vaginae",
                                       "Gardnerella",
                                       "Gardnerella uncultured",
                                       "Megasphaera unidentified",
                                       "Sneathia sanguinegens",
                                       "Streptococcus agalactiae",
                                       "Ureaplasma",
                                       "Other"))

plotdata_16s <- plotdata_16s %>% mutate(order =
                                          case_when(
                                            variable == "lantri" ~ 1,
                                            variable == "lcrisp" ~ 2,
                                            variable == "liners" ~ 3,
                                            variable == "ljens" ~ 4,
                                            variable == "lpara" ~ 5,
                                            variable == "lacto" ~ 6,
                                            variable == "lmet" ~ 7,
                                            variable == "bif" ~ 8,
                                            variable == "avag" ~ 9,
                                            variable == "gard" ~ 10,
                                            variable == "gard_un" ~ 11,
                                            variable == "mega" ~ 12,
                                            variable == "ssang" ~ 13,
                                            variable == "strepa" ~ 14,
                                            variable == "urea" ~ 15,
                                            variable == "other" ~ 16
                                  ))

plotdata_16s <- plotdata_16s[order(plotdata_16s$order),]
plotdata_16s$sample <- factor(plotdata_16s$sample,
                              levels = leaf$label)
plotdata_16s <- plotdata_16s[order(plotdata_16s$sample),]
# plotdata_16s should have 1152 obs of 6 vars

# merge in case status, add variable for y position of case status tiles in plot
plot_case <- data[, colnames(data) %in% c("sample", "C_CASE2")]
plot_case$sample <- factor(plot_case$sample,
                           levels = leaf$label)
plotdata_16s <- left_join(plotdata_16s, plot_case)
plotdata_16s$case_y <- -0.04
# plotdata_16s should have 1152 obs of 8 vars



### generate 16S composition stacked bar plot
# will pull legend for later patchworked final plot
stacked_bar_16s <- ggplot(data = plotdata_16s, aes(x = factor(sample), y = value)) +
  geom_bar(stat = "identity",
           aes(fill = factor(order, levels = c("16",
                                               "15",
                                               "14",
                                               "13",
                                               "12",
                                               "11",
                                               "10",
                                               "9",
                                               "8",
                                               "7",
                                               "6",
                                               "5",
                                               "4",
                                               "3",
                                               "2",
                                               "1")))) +
  scale_fill_manual(values = c("gray48",
                               "#006837",
                               "#31a354",
                               "#c2e699",
                               "#ffffcc",
                               "#fbb4b9",
                               "#f768a1",
                               "#c51b8a",
                               "#7a0177",
                               
                               "#c6dbef",
                               "#9ecae1",
                               "#6baed6",
                               "#2171b5",
                               "#08519c",
                               "#08306b",
                               "#031737"),
                    labels = c("Other",
                               expression(paste(italic("Ureaplasma"), "")),
                               expression(paste(italic("Streptococcus agalactiae"), "")),
                               expression(paste(italic("Sneathia sanguinegens"), "")),
                               expression(paste(italic("Megasphaera"), " unidentified")),
                               expression(paste(italic("Gardnerella"), " uncultured")),
                               expression(paste(italic("Gardnerella"), "")),
                               expression(paste(italic("Fannyhessea vaginae"), "")),
                               expression(paste(italic("Bifidobacteriaceae"), "")),
                               expression(paste(italic("Lactobacillus"), " metagenome")),
                               expression(paste(italic("Lactobacillus"), "")),
                               expression(paste(italic("Lactobacillus paragasseri"), "")),
                               expression(paste(italic("Lactobacillus jensenii"), "")),
                               expression(paste(italic("Lactobacillus iners"), "")),
                               expression(paste(italic("Lactobacillus cripatus")," ")),
                               expression(paste(italic("Lactobacillus antri"), ""))),
                    "Bacterial taxa") +
  guides(fill = guide_legend(reverse = T,
                             override.aes = list(size = 4))) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x="", y = "Taxon relative abundance") +
  geom_rect(aes(xmin = 0.5, xmax = 72.5,
                ymin = 1.0275, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 0.5, xmax = 0.35, 
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 24.43, xmax = 24.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 55.43, xmax = 55.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 72.35, xmax = 72.5,
                ymin = 1.01, ymax = 1.03),
            fill = "black") +
  annotate(geom = "text", x = 64, y = 1.06,
           label = expression(paste(italic("L. crispatus"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 40, y = 1.06,
           label = expression(paste(italic("L. iners"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 12.5, y = 1.06,
           label = "Mixed", 
           color = "black", size = 2.5, angle = 90) +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stacked_bar_16s

# pull legend for patchworking
legend_16 <- as_ggplot(get_legend(stacked_bar_16s))



### generate 16S composition stacked bar plot with tiles for case status, no legend
# will use plot for later patchworked final plot
stacked_bar_16s_case_nolegend <- ggplot(data = plotdata_16s, aes(x = factor(sample), y = value)) +
  geom_bar(stat = "identity",
           aes(fill = factor(order, levels = c("16",
                                               "15",
                                               "14",
                                               "13",
                                               "12",
                                               "11",
                                               "10",
                                               "9",
                                               "8",
                                               "7",
                                               "6",
                                               "5",
                                               "4",
                                               "3",
                                               "2",
                                               "1")))) +
  scale_fill_manual(values = c("gray48",
                               "#006837",
                               "#31a354",
                               "#c2e699",
                               "#ffffcc",
                               "#fbb4b9",
                               "#f768a1",
                               "#c51b8a",
                               "#7a0177",
                               
                               "#c6dbef",
                               "#9ecae1",
                               "#6baed6",
                               "#2171b5",
                               "#08519c",
                               "#08306b",
                               "#031737"),
                    labels = c("Other",
                               expression(paste(italic("Ureaplasma"), "")),
                               expression(paste(italic("Streptococcus agalactiae"), "")),
                               expression(paste(italic("Sneathia sanguinegens"), "")),
                               expression(paste(italic("Megasphaera"), " unidentified")),
                               expression(paste(italic("Gardnerella"), " uncultured")),
                               expression(paste(italic("Gardnerella"), "")),
                               expression(paste(italic("Fannyhessea vaginae"), "")),
                               expression(paste(italic("Bifidobacteriaceae"), "")),
                               expression(paste(italic("Lactobacillus"), " metagenome")),
                               expression(paste(italic("Lactobacillus"), "")),
                               expression(paste(italic("Lactobacillus paragasseri"), "")),
                               expression(paste(italic("Lactobacillus jensenii"), "")),
                               expression(paste(italic("Lactobacillus iners"), "")),
                               expression(paste(italic("Lactobacillus cripatus")," ")),
                               expression(paste(italic("Lactobacillus antri"), ""))),
                    "Bacterial taxa") +
  guides(fill = guide_legend(reverse = T,
                             override.aes = list(size = 4),
                             order = 3)) +
  new_scale("fill") +
  geom_tile(data = plotdata_16s, aes(x = sample, y = case_y,
                                     fill = factor(C_CASE2, levels = c(1, 0)),
                                     height = 0.05, width = 0.875)) +
  scale_fill_manual(values = c("red3", "black"),
                    labels = c("Preterm case", "Term control"),
                    "Birth outcome") +
  guides(fill = guide_legend(reverse = F, 
                             override.aes = list(size = 4),
                             order = 2)) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x = "", y = "Taxon relative abundance") +
  geom_rect(aes(xmin = 0.5, xmax = 72.5,
                ymin = 1.0275, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 0.5, xmax = 0.35, 
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 24.43, xmax = 24.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 55.43, xmax = 55.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 72.35, xmax = 72.5,
                ymin = 1.01, ymax = 1.03),
            fill = "black") +
  annotate(geom = "text", x = 64, y = 1.06,
           label = expression(paste(italic("L. crispatus"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 40, y = 1.06,
           label = expression(paste(italic("L. iners"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 12.5, y = 1.06,
           label = "Mixed", 
           color = "black", size = 2.5, angle = 90) +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stacked_bar_16s_case_nolegend



### generate plot of tiles for case status
# will pull legend for later patchworked final plot
tiles <- ggplot(data = plotdata_16s, aes(x = factor(sample), y = value)) +
  geom_tile(data = plotdata_16s, aes(x = sample, y = case_y,
                                     fill = factor(C_CASE2, levels = c(1, 0)),
                                     height = 0.05, width = 0.875)) +
  scale_fill_manual(values = c("red3", "black"),
                    labels = c("Preterm case", "Term control"),
                    "Birth outcome") +
  guides(fill = guide_legend(reverse = F, 
                             override.aes = list(size = 4),
                             order = 2)) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x = "", y = "") +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# tiles

# pull legend for patchworking
legend_tile <- as_ggplot(get_legend(tiles))



### generate plot of dendrogram 
# to use in final patchworked plot
dendrogram <- ggplot(segment(dendro_data(dendro, type = "rectangle"))) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lwd = .35) +
  coord_flip() + 
  scale_y_reverse() +
  theme_void() +
  scale_x_reverse(expand = expansion(mult = 0.01))
# dendrogram



### generate plot spacer 
# to use in final patchworked plot
plot_spacer <- plot_spacer()



### prepare long WMGS dataset for plotting

# subset to observed metagenome level 1 KO relative abundances, sample, case status
plotdata_meta <- data[, c(1, 4305:4312, 4421)]

# generate long dataframe for observed metagenome KO (highest level) stacked bar plot
plotdata_meta <- melt(data = plotdata_meta, id.vars = c("sample", "C_CASE2"),
                 measure.vars = c("met_ra", "gen_ra", "env_ra", "cel_ra", 
                                  "org_ra", "hum_ra", "bri_ra", "not_ra"))
plotdata_meta$sample <- factor(plotdata_meta$sample,
                          levels = leaf$label)
plotdata_meta <- plotdata_meta[order(plotdata_meta$sample),]
# plotdata_meta should have 576 obs of 4 vars



### generate observed metagenome KO level 1 stacked bar plot
# will pull legend for final patchworked plot
stacked_bar_meta <- ggplot(data = plotdata_meta, aes(x = factor(sample), y = value)) +
  geom_bar(stat = "identity", 
           aes(fill = variable)) +
  scale_fill_manual("Functional category",
                    breaks=c("met_ra", 
                             "gen_ra", 
                             "env_ra", 
                             "cel_ra", 
                             "org_ra", 
                             "hum_ra", 
                             "bri_ra", 
                             "not_ra"),
                    values=c("#CC0000",
                             "#FF7733",
                             "#FFD500",
                             "#41ab5d",
                             "#00FFFF",
                             "#3377FF",
                             "#FF66B3",
                             "#8c510a"),
                    labels=c("Metabolism",
                             "Genetic information processing",
                             "Environmental information processing",
                             "Cellular processes",
                             "Organismal systems",
                             "Human diseases",
                             "BRITE hierarchies",
                             "Not included in pathway or BRITE")) +
  guides(fill = guide_legend(reverse = F, override.aes = list(size = 4))) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x = "", y = "KO functional category relative abundance") +
  geom_rect(aes(xmin = 0.5, xmax = 72.5,
                ymin = 1.0275, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 0.5, xmax = 0.35, 
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 24.43, xmax = 24.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 55.43, xmax = 55.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 72.35, xmax = 72.5,
                ymin = 1.01, ymax = 1.03),
            fill = "black") +
  annotate(geom = "text", x = 64, y = 1.06,
           label = expression(paste(italic("L. crispatus"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 40, y = 1.06,
           label = expression(paste(italic("L. iners"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 12.5, y = 1.06,
           label = "Mixed", 
           color = "black", size = 2.5, angle = 90) +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stacked_bar_meta

# pull legend for patchworking
legend_meta <- as_ggplot(get_legend(stacked_bar_meta))



### generate observed metagenome KO level 1 stacked bar plot without legend
# will use in for final patchworked plot
stacked_bar_meta_nolegend <- ggplot(data = plotdata_meta, aes(x = factor(sample), y = value)) +
  geom_bar(stat = "identity", 
           aes(fill = variable)) +
  scale_fill_manual("Functional category",
                    breaks=c("met_ra", 
                             "gen_ra", 
                             "env_ra", 
                             "cel_ra", 
                             "org_ra", 
                             "hum_ra", 
                             "bri_ra", 
                             "not_ra"),
                    values=c("#CC0000",
                             "#FF7733",
                             "#FFD500",
                             "#41ab5d",
                             "#00FFFF",
                             "#3377FF",
                             "#FF66B3",
                             "#8c510a"),
                    labels=c("Metabolism",
                             "Genetic information processing",
                             "Environmental information processing",
                             "Cellular processes",
                             "Organismal systems",
                             "Human diseases",
                             "BRITE hierarchies",
                             "Not included in pathway or BRITE")) +
  guides(fill = guide_legend(reverse = F, override.aes = list(size = 4))) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x = "", y = "KO functional category relative abundance") +
  geom_rect(aes(xmin = 0.5, xmax = 72.5,
                ymin = 1.0275, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 0.5, xmax = 0.35, 
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 24.43, xmax = 24.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 55.43, xmax = 55.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 72.35, xmax = 72.5,
                ymin = 1.01, ymax = 1.03),
            fill = "black") +
  annotate(geom = "text", x = 64, y = 1.06,
           label = expression(paste(italic("L. crispatus"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 40, y = 1.06,
           label = expression(paste(italic("L. iners"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 12.5, y = 1.06,
           label = "Mixed", 
           color = "black", size = 2.5, angle = 90) +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stacked_bar_meta_nolegend


### generate figure 1 patchworked plot
(stacked_bars_patch<-(dendrogram | 
                        stacked_bar_16s_case_nolegend | 
                        stacked_bar_meta_nolegend | 
                        legend_tile | 
                        legend_meta | 
                        legend_16 | 
                        plot_spacer)) + 
  plot_layout(design = c(area(1, 1, 20, 10), #dendrogram area
                         area(1, 9, 20, 58), #16S stacked bar and case status tile area
                         area(1, 58, 20, 107), #metagenome stacked bar area
                         area(0, 114, 2), #case status legend area
                         area(4.5, 123, 8), #metagenome legend area
                         area(11, 119, 20), #16S legend area
                         area(1, 135, 20, 138))) #spacer area






######################################################################################
### GENERATE SUPPLEMENTAL FIGURE 2 OBSERVED METAGENOME KO LEVEL 2 STACKED BAR PLOT ###
######################################################################################



### generate long dataset for plotting

# subset to observed metagenome level 2 KO relative abundances, sample, case status
plotdata_meta2 <- data[, c(1, 4367:4419, 4421)]

# generate long datasdataframeet for plotting
plotdata_meta2 <- melt(data = plotdata_meta2, id.vars = c("sample", "C_CASE2"),
                       measure.vars = colnames(plotdata_meta2)[colnames(plotdata_meta2) %!in% 
                                                                 c("sample","C_CASE2")])
plotdata_meta2$sample <- factor(plotdata_meta2$sample,
                               levels = leaf$label)
plotdata_meta2 <- plotdata_meta2[order(plotdata_meta2$sample),]

# add y values for tile location for case status
plotdata_meta2$case_y <- -0.04
# plotdata_meta2 should have 3816 obs of 5 vars

### plot ordered by 16S hierarchical clustering dendrogram with tiles for race and birth
stacked_bar_meta2 <- ggplot(data = plotdata_meta2, aes(x = factor(sample), y = value)) +
  geom_bar(stat = "identity",
           aes(fill = variable)) +
  scale_fill_manual("Functional category",
                    values=c("#521414", 
                             "#660000", 
                             "#8A0F0F", 
                             "#CC0000", 
                             "#FF3333", 
                             "#EB4747", 
                             "#FF6666", 
                             "#E87D7D", 
                             "#FF9999", 
                             "#FFCCCC", 
                             "#F7D4D4",
                             
                             "#CC4400",
                             "#FF5500",
                             "#FF7733",
                             "#FF9966",
                             "#FFBB99",
                             
                             "#FFD500",
                             "#FFEE99",
                             
                             "#005a32",
                             "#238b45",
                             "#41ab5d",
                             "#74c476",
                             "#a1d99b",
                             "#c7e9c0",
                             
                             "#1F4747",
                             "#006666",
                             "#178282",
                             "#1FADAD",
                             "#00CCCC",
                             "#00FFFF",
                             "#66FFFF",
                             "#99FFFF",
                             "#CCFFFF",
                             
                             "#002266",
                             "#003399",
                             "#0055FF",
                             "#3377FF",
                             "#6699FF",
                             "#99BBFF",
                             "#CCDDFF",
                             
                             "#330066",
                             "#6600CC",
                             "#9933FF",
                             "#CC99FF",
                             "#E5CCFF",
                             
                             "#FF3399",
                             "#FF66B3",
                             "#FF99CC",
                             "#FFCCE6",
                             
                             "#543005",
                             "#8c510a",
                             "#bf812d",
                             "#dfc27d"),
                    labels=c("Carbohydrate metabolism",
                             "Energy metabolism",
                             "Lipid metabolism",
                             "Nucleotide metabolism",
                             "Amino acid metabolism",
                             "Metabolism of other amino acids",
                             "Glycan biosynthesis and metabolism",
                             "Metabolism of cofactors and vitamins",
                             "Metabolism of terpenoids and polyketides",
                             "Biosynthesis of other secondary metabolites",
                             "Xenobiotics degradation and metabolism",
                             
                             "Transcription",
                             "Translation",
                             "Folding, sorting and degradation",
                             "Replication and repair",
                             "Information processing in viruses",
                             
                             "Membrane transport",
                             "Signal transduction",
                             
                             "Transport and catabolism",
                             "Cell motility",
                             "Cell growth and death",
                             "Cellular community - eukaryotes",
                             "Cellular community - prokaryotes",
                             "Aging",
                             
                             "Immune system",
                             "Endocrine system",
                             "Circulatory system",
                             "Digestive system",
                             "Excretory system",
                             "Nervous system",
                             "Sensory system",
                             "Development and regulation",
                             "Environmental adaptation",
                             
                             "Cancer: overview",
                             "Cancer: specific types",
                             "Immune disease",
                             "Neurodegenerative disease",
                             "Substance dependence",
                             "Cardiovascular disease",
                             "Endocrine and metabolic disease",
                             
                             "Infectious disease: bacterial",
                             "Infectious disease: viral",
                             "Infectious disease: parasitic",
                             "Drug resistance: antimicrobial",
                             "Drug resistance: antineoplastic",
                             
                             "Protein families: metabolism",
                             "Protein families: genetic information processing",
                             "Protein families: signaling and cellular processes",
                             "Viral protein families",
                             
                             "Unclassified: metabolism",
                             "Unclassified: genetic information processing",
                             "Unclassified: signaling and cellular processes",
                             "Poorly characterized")) +
  guides(fill = guide_legend(reverse = F,override.aes = list(size = 1), 
                             order = 3, ncol = 2)) +
  new_scale("fill") +
  geom_tile(data = plotdata_meta2, aes(x = sample, y = case_y,
                                       fill = factor(C_CASE2, levels = c(1,0)),
                                       height = 0.05, width = 0.875)) +
  scale_fill_manual(values = c("red3", "black"),
                    labels = c("Preterm case", "Term control"),
                    "Birth outcome") +
  guides(fill = guide_legend(reverse = F, override.aes = list(size = 1), order = 2)) +
  scale_x_discrete(limits = rev, labels = NULL) +
  labs(x = "", y = "KO functional category relative abundance") +
  geom_rect(aes(xmin = 0.5, xmax = 72.5,
                ymin = 1.0275, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 0.5, xmax = 0.35, 
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 24.43, xmax = 24.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 55.43, xmax = 55.58,
                ymin = 1.01, ymax = 1.03),
            fill = "black") + 
  geom_rect(aes(xmin = 72.35, xmax = 72.5,
                ymin = 1.01, ymax = 1.03),
            fill = "black") +
  annotate(geom = "text", x = 64, y = 1.06,
           label = expression(paste(italic("L. crispatus"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 40, y = 1.06,
           label = expression(paste(italic("L. iners"), " dominated")),
           color = "black", size = 2.5, angle = 90) +
  annotate(geom = "text", x = 12.5, y = 1.06,
           label = "Mixed", 
           color = "black", size = 2.5, angle = 90) +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text = element_text(size = 5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stacked_bar_meta2



### patchwork final plot
(dendrogram | plot_spacer | stacked_bar_meta2) +
  plot_layout(widths = c(1.5, -1.3, 9.5))






###################################################
### OBSERVED METAGENOME HIERARCHICAL CLUSTERING ###
###################################################



### generate transposed counts dataframe for observed metagenome data 

# generate matrix of KO counts
counts_meta <- data.frame(t(data[, 233:2263]))
# counts_meta should have 2031 obs of 72 vars



### generate phyloseq object for observed metagenome data

# generate KO table (equivalent of OTU table)
phylo_ko <- otu_table(counts_meta, taxa_are_rows = T)

# generate phyloseq object
phylo_meta <- phyloseq(phylo_ko)



### hierarchical clustering

# run clValid to identify optimal number of clusters
valid_meta <- clValid(counts_meta, 
                      nClust = 2:10, 
                      clMethods = "hierarchical", 
                      validation = "internal", 
                      maxitems = 2031, 
                      metric = "euclidean", 
                      method = "ward",
                      neighbSize = 10)
summary(valid_meta)
optimalScores(valid_meta)
# using 2 clusters based on internal validation statistics

# generate JSD matrix
jsd_meta <- phyloseq::distance(phylo_ko, method = "jsd")
# jsd_meta should have 2556 elements

# hierarchical clustering to generate clusters and dendrogram sorted by average distance
hclust_meta <- dendsort(hclust(jsd_meta,
                               method = "ward.D",
                               members = NULL),
                        type = "average")

# plot and save dendrogram, add rectangles around clusters, save order of leaf labels
dendro_meta <- as.dendrogram(hclust_meta)
plot(hclust_meta)
rect.hclust(hclust_meta, 
            k = 2,
            border = "red")
leaf_meta <- dendro_data(dendro_meta)$labels[3]

# generate cluster variable, add to data, generate factor version
cluster_meta <- cutree(hclust_meta,
                       k = 2)
cluster_meta <- cluster_meta[order(names(cluster_meta))] 
data <- data[order(data$sample),] 
data <- cbind(data, cluster_meta)
data$cluster_meta_f <- factor(data$cluster_meta)
# data should have 72 obs of 4433 vars






####################################################
### ESTIMATE OBSERVED METAGENOME ALPHA DIVERSITY ###
####################################################



### estimate alpha diversity
# estimate_richness only takes integers as input
# metagenome count data are non-integers as an artifact of processing
# coercing metagenome KO counts to integer before estimating alpha diversity
alpha_meta <- estimate_richness(transform_sample_counts(phylo_meta, as.integer),
                                split = T,
                                measures = c("Observed",
                                             "Shannon",
                                             "Simpson",
                                             "InvSimpson"))
# alpha_meta should have 72 obs of 4 vars

# add alpha_meta to data
rownames(alpha_meta) <- data$sample
alpha_meta <- alpha_meta[order(rownames(alpha_meta)),] 
data <- cbind(data, alpha_meta)
# data should have 72 obs of 4437 vars






##################################################
### PREDICTED METAGENOME PREP AND DESCRIPTIVES ###
##################################################



### generate vectors of predicted gene families that were observed ("true positives") and were not observed 
### ("false positives")
### generate vectors of observed gene families that were not predicted ("false negatives")
### generate vector of gene families that were observed and predicted by PICRUSt2 and Tax4Fun2

obs_feat <- colnames(data[, 2265:4295])
# obs_feat should have 2031 elements

pi_feat <- colnames(pi[, 1:(ncol(pi) - 2)])
# pi_feat should have 5882 elements
pi_fp <- pi_feat[pi_feat %!in% obs_feat]
# pi_fp should have 4379 elements
pi_tp <- obs_feat[obs_feat %in% pi_feat]
# pi_tp should have 1503 elements
pi_fn <- obs_feat[obs_feat %!in% pi_feat]
# pi_fn should have 528 elements

tax_feat <- colnames(tax[, 1:(ncol(tax) - 2)])
# tax_feat should have 7049 elements
tax_fp <- tax_feat[tax_feat %!in% obs_feat]
# tax_fp should have 5543 elements
tax_tp <- obs_feat[obs_feat %in% tax_feat]
# tax_tp should have 1506 elements
tax_fn <- obs_feat[obs_feat %!in% tax_feat]
# tax_fn should have 525 elements

tp <- pi_tp[pi_tp %in% tax_tp]
# tp should have 1490 elements



### pull weighted NSTI values from pi dataframe, percent of reads discarded from tax dataframe
### merge cluster_f into both, save in own objects with sample IDs

wnsti <- pi[, colnames(pi) %in% c("sample", "wnsti")]
wnsti <- left_join(wnsti, phylo_sample[, colnames(phylo_sample) %in% c("sample", "cluster_f")])
# wnsti should have 72 obs of 3 vars

discard <- tax[, colnames(tax) %in% c("sample", "reads")]
discard <- left_join(discard, phylo_sample[, colnames(phylo_sample) %in% c("sample", "cluster_f")])
# discard should have 72 obs of 3 vars



# summarize PICRUSt2 wNSTI, Tax4Fun2 % reads discarded, overall and by cluster

summary(wnsti$wnsti)
summary(wnsti$wnsti[wnsti$cluster_f == "crisp"])
summary(wnsti$wnsti[wnsti$cluster_f == "iners"])
summary(wnsti$wnsti[wnsti$cluster_f == "mixed"])

summary(discard$reads)
summary(discard$reads[discard$cluster_f == "crisp"])
summary(discard$reads[discard$cluster_f == "iners"])
summary(discard$reads[discard$cluster_f == "mixed"])



### generate dataframe of observed gene family relative abundances
### restrict this dataframe, pi, tax to gene families observed and predicted by both PICRUSt2 and Tax4Fun2 (tp)
### merge in case status, cluster

obs <- data[, c(2265:4295, 1)]
obs <- obs [, colnames(obs) %in% c("sample", tp)]
# obs should have 72 obs of 1491 vars

pi <- pi[, colnames(pi) %in% c("sample", tp)]
# pi should have 72 obs of 1491 vars

tax <- tax[, colnames(tax) %in% c("sample", tp)]
# tax should have 72 obs of 1491 vars

# order variables in obs, pi, tax in same order
obs <- obs[, order(colnames(obs))]
pi <- pi[, order(colnames(pi))]
tax <- tax[, order(colnames(tax))]

# order observations in obs, pi, tax in same order
obs <- obs[order(obs$sample),]
pi <- pi[order(pi$sample),]
tax <- tax[order(tax$sample),]

# merge in case status and cluster
obs <- left_join(obs, data[, colnames(data) %in% c("sample", "C_CASE2", "cluster")])
pi <- left_join(pi, data[, colnames(data) %in% c("sample", "C_CASE2", "cluster")])
tax <- left_join(tax, data[, colnames(data) %in% c("sample", "C_CASE2", "cluster")])
# obs, pi, and tax should all have 72 obs of 1493 vars






#######################################################################
### METAGENOME INFERENCE COMPARISON, UNPERMUTED, CLUSTER STRATIFIED ###
#######################################################################



### generate versions of obs, pi, tax restricted by cluster

# restrict by cluster
obs_1 <- obs[obs$cluster == 1,]
obs_2 <- obs[obs$cluster == 2,]
obs_3 <- obs[obs$cluster == 3,]

pi_1 <- pi[pi$cluster == 1,]
pi_2 <- pi[pi$cluster == 2,]
pi_3 <- pi[pi$cluster == 3,]

tax_1 <- tax[tax$cluster == 1,]
tax_2 <- tax[tax$cluster == 2,]
tax_3 <- tax[tax$cluster == 3,]
# all _1 should have 24 obs of 1493 vars
# all _2 should have 31 obs of 1493 vars
# all _3 should have 17 obs of 1493 vars

# generate vector of column names that sum to 0 in obs_1, pi_1, tax_1
obs_1_miss <- colnames(obs_1[, colSums(obs_1[, 1:1490]) == 0])
# obs_1_miss should have 94 elements
pi_1_miss <- colnames(pi_1[, colSums(pi_1[, 1:1490]) == 0])
# pi_1_miss should have 4 elements
tax_1_miss <- colnames(tax_1[, colSums(tax_1[, 1:1490]) == 0])
# tax_1_miss should have 3 elements

# generate vector of all column names that sum to 0 in either obs_1 or pi_1
# drop these columns
obs_pi_1_miss <- unique(c(obs_1_miss, pi_1_miss))
# obs_pi_1_miss should have 98 elements
obs_pi_1 <- obs_1[, colnames(obs_1) %!in% obs_pi_1_miss]
pi_1 <- pi_1[, colnames(pi_1) %!in% obs_pi_1_miss]
# both should have 24 obs of 1395 vars

# generate vector of all column names that sum to 0 in either obs_1 or tax_1
# drop these columns
obs_tax_1_miss <- unique(c(obs_1_miss, tax_1_miss))
# obs_tax_1_miss should have 96 elements
obs_tax_1 <- obs_1[, colnames(obs_1) %!in% obs_tax_1_miss]
tax_1 <- tax_1[, colnames(tax_1) %!in% obs_tax_1_miss]
# both should have 24 obs of 1397 vars



# generate vector of column names that sum to 0 in obs_2, pi_2, tax_2
check <- obs_2[, 1:1490]
check[nrow(check)+1,] <- colSums(check)
check <- as.numeric(t(check[nrow(check),]))
names(check) <- colnames(obs_2[, 1:1490])
obs_2_miss <- names(check[check == 0])
# obs_2_miss should have 356 elements
pi_2_miss <- colnames(pi_2[, colSums(pi_2[, 1:1490]) == 0])
# pi_2_miss should have 10 elements
tax_2_miss <- colnames(tax_2[, colSums(tax_2[, 1:1490]) == 0])
# tax_2_miss should have 7 elements

# generate vector of all column names that sum to 0 in either obs_2 or pi_2
# drop these columns
obs_pi_2_miss <- unique(c(obs_2_miss, pi_2_miss))
# obs_tax_2_miss should have 361 elements
obs_pi_2 <- obs_2[, colnames(obs_2) %!in% obs_pi_2_miss]
pi_2 <- pi_2[, colnames(pi_2) %!in% obs_pi_2_miss]
# both should have 31 obs of 1132 vars

# generate vector of all column names that sum to 0 in either obs_2 or tax_2
# drop these columns
obs_tax_2_miss <- unique(c(obs_2_miss, tax_2_miss))
# obs_tax_2_miss should have 357 elements
obs_tax_2 <- obs_2[, colnames(obs_2) %!in% obs_tax_2_miss]
tax_2 <- tax_2[, colnames(tax_2) %!in% obs_tax_2_miss]
# both should have 31 obs of 1136 vars



# generate vector of column names that sum to 0 in obs_3, pi_3, tax_3
check <- obs_3[, 1:1490]
check[nrow(check)+1,] <- colSums(check)
check <- as.numeric(t(check[nrow(check),]))
names(check) <- colnames(obs_3[, 1:1490])
obs_3_miss <- names(check[check == 0])
# obs_3_miss should have 468 elements
pi_3_miss <- colnames(pi_3[, colSums(pi_3[, 1:1490]) == 0])
# pi_3_miss should have 18 elements
tax_3_miss <- colnames(tax_3[, colSums(tax_3[, 1:1490]) == 0])
# tax_3_miss should have 14 elements

# generate vector of all column names that sum to 0 in either obs_3 or pi_3
# drop these columns
obs_pi_3_miss <- unique(c(obs_3_miss, pi_3_miss))
# obs_tax_3_miss should have 474 elements
obs_pi_3 <- obs_3[, colnames(obs_3) %!in% obs_pi_3_miss]
pi_3 <- pi_3[, colnames(pi_3) %!in% obs_pi_3_miss]
# both should have 17 obs of 1019 vars

# generate vector of all column names that sum to 0 in either obs_3 or tax_3
# drop these columns
obs_tax_3_miss <- unique(c(obs_3_miss, tax_3_miss))
# obs_tax_3_miss should have 470 elements
obs_tax_3 <- obs_3[, colnames(obs_3) %!in% obs_tax_3_miss]
tax_3 <- tax_3[, colnames(tax_3) %!in% obs_tax_3_miss]
# both should have 17 obs of 1023 vars



### estimate spearman correlation coefficients between observed and predicted KO relative abundances
### overall and restricted by cluster, unpermuted

# for PICRUSt2 predictions
obs_pi_corr <- vector()
obs_pi_corr_1 <- vector()
obs_pi_corr_2 <- vector()
obs_pi_corr_3 <- vector()

for(each in tp){
  corr <- cor(obs[, colnames(obs) == each],
              pi[, colnames(pi) == each],
              method = "spearman")
  obs_pi_corr <- append(obs_pi_corr, corr)
}

for(each in tp[tp %!in% obs_pi_1_miss]){
  corr1 <- cor(obs[obs$cluster == 1, colnames(obs) == each],
               pi[pi$cluster == 1, colnames(pi) == each],
               method = "spearman")
  obs_pi_corr_1 <- append(obs_pi_corr_1, corr1)
}

for(each in tp[tp %!in% obs_pi_2_miss]){
  corr2 <- cor(obs[obs$cluster == 2, colnames(obs) == each],
               pi[pi$cluster == 2, colnames(pi) == each],
               method = "spearman")
  obs_pi_corr_2 <- append(obs_pi_corr_2, corr2)
}

for(each in tp[tp %!in% obs_pi_3_miss]){
  corr3 <- cor(obs[obs$cluster == 3, colnames(obs) == each],
               pi[pi$cluster == 3, colnames(pi) == each],
               method = "spearman")
  obs_pi_corr_3 <- append(obs_pi_corr_3, corr3)
}

summary(obs_pi_corr)
#obs_pi_corr should have 1490 elements, median 0.20
summary(obs_pi_corr_1)
#obs_pi_corr_1 should have 1392 elements, median 0.21
summary(obs_pi_corr_2)
#obs_pi_corr_2 should have 1129 elements, median 0.06
summary(obs_pi_corr_3)
#obs_pi_corr_3 should have 1016 elements, median 0.24



# for Tax4Fun2 predictions
obs_tax_corr <- vector()
obs_tax_corr_1 <- vector()
obs_tax_corr_2 <- vector()
obs_tax_corr_3 <- vector()

for(each in tp){
  corr <- cor(obs[, colnames(obs) == each],
              tax[, colnames(tax) == each],
              method = "spearman")
  obs_tax_corr <- append(obs_tax_corr, corr)
}

for(each in tp[tp %!in% obs_tax_1_miss]){
  corr1 <- cor(obs[obs$cluster == 1, colnames(obs) == each],
               tax[tax$cluster == 1, colnames(tax) == each],
               method = "spearman")
  obs_tax_corr_1 <- append(obs_tax_corr_1, corr1)
}

for(each in tp[tp %!in% obs_tax_2_miss]){
  corr2 <- cor(obs[obs$cluster == 2, colnames(obs) == each],
               tax[tax$cluster == 2, colnames(tax) == each],
               method = "spearman")
  obs_tax_corr_2 <- append(obs_tax_corr_2, corr2)
}

for(each in tp[tp %!in% obs_tax_3_miss]){
  corr3 <- cor(obs[obs$cluster == 3, colnames(obs) == each],
               tax[tax$cluster == 3, colnames(tax) == each],
               method = "spearman")
  obs_tax_corr_3 <- append(obs_tax_corr_3, corr3)
}

summary(obs_tax_corr)
#obs_tax_corr should have 1490 elements, median 0.22
summary(obs_tax_corr_1)
#obs_tax_corr_1 should have 1394 elements, median 0.24
summary(obs_tax_corr_2)
#obs_tax_corr_2 should have 1133 elements, median 0.11
summary(obs_tax_corr_3)
#obs_tax_corr_3 should have 1020 elements, median 0.25



### inference using wilcoxon test comparing observed and predicted KO relative abundances between
### preterm birth cases and term birth controls, overall and stratified by cluster

# inference for observed KO relative abundances
obs_p <- vector()
obs_pi_p_1 <- vector()
obs_pi_p_1_gene <- vector()
obs_pi_p_2 <- vector()
obs_pi_p_2_gene <- vector()
obs_pi_p_3 <- vector()
obs_pi_p_3_gene <- vector()
obs_tax_p_1 <- vector()
obs_tax_p_1_gene <- vector()
obs_tax_p_2 <- vector()
obs_tax_p_2_gene <- vector()
obs_tax_p_3 <- vector()
obs_tax_p_3_gene <- vector()
for(each in tp){
  p <- wilcox.test(obs[obs$C_CASE2 == 1, colnames(obs) == each], 
                   obs[obs$C_CASE2 == 0, colnames(obs) == each])[[3]][1]
  p_t <- log10(p)*
    sign(mean(obs[obs$C_CASE2 == 1, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0, colnames(obs) == each]))
  obs_p <- append(obs_p, p_t)
}
for(each in tp[tp %!in% obs_pi_1_miss]){
  p1 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 1, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 1, colnames(obs) == each])[[3]][1]
  p_t1 <- log10(p1)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 1, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 1, colnames(obs) == each]))
  obs_pi_p_1 <- append(obs_pi_p_1, p_t1)
  obs_pi_p_1_gene <- append(obs_pi_p_1_gene, each)
}
for(each in tp[tp %!in% obs_pi_2_miss]){  
  p2 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 2, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 2, colnames(obs) == each])[[3]][1]
  p_t2 <- log10(p2)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 2, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 2, colnames(obs) == each]))
  obs_pi_p_2 <- append(obs_pi_p_2, p_t2)
  obs_pi_p_2_gene <- append(obs_pi_p_2_gene, each)
}
for(each in tp[tp %!in% obs_pi_3_miss]){  
  p3 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 3, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 3, colnames(obs) == each])[[3]][1]
  p_t3 <- log10(p3)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 3, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 3, colnames(obs) == each]))
  obs_pi_p_3 <- append(obs_pi_p_3, p_t3)
  obs_pi_p_3_gene <- append(obs_pi_p_3_gene, each)
}
for(each in tp[tp %!in% obs_tax_1_miss]){
  p1 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 1, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 1, colnames(obs) == each])[[3]][1]
  p_t1 <- log10(p1)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 1, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 1, colnames(obs) == each]))
  obs_tax_p_1 <- append(obs_tax_p_1, p_t1)
  obs_tax_p_1_gene <- append(obs_tax_p_1_gene, each)
}
for(each in tp[tp %!in% obs_tax_2_miss]){  
  p2 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 2, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 2, colnames(obs) == each])[[3]][1]
  p_t2 <- log10(p2)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 2, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 2, colnames(obs) == each]))
  obs_tax_p_2 <- append(obs_tax_p_2, p_t2)
  obs_tax_p_2_gene <- append(obs_tax_p_2_gene, each)
}
for(each in tp[tp %!in% obs_tax_3_miss]){  
  p3 <- wilcox.test(obs[obs$C_CASE2 == 1 & obs$cluster == 3, colnames(obs) == each], 
                    obs[obs$C_CASE2 == 0 & obs$cluster == 3, colnames(obs) == each])[[3]][1]
  p_t3 <- log10(p3)*
    sign(mean(obs[obs$C_CASE2 == 1  & obs$cluster == 3, colnames(obs) == each])-
           mean(obs[obs$C_CASE2 == 0 & obs$cluster == 3, colnames(obs) == each]))
  obs_tax_p_3 <- append(obs_tax_p_3, p_t3)
  obs_tax_p_3_gene <- append(obs_tax_p_3_gene, each)
}
summary(obs_p)
# obs_p should have 1490 elements
summary(obs_pi_p_1)
# obs_pi_p_1 should have 1392 elements
summary(obs_pi_p_2)
# obs_pi_p_2 should have 1129 elements
summary(obs_pi_p_3)
# obs_pi_p_3 should have 1016 elements
summary(obs_tax_p_1)
# obs_tax_p_1 should have 1394 elements
summary(obs_tax_p_2)
# obs_tax_p_2 should have 1133 elements
summary(obs_tax_p_3)
# obs_tax_p_3 should have 1020 elements



# inference for PICRUSt2 KO relative abundances
pi_p <- vector()
pi_p_1 <- vector()
pi_p_1_gene <- vector()
pi_p_2 <- vector()
pi_p_2_gene <- vector()
pi_p_3 <- vector()
pi_p_3_gene <- vector()
for(each in tp){
  p <- wilcox.test(pi[pi$C_CASE2 == 1, colnames(pi) == each], 
                 pi[pi$C_CASE2 == 0, colnames(pi) == each])[[3]][1]
  p_t <- log10(p)*
    sign(mean(pi[pi$C_CASE2 == 1, colnames(pi) == each])-
           mean(pi[pi$C_CASE2 == 0, colnames(pi) == each]))
  pi_p <- append(pi_p, p_t)
}
for(each in tp[tp %!in% obs_pi_1_miss]){
  p1 <- wilcox.test(pi[pi$C_CASE2 == 1 & pi$cluster == 1, colnames(pi) == each], 
                  pi[pi$C_CASE2 == 0 & pi$cluster == 1, colnames(pi) == each])[[3]][1]
  p_t1 <- log10(p1)*
    sign(mean(pi[pi$C_CASE2 == 1  & pi$cluster == 1, colnames(pi) == each])-
           mean(pi[pi$C_CASE2 == 0 & pi$cluster == 1, colnames(pi) == each]))
  pi_p_1 <- append(pi_p_1, p_t1)
  pi_p_1_gene <- append(pi_p_1_gene, each)
}
for(each in tp[tp %!in% obs_pi_2_miss]){
  p2 <- wilcox.test(pi[pi$C_CASE2 == 1 & pi$cluster == 2, colnames(pi) == each], 
                  pi[pi$C_CASE2 == 0 & pi$cluster == 2, colnames(pi) == each])[[3]][1]
  p_t2 <- log10(p2)*
    sign(mean(pi[pi$C_CASE2 == 1  & pi$cluster == 2, colnames(pi) == each])-
           mean(pi[pi$C_CASE2 == 0 & pi$cluster == 2, colnames(pi) == each]))
  pi_p_2 <- append(pi_p_2, p_t2)
  pi_p_2_gene <- append(pi_p_2_gene, each)
}
for(each in tp[tp %!in% obs_pi_3_miss]){
  p3 <- wilcox.test(pi[pi$C_CASE2 == 1 & pi$cluster == 3, colnames(pi) == each], 
                  pi[pi$C_CASE2 == 0 & pi$cluster == 3, colnames(pi) == each])[[3]][1]
  p_t3 <- log10(p3)*
    sign(mean(pi[pi$C_CASE2 == 1  & pi$cluster == 3, colnames(pi) == each])-
           mean(pi[pi$C_CASE2 == 0 & pi$cluster == 3, colnames(pi) == each]))
  pi_p_3 <- append(pi_p_3, p_t3)
  pi_p_3_gene <- append(pi_p_3_gene, each)
}
summary(pi_p)
# pi_p should have 1490 elements
summary(pi_p_1)
# pi_p_1 should have 1392 elements
summary(pi_p_2)
# pi_p_2 should have 1129 elements
summary(pi_p_3)
# pi_p_3 should have 1016 elements



# inference for Tax4Fun2 KO relative abundances
tax_p <- vector()
tax_p_1 <- vector()
tax_p_1_gene <- vector()
tax_p_2 <- vector()
tax_p_2_gene <- vector()
tax_p_3 <- vector()
tax_p_3_gene <- vector()
for(each in tp){
  p <- wilcox.test(tax[tax$C_CASE2 == 1, colnames(tax) == each], 
                   tax[tax$C_CASE2 == 0, colnames(tax) == each])[[3]][1]
  p_t <- log10(p)*
    sign(mean(tax[tax$C_CASE2 == 1, colnames(tax) == each])-
           mean(tax[tax$C_CASE2 == 0, colnames(tax) == each]))
  tax_p <- append(tax_p, p_t)
}
for(each in tp[tp %!in% obs_tax_1_miss]){
  p1 <- wilcox.test(tax[tax$C_CASE2 == 1 & tax$cluster == 1, colnames(tax) == each], 
                    tax[tax$C_CASE2 == 0 & tax$cluster == 1, colnames(tax) == each])[[3]][1]
  p_t1 <- log10(p1)*
    sign(mean(tax[tax$C_CASE2 == 1  & tax$cluster == 1, colnames(tax) == each])-
           mean(tax[tax$C_CASE2 == 0 & tax$cluster == 1, colnames(tax) == each]))
  tax_p_1 <- append(tax_p_1, p_t1)
  tax_p_1_gene <- append(tax_p_1_gene, each)
}
for(each in tp[tp %!in% obs_tax_2_miss]){
  p2 <- wilcox.test(tax[tax$C_CASE2 == 1 & tax$cluster == 2, colnames(tax) == each], 
                    tax[tax$C_CASE2 == 0 & tax$cluster == 2, colnames(tax) == each])[[3]][1]
  p_t2 <- log10(p2)*
    sign(mean(tax[tax$C_CASE2 == 1  & tax$cluster == 2, colnames(tax) == each])-
           mean(tax[tax$C_CASE2 == 0 & tax$cluster == 2, colnames(tax) == each]))
  tax_p_2 <- append(tax_p_2, p_t2)
  tax_p_2_gene <- append(tax_p_2_gene, each)
}
for(each in tp[tp %!in% obs_tax_3_miss]){
  p3 <- wilcox.test(tax[tax$C_CASE2 == 1 & tax$cluster == 3, colnames(tax) == each], 
                    tax[tax$C_CASE2 == 0 & tax$cluster == 3, colnames(tax) == each])[[3]][1]
  p_t3 <- log10(p3)*
    sign(mean(tax[tax$C_CASE2 == 1  & tax$cluster == 3, colnames(tax) == each])-
           mean(tax[tax$C_CASE2 == 0 & tax$cluster == 3, colnames(tax) == each]))
  tax_p_3 <- append(tax_p_3, p_t3)
  tax_p_3_gene <- append(tax_p_3_gene, each)
}
summary(tax_p)
# tax_p should have 1490 elements
summary(tax_p_1)
# tax_p_1 should have 1394 elements
summary(tax_p_2)
# tax_p_2 should have 1133 elements
summary(tax_p_3)
# tax_p_3 should have 1020 elements



### estimate spearman correlation coefficients between transformed p values
(obs_pi_p_corr <- cor(obs_p, pi_p, method = "spearman"))
(obs_pi_p_corr_1 <- cor(obs_pi_p_1, pi_p_1, method = "spearman"))
(obs_pi_p_corr_2 <- cor(obs_pi_p_2, pi_p_2, method = "spearman"))
(obs_pi_p_corr_3 <- cor(obs_pi_p_3, pi_p_3, method = "spearman"))

(obs_tax_p_corr <- cor(obs_p, tax_p, method = "spearman"))
(obs_tax_p_corr_1 <- cor(obs_tax_p_1, tax_p_1, method = "spearman"))
(obs_tax_p_corr_2 <- cor(obs_tax_p_2, tax_p_2, method = "spearman"))
(obs_tax_p_corr_3 <- cor(obs_tax_p_3, tax_p_3, method = "spearman"))



### generate dataframe of transformed p values for plotting

obs_pi_p_1_df <- data.frame(obs_pi_p_1)
obs_pi_p_1_df$gene <- obs_pi_p_1_gene
colnames(obs_pi_p_1_df) <- c("p_obs", "gene")
obs_pi_p_1_df$cluster <- 1

obs_pi_p_2_df <- data.frame(obs_pi_p_2)
obs_pi_p_2_df$gene <- obs_pi_p_2_gene
colnames(obs_pi_p_2_df) <- c("p_obs", "gene")
obs_pi_p_2_df$cluster <- 2

obs_pi_p_3_df <- data.frame(obs_pi_p_3)
obs_pi_p_3_df$gene <- obs_pi_p_3_gene
colnames(obs_pi_p_3_df) <- c("p_obs", "gene")
obs_pi_p_3_df$cluster <- 3

obs_tax_p_1_df <- data.frame(obs_tax_p_1)
obs_tax_p_1_df$gene <- obs_tax_p_1_gene
colnames(obs_tax_p_1_df) <- c("p_obs", "gene")
obs_tax_p_1_df$cluster <- 1

obs_tax_p_2_df <- data.frame(obs_tax_p_2)
obs_tax_p_2_df$gene <- obs_tax_p_2_gene
colnames(obs_tax_p_2_df) <- c("p_obs", "gene")
obs_tax_p_2_df$cluster <- 2

obs_tax_p_3_df <- data.frame(obs_tax_p_3)
obs_tax_p_3_df$gene <- obs_tax_p_3_gene
colnames(obs_tax_p_3_df) <- c("p_obs", "gene")
obs_tax_p_3_df$cluster <- 3

pi_p_1_df <- data.frame(pi_p_1)
pi_p_1_df$gene <- pi_p_1_gene
colnames(pi_p_1_df) <- c("p_pred", "gene")
pi_p_1_df$cluster <- 1
pi_p_1_df$method <- "PICRUSt2"

pi_p_2_df <- data.frame(pi_p_2)
pi_p_2_df$gene <- pi_p_2_gene
colnames(pi_p_2_df) <- c("p_pred", "gene")
pi_p_2_df$cluster <- 2
pi_p_2_df$method <- "PICRUSt2"

pi_p_3_df <- data.frame(pi_p_3)
pi_p_3_df$gene <- pi_p_3_gene
colnames(pi_p_3_df) <- c("p_pred", "gene")
pi_p_3_df$cluster <- 3
pi_p_3_df$method <- "PICRUSt2"

tax_p_1_df <- data.frame(tax_p_1)
tax_p_1_df$gene <- tax_p_1_gene
colnames(tax_p_1_df) <- c("p_pred", "gene")
tax_p_1_df$cluster <- 1
tax_p_1_df$method <- "Tax4Fun2"

tax_p_2_df <- data.frame(tax_p_2)
tax_p_2_df$gene <- tax_p_2_gene
colnames(tax_p_2_df) <- c("p_pred", "gene")
tax_p_2_df$cluster <- 2
tax_p_2_df$method <- "Tax4Fun2"

tax_p_3_df <- data.frame(tax_p_3)
tax_p_3_df$gene <- tax_p_3_gene
colnames(tax_p_3_df) <- c("p_pred", "gene")
tax_p_3_df$cluster <- 3
tax_p_3_df$method <- "Tax4Fun2"

obs_p <- rbind(obs_pi_p_1_df, obs_pi_p_2_df, obs_pi_p_3_df, 
               obs_tax_p_1_df, obs_tax_p_2_df, obs_tax_p_3_df) 
# obs_p should have 7084 obs of 3 vars
pred_p <- rbind(pi_p_1_df, pi_p_2_df, pi_p_3_df, 
                tax_p_1_df, tax_p_2_df, tax_p_3_df)
# pred_p should have 7084 obs of 4 vars

plot_pvp_cluster <- left_join(obs_p, pred_p)
plot_pvp_cluster$cluster <- factor(plot_pvp_cluster$cluster,
                                   levels = c(3, 2, 1))
plot_pvp_cluster$method <- factor(plot_pvp_cluster$method)
# plot_pvp_cluster should have 14152 obs of 5 vars






#####################################################################
### METAGENOME INFERENCE COMPARISON, PERMUTED, CLUSTER STRATIFIED ###
#####################################################################



### generate permuted datasets overall and stratified by cluster

set.seed(1)



# for obs and pi
# generate lists with 100 permuted versions of obs, pi, both overall
# generate lists with 100 permuted versions of obs, pi, both restricted to samples in cluster 1, genes obs/pred by both obs and pi in cluster 1
# generate lists with 100 permuted versions of obs, pi, both restricted to samples in cluster 2, genes obs/pred by both obs and pi in cluster 2
# generate lists with 100 permuted versions of obs, pi, both restricted to samples in cluster 3, genes obs/pred by both obs and pi in cluster 3
obs_pi_perms <- list()
pi_obs_perms <- list()
obs_pi_perms_1 <- list()
pi_obs_perms_1 <- list()
obs_pi_perms_2 <- list()
pi_obs_perms_2 <- list()
obs_pi_perms_3 <- list()
pi_obs_perms_3 <- list()
for(x in 1:100){
  obs_pi_perm <- obs
  pi_obs_perm <- pi
  obs_pi_perm_1 <- obs_1
  pi_obs_perm_1 <- pi_1
  obs_pi_perm_2 <- obs_2
  pi_obs_perm_2 <- pi_2
  obs_pi_perm_3 <- obs_3
  pi_obs_perm_3 <- pi_3
  for(each in tp){
    obs_pi_perm[, colnames(obs_pi_perm) == each] <- sample(obs_pi_perm[, colnames(obs_pi_perm) == each])
    pi_obs_perm[, colnames(pi_obs_perm) == each] <- sample(pi_obs_perm[, colnames(pi_obs_perm) == each])
  }
  for(each in tp[tp %!in% obs_pi_1_miss]){
    obs_pi_perm_1[, colnames(obs_pi_perm_1) == each] <- sample(obs_pi_perm_1[, colnames(obs_pi_perm_1) == each])
    pi_obs_perm_1[, colnames(pi_obs_perm_1) == each] <- sample(pi_obs_perm_1[, colnames(pi_obs_perm_1) == each])
  }
  for(each in tp[tp %!in% obs_pi_2_miss]){
    obs_pi_perm_2[, colnames(obs_pi_perm_2) == each] <- sample(obs_pi_perm_2[, colnames(obs_pi_perm_2) == each])
    pi_obs_perm_2[, colnames(pi_obs_perm_2) == each] <- sample(pi_obs_perm_2[, colnames(pi_obs_perm_2) == each])
  }
  for(each in tp[tp %!in% obs_pi_3_miss]){
    obs_pi_perm_3[, colnames(obs_pi_perm_3) == each] <- sample(obs_pi_perm_3[, colnames(obs_pi_perm_3) == each])
    pi_obs_perm_3[, colnames(pi_obs_perm_3) == each] <- sample(pi_obs_perm_3[, colnames(pi_obs_perm_3) == each])
  }
  obs_pi_perms <- c(obs_pi_perms, list(obs_pi_perm))
  pi_obs_perms <- c(pi_obs_perms, list(pi_obs_perm))
  obs_pi_perms_1 <- c(obs_pi_perms_1, list(obs_pi_perm_1))
  pi_obs_perms_1 <- c(pi_obs_perms_1, list(pi_obs_perm_1))
  obs_pi_perms_2 <- c(obs_pi_perms_2, list(obs_pi_perm_2))
  pi_obs_perms_2 <- c(pi_obs_perms_2, list(pi_obs_perm_2))
  obs_pi_perms_3 <- c(obs_pi_perms_3, list(obs_pi_perm_3))
  pi_obs_perms_3 <- c(pi_obs_perms_3, list(pi_obs_perm_3))
}



# for obs and tax
# generate lists with 100 permuted versions of obs, tax, both overall
# generate lists with 100 permuted versions of obs, tax, both restricted to samples in cluster 1, genes obs/pred by both obs and tax in cluster 1
# generate lists with 100 permuted versions of obs, tax, both restricted to samples in cluster 2, genes obs/pred by both obs and tax in cluster 2
# generate lists with 100 permuted versions of obs, tax, both restricted to samples in cluster 3, genes obs/pred by both obs and tax in cluster 3
obs_tax_perms <- list()
tax_obs_perms <- list()
obs_tax_perms_1 <- list()
tax_obs_perms_1 <- list()
obs_tax_perms_2 <- list()
tax_obs_perms_2 <- list()
obs_tax_perms_3 <- list()
tax_obs_perms_3 <- list()
for(x in 1:100){
  obs_tax_perm <- obs
  tax_obs_perm <- tax
  obs_tax_perm_1 <- obs_1
  tax_obs_perm_1 <- tax_1
  obs_tax_perm_2 <- obs_2
  tax_obs_perm_2 <- tax_2
  obs_tax_perm_3 <- obs_3
  tax_obs_perm_3 <- tax_3
  for(each in tp){
    obs_tax_perm[, colnames(obs_tax_perm) == each] <- sample(obs_tax_perm[, colnames(obs_tax_perm) == each])
    tax_obs_perm[, colnames(tax_obs_perm) == each] <- sample(tax_obs_perm[, colnames(tax_obs_perm) == each])
  }
  for(each in tp[tp %!in% obs_tax_1_miss]){
    obs_tax_perm_1[, colnames(obs_tax_perm_1) == each] <- sample(obs_tax_perm_1[, colnames(obs_tax_perm_1) == each])
    tax_obs_perm_1[, colnames(tax_obs_perm_1) == each] <- sample(tax_obs_perm_1[, colnames(tax_obs_perm_1) == each])
  }
  for(each in tp[tp %!in% obs_tax_2_miss]){
    obs_tax_perm_2[, colnames(obs_tax_perm_2) == each] <- sample(obs_tax_perm_2[, colnames(obs_tax_perm_2) == each])
    tax_obs_perm_2[, colnames(tax_obs_perm_2) == each] <- sample(tax_obs_perm_2[, colnames(tax_obs_perm_2) == each])
  }
  for(each in tp[tp %!in% obs_tax_3_miss]){
    obs_tax_perm_3[, colnames(obs_tax_perm_3) == each] <- sample(obs_tax_perm_3[, colnames(obs_tax_perm_3) == each])
    tax_obs_perm_3[, colnames(tax_obs_perm_3) == each] <- sample(tax_obs_perm_3[, colnames(tax_obs_perm_3) == each])
  }
  obs_tax_perms <- c(obs_tax_perms, list(obs_tax_perm))
  tax_obs_perms <- c(tax_obs_perms, list(tax_obs_perm))
  obs_tax_perms_1 <- c(obs_tax_perms_1, list(obs_tax_perm_1))
  tax_obs_perms_1 <- c(tax_obs_perms_1, list(tax_obs_perm_1))
  obs_tax_perms_2 <- c(obs_tax_perms_2, list(obs_tax_perm_2))
  tax_obs_perms_2 <- c(tax_obs_perms_2, list(tax_obs_perm_2))
  obs_tax_perms_3 <- c(obs_tax_perms_3, list(obs_tax_perm_3))
  tax_obs_perms_3 <- c(tax_obs_perms_3, list(tax_obs_perm_3))
}



### generate correlations using permuted datasets

# generate vectors of spearman correlation coefficients between permuted gene relative abundance 
# across to 100 permuted datasets
obs_pi_perms_corr <- vector()
obs_pi_perms_corr_1 <- vector()
obs_pi_perms_corr_2 <- vector()
obs_pi_perms_corr_3 <- vector()
obs_tax_perms_corr <- vector()
obs_tax_perms_corr_1 <- vector()
obs_tax_perms_corr_2 <- vector()
obs_tax_perms_corr_3 <- vector()
for(x in 1:100){
  for(each in tp){
    corr <- cor(obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) == each], 
              pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) == each], method = "spearman")
    obs_pi_perms_corr <- append(obs_pi_perms_corr, corr)
    corr <- cor(obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) == each], 
              obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) == each], method = "spearman")
    obs_tax_perms_corr <- append(obs_tax_perms_corr, corr)
  }
  for(each in tp[tp %!in% obs_pi_1_miss]){
    corr <- cor(obs_pi_perms_1[[x]][, colnames(obs_pi_perms_1[[x]]) == each], 
              pi_obs_perms_1[[x]][, colnames(pi_obs_perms_1[[x]]) == each], method = "spearman")
    obs_pi_perms_corr_1 <- append(obs_pi_perms_corr_1, corr)
  }
  for(each in tp[tp %!in% obs_pi_2_miss]){
    corr <- cor(obs_pi_perms_2[[x]][, colnames(obs_pi_perms_2[[x]]) == each], 
              pi_obs_perms_2[[x]][, colnames(pi_obs_perms_2[[x]]) == each], method = "spearman")
    obs_pi_perms_corr_2 <- append(obs_pi_perms_corr_2, corr)
  }
  for(each in tp[tp %!in% obs_pi_3_miss]){
    corr <- cor(obs_pi_perms_3[[x]][, colnames(obs_pi_perms_3[[x]]) == each], 
              pi_obs_perms_3[[x]][, colnames(pi_obs_perms_3[[x]]) == each], method = "spearman")
    obs_pi_perms_corr_3 <- append(obs_pi_perms_corr_3, corr)
  }
  for(each in tp[tp %!in% obs_tax_1_miss]){
    corr <- cor(obs_tax_perms_1[[x]][, colnames(obs_tax_perms_1[[x]]) == each], 
              obs_pi_perms_1[[x]][, colnames(obs_pi_perms_1[[x]]) == each], method = "spearman")
    obs_tax_perms_corr_1 <- append(obs_tax_perms_corr_1, corr)
  }
  for(each in tp[tp %!in% obs_tax_2_miss]){
    corr <- cor(obs_tax_perms_2[[x]][, colnames(obs_tax_perms_2[[x]]) == each], 
              obs_pi_perms_2[[x]][, colnames(obs_pi_perms_2[[x]]) == each], method = "spearman")
    obs_tax_perms_corr_2 <- append(obs_tax_perms_corr_2, corr)
  }
  for(each in tp[tp %!in% obs_tax_3_miss]){
    corr <- cor(obs_tax_perms_3[[x]][, colnames(obs_tax_perms_3[[x]]) == each], 
              obs_pi_perms_3[[x]][, colnames(obs_pi_perms_3[[x]]) == each], method = "spearman")
    obs_tax_perms_corr_3 <- append(obs_tax_perms_corr_3, corr)
  }
}
summary(obs_pi_perms_corr)
summary(obs_pi_perms_corr_1)
summary(obs_pi_perms_corr_2)
summary(obs_pi_perms_corr_3)
summary(obs_tax_perms_corr)
summary(obs_tax_perms_corr_1)
summary(obs_tax_perms_corr_2)
summary(obs_tax_perms_corr_3)



### generate lists of p value for wilcoxon test for each permuted gene relative abundance across 
### 100 permuted datasets, overall and stratified by cluster restricting by sample and gene
### comparing preterm birth cases to term birth controls
obs_pi_perms_p <- list()
obs_pi_perms_p_1 <- list()
obs_pi_perms_p_2 <- list()
obs_pi_perms_p_3 <- list()
obs_tax_perms_p <- list()
obs_tax_perms_p_1 <- list()
obs_tax_perms_p_2 <- list()
obs_tax_perms_p_3 <- list()
pi_obs_perms_p <- list()
pi_obs_perms_p_1 <- list()
pi_obs_perms_p_2 <- list()
pi_obs_perms_p_3 <- list()
tax_obs_perms_p <- list()
tax_obs_perms_p_1 <- list()
tax_obs_perms_p_2 <- list()
tax_obs_perms_p_3 <- list()
for(x in 1:100){
  p_obs_pi <- vector()
  p_obs_pi_1 <- vector()
  p_obs_pi_2 <- vector()
  p_obs_pi_3 <- vector()
  p_obs_tax <- vector()
  p_obs_tax_1 <- vector()
  p_obs_tax_2 <- vector()
  p_obs_tax_3 <- vector()
  p_pi_obs <- vector()
  p_pi_obs_1 <- vector()
  p_pi_obs_2 <- vector()
  p_pi_obs_3 <- vector()
  p_tax_obs <- vector()
  p_tax_obs_1 <- vector()
  p_tax_obs_2 <- vector()
  p_tax_obs_3 <- vector()
  for(each in tp){
    p <- wilcox.test(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each], 
                   obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each])-
             mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each]))
    p_obs_pi <- append(p_obs_pi, p_t)
    p <- wilcox.test(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 1, colnames(obs_tax_perms[[x]]) == each], 
                   obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 0, colnames(obs_tax_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 1, colnames(obs_tax_perms[[x]]) == each])-
             mean(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 0, colnames(obs_tax_perms[[x]]) == each]))
    p_obs_tax <- append(p_obs_tax, p_t)
    p <- wilcox.test(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 1, colnames(pi_obs_perms[[x]]) == each], 
                   pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 0, colnames(pi_obs_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 1, colnames(pi_obs_perms[[x]]) == each])-
             mean(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 0, colnames(pi_obs_perms[[x]]) == each]))
    p_pi_obs <- append(p_pi_obs, p_t)
    p <- wilcox.test(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each], 
                   obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each])-
             mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each]))
    p_tax_obs <- append(p_tax_obs, p_t)
  }
  for(each in tp[tp %!in% obs_pi_1_miss]){
    p <- wilcox.test(obs_pi_perms_1[[x]][obs_pi_perms_1[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_1[[x]]) == each], 
                   obs_pi_perms_1[[x]][obs_pi_perms_1[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_1[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_pi_perms_1[[x]][obs_pi_perms_1[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_1[[x]]) == each])-
             mean(obs_pi_perms_1[[x]][obs_pi_perms_1[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_1[[x]]) == each]))
    p_obs_pi_1 <- append(p_obs_pi_1, p_t)
    p <- wilcox.test(pi_obs_perms_1[[x]][pi_obs_perms_1[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_1[[x]]) == each], 
                   pi_obs_perms_1[[x]][pi_obs_perms_1[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_1[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(pi_obs_perms_1[[x]][pi_obs_perms_1[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_1[[x]]) == each])-
             mean(pi_obs_perms_1[[x]][pi_obs_perms_1[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_1[[x]]) == each]))
    p_pi_obs_1 <- append(p_pi_obs_1, p_t)
  }
  for(each in tp[tp %!in% obs_pi_2_miss]){
    p <- wilcox.test(obs_pi_perms_2[[x]][obs_pi_perms_2[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_2[[x]]) == each], 
                   obs_pi_perms_2[[x]][obs_pi_perms_2[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_2[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_pi_perms_2[[x]][obs_pi_perms_2[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_2[[x]]) == each])-
             mean(obs_pi_perms_2[[x]][obs_pi_perms_2[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_2[[x]]) == each]))
    p_obs_pi_2 <- append(p_obs_pi_2, p_t)
    p <- wilcox.test(pi_obs_perms_2[[x]][pi_obs_perms_2[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_2[[x]]) == each], 
                   pi_obs_perms_2[[x]][pi_obs_perms_2[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_2[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(pi_obs_perms_2[[x]][pi_obs_perms_2[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_2[[x]]) == each])-
             mean(pi_obs_perms_2[[x]][pi_obs_perms_2[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_2[[x]]) == each]))
    p_pi_obs_2 <- append(p_pi_obs_2, p_t)
  }
  for(each in tp[tp %!in% obs_pi_3_miss]){
    p <- wilcox.test(obs_pi_perms_3[[x]][obs_pi_perms_3[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_3[[x]]) == each], 
                   obs_pi_perms_3[[x]][obs_pi_perms_3[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_3[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_pi_perms_3[[x]][obs_pi_perms_3[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_3[[x]]) == each])-
             mean(obs_pi_perms_3[[x]][obs_pi_perms_3[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_3[[x]]) == each]))
    p_obs_pi_3 <- append(p_obs_pi_3, p_t)
    p <- wilcox.test(pi_obs_perms_3[[x]][pi_obs_perms_3[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_3[[x]]) == each], 
                   pi_obs_perms_3[[x]][pi_obs_perms_3[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_3[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(pi_obs_perms_3[[x]][pi_obs_perms_3[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_3[[x]]) == each])-
             mean(pi_obs_perms_3[[x]][pi_obs_perms_3[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_3[[x]]) == each]))
    p_pi_obs_3 <- append(p_pi_obs_3, p_t)
  }
  for(each in tp[tp %!in% obs_tax_1_miss]){
    p <- wilcox.test(obs_tax_perms_1[[x]][obs_tax_perms_1[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_1[[x]]) == each], 
                   obs_tax_perms_1[[x]][obs_tax_perms_1[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_1[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_tax_perms_1[[x]][obs_tax_perms_1[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_1[[x]]) == each])-
             mean(obs_tax_perms_1[[x]][obs_tax_perms_1[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_1[[x]]) == each]))
    p_obs_tax_1 <- append(p_obs_tax_1, p_t)
    p <- wilcox.test(tax_obs_perms_1[[x]][tax_obs_perms_1[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_1[[x]]) == each], 
                   tax_obs_perms_1[[x]][tax_obs_perms_1[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_1[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(tax_obs_perms_1[[x]][tax_obs_perms_1[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_1[[x]]) == each])-
             mean(tax_obs_perms_1[[x]][tax_obs_perms_1[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_1[[x]]) == each]))
    p_tax_obs_1 <- append(p_tax_obs_1, p_t)
  }
  for(each in tp[tp %!in% obs_tax_2_miss]){
    p <- wilcox.test(obs_tax_perms_2[[x]][obs_tax_perms_2[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_2[[x]]) == each], 
                   obs_tax_perms_2[[x]][obs_tax_perms_2[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_2[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_tax_perms_2[[x]][obs_tax_perms_2[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_2[[x]]) == each])-
             mean(obs_tax_perms_2[[x]][obs_tax_perms_2[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_2[[x]]) == each]))
    p_obs_tax_2 <- append(p_obs_tax_2, p_t)
    p <- wilcox.test(tax_obs_perms_2[[x]][tax_obs_perms_2[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_2[[x]]) == each], 
                   tax_obs_perms_2[[x]][tax_obs_perms_2[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_2[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(tax_obs_perms_2[[x]][tax_obs_perms_2[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_2[[x]]) == each])-
             mean(tax_obs_perms_2[[x]][tax_obs_perms_2[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_2[[x]]) == each]))
    p_tax_obs_2 <- append(p_tax_obs_2, p_t)
  }
  for(each in tp[tp %!in% obs_tax_3_miss]){
    p <- wilcox.test(obs_tax_perms_3[[x]][obs_tax_perms_3[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_3[[x]]) == each], 
                   obs_tax_perms_3[[x]][obs_tax_perms_3[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_3[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(obs_tax_perms_3[[x]][obs_tax_perms_3[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_3[[x]]) == each])-
             mean(obs_tax_perms_3[[x]][obs_tax_perms_3[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_3[[x]]) == each]))
    p_obs_tax_3 <- append(p_obs_tax_3, p_t)
    p <- wilcox.test(tax_obs_perms_3[[x]][tax_obs_perms_3[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_3[[x]]) == each], 
                   tax_obs_perms_3[[x]][tax_obs_perms_3[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_3[[x]]) == each])[[3]][1]
    p_t <- log10(p)*
      sign(mean(tax_obs_perms_3[[x]][tax_obs_perms_3[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_3[[x]]) == each])-
             mean(tax_obs_perms_3[[x]][tax_obs_perms_3[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_3[[x]]) == each]))
    p_tax_obs_3 <- append(p_tax_obs_3, p_t)
  }
  obs_pi_perms_p <- c(obs_pi_perms_p, list(p_obs_pi))
  obs_pi_perms_p_1 <- c(obs_pi_perms_p_1, list(p_obs_pi_1))
  obs_pi_perms_p_2 <- c(obs_pi_perms_p_2, list(p_obs_pi_2))
  obs_pi_perms_p_3 <- c(obs_pi_perms_p_3, list(p_obs_pi_3))
  obs_tax_perms_p <- c(obs_tax_perms_p, list(p_obs_tax))
  obs_tax_perms_p_1 <- c(obs_tax_perms_p_1, list(p_obs_tax_1))
  obs_tax_perms_p_2 <- c(obs_tax_perms_p_2, list(p_obs_tax_2))
  obs_tax_perms_p_3 <- c(obs_tax_perms_p_3, list(p_obs_tax_3))
  pi_obs_perms_p <- c(pi_obs_perms_p, list(p_pi_obs))
  pi_obs_perms_p_1 <- c(pi_obs_perms_p_1, list(p_pi_obs_1))
  pi_obs_perms_p_2 <- c(pi_obs_perms_p_2, list(p_pi_obs_2))
  pi_obs_perms_p_3 <- c(pi_obs_perms_p_3, list(p_pi_obs_3))
  tax_obs_perms_p <- c(tax_obs_perms_p, list(p_tax_obs))
  tax_obs_perms_p_1 <- c(tax_obs_perms_p_1, list(p_tax_obs_1))
  tax_obs_perms_p_2 <- c(tax_obs_perms_p_2, list(p_tax_obs_2))
  tax_obs_perms_p_3 <- c(tax_obs_perms_p_3, list(p_tax_obs_3))
}



### generate vectors of spearman correlation coefficients between transformed p values from 100 
### permuted datasets overall and stratified by cluster 
### comparing preterm birth cases to term birth controls
obs_pi_perms_p_corr <- vector()
obs_pi_perms_p_corr_1 <- vector()
obs_pi_perms_p_corr_2 <- vector()
obs_pi_perms_p_corr_3 <- vector()
obs_tax_perms_p_corr <- vector()
obs_tax_perms_p_corr_1 <- vector()
obs_tax_perms_p_corr_2 <- vector()
obs_tax_perms_p_corr_3 <- vector()
for(x in 1:100){
  corr_pi <- cor(obs_pi_perms_p[[x]], pi_obs_perms_p[[x]], method = "spearman")
  corr_pi_1 <- cor(obs_pi_perms_p_1[[x]], pi_obs_perms_p_1[[x]], method = "spearman")
  corr_pi_2 <- cor(obs_pi_perms_p_2[[x]], pi_obs_perms_p_2[[x]], method = "spearman")
  corr_pi_3 <- cor(obs_pi_perms_p_3[[x]], pi_obs_perms_p_3[[x]], method = "spearman")
  corr_tax <- cor(obs_tax_perms_p[[x]], tax_obs_perms_p[[x]], method = "spearman")
  corr_tax_1 <- cor(obs_tax_perms_p_1[[x]], tax_obs_perms_p_1[[x]], method = "spearman")
  corr_tax_2 <- cor(obs_tax_perms_p_2[[x]], tax_obs_perms_p_2[[x]], method = "spearman")
  corr_tax_3 <- cor(obs_tax_perms_p_3[[x]], tax_obs_perms_p_3[[x]], method = "spearman")
  obs_pi_perms_p_corr <- append(obs_pi_perms_p_corr, corr_pi)
  obs_pi_perms_p_corr_1 <- append(obs_pi_perms_p_corr_1, corr_pi_1)
  obs_pi_perms_p_corr_2 <- append(obs_pi_perms_p_corr_2, corr_pi_2)
  obs_pi_perms_p_corr_3 <- append(obs_pi_perms_p_corr_3, corr_pi_3)
  obs_tax_perms_p_corr <- append(obs_tax_perms_p_corr, corr_tax)
  obs_tax_perms_p_corr_1 <- append(obs_tax_perms_p_corr_1, corr_tax_1)
  obs_tax_perms_p_corr_2 <- append(obs_tax_perms_p_corr_2, corr_tax_2)
  obs_tax_perms_p_corr_3 <- append(obs_tax_perms_p_corr_3, corr_tax_3)
}
summary(obs_pi_perms_p_corr)
summary(obs_pi_perms_p_corr_1)
summary(obs_pi_perms_p_corr_2)
summary(obs_pi_perms_p_corr_3)
summary(obs_tax_perms_p_corr)
summary(obs_tax_perms_p_corr_1)
summary(obs_tax_perms_p_corr_2)
summary(obs_tax_perms_p_corr_3)


### generate dataframe for plotting relative abundance correlations from cluster stratified analysis

plot_pi_un <- data.frame(matrix(ncol = 0, nrow = length(obs_pi_corr) + 
                                length(obs_pi_corr_1) + 
                                length(obs_pi_corr_2) + 
                                length(obs_pi_corr_3)))
plot_pi_un$corr <- c(obs_pi_corr, obs_pi_corr_1, obs_pi_corr_2, obs_pi_corr_3)
plot_pi_un$method <- "PICRUSt2"
plot_pi_un$perm <- "Unpermuted"
plot_pi_un$method_perm <- "PICRUSt2,  Unpermuted"
plot_pi_un$cluster <- c(rep("All", length(obs_pi_corr)), 
                  rep("Mixed", length(obs_pi_corr_1)), 
                  rep("L. iners dominated", length(obs_pi_corr_2)), 
                  rep("L. crispatus dominated", length(obs_pi_corr_3)))

plot_tax_un <- data.frame(matrix(ncol = 0, nrow = length(obs_tax_corr) + 
                                    length(obs_tax_corr_1) + 
                                    length(obs_tax_corr_2) + 
                                    length(obs_tax_corr_3)))
plot_tax_un$corr <- c(obs_tax_corr, obs_tax_corr_1, obs_tax_corr_2, obs_tax_corr_3)
plot_tax_un$method <- "Tax4Fun2"
plot_tax_un$perm <- "Unpermuted"
plot_tax_un$method_perm <- "Tax4Fun2,  Unpermuted"
plot_tax_un$cluster <- c(rep("All", length(obs_tax_corr)), 
                      rep("Mixed", length(obs_tax_corr_1)), 
                      rep("L. iners dominated", length(obs_tax_corr_2)), 
                      rep("L. crispatus dominated", length(obs_tax_corr_3)))

plot_pi_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_pi_perms_corr) + 
                                  length(obs_pi_perms_corr_1) + 
                                  length(obs_pi_perms_corr_2) + 
                                  length(obs_pi_perms_corr_3)))
plot_pi_perm$corr <- c(obs_pi_perms_corr, obs_pi_perms_corr_1, obs_pi_perms_corr_2, obs_pi_perms_corr_3)
plot_pi_perm$method <- "PICRUSt2"
plot_pi_perm$perm <- "Permuted"
plot_pi_perm$method_perm <- "PICRUSt2,  Permuted"
plot_pi_perm$cluster <- c(rep("All", length(obs_pi_perms_corr)), 
                    rep("Mixed", length(obs_pi_perms_corr_1)), 
                    rep("L. iners dominated", length(obs_pi_perms_corr_2)), 
                    rep("L. crispatus dominated", length(obs_pi_perms_corr_3)))

plot_tax_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_tax_perms_corr) + 
                                      length(obs_tax_perms_corr_1) + 
                                      length(obs_tax_perms_corr_2) + 
                                      length(obs_tax_perms_corr_3)))
plot_tax_perm$corr <- c(obs_tax_perms_corr, obs_tax_perms_corr_1, obs_tax_perms_corr_2, obs_tax_perms_corr_3)
plot_tax_perm$method <- "Tax4Fun2"
plot_tax_perm$perm <- "Permuted"
plot_tax_perm$method_perm <- "Tax4Fun2,  Permuted"
plot_tax_perm$cluster <- c(rep("All", length(obs_tax_perms_corr)), 
                        rep("Mixed", length(obs_tax_perms_corr_1)), 
                        rep("L. iners dominated", length(obs_tax_perms_corr_2)), 
                        rep("L. crispatus dominated", length(obs_tax_perms_corr_3)))

plot_ra_cluster <- as.data.frame(rbind(plot_pi_un, plot_tax_un, plot_pi_perm, plot_tax_perm))
plot_ra_cluster$cluster <- factor(plot_ra_cluster$cluster,
                                  levels = c("All",
                                             "L. crispatus dominated",
                                             "L. iners dominated",
                                             "Mixed"))
# plot_ra_cluster should have 1016464 obs of 5 vars



### generate dataframe for plotting transformed p value (comparing preterm birth cases to term 
### birth controls) correlations from cluster stratified analysis

plot_pi_un <- data.frame(matrix(ncol = 0, nrow = 4))
plot_pi_un$corr <- c(obs_pi_p_corr, 
                   obs_pi_p_corr_1, 
                   obs_pi_p_corr_2, 
                   obs_pi_p_corr_3)
plot_pi_un$method <- "PICRUSt2"
plot_pi_un$perm <- "Unpermuted"
plot_pi_un$method_perm <- "PICRUSt2,  Unpermuted"
plot_pi_un$cluster <- c("All", "Mixed", "L. iners dominated", "L. crispatus dominated")

plot_tax_un <- data.frame(matrix(ncol = 0, nrow = 4))
plot_tax_un$corr <- c(obs_tax_p_corr, 
                       obs_tax_p_corr_1, 
                       obs_tax_p_corr_2, 
                       obs_tax_p_corr_3)
plot_tax_un$method <- "Tax4Fun2"
plot_tax_un$perm <- "Unpermuted"
plot_tax_un$method_perm <- "Tax4Fun2,  Unpermuted"
plot_tax_un$cluster <- c("All", "Mixed", "L. iners dominated", "L. crispatus dominated")

plot_pi_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_pi_perms_p_corr)+
                                  length(obs_pi_perms_p_corr_1)+
                                  length(obs_pi_perms_p_corr_2)+
                                  length(obs_pi_perms_p_corr_3)))
plot_pi_perm$corr <- c(obs_pi_perms_p_corr, 
                     obs_pi_perms_p_corr_1, 
                     obs_pi_perms_p_corr_2, 
                     obs_pi_perms_p_corr_3)
plot_pi_perm$method <- "PICRUSt2"
plot_pi_perm$perm <- "Permuted"
plot_pi_perm$method_perm <- "PICRUSt2,  Permuted"
plot_pi_perm$cluster <- c(rep("All", length(obs_pi_perms_p_corr)), 
                    rep("Mixed", length(obs_pi_perms_p_corr_1)), 
                    rep("L. iners dominated", length(obs_pi_perms_p_corr_2)), 
                    rep("L. crispatus dominated", length(obs_pi_perms_p_corr_3)))

plot_tax_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_tax_perms_p_corr)+
                                      length(obs_tax_perms_p_corr_1)+
                                      length(obs_tax_perms_p_corr_2)+
                                      length(obs_tax_perms_p_corr_3)))
plot_tax_perm$corr <- c(obs_tax_perms_p_corr, 
                         obs_tax_perms_p_corr_1, 
                         obs_tax_perms_p_corr_2, 
                         obs_tax_perms_p_corr_3)
plot_tax_perm$method <- "Tax4Fun2"
plot_tax_perm$perm <- "Permuted"
plot_tax_perm$method_perm <- "Tax4Fun2,  Permuted"
plot_tax_perm$cluster <- c(rep("All", length(obs_tax_perms_p_corr)), 
                        rep("Mixed", length(obs_tax_perms_p_corr_1)), 
                        rep("L. iners dominated", length(obs_tax_perms_p_corr_2)), 
                        rep("L. crispatus dominated", length(obs_tax_perms_p_corr_3)))

plot_pcorr_cluster <- as.data.frame(rbind(plot_pi_un, plot_tax_un, plot_pi_perm, plot_tax_perm))
plot_pcorr_cluster$cluster <- factor(plot_pcorr_cluster$cluster, 
                                     levels = c("All",
                                                "L. crispatus dominated",
                                                "L. iners dominated",
                                                "Mixed"))
# plot_pcorr_cluster should have 808 obs of 5 vars






######################################################################################
### METAGENOME INFERENCE COMPARISON, UNPERMUTED, KO FUNCTIONAL CATEGORY STRATIFIED ###
######################################################################################



### generate dataframe with "true positive" genes and their level 1 KO functional categories

ko_cat <- as.data.frame(substr(tp, 1, 6))
colnames(ko_cat) <- "ko"

ko_trim <- ko[, colnames(ko) %in% c("level1", "level4")]
colnames(ko_trim) <- c("level1", "ko")

ko_cat <- left_join(ko_cat, ko_trim)

# make level1 factor variable
ko_cat <- ko_cat %>% mutate(cat = case_when(level1 == "09100 Metabolism" ~ "Metabolism", 
                                            level1 == "09130 Environmental Information Processing" ~ "Environmental information processing", 
                                            level1 == "09150 Organismal Systems" ~ "Organismal systems", 
                                            level1 == "09160 Human Diseases" ~ "Human diseases", 
                                            level1 == "09180 Brite Hierarchies" ~ "BRITE hierarchies", 
                                            level1 == "09190 Not Included in Pathway or Brite" ~ "Not included in pathway or BRITE", 
                                            level1 == "09120 Genetic Information Processing" ~ "Genetic information processing", 
                                            level1 == "09140 Cellular Processes" ~ "Cellular processes"))
ko_cat$cat <- factor(ko_cat$cat, levels = c("Metabolism", 
                                            "Genetic information processing", 
                                            "Environmental information processing", 
                                            "Cellular processes", 
                                            "Organismal systems", 
                                            "Human diseases", 
                                            "BRITE hierarchies", 
                                            "Not included in pathway or BRITE"))

# add "_ra" to end of ko names to match other dataframes
ko_cat$ko <- paste0(ko_cat$ko, "_ra")

# subset ko_cat to unique observations, remove observations with cat == NA
ko_cat <- unique(ko_cat)
ko_cat <- ko_cat[is.na(ko_cat$cat)!= T, ]
# ko_cat should have 2168 obs of 3 vars



### restrict obs, pi, tax to KO realtive abundances for KOs with level-1 KO functional category,
### sample variable, and case variable

obs <- obs[, colnames(obs) %in% c(ko_cat$ko, "sample", "C_CASE2")]
obs <- obs[, order(colnames(obs))]
obs <- obs[order(rownames(obs)), ]
# obs should have 72 obs of 1489 vars

pi <- pi[, colnames(pi) %in% c(ko_cat$ko, "sample", "C_CASE2")]
pi <- pi[, order(colnames(pi))]
pi <- pi[order(rownames(pi)), ]
# pi should have 72 obs of 1489 vars

tax <- tax[, colnames(tax) %in% c(ko_cat$ko, "sample", "C_CASE2")]
tax <- tax[, order(colnames(tax))]
tax <- tax[order(rownames(tax)), ]
# tax should have 72 obs of 1489 vars



### estimate spearman correlation coefficients between observed and shared KO relative abundances 
### overall and restricted by level-1 KO functional category

obs_pi_corr <- vector()
obs_pi_corr_cat <- vector()
for(each in tp){
  cat <- ko_cat$cat[ko_cat$ko == each]
  obs_pi_corr_cat <- append(obs_pi_corr_cat, cat)
  corr <- cor(obs[, colnames(obs) == each], pi[, colnames(pi) == each], method = "spearman")
  obs_pi_corr <- append(obs_pi_corr, rep(corr, length(cat)))
}
names(obs_pi_corr) <- obs_pi_corr_cat
summary(obs_pi_corr)
# obs_pi_corr should have 2168 entries
summary(obs_pi_corr[names(obs_pi_corr) == 1]) # Metabolism
summary(obs_pi_corr[names(obs_pi_corr) == 2]) # Genetic information processing
summary(obs_pi_corr[names(obs_pi_corr) == 3]) # Environmental information processing
summary(obs_pi_corr[names(obs_pi_corr) == 4]) # Cellular processes
summary(obs_pi_corr[names(obs_pi_corr) == 5]) # Organismal systems
summary(obs_pi_corr[names(obs_pi_corr) == 6]) # Human diseases
summary(obs_pi_corr[names(obs_pi_corr) == 7]) # BRITE hierarchies
summary(obs_pi_corr[names(obs_pi_corr) == 8]) # Not included in pathway or BRITE

obs_tax_corr <- vector()
obs_tax_corr_cat <- vector()
for(each in tp){
  cat <- ko_cat$cat[ko_cat$ko == each]
  obs_tax_corr_cat <- append(obs_tax_corr_cat, cat)
  corr <- cor(obs[, colnames(obs) == each], tax[, colnames(tax) == each], method = "spearman")
  obs_tax_corr <- append(obs_tax_corr, rep(corr, length(cat)))
}
names(obs_tax_corr) <- obs_tax_corr_cat
summary(obs_tax_corr)
# obs_tax_corr should have 2168 elements
summary(obs_tax_corr[names(obs_tax_corr) == 1]) # Metabolism
summary(obs_tax_corr[names(obs_tax_corr) == 2]) # Genetic information processing
summary(obs_tax_corr[names(obs_tax_corr) == 3]) # Environmental information processing
summary(obs_tax_corr[names(obs_tax_corr) == 4]) # Cellular processes
summary(obs_tax_corr[names(obs_tax_corr) == 5]) # Organismal systems
summary(obs_tax_corr[names(obs_tax_corr) == 6]) # Human diseases
summary(obs_tax_corr[names(obs_tax_corr) == 7]) # BRITE hierarchies
summary(obs_tax_corr[names(obs_tax_corr) == 8]) # Not included in pathway or BRITE



### inference using wilcoxon tests comparing KO relative abundances between preterm birth cases 
### and term birth controls

obs_p <- vector()
obs_p_cat <- vector()
gene <- vector()
for(each in tp){
  for(level1 in ko_cat$cat[ko_cat$ko == each]){
    p <- wilcox.test(obs[obs$C_CASE2 == 1, colnames(obs) == each],
                   obs[obs$C_CASE2 == 0, colnames(obs) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(obs[obs$C_CASE2 == 1, colnames(obs) == each])-
                         mean(obs[obs$C_CASE2 == 0, colnames(obs) == each]))
    obs_p <- append(obs_p, p_t)
    obs_p_cat <- append(obs_p_cat, level1)
    gene <- append(gene, each)
  }
}
obs_inf <- data.frame(matrix(ncol = 0, nrow = length(obs_p)))
obs_inf$p <- obs_p
obs_inf$cat <- obs_p_cat
obs_inf$gene <- gene
# obs_inf should have 2168 obs of 3 vars
summary(obs_inf$p)

pi_p <- vector()
pi_p_cat <- vector()
gene <- vector()
for(each in tp){
  for(level1 in ko_cat$cat[ko_cat$ko == each]){
    p <- wilcox.test(pi[pi$C_CASE2 == 1, colnames(pi) == each],
                   pi[pi$C_CASE2 == 0, colnames(pi) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(pi[pi$C_CASE2 == 1, colnames(pi) == each])-
                         mean(pi[pi$C_CASE2 == 0, colnames(pi) == each]))
    pi_p <- append(pi_p, p_t)
    pi_p_cat <- append(pi_p_cat, level1)
    gene <- append(gene, each)
  }
}
pi_inf <- data.frame(matrix(ncol = 0, nrow = length(pi_p)))
pi_inf$p <- pi_p
pi_inf$cat <- pi_p_cat
pi_inf$gene <- gene
# pi_inf should have 2168 obs of 3 vars
summary(pi_inf$p)

tax_p <- vector()
tax_p_cat <- vector()
gene <- vector()
for(each in tp){
  for(level1 in ko_cat$cat[ko_cat$ko == each]){
    p <- wilcox.test(tax[tax$C_CASE2 == 1, colnames(tax) == each],
                   tax[tax$C_CASE2 == 0, colnames(tax) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(tax[tax$C_CASE2 == 1, colnames(tax) == each])-
                         mean(tax[tax$C_CASE2 == 0, colnames(tax) == each]))
    tax_p <- append(tax_p, p_t)
    tax_p_cat <- append(tax_p_cat, level1)
    gene <- append(gene, each)
  }
}
tax_inf <- data.frame(matrix(ncol = 0, nrow = length(tax_p)))
tax_inf$p <- tax_p
tax_inf$cat <- tax_p_cat
tax_inf$gene <- gene
# tax_inf should have 2168 obs of 3 vars
summary(tax_inf$p)



# calculate spearman correlation coefficients between transformed p values
obs_pi_p_corr <- cor(obs_inf$p, pi_inf$p, method = "spearman")
obs_pi_p_corr_met <- cor(obs_inf$p[obs_inf$cat == "Metabolism"],
                                pi_inf$p[pi_inf$cat == "Metabolism"],method = "spearman")
obs_pi_p_corr_gen <- cor(obs_inf$p[obs_inf$cat == "Genetic information processing"],
                                pi_inf$p[pi_inf$cat == "Genetic information processing"],method = "spearman")
obs_pi_p_corr_env <- cor(obs_inf$p[obs_inf$cat == "Environmental information processing"],
                                pi_inf$p[pi_inf$cat == "Environmental information processing"],method = "spearman")
obs_pi_p_corr_cel <- cor(obs_inf$p[obs_inf$cat == "Cellular processes"],
                                pi_inf$p[pi_inf$cat == "Cellular processes"],method = "spearman")
obs_pi_p_corr_org <- cor(obs_inf$p[obs_inf$cat == "Organismal systems"],
                                pi_inf$p[pi_inf$cat == "Organismal systems"],method = "spearman")
obs_pi_p_corr_hum <- cor(obs_inf$p[obs_inf$cat == "Human diseases"],
                                pi_inf$p[pi_inf$cat == "Human diseases"],method = "spearman")
obs_pi_p_corr_bri <- cor(obs_inf$p[obs_inf$cat == "BRITE hierarchies"],
                                pi_inf$p[pi_inf$cat == "BRITE hierarchies"],method = "spearman")
obs_pi_p_corr_not <- cor(obs_inf$p[obs_inf$cat == "Not included in pathway or BRITE"],
                                pi_inf$p[pi_inf$cat == "Not included in pathway or BRITE"],method = "spearman")

obs_tax_p_corr <- cor(obs_inf$p, tax_inf$p, method = "spearman")
obs_tax_p_corr_met <- cor(obs_inf$p[obs_inf$cat == "Metabolism"],
                                    tax_inf$p[tax_inf$cat == "Metabolism"],method = "spearman")
obs_tax_p_corr_gen <- cor(obs_inf$p[obs_inf$cat == "Genetic information processing"],
                                    tax_inf$p[tax_inf$cat == "Genetic information processing"],method = "spearman")
obs_tax_p_corr_env <- cor(obs_inf$p[obs_inf$cat == "Environmental information processing"],
                                    tax_inf$p[tax_inf$cat == "Environmental information processing"],method = "spearman")
obs_tax_p_corr_cel <- cor(obs_inf$p[obs_inf$cat == "Cellular processes"],
                                    tax_inf$p[tax_inf$cat == "Cellular processes"],method = "spearman")
obs_tax_p_corr_org <- cor(obs_inf$p[obs_inf$cat == "Organismal systems"],
                                    tax_inf$p[tax_inf$cat == "Organismal systems"],method = "spearman")
obs_tax_p_corr_hum <- cor(obs_inf$p[obs_inf$cat == "Human diseases"],
                                    tax_inf$p[tax_inf$cat == "Human diseases"],method = "spearman")
obs_tax_p_corr_bri <- cor(obs_inf$p[obs_inf$cat == "BRITE hierarchies"],
                                    tax_inf$p[tax_inf$cat == "BRITE hierarchies"],method = "spearman")
obs_tax_p_corr_not <- cor(obs_inf$p[obs_inf$cat == "Not included in pathway or BRITE"],
                                    tax_inf$p[tax_inf$cat == "Not included in pathway or BRITE"],method = "spearman")



### generate dataframe of transformed p values for plotting

colnames(obs_inf) <- c("p_obs", "cat", "gene")
colnames(pi_inf) <- c("p_pi", "cat", "gene")
colnames(tax_inf) <- c("p_tax", "cat", "gene")
plot_pvp_ko <- left_join(obs_inf, pi_inf)
plot_pvp_ko <- left_join(plot_pvp_ko, tax_inf)
# plot_pvp_ko should have 2168 obs of 5 vars






####################################################################################
### METAGENOME INFERENCE COMPARISON, PERMUTED, KO FUNCTIONAL CATEGORY STRATIFIED ###
####################################################################################



### generate permuted datasets overall

set.seed(1)

# generate lists with 100 permuted versions of obs, pi
obs_pi_perms <- list()
pi_obs_perms <- list()
for(x in 1:100){
  obs_pi_perm <- obs
  pi_obs_perm <- pi
  for(each in ko_cat$ko){
    obs_pi_perm[, colnames(obs_pi_perm) == each] <- sample(obs_pi_perm[, colnames(obs_pi_perm) == each])
    pi_obs_perm[, colnames(pi_obs_perm) == each] <- sample(pi_obs_perm[, colnames(pi_obs_perm) == each])
  }
  obs_pi_perms <- c(obs_pi_perms, list(obs_pi_perm))
  pi_obs_perms <- c(pi_obs_perms, list(pi_obs_perm))
}

# generate lists with 100 permuted versions of obs, tax
obs_tax_perms <- list()
tax_obs_perms <- list()
for(x in 1:100){
  obs_tax_perm <- obs
  tax_obs_perm <- tax
  for(each in ko_cat$ko){
    obs_tax_perm[, colnames(obs_tax_perm) == each] <- sample(obs_tax_perm[, colnames(obs_tax_perm) == each])
    tax_obs_perm[, colnames(tax_obs_perm) == each] <- sample(tax_obs_perm[, colnames(tax_obs_perm) == each])
  }
  obs_tax_perms <- c(obs_tax_perms, list(obs_tax_perm))
  tax_obs_perms <- c(tax_obs_perms, list(tax_obs_perm))
}



### generate vector of correlation coefficients between permuted gene relative abundance across  
### 100 permuted datasets, overall

# for obs and pi
obs_pi_perms_corr <- vector()
for(x in 1:100){
  for(each in ko_cat$ko){
    corr <- cor(obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) == each], 
                pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) == each], 
                method = "spearman")
    obs_pi_perms_corr <- append(obs_pi_perms_corr, corr)
  }
}
summary(obs_pi_perms_corr)



# for obs and tax
obs_tax_perms_corr <- vector()
for(x in 1:100){
  for(each in ko_cat$ko){
    corr <- cor(obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) == each], 
                tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) == each], 
                method = "spearman")
    obs_tax_perms_corr <- append(obs_tax_perms_corr, corr)
  }
}
summary(obs_tax_perms_corr)



### generate lists of p value for wilcoxon test for each permuted gene relative abundance across 
### 100 permuted datasets, comparing preterm birth cases to term birth controls, overall

# for obs
obs_pi_perms_p <- list()
for(x in 1:100){
  p_vec <- vector()
  for(each in ko_cat$ko){
    p <- wilcox.test(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each], 
                   obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 1, colnames(obs_pi_perms[[x]]) == each])-
                         mean(obs_pi_perms[[x]][obs_pi_perms[[x]]$C_CASE2 == 0, colnames(obs_pi_perms[[x]]) == each]))
    p_vec <- append(p_vec, p_t)
  }
  obs_pi_perms_p <- c(obs_pi_perms_p, list(p_vec))
}

obs_tax_perms_p <- list()
for(x in 1:100){
  p_vec <- vector()
  for(each in ko_cat$ko){
    p <- wilcox.test(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 1, colnames(obs_tax_perms[[x]]) == each], 
                     obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 0, colnames(obs_tax_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 1, colnames(obs_tax_perms[[x]]) == each])-
                           mean(obs_tax_perms[[x]][obs_tax_perms[[x]]$C_CASE2 == 0, colnames(obs_tax_perms[[x]]) == each]))
    p_vec <- append(p_vec, p_t)
  }
  obs_tax_perms_p <- c(obs_tax_perms_p, list(p_vec))
}

# for pi
pi_obs_perms_p <- list()
for(x in 1:100){
  p_vec <- vector()
  for(each in ko_cat$ko){
    p <- wilcox.test(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 1, colnames(pi_obs_perms[[x]]) == each], 
                     pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 0, colnames(pi_obs_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 1, colnames(pi_obs_perms[[x]]) == each])-
                           mean(pi_obs_perms[[x]][pi_obs_perms[[x]]$C_CASE2 == 0, colnames(pi_obs_perms[[x]]) == each]))
    p_vec <- append(p_vec, p_t)
  }
  pi_obs_perms_p <- c(pi_obs_perms_p, list(p_vec))
}

# for tax
tax_obs_perms_p <- list()
for(x in 1:100){
  p_vec <- vector()
  for(each in ko_cat$ko){
    p <- wilcox.test(tax_obs_perms[[x]][tax_obs_perms[[x]]$C_CASE2 == 1, colnames(tax_obs_perms[[x]]) == each], 
                     tax_obs_perms[[x]][tax_obs_perms[[x]]$C_CASE2 == 0, colnames(tax_obs_perms[[x]]) == each])[[3]][1]
    p_t <- log10(p)*sign(mean(tax_obs_perms[[x]][tax_obs_perms[[x]]$C_CASE2 == 1, colnames(tax_obs_perms[[x]]) == each])-
                           mean(tax_obs_perms[[x]][tax_obs_perms[[x]]$C_CASE2 == 0, colnames(tax_obs_perms[[x]]) == each]))
    p_vec <- append(p_vec, p_t)
  }
  tax_obs_perms_p <- c(tax_obs_perms_p, list(p_vec))
}



### generate vector of correlation coefficients between transformed p values (comparing preterm
### birth cases to term birth controls) across 100 permuted datasets, overall

# for obs and pi
obs_pi_perms_p_corr <- vector()
for(x in 1:100){
  corr <- cor(obs_pi_perms_p[[x]], pi_obs_perms_p[[x]], method = "spearman")
  obs_pi_perms_p_corr <- append(obs_pi_perms_p_corr, corr)
}
summary(obs_pi_perms_p_corr)

# for obs and tax
obs_tax_perms_p_corr <- vector()
for(x in 1:100){
  corr <- cor(obs_tax_perms_p[[x]], tax_obs_perms_p[[x]], method = "spearman")
  obs_tax_perms_p_corr <- append(obs_tax_perms_p_corr, corr)
}
summary(obs_tax_perms_p_corr)



### generate lists with 100 permuted versions of obs, pi, tax restricted by level-1
### KO functional category
obs_pi_perms_met <- list()
obs_pi_perms_gen <- list()
obs_pi_perms_env <- list()
obs_pi_perms_cel <- list()
obs_pi_perms_org <- list()
obs_pi_perms_hum <- list()
obs_pi_perms_bri <- list()
obs_pi_perms_not <- list()

obs_tax_perms_met <- list()
obs_tax_perms_gen <- list()
obs_tax_perms_env <- list()
obs_tax_perms_cel <- list()
obs_tax_perms_org <- list()
obs_tax_perms_hum <- list()
obs_tax_perms_bri <- list()
obs_tax_perms_not <- list()

pi_obs_perms_met <- list()
pi_obs_perms_gen <- list()
pi_obs_perms_env <- list()
pi_obs_perms_cel <- list()
pi_obs_perms_org <- list()
pi_obs_perms_hum <- list()
pi_obs_perms_bri <- list()
pi_obs_perms_not <- list()

tax_obs_perms_met <- list()
tax_obs_perms_gen <- list()
tax_obs_perms_env <- list()
tax_obs_perms_cel <- list()
tax_obs_perms_org <- list()
tax_obs_perms_hum <- list()
tax_obs_perms_bri <- list()
tax_obs_perms_not <- list()

for(x in 1:100){
  obs_pi_perms_met[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Metabolism"], "C_CASE2", "sample")]
  obs_pi_perms_gen[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Genetic information processing"], "C_CASE2", "sample")]
  obs_pi_perms_env[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Environmental information processing"], "C_CASE2", "sample")]
  obs_pi_perms_cel[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Cellular processes"], "C_CASE2", "sample")]
  obs_pi_perms_org[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Organismal systems"], "C_CASE2", "sample")]
  obs_pi_perms_hum[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Human diseases"], "C_CASE2", "sample")]
  obs_pi_perms_bri[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "BRITE hierarchies"], "C_CASE2", "sample")]
  obs_pi_perms_not[[x]] <- obs_pi_perms[[x]][, colnames(obs_pi_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Not included in pathway or BRITE"], "C_CASE2", "sample")]
  
  obs_tax_perms_met[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Metabolism"], "C_CASE2", "sample")]
  obs_tax_perms_gen[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Genetic information processing"], "C_CASE2", "sample")]
  obs_tax_perms_env[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Environmental information processing"], "C_CASE2", "sample")]
  obs_tax_perms_cel[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Cellular processes"], "C_CASE2", "sample")]
  obs_tax_perms_org[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Organismal systems"], "C_CASE2", "sample")]
  obs_tax_perms_hum[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Human diseases"], "C_CASE2", "sample")]
  obs_tax_perms_bri[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "BRITE hierarchies"], "C_CASE2", "sample")]
  obs_tax_perms_not[[x]] <- obs_tax_perms[[x]][, colnames(obs_tax_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Not included in pathway or BRITE"], "C_CASE2", "sample")]
  
  pi_obs_perms_met[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Metabolism"], "C_CASE2", "sample")]
  pi_obs_perms_gen[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Genetic information processing"], "C_CASE2", "sample")]
  pi_obs_perms_env[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Environmental information processing"], "C_CASE2", "sample")]
  pi_obs_perms_cel[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Cellular processes"], "C_CASE2", "sample")]
  pi_obs_perms_org[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Organismal systems"], "C_CASE2", "sample")]
  pi_obs_perms_hum[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Human diseases"], "C_CASE2", "sample")]
  pi_obs_perms_bri[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "BRITE hierarchies"], "C_CASE2", "sample")]
  pi_obs_perms_not[[x]] <- pi_obs_perms[[x]][, colnames(pi_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Not included in pathway or BRITE"], "C_CASE2", "sample")]
  
  tax_obs_perms_met[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Metabolism"], "C_CASE2", "sample")]
  tax_obs_perms_gen[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Genetic information processing"], "C_CASE2", "sample")]
  tax_obs_perms_env[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Environmental information processing"], "C_CASE2", "sample")]
  tax_obs_perms_cel[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Cellular processes"], "C_CASE2", "sample")]
  tax_obs_perms_org[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Organismal systems"], "C_CASE2", "sample")]
  tax_obs_perms_hum[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Human diseases"], "C_CASE2", "sample")]
  tax_obs_perms_bri[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "BRITE hierarchies"], "C_CASE2", "sample")]
  tax_obs_perms_not[[x]] <- tax_obs_perms[[x]][, colnames(tax_obs_perms[[x]]) %in% c(ko_cat$ko[ko_cat$cat == "Not included in pathway or BRITE"], "C_CASE2", "sample")]
}



### generate vectors of correlation coefficients between permuted gene relative abundance 
### across to 100 permuted datasets,  restricted by KO category
obs_pi_perms_corr_met <- vector()
obs_pi_perms_corr_gen <- vector()
obs_pi_perms_corr_env <- vector()
obs_pi_perms_corr_cel <- vector()
obs_pi_perms_corr_org <- vector()
obs_pi_perms_corr_hum <- vector()
obs_pi_perms_corr_bri <- vector()
obs_pi_perms_corr_not <- vector()

obs_tax_perms_corr_met <- vector()
obs_tax_perms_corr_gen <- vector()
obs_tax_perms_corr_env <- vector()
obs_tax_perms_corr_cel <- vector()
obs_tax_perms_corr_org <- vector()
obs_tax_perms_corr_hum <- vector()
obs_tax_perms_corr_bri <- vector()
obs_tax_perms_corr_not <- vector()

for(x in 1:100){
  for(each in tp){
    corr_met <- cor(obs_pi_perms_met[[x]][, colnames(obs_pi_perms_met[[x]]) == each], 
                    pi_obs_perms_met[[x]][, colnames(pi_obs_perms_met[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_met <- append(obs_pi_perms_corr_met, corr_met)
    corr_gen <- cor(obs_pi_perms_gen[[x]][, colnames(obs_pi_perms_gen[[x]]) == each], 
                    pi_obs_perms_gen[[x]][, colnames(pi_obs_perms_gen[[x]]) == each],
                    method = "spearman")
    obs_pi_perms_corr_gen <- append(obs_pi_perms_corr_gen, corr_gen)
    corr_env <- cor(obs_pi_perms_env[[x]][, colnames(obs_pi_perms_env[[x]]) == each], 
                    pi_obs_perms_env[[x]][, colnames(pi_obs_perms_env[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_env <- append(obs_pi_perms_corr_env, corr_env)
    corr_cel <- cor(obs_pi_perms_cel[[x]][, colnames(obs_pi_perms_cel[[x]]) == each], 
                    pi_obs_perms_cel[[x]][, colnames(pi_obs_perms_cel[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_cel <- append(obs_pi_perms_corr_cel, corr_cel)
    corr_org <- cor(obs_pi_perms_org[[x]][, colnames(obs_pi_perms_org[[x]]) == each], 
                    pi_obs_perms_org[[x]][, colnames(pi_obs_perms_org[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_org <- append(obs_pi_perms_corr_org, corr_org)
    corr_hum <- cor(obs_pi_perms_hum[[x]][, colnames(obs_pi_perms_hum[[x]]) == each], 
                    pi_obs_perms_hum[[x]][, colnames(pi_obs_perms_hum[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_hum <- append(obs_pi_perms_corr_hum, corr_hum)
    corr_bri <- cor(obs_pi_perms_bri[[x]][, colnames(obs_pi_perms_bri[[x]]) == each], 
                    pi_obs_perms_bri[[x]][, colnames(pi_obs_perms_bri[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_bri <- append(obs_pi_perms_corr_bri, corr_bri)
    corr_not <- cor(obs_pi_perms_not[[x]][, colnames(obs_pi_perms_not[[x]]) == each], 
                    pi_obs_perms_not[[x]][, colnames(pi_obs_perms_not[[x]]) == each], 
                    method = "spearman")
    obs_pi_perms_corr_not <- append(obs_pi_perms_corr_not, corr_not)
    
    corr_met <- cor(obs_tax_perms_met[[x]][, colnames(obs_tax_perms_met[[x]]) == each], 
                    tax_obs_perms_met[[x]][, colnames(tax_obs_perms_met[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_met <- append(obs_tax_perms_corr_met, corr_met)
    corr_gen <- cor(obs_tax_perms_gen[[x]][, colnames(obs_tax_perms_gen[[x]]) == each], 
                    tax_obs_perms_gen[[x]][, colnames(tax_obs_perms_gen[[x]]) == each],
                    method = "spearman")
    obs_tax_perms_corr_gen <- append(obs_tax_perms_corr_gen, corr_gen)
    corr_env <- cor(obs_tax_perms_env[[x]][, colnames(obs_tax_perms_env[[x]]) == each], 
                    tax_obs_perms_env[[x]][, colnames(tax_obs_perms_env[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_env <- append(obs_tax_perms_corr_env, corr_env)
    corr_cel <- cor(obs_tax_perms_cel[[x]][, colnames(obs_tax_perms_cel[[x]]) == each], 
                    tax_obs_perms_cel[[x]][, colnames(tax_obs_perms_cel[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_cel <- append(obs_tax_perms_corr_cel, corr_cel)
    corr_org <- cor(obs_tax_perms_org[[x]][, colnames(obs_tax_perms_org[[x]]) == each], 
                    tax_obs_perms_org[[x]][, colnames(tax_obs_perms_org[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_org <- append(obs_tax_perms_corr_org, corr_org)
    corr_hum <- cor(obs_tax_perms_hum[[x]][, colnames(obs_tax_perms_hum[[x]]) == each], 
                    tax_obs_perms_hum[[x]][, colnames(tax_obs_perms_hum[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_hum <- append(obs_tax_perms_corr_hum, corr_hum)
    corr_bri <- cor(obs_tax_perms_bri[[x]][, colnames(obs_tax_perms_bri[[x]]) == each], 
                    tax_obs_perms_bri[[x]][, colnames(tax_obs_perms_bri[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_bri <- append(obs_tax_perms_corr_bri, corr_bri)
    corr_not <- cor(obs_tax_perms_not[[x]][, colnames(obs_tax_perms_not[[x]]) == each], 
                    tax_obs_perms_not[[x]][, colnames(tax_obs_perms_not[[x]]) == each], 
                    method = "spearman")
    obs_tax_perms_corr_not <- append(obs_tax_perms_corr_not, corr_not)
  }
}
summary(obs_pi_perms_corr_met)
summary(obs_pi_perms_corr_gen)
summary(obs_pi_perms_corr_env)
summary(obs_pi_perms_corr_cel)
summary(obs_pi_perms_corr_org)
summary(obs_pi_perms_corr_hum)
summary(obs_pi_perms_corr_bri)
summary(obs_pi_perms_corr_not)

summary(obs_tax_perms_corr_met)
summary(obs_tax_perms_corr_gen)
summary(obs_tax_perms_corr_env)
summary(obs_tax_perms_corr_cel)
summary(obs_tax_perms_corr_org)
summary(obs_tax_perms_corr_hum)
summary(obs_tax_perms_corr_bri)
summary(obs_tax_perms_corr_not)



### generate list of p values for wilcoxon test comparing preterm birth cases to 
### term birth controls, across 100 permuted datasets and restricted by KO functional category
obs_pi_perms_p_met <- list()
obs_pi_perms_p_gen <- list()
obs_pi_perms_p_env <- list()
obs_pi_perms_p_cel <- list()
obs_pi_perms_p_org <- list()
obs_pi_perms_p_hum <- list()
obs_pi_perms_p_bri <- list()
obs_pi_perms_p_not <- list()

obs_tax_perms_p_met <- list()
obs_tax_perms_p_gen <- list()
obs_tax_perms_p_env <- list()
obs_tax_perms_p_cel <- list()
obs_tax_perms_p_org <- list()
obs_tax_perms_p_hum <- list()
obs_tax_perms_p_bri <- list()
obs_tax_perms_p_not <- list()

pi_obs_perms_p_met <- list()
pi_obs_perms_p_gen <- list()
pi_obs_perms_p_env <- list()
pi_obs_perms_p_cel <- list()
pi_obs_perms_p_org <- list()
pi_obs_perms_p_hum <- list()
pi_obs_perms_p_bri <- list()
pi_obs_perms_p_not <- list()

tax_obs_perms_p_met <- list()
tax_obs_perms_p_gen <- list()
tax_obs_perms_p_env <- list()
tax_obs_perms_p_cel <- list()
tax_obs_perms_p_org <- list()
tax_obs_perms_p_hum <- list()
tax_obs_perms_p_bri <- list()
tax_obs_perms_p_not <- list()

for(x in 1:100){
  p_vec_met <- vector()
  p_vec_gen <- vector()
  p_vec_env <- vector()
  p_vec_cel <- vector()
  p_vec_org <- vector()
  p_vec_hum <- vector()
  p_vec_bri <- vector()
  p_vec_not <- vector()
  for(each in colnames(obs_pi_perms_met[[1]])[2:(length(colnames(obs_pi_perms_met[[1]])) - 1)]){
    p_met <- wilcox.test(obs_pi_perms_met[[x]][obs_pi_perms_met[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_met[[x]]) == each], obs_pi_perms_met[[x]][obs_pi_perms_met[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_met[[x]]) == each])[[3]][1]
    p_t_met <- log10(p_met) * sign(mean(obs_pi_perms_met[[x]][obs_pi_perms_met[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_met[[x]]) == each])- mean(obs_pi_perms_met[[x]][obs_pi_perms_met[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_met[[x]]) == each]))
    p_vec_met <- append(p_vec_met, p_t_met)
  }
  for(each in colnames(obs_pi_perms_gen[[1]])[2:(length(colnames(obs_pi_perms_gen[[1]])) - 1)]){
    p_gen <- wilcox.test(obs_pi_perms_gen[[x]][obs_pi_perms_gen[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_gen[[x]]) == each], obs_pi_perms_gen[[x]][obs_pi_perms_gen[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_gen[[x]]) == each])[[3]][1]
    p_t_gen <- log10(p_gen) * sign(mean(obs_pi_perms_gen[[x]][obs_pi_perms_gen[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_gen[[x]]) == each])- mean(obs_pi_perms_gen[[x]][obs_pi_perms_gen[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_gen[[x]]) == each]))
    p_vec_gen <- append(p_vec_gen, p_t_gen)
  }
  for(each in colnames(obs_pi_perms_env[[1]])[2:(length(colnames(obs_pi_perms_env[[1]])) - 1)]){
    p_env <- wilcox.test(obs_pi_perms_env[[x]][obs_pi_perms_env[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_env[[x]]) == each], obs_pi_perms_env[[x]][obs_pi_perms_env[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_env[[x]]) == each])[[3]][1]
    p_t_env <- log10(p_env) * sign(mean(obs_pi_perms_env[[x]][obs_pi_perms_env[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_env[[x]]) == each])- mean(obs_pi_perms_env[[x]][obs_pi_perms_env[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_env[[x]]) == each]))
    p_vec_env <- append(p_vec_env, p_t_env)
  }
  for(each in colnames(obs_pi_perms_cel[[1]])[2:(length(colnames(obs_pi_perms_cel[[1]])) - 1)]){
    p_cel <- wilcox.test(obs_pi_perms_cel[[x]][obs_pi_perms_cel[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_cel[[x]]) == each], obs_pi_perms_cel[[x]][obs_pi_perms_cel[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_cel[[x]]) == each])[[3]][1]
    p_t_cel <- log10(p_cel) * sign(mean(obs_pi_perms_cel[[x]][obs_pi_perms_cel[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_cel[[x]]) == each])- mean(obs_pi_perms_cel[[x]][obs_pi_perms_cel[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_cel[[x]]) == each]))
    p_vec_cel <- append(p_vec_cel, p_t_cel)
  }
  for(each in colnames(obs_pi_perms_org[[1]])[2:(length(colnames(obs_pi_perms_org[[1]])) - 1)]){
    p_org <- wilcox.test(obs_pi_perms_org[[x]][obs_pi_perms_org[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_org[[x]]) == each], obs_pi_perms_org[[x]][obs_pi_perms_org[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_org[[x]]) == each])[[3]][1]
    p_t_org <- log10(p_org) * sign(mean(obs_pi_perms_org[[x]][obs_pi_perms_org[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_org[[x]]) == each])- mean(obs_pi_perms_org[[x]][obs_pi_perms_org[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_org[[x]]) == each]))
    p_vec_org <- append(p_vec_org, p_t_org) 
  }
  for(each in colnames(obs_pi_perms_hum[[1]])[2:(length(colnames(obs_pi_perms_hum[[1]])) - 1)]){
    p_hum <- wilcox.test(obs_pi_perms_hum[[x]][obs_pi_perms_hum[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_hum[[x]]) == each], obs_pi_perms_hum[[x]][obs_pi_perms_hum[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_hum[[x]]) == each])[[3]][1]
    p_t_hum <- log10(p_hum) * sign(mean(obs_pi_perms_hum[[x]][obs_pi_perms_hum[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_hum[[x]]) == each])- mean(obs_pi_perms_hum[[x]][obs_pi_perms_hum[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_hum[[x]]) == each]))
    p_vec_hum <- append(p_vec_hum, p_t_hum)
  }
  for(each in colnames(obs_pi_perms_bri[[1]])[2:(length(colnames(obs_pi_perms_bri[[1]])) - 1)]){
    p_bri <- wilcox.test(obs_pi_perms_bri[[x]][obs_pi_perms_bri[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_bri[[x]]) == each], obs_pi_perms_bri[[x]][obs_pi_perms_bri[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_bri[[x]]) == each])[[3]][1]
    p_t_bri <- log10(p_bri) * sign(mean(obs_pi_perms_bri[[x]][obs_pi_perms_bri[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_bri[[x]]) == each])- mean(obs_pi_perms_bri[[x]][obs_pi_perms_bri[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_bri[[x]]) == each]))
    p_vec_bri <- append(p_vec_bri, p_t_bri)
  }
  for(each in colnames(obs_pi_perms_not[[1]])[2:(length(colnames(obs_pi_perms_not[[1]])) - 1)]){
    p_not <- wilcox.test(obs_pi_perms_not[[x]][obs_pi_perms_not[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_not[[x]]) == each], obs_pi_perms_not[[x]][obs_pi_perms_not[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_not[[x]]) == each])[[3]][1]
    p_t_not <- log10(p_not) * sign(mean(obs_pi_perms_not[[x]][obs_pi_perms_not[[x]]$C_CASE2 == 1, colnames(obs_pi_perms_not[[x]]) == each])- mean(obs_pi_perms_not[[x]][obs_pi_perms_not[[x]]$C_CASE2 == 0, colnames(obs_pi_perms_not[[x]]) == each]))
    p_vec_not <- append(p_vec_not, p_t_not)
  }
  obs_pi_perms_p_met <- c(obs_pi_perms_p_met, list(p_vec_met))
  obs_pi_perms_p_gen <- c(obs_pi_perms_p_gen, list(p_vec_gen))
  obs_pi_perms_p_env <- c(obs_pi_perms_p_env, list(p_vec_env))
  obs_pi_perms_p_cel <- c(obs_pi_perms_p_cel, list(p_vec_cel)) 
  obs_pi_perms_p_org <- c(obs_pi_perms_p_org, list(p_vec_org))
  obs_pi_perms_p_hum <- c(obs_pi_perms_p_hum, list(p_vec_hum))
  obs_pi_perms_p_bri <- c(obs_pi_perms_p_bri, list(p_vec_bri))
  obs_pi_perms_p_not <- c(obs_pi_perms_p_not, list(p_vec_not))
}

for(x in 1:100){
  p_vec_met <- vector()
  p_vec_gen <- vector()
  p_vec_env <- vector()
  p_vec_cel <- vector()
  p_vec_org <- vector()
  p_vec_hum <- vector()
  p_vec_bri <- vector()
  p_vec_not <- vector()
  for(each in colnames(obs_tax_perms_met[[1]])[2:(length(colnames(obs_tax_perms_met[[1]])) - 1)]){
    p_met <- wilcox.test(obs_tax_perms_met[[x]][obs_tax_perms_met[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_met[[x]]) == each], obs_tax_perms_met[[x]][obs_tax_perms_met[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_met[[x]]) == each])[[3]][1]
    p_t_met <- log10(p_met) * sign(mean(obs_tax_perms_met[[x]][obs_tax_perms_met[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_met[[x]]) == each])- mean(obs_tax_perms_met[[x]][obs_tax_perms_met[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_met[[x]]) == each]))
    p_vec_met <- append(p_vec_met, p_t_met)
  }
  for(each in colnames(obs_tax_perms_gen[[1]])[2:(length(colnames(obs_tax_perms_gen[[1]])) - 1)]){
    p_gen <- wilcox.test(obs_tax_perms_gen[[x]][obs_tax_perms_gen[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_gen[[x]]) == each], obs_tax_perms_gen[[x]][obs_tax_perms_gen[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_gen[[x]]) == each])[[3]][1]
    p_t_gen <- log10(p_gen) * sign(mean(obs_tax_perms_gen[[x]][obs_tax_perms_gen[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_gen[[x]]) == each])- mean(obs_tax_perms_gen[[x]][obs_tax_perms_gen[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_gen[[x]]) == each]))
    p_vec_gen <- append(p_vec_gen, p_t_gen)
  }
  for(each in colnames(obs_tax_perms_env[[1]])[2:(length(colnames(obs_tax_perms_env[[1]])) - 1)]){
    p_env <- wilcox.test(obs_tax_perms_env[[x]][obs_tax_perms_env[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_env[[x]]) == each], obs_tax_perms_env[[x]][obs_tax_perms_env[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_env[[x]]) == each])[[3]][1]
    p_t_env <- log10(p_env) * sign(mean(obs_tax_perms_env[[x]][obs_tax_perms_env[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_env[[x]]) == each])- mean(obs_tax_perms_env[[x]][obs_tax_perms_env[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_env[[x]]) == each]))
    p_vec_env <- append(p_vec_env, p_t_env)
  }
  for(each in colnames(obs_tax_perms_cel[[1]])[2:(length(colnames(obs_tax_perms_cel[[1]])) - 1)]){
    p_cel <- wilcox.test(obs_tax_perms_cel[[x]][obs_tax_perms_cel[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_cel[[x]]) == each], obs_tax_perms_cel[[x]][obs_tax_perms_cel[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_cel[[x]]) == each])[[3]][1]
    p_t_cel <- log10(p_cel) * sign(mean(obs_tax_perms_cel[[x]][obs_tax_perms_cel[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_cel[[x]]) == each])- mean(obs_tax_perms_cel[[x]][obs_tax_perms_cel[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_cel[[x]]) == each]))
    p_vec_cel <- append(p_vec_cel, p_t_cel)
  }
  for(each in colnames(obs_tax_perms_org[[1]])[2:(length(colnames(obs_tax_perms_org[[1]])) - 1)]){
    p_org <- wilcox.test(obs_tax_perms_org[[x]][obs_tax_perms_org[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_org[[x]]) == each], obs_tax_perms_org[[x]][obs_tax_perms_org[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_org[[x]]) == each])[[3]][1]
    p_t_org <- log10(p_org) * sign(mean(obs_tax_perms_org[[x]][obs_tax_perms_org[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_org[[x]]) == each])- mean(obs_tax_perms_org[[x]][obs_tax_perms_org[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_org[[x]]) == each]))
    p_vec_org <- append(p_vec_org, p_t_org) 
  }
  for(each in colnames(obs_tax_perms_hum[[1]])[2:(length(colnames(obs_tax_perms_hum[[1]])) - 1)]){
    p_hum <- wilcox.test(obs_tax_perms_hum[[x]][obs_tax_perms_hum[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_hum[[x]]) == each], obs_tax_perms_hum[[x]][obs_tax_perms_hum[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_hum[[x]]) == each])[[3]][1]
    p_t_hum <- log10(p_hum) * sign(mean(obs_tax_perms_hum[[x]][obs_tax_perms_hum[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_hum[[x]]) == each])- mean(obs_tax_perms_hum[[x]][obs_tax_perms_hum[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_hum[[x]]) == each]))
    p_vec_hum <- append(p_vec_hum, p_t_hum)
  }
  for(each in colnames(obs_tax_perms_bri[[1]])[2:(length(colnames(obs_tax_perms_bri[[1]])) - 1)]){
    p_bri <- wilcox.test(obs_tax_perms_bri[[x]][obs_tax_perms_bri[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_bri[[x]]) == each], obs_tax_perms_bri[[x]][obs_tax_perms_bri[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_bri[[x]]) == each])[[3]][1]
    p_t_bri <- log10(p_bri) * sign(mean(obs_tax_perms_bri[[x]][obs_tax_perms_bri[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_bri[[x]]) == each])- mean(obs_tax_perms_bri[[x]][obs_tax_perms_bri[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_bri[[x]]) == each]))
    p_vec_bri <- append(p_vec_bri, p_t_bri)
  }
  for(each in colnames(obs_tax_perms_not[[1]])[2:(length(colnames(obs_tax_perms_not[[1]])) - 1)]){
    p_not <- wilcox.test(obs_tax_perms_not[[x]][obs_tax_perms_not[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_not[[x]]) == each], obs_tax_perms_not[[x]][obs_tax_perms_not[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_not[[x]]) == each])[[3]][1]
    p_t_not <- log10(p_not) * sign(mean(obs_tax_perms_not[[x]][obs_tax_perms_not[[x]]$C_CASE2 == 1, colnames(obs_tax_perms_not[[x]]) == each])- mean(obs_tax_perms_not[[x]][obs_tax_perms_not[[x]]$C_CASE2 == 0, colnames(obs_tax_perms_not[[x]]) == each]))
    p_vec_not <- append(p_vec_not, p_t_not)
  }
  obs_tax_perms_p_met <- c(obs_tax_perms_p_met, list(p_vec_met))
  obs_tax_perms_p_gen <- c(obs_tax_perms_p_gen, list(p_vec_gen))
  obs_tax_perms_p_env <- c(obs_tax_perms_p_env, list(p_vec_env))
  obs_tax_perms_p_cel <- c(obs_tax_perms_p_cel, list(p_vec_cel)) 
  obs_tax_perms_p_org <- c(obs_tax_perms_p_org, list(p_vec_org))
  obs_tax_perms_p_hum <- c(obs_tax_perms_p_hum, list(p_vec_hum))
  obs_tax_perms_p_bri <- c(obs_tax_perms_p_bri, list(p_vec_bri))
  obs_tax_perms_p_not <- c(obs_tax_perms_p_not, list(p_vec_not))
}

for(x in 1:100){
  p_vec_met <- vector()
  p_vec_gen <- vector()
  p_vec_env <- vector()
  p_vec_cel <- vector()
  p_vec_org <- vector()
  p_vec_hum <- vector()
  p_vec_bri <- vector()
  p_vec_not <- vector()
  for(each in colnames(pi_obs_perms_met[[1]])[2:(length(colnames(pi_obs_perms_met[[1]])) - 1)]){
    p_met <- wilcox.test(pi_obs_perms_met[[x]][pi_obs_perms_met[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_met[[x]]) == each], pi_obs_perms_met[[x]][pi_obs_perms_met[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_met[[x]]) == each])[[3]][1]
    p_t_met <- log10(p_met) * sign(mean(pi_obs_perms_met[[x]][pi_obs_perms_met[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_met[[x]]) == each])- mean(pi_obs_perms_met[[x]][pi_obs_perms_met[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_met[[x]]) == each]))
    p_vec_met <- append(p_vec_met, p_t_met)
  }
  for(each in colnames(pi_obs_perms_gen[[1]])[2:(length(colnames(pi_obs_perms_gen[[1]])) - 1)]){
    p_gen <- wilcox.test(pi_obs_perms_gen[[x]][pi_obs_perms_gen[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_gen[[x]]) == each], pi_obs_perms_gen[[x]][pi_obs_perms_gen[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_gen[[x]]) == each])[[3]][1]
    p_t_gen <- log10(p_gen) * sign(mean(pi_obs_perms_gen[[x]][pi_obs_perms_gen[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_gen[[x]]) == each])- mean(pi_obs_perms_gen[[x]][pi_obs_perms_gen[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_gen[[x]]) == each]))
    p_vec_gen <- append(p_vec_gen, p_t_gen)
  }
  for(each in colnames(pi_obs_perms_env[[1]])[2:(length(colnames(pi_obs_perms_env[[1]])) - 1)]){
    p_env <- wilcox.test(pi_obs_perms_env[[x]][pi_obs_perms_env[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_env[[x]]) == each], pi_obs_perms_env[[x]][pi_obs_perms_env[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_env[[x]]) == each])[[3]][1]
    p_t_env <- log10(p_env) * sign(mean(pi_obs_perms_env[[x]][pi_obs_perms_env[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_env[[x]]) == each])- mean(pi_obs_perms_env[[x]][pi_obs_perms_env[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_env[[x]]) == each]))
    p_vec_env <- append(p_vec_env, p_t_env)
  }
  for(each in colnames(pi_obs_perms_cel[[1]])[2:(length(colnames(pi_obs_perms_cel[[1]])) - 1)]){
    p_cel <- wilcox.test(pi_obs_perms_cel[[x]][pi_obs_perms_cel[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_cel[[x]]) == each], pi_obs_perms_cel[[x]][pi_obs_perms_cel[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_cel[[x]]) == each])[[3]][1]
    p_t_cel <- log10(p_cel) * sign(mean(pi_obs_perms_cel[[x]][pi_obs_perms_cel[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_cel[[x]]) == each])- mean(pi_obs_perms_cel[[x]][pi_obs_perms_cel[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_cel[[x]]) == each]))
    p_vec_cel <- append(p_vec_cel, p_t_cel)
  }
  for(each in colnames(pi_obs_perms_org[[1]])[2:(length(colnames(pi_obs_perms_org[[1]])) - 1)]){
    p_org <- wilcox.test(pi_obs_perms_org[[x]][pi_obs_perms_org[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_org[[x]]) == each], pi_obs_perms_org[[x]][pi_obs_perms_org[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_org[[x]]) == each])[[3]][1]
    p_t_org <- log10(p_org) * sign(mean(pi_obs_perms_org[[x]][pi_obs_perms_org[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_org[[x]]) == each])- mean(pi_obs_perms_org[[x]][pi_obs_perms_org[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_org[[x]]) == each]))
    p_vec_org <- append(p_vec_org, p_t_org) 
  }
  for(each in colnames(pi_obs_perms_hum[[1]])[2:(length(colnames(pi_obs_perms_hum[[1]])) - 1)]){
    p_hum <- wilcox.test(pi_obs_perms_hum[[x]][pi_obs_perms_hum[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_hum[[x]]) == each], pi_obs_perms_hum[[x]][pi_obs_perms_hum[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_hum[[x]]) == each])[[3]][1]
    p_t_hum <- log10(p_hum) * sign(mean(pi_obs_perms_hum[[x]][pi_obs_perms_hum[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_hum[[x]]) == each])- mean(pi_obs_perms_hum[[x]][pi_obs_perms_hum[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_hum[[x]]) == each]))
    p_vec_hum <- append(p_vec_hum, p_t_hum)
  }
  for(each in colnames(pi_obs_perms_bri[[1]])[2:(length(colnames(pi_obs_perms_bri[[1]])) - 1)]){
    p_bri <- wilcox.test(pi_obs_perms_bri[[x]][pi_obs_perms_bri[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_bri[[x]]) == each], pi_obs_perms_bri[[x]][pi_obs_perms_bri[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_bri[[x]]) == each])[[3]][1]
    p_t_bri <- log10(p_bri) * sign(mean(pi_obs_perms_bri[[x]][pi_obs_perms_bri[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_bri[[x]]) == each])- mean(pi_obs_perms_bri[[x]][pi_obs_perms_bri[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_bri[[x]]) == each]))
    p_vec_bri <- append(p_vec_bri, p_t_bri)
  }
  for(each in colnames(pi_obs_perms_not[[1]])[2:(length(colnames(pi_obs_perms_not[[1]])) - 1)]){
    p_not <- wilcox.test(pi_obs_perms_not[[x]][pi_obs_perms_not[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_not[[x]]) == each], pi_obs_perms_not[[x]][pi_obs_perms_not[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_not[[x]]) == each])[[3]][1]
    p_t_not <- log10(p_not) * sign(mean(pi_obs_perms_not[[x]][pi_obs_perms_not[[x]]$C_CASE2 == 1, colnames(pi_obs_perms_not[[x]]) == each])- mean(pi_obs_perms_not[[x]][pi_obs_perms_not[[x]]$C_CASE2 == 0, colnames(pi_obs_perms_not[[x]]) == each]))
    p_vec_not <- append(p_vec_not, p_t_not)
  }
  pi_obs_perms_p_met <- c(pi_obs_perms_p_met, list(p_vec_met))
  pi_obs_perms_p_gen <- c(pi_obs_perms_p_gen, list(p_vec_gen))
  pi_obs_perms_p_env <- c(pi_obs_perms_p_env, list(p_vec_env))
  pi_obs_perms_p_cel <- c(pi_obs_perms_p_cel, list(p_vec_cel)) 
  pi_obs_perms_p_org <- c(pi_obs_perms_p_org, list(p_vec_org))
  pi_obs_perms_p_hum <- c(pi_obs_perms_p_hum, list(p_vec_hum))
  pi_obs_perms_p_bri <- c(pi_obs_perms_p_bri, list(p_vec_bri))
  pi_obs_perms_p_not <- c(pi_obs_perms_p_not, list(p_vec_not))
}

for(x in 1:100){
  p_vec_met <- vector()
  p_vec_gen <- vector()
  p_vec_env <- vector()
  p_vec_cel <- vector()
  p_vec_org <- vector()
  p_vec_hum <- vector()
  p_vec_bri <- vector()
  p_vec_not <- vector()
  for(each in colnames(tax_obs_perms_met[[1]])[2:(length(colnames(tax_obs_perms_met[[1]])) - 1)]){
    p_met <- wilcox.test(tax_obs_perms_met[[x]][tax_obs_perms_met[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_met[[x]]) == each], tax_obs_perms_met[[x]][tax_obs_perms_met[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_met[[x]]) == each])[[3]][1]
    p_t_met <- log10(p_met) * sign(mean(tax_obs_perms_met[[x]][tax_obs_perms_met[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_met[[x]]) == each])- mean(tax_obs_perms_met[[x]][tax_obs_perms_met[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_met[[x]]) == each]))
    p_vec_met <- append(p_vec_met, p_t_met)
  }
  for(each in colnames(tax_obs_perms_gen[[1]])[2:(length(colnames(tax_obs_perms_gen[[1]])) - 1)]){
    p_gen <- wilcox.test(tax_obs_perms_gen[[x]][tax_obs_perms_gen[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_gen[[x]]) == each], tax_obs_perms_gen[[x]][tax_obs_perms_gen[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_gen[[x]]) == each])[[3]][1]
    p_t_gen <- log10(p_gen) * sign(mean(tax_obs_perms_gen[[x]][tax_obs_perms_gen[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_gen[[x]]) == each])- mean(tax_obs_perms_gen[[x]][tax_obs_perms_gen[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_gen[[x]]) == each]))
    p_vec_gen <- append(p_vec_gen, p_t_gen)
  }
  for(each in colnames(tax_obs_perms_env[[1]])[2:(length(colnames(tax_obs_perms_env[[1]])) - 1)]){
    p_env <- wilcox.test(tax_obs_perms_env[[x]][tax_obs_perms_env[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_env[[x]]) == each], tax_obs_perms_env[[x]][tax_obs_perms_env[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_env[[x]]) == each])[[3]][1]
    p_t_env <- log10(p_env) * sign(mean(tax_obs_perms_env[[x]][tax_obs_perms_env[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_env[[x]]) == each])- mean(tax_obs_perms_env[[x]][tax_obs_perms_env[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_env[[x]]) == each]))
    p_vec_env <- append(p_vec_env, p_t_env)
  }
  for(each in colnames(tax_obs_perms_cel[[1]])[2:(length(colnames(tax_obs_perms_cel[[1]])) - 1)]){
    p_cel <- wilcox.test(tax_obs_perms_cel[[x]][tax_obs_perms_cel[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_cel[[x]]) == each], tax_obs_perms_cel[[x]][tax_obs_perms_cel[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_cel[[x]]) == each])[[3]][1]
    p_t_cel <- log10(p_cel) * sign(mean(tax_obs_perms_cel[[x]][tax_obs_perms_cel[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_cel[[x]]) == each])- mean(tax_obs_perms_cel[[x]][tax_obs_perms_cel[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_cel[[x]]) == each]))
    p_vec_cel <- append(p_vec_cel, p_t_cel)
  }
  for(each in colnames(tax_obs_perms_org[[1]])[2:(length(colnames(tax_obs_perms_org[[1]])) - 1)]){
    p_org <- wilcox.test(tax_obs_perms_org[[x]][tax_obs_perms_org[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_org[[x]]) == each], tax_obs_perms_org[[x]][tax_obs_perms_org[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_org[[x]]) == each])[[3]][1]
    p_t_org <- log10(p_org) * sign(mean(tax_obs_perms_org[[x]][tax_obs_perms_org[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_org[[x]]) == each])- mean(tax_obs_perms_org[[x]][tax_obs_perms_org[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_org[[x]]) == each]))
    p_vec_org <- append(p_vec_org, p_t_org) 
  }
  for(each in colnames(tax_obs_perms_hum[[1]])[2:(length(colnames(tax_obs_perms_hum[[1]])) - 1)]){
    p_hum <- wilcox.test(tax_obs_perms_hum[[x]][tax_obs_perms_hum[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_hum[[x]]) == each], tax_obs_perms_hum[[x]][tax_obs_perms_hum[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_hum[[x]]) == each])[[3]][1]
    p_t_hum <- log10(p_hum) * sign(mean(tax_obs_perms_hum[[x]][tax_obs_perms_hum[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_hum[[x]]) == each])- mean(tax_obs_perms_hum[[x]][tax_obs_perms_hum[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_hum[[x]]) == each]))
    p_vec_hum <- append(p_vec_hum, p_t_hum)
  }
  for(each in colnames(tax_obs_perms_bri[[1]])[2:(length(colnames(tax_obs_perms_bri[[1]])) - 1)]){
    p_bri <- wilcox.test(tax_obs_perms_bri[[x]][tax_obs_perms_bri[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_bri[[x]]) == each], tax_obs_perms_bri[[x]][tax_obs_perms_bri[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_bri[[x]]) == each])[[3]][1]
    p_t_bri <- log10(p_bri) * sign(mean(tax_obs_perms_bri[[x]][tax_obs_perms_bri[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_bri[[x]]) == each])- mean(tax_obs_perms_bri[[x]][tax_obs_perms_bri[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_bri[[x]]) == each]))
    p_vec_bri <- append(p_vec_bri, p_t_bri)
  }
  for(each in colnames(tax_obs_perms_not[[1]])[2:(length(colnames(tax_obs_perms_not[[1]])) - 1)]){
    p_not <- wilcox.test(tax_obs_perms_not[[x]][tax_obs_perms_not[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_not[[x]]) == each], tax_obs_perms_not[[x]][tax_obs_perms_not[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_not[[x]]) == each])[[3]][1]
    p_t_not <- log10(p_not) * sign(mean(tax_obs_perms_not[[x]][tax_obs_perms_not[[x]]$C_CASE2 == 1, colnames(tax_obs_perms_not[[x]]) == each])- mean(tax_obs_perms_not[[x]][tax_obs_perms_not[[x]]$C_CASE2 == 0, colnames(tax_obs_perms_not[[x]]) == each]))
    p_vec_not <- append(p_vec_not, p_t_not)
  }
  tax_obs_perms_p_met <- c(tax_obs_perms_p_met, list(p_vec_met))
  tax_obs_perms_p_gen <- c(tax_obs_perms_p_gen, list(p_vec_gen))
  tax_obs_perms_p_env <- c(tax_obs_perms_p_env, list(p_vec_env))
  tax_obs_perms_p_cel <- c(tax_obs_perms_p_cel, list(p_vec_cel)) 
  tax_obs_perms_p_org <- c(tax_obs_perms_p_org, list(p_vec_org))
  tax_obs_perms_p_hum <- c(tax_obs_perms_p_hum, list(p_vec_hum))
  tax_obs_perms_p_bri <- c(tax_obs_perms_p_bri, list(p_vec_bri))
  tax_obs_perms_p_not <- c(tax_obs_perms_p_not, list(p_vec_not))
}



### generate vectors of correlation coefficients between transforned p values from tests 
### comparing preterm birth cases to term birth controls, across 100 permuted dataframes,
### restricted by level-1 KO functional category

obs_pi_perms_p_corr_met <- vector()
obs_pi_perms_p_corr_gen <- vector()
obs_pi_perms_p_corr_env <- vector()
obs_pi_perms_p_corr_cel <- vector()
obs_pi_perms_p_corr_org <- vector()
obs_pi_perms_p_corr_hum <- vector()
obs_pi_perms_p_corr_bri <- vector()
obs_pi_perms_p_corr_not <- vector()

obs_tax_perms_p_corr_met <- vector()
obs_tax_perms_p_corr_gen <- vector()
obs_tax_perms_p_corr_env <- vector()
obs_tax_perms_p_corr_cel <- vector()
obs_tax_perms_p_corr_org <- vector()
obs_tax_perms_p_corr_hum <- vector()
obs_tax_perms_p_corr_bri <- vector()
obs_tax_perms_p_corr_not <- vector()
for(x in 1:100){
  corr_met <- cor(obs_pi_perms_p_met[[x]], pi_obs_perms_p_met[[x]], method = "spearman")
  obs_pi_perms_p_corr_met <- append(obs_pi_perms_p_corr_met, corr_met)
  corr_gen <- cor(obs_pi_perms_p_gen[[x]], pi_obs_perms_p_gen[[x]], method = "spearman")
  obs_pi_perms_p_corr_gen <- append(obs_pi_perms_p_corr_gen, corr_gen)
  corr_env <- cor(obs_pi_perms_p_env[[x]], pi_obs_perms_p_env[[x]], method = "spearman")
  obs_pi_perms_p_corr_env <- append(obs_pi_perms_p_corr_env, corr_env)
  corr_cel <- cor(obs_pi_perms_p_cel[[x]], pi_obs_perms_p_cel[[x]], method = "spearman")
  obs_pi_perms_p_corr_cel <- append(obs_pi_perms_p_corr_cel, corr_cel)
  corr_org <- cor(obs_pi_perms_p_org[[x]], pi_obs_perms_p_org[[x]], method = "spearman")
  obs_pi_perms_p_corr_org <- append(obs_pi_perms_p_corr_org, corr_org)
  corr_hum <- cor(obs_pi_perms_p_hum[[x]], pi_obs_perms_p_hum[[x]], method = "spearman")
  obs_pi_perms_p_corr_hum <- append(obs_pi_perms_p_corr_hum, corr_hum)
  corr_bri <- cor(obs_pi_perms_p_bri[[x]], pi_obs_perms_p_bri[[x]], method = "spearman")
  obs_pi_perms_p_corr_bri <- append(obs_pi_perms_p_corr_bri, corr_bri)
  corr_not <- cor(obs_pi_perms_p_not[[x]], pi_obs_perms_p_not[[x]], method = "spearman")
  obs_pi_perms_p_corr_not <- append(obs_pi_perms_p_corr_not, corr_not)
  
  corr_met <- cor(obs_tax_perms_p_met[[x]], tax_obs_perms_p_met[[x]], method = "spearman")
  obs_tax_perms_p_corr_met <- append(obs_tax_perms_p_corr_met, corr_met)
  corr_gen <- cor(obs_tax_perms_p_gen[[x]], tax_obs_perms_p_gen[[x]], method = "spearman")
  obs_tax_perms_p_corr_gen <- append(obs_tax_perms_p_corr_gen, corr_gen)
  corr_env <- cor(obs_tax_perms_p_env[[x]], tax_obs_perms_p_env[[x]], method = "spearman")
  obs_tax_perms_p_corr_env <- append(obs_tax_perms_p_corr_env, corr_env)
  corr_cel <- cor(obs_tax_perms_p_cel[[x]], tax_obs_perms_p_cel[[x]], method = "spearman")
  obs_tax_perms_p_corr_cel <- append(obs_tax_perms_p_corr_cel, corr_cel)
  corr_org <- cor(obs_tax_perms_p_org[[x]], tax_obs_perms_p_org[[x]], method = "spearman")
  obs_tax_perms_p_corr_org <- append(obs_tax_perms_p_corr_org, corr_org)
  corr_hum <- cor(obs_tax_perms_p_hum[[x]], tax_obs_perms_p_hum[[x]], method = "spearman")
  obs_tax_perms_p_corr_hum <- append(obs_tax_perms_p_corr_hum, corr_hum)
  corr_bri <- cor(obs_tax_perms_p_bri[[x]], tax_obs_perms_p_bri[[x]], method = "spearman")
  obs_tax_perms_p_corr_bri <- append(obs_tax_perms_p_corr_bri, corr_bri)
  corr_not <- cor(obs_tax_perms_p_not[[x]], tax_obs_perms_p_not[[x]], method = "spearman")
  obs_tax_perms_p_corr_not <- append(obs_tax_perms_p_corr_not, corr_not)
}
summary(obs_pi_perms_p_corr_met)
summary(obs_pi_perms_p_corr_gen)
summary(obs_pi_perms_p_corr_env)
summary(obs_pi_perms_p_corr_cel)
summary(obs_pi_perms_p_corr_org)
summary(obs_pi_perms_p_corr_hum)
summary(obs_pi_perms_p_corr_bri)
summary(obs_pi_perms_p_corr_not)

summary(obs_tax_perms_p_corr_met)
summary(obs_tax_perms_p_corr_gen)
summary(obs_tax_perms_p_corr_env)
summary(obs_tax_perms_p_corr_cel)
summary(obs_tax_perms_p_corr_org)
summary(obs_tax_perms_p_corr_hum)
summary(obs_tax_perms_p_corr_bri)
summary(obs_tax_perms_p_corr_not)



### generate dataframe for plotting relative abundance correlation coefficients, 
### unpermuted and permuted, from KO functional category stratified analysis

plot_pi_un <- data.frame(matrix(ncol = 0, nrow = length(obs_pi_corr)*2))
plot_pi_un$corr <- rep(obs_pi_corr, 2)
plot_pi_un$method <- "PICRUSt2"
plot_pi_un$perm <- "Unpermuted"
plot_pi_un$method_perm <- "PICRUSt2,  Unpermuted"
plot_pi_un$cat <- c(obs_pi_corr_cat, rep("All", length(obs_pi_corr)))

plot_tax_un <- data.frame(matrix(ncol = 0, nrow = length(obs_tax_corr)*2))
plot_tax_un$corr <- rep(obs_tax_corr, 2)
plot_tax_un$method <- "Tax4Fun2"
plot_tax_un$perm <- "Unpermuted"
plot_tax_un$method_perm <- "Tax4Fun2,  Unpermuted"
plot_tax_un$cat <- c(obs_tax_corr_cat, rep("All", length(obs_tax_corr)))

plot_pi_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_pi_perms_corr)*2))
plot_pi_perm$corr <- rep(c(obs_pi_perms_corr_met, 
                         obs_pi_perms_corr_gen, 
                         obs_pi_perms_corr_env, 
                         obs_pi_perms_corr_cel, 
                         obs_pi_perms_corr_org, 
                         obs_pi_perms_corr_hum, 
                         obs_pi_perms_corr_bri, 
                         obs_pi_perms_corr_not), 
                       2)
plot_pi_perm$method <- "PICRUSt2"
plot_pi_perm$perm <- "Permuted"
plot_pi_perm$method_perm <- "PICRUSt2,  Permuted"
plot_pi_perm$cat <- c(rep("Metabolism", length(obs_pi_perms_corr_met)), 
                    rep("Genetic information processing", length(obs_pi_perms_corr_gen)), 
                    rep("Environmental information processing", length(obs_pi_perms_corr_env)), 
                    rep("Cellular processes", length(obs_pi_perms_corr_cel)), 
                    rep("Organismal systems", length(obs_pi_perms_corr_org)), 
                    rep("Human diseases", length(obs_pi_perms_corr_hum)), 
                    rep("BRITE hierarchies", length(obs_pi_perms_corr_bri)), 
                    rep("Not included in pathway or BRITE", length(obs_pi_perms_corr_not)), 
                    rep("All", length(obs_pi_perms_corr)))

plot_tax_perm <- data.frame(matrix(ncol = 0, nrow = length(obs_tax_perms_corr)*2))
plot_tax_perm$corr <- rep(c(obs_tax_perms_corr_met, 
                             obs_tax_perms_corr_gen, 
                             obs_tax_perms_corr_env, 
                             obs_tax_perms_corr_cel, 
                             obs_tax_perms_corr_org, 
                             obs_tax_perms_corr_hum, 
                             obs_tax_perms_corr_bri, 
                             obs_tax_perms_corr_not), 
                           2)
plot_tax_perm$method <- "Tax4Fun2"
plot_tax_perm$perm <- "Permuted"
plot_tax_perm$method_perm <- "Tax4Fun2,  Permuted"
plot_tax_perm$cat <- c(rep("Metabolism", length(obs_tax_perms_corr_met)), 
                        rep("Genetic information processing", length(obs_tax_perms_corr_gen)), 
                        rep("Environmental information processing", length(obs_tax_perms_corr_env)), 
                        rep("Cellular processes", length(obs_tax_perms_corr_cel)), 
                        rep("Organismal systems", length(obs_tax_perms_corr_org)), 
                        rep("Human diseases", length(obs_tax_perms_corr_hum)), 
                        rep("BRITE hierarchies", length(obs_tax_perms_corr_bri)), 
                        rep("Not included in pathway or BRITE", length(obs_tax_perms_corr_not)), 
                        rep("All", length(obs_tax_perms_corr)))

plot_ra_ko <- as.data.frame(rbind(plot_pi_un, plot_tax_un, plot_pi_perm, plot_tax_perm))
# plot_ra_ko should have 875872 obs of 5 vars

# make KO functional category variable into factor
# some values are numeric from obs_pi_corr_cat and obs_tax_corr_cat
# so need to make those into corresponding character values first
plot_ra_ko <- plot_ra_ko %>% mutate(cat = case_when(cat == "1" ~ "Metabolism",
                                                        cat == "2" ~ "Genetic information processing",
                                                        cat == "3" ~ "Environmental information processing",
                                                        cat == "4" ~ "Cellular processes",
                                                        cat == "5" ~ "Organismal systems",
                                                        cat == "6" ~ "Human diseases",
                                                        cat == "7" ~ "BRITE hierarchies",
                                                        cat == "8" ~ "Not included in pathway or BRITE",
                                                        cat == "Metabolism" ~ "Metabolism",
                                                        cat == "Genetic information processing" ~ "Genetic information processing",
                                                        cat == "Environmental information processing" ~ "Environmental information processing",
                                                        cat == "Cellular processes" ~ "Cellular processes",
                                                        cat == "Organismal systems" ~ "Organismal systems",
                                                        cat == "Human diseases" ~ "Human diseases",
                                                        cat == "BRITE hierarchies" ~ "BRITE hierarchies",
                                                        cat == "Not included in pathway or BRITE" ~ "Not included in pathway or BRITE",
                                                        cat == "All" ~ "All"))
plot_ra_ko$cat <- factor(plot_ra_ko$cat,
                            levels = c("All",
                                       "Metabolism",
                                       "Genetic information processing",
                                       "Environmental information processing",
                                       "Cellular processes",
                                       "Organismal systems",
                                       "Human diseases",
                                       "BRITE hierarchies",
                                       "Not included in pathway or BRITE"))



### generate dataframe for plotting p value correlation coefficients, 
### unpermuted and permuted, from KO functional category stratified analysis

plot_pi_un <- data.frame(matrix(ncol = 0, nrow = length(unique(ko_cat$cat)) + 1))
plot_pi_un$corr <- c(obs_pi_p_corr, 
                   obs_pi_p_corr_met, 
                   obs_pi_p_corr_gen, 
                   obs_pi_p_corr_env, 
                   obs_pi_p_corr_cel, 
                   obs_pi_p_corr_org, 
                   obs_pi_p_corr_hum, 
                   obs_pi_p_corr_bri, 
                   obs_pi_p_corr_not)
plot_pi_un$method <- "PICRUSt2"
plot_pi_un$perm <- "Unpermuted"
plot_pi_un$method_perm <- "PICRUSt2,  Unpermuted"
plot_pi_un$cat <- c("All", "Metabolism", "Genetic information processing", "Environmental information processing", 
                  "Cellular processes", "Organismal systems", "Human diseases", "BRITE hierarchies", 
                  "Not included in pathway or BRITE")

plot_tax_un <- data.frame(matrix(ncol = 0, nrow = length(unique(ko_cat$cat)) + 1))
plot_tax_un$corr <- c(obs_tax_p_corr, 
                       obs_tax_p_corr_met, 
                       obs_tax_p_corr_gen, 
                       obs_tax_p_corr_env, 
                       obs_tax_p_corr_cel, 
                       obs_tax_p_corr_org, 
                       obs_tax_p_corr_hum, 
                       obs_tax_p_corr_bri, 
                       obs_tax_p_corr_not)
plot_tax_un$method <- "Tax4Fun2"
plot_tax_un$perm <- "Unpermuted"
plot_tax_un$method_perm <- "Tax4Fun2,  Unpermuted"
plot_tax_un$cat <- c("All", "Metabolism", "Genetic information processing", "Environmental information processing", 
                      "Cellular processes", "Organismal systems", "Human diseases", "BRITE hierarchies", 
                      "Not included in pathway or BRITE")

plot_pi_perm <- data.frame(matrix(ncol = 0, nrow = (length(unique(ko_cat$cat)) + 1) * 100))
plot_pi_perm$corr <- c(obs_pi_perms_p_corr, 
                     obs_pi_perms_p_corr_met, 
                     obs_pi_perms_p_corr_gen, 
                     obs_pi_perms_p_corr_env, 
                     obs_pi_perms_p_corr_cel, 
                     obs_pi_perms_p_corr_org, 
                     obs_pi_perms_p_corr_hum, 
                     obs_pi_perms_p_corr_bri, 
                     obs_pi_perms_p_corr_not)
plot_pi_perm$method <- "PICRUSt2"
plot_pi_perm$perm <- "Permuted"
plot_pi_perm$method_perm <- "PICRUSt2,  Permuted"
plot_pi_perm$cat <- c(rep("All", length(obs_pi_perms_p_corr)), 
                    rep("Metabolism", length(obs_pi_perms_p_corr_met)), 
                    rep("Genetic information processing", length(obs_pi_perms_p_corr_gen)), 
                    rep("Environmental information processing", length(obs_pi_perms_p_corr_env)), 
                    rep("Cellular processes", length(obs_pi_perms_p_corr_cel)), 
                    rep("Organismal systems", length(obs_pi_perms_p_corr_org)), 
                    rep("Human diseases", length(obs_pi_perms_p_corr_hum)), 
                    rep("BRITE hierarchies", length(obs_pi_perms_p_corr_bri)), 
                    rep("Not included in pathway or BRITE", length(obs_pi_perms_p_corr_not)))

plot_tax_perm <- data.frame(matrix(ncol = 0, nrow = (length(unique(ko_cat$cat)) + 1) * 100))
plot_tax_perm$corr <- c(obs_tax_perms_p_corr, 
                         obs_tax_perms_p_corr_met, 
                         obs_tax_perms_p_corr_gen, 
                         obs_tax_perms_p_corr_env, 
                         obs_tax_perms_p_corr_cel, 
                         obs_tax_perms_p_corr_org, 
                         obs_tax_perms_p_corr_hum, 
                         obs_tax_perms_p_corr_bri, 
                         obs_tax_perms_p_corr_not)
plot_tax_perm$method <- "Tax4Fun2"
plot_tax_perm$perm <- "Permuted"
plot_tax_perm$method_perm <- "Tax4Fun2,  Permuted"
plot_tax_perm$cat <- c(rep("All", length(obs_tax_perms_p_corr)), 
                        rep("Metabolism", length(obs_tax_perms_p_corr_met)), 
                        rep("Genetic information processing", length(obs_tax_perms_p_corr_gen)), 
                        rep("Environmental information processing", length(obs_tax_perms_p_corr_env)), 
                        rep("Cellular processes", length(obs_tax_perms_p_corr_cel)), 
                        rep("Organismal systems", length(obs_tax_perms_p_corr_org)), 
                        rep("Human diseases", length(obs_tax_perms_p_corr_hum)), 
                        rep("BRITE hierarchies", length(obs_tax_perms_p_corr_bri)), 
                        rep("Not included in pathway or BRITE", length(obs_tax_perms_p_corr_not)))

plot_pcorr_ko <- as.data.frame(rbind(plot_pi_un, plot_tax_un, plot_pi_perm, plot_tax_perm))
# plot_pcorr_ko should have 1818 obs of 5 vars

# make KO functional category variable into factor
plot_pcorr_ko$cat <- factor(plot_pcorr_ko$cat,
                            levels = c("All",
                                       "Metabolism",
                                       "Genetic information processing",
                                       "Environmental information processing",
                                       "Cellular processes",
                                       "Organismal systems",
                                       "Human diseases",
                                       "BRITE hierarchies",
                                       "Not included in pathway or BRITE"))







#########################################################
### GENERATE FIGURE 4 CLUSTER STRATIFIED RESULTS PLOT ###
#########################################################



### generate plot of relative abundance correlations from cluster stratified analysis

plot_ra_cluster_nolegend <- ggplot(data = plot_ra_cluster, 
                                   aes(x = as.factor(method_perm), y = corr, color = cluster)) + 
  geom_boxplot(outlier.shape = NA, fill = "white", 
               lwd = 0.35, show.legend = F) + 
  stat_boxplot(geom = "errorbar", width = 0.5, lwd = 0.35,
               position = position_dodge(width = 0.75), show.legend = F) + 
  geom_jitter(data = plot_ra_cluster, 
              aes(x = as.factor(method_perm), y = corr, color = as.factor(cluster), alpha = 0.001), 
              shape = 20, cex = 1,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1, jitter.height = 0), ) + 
  xlab("") + 
  ylab("Spearman rho") + 
  scale_color_manual(values = c("black", "#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c("All", 
                              expression(paste(italic("L. crispatus"), " dominated")), 
                              expression(paste(italic("L. iners"), " dominated")), 
                              "Mixed"), 
                     "Cluster") + 
  guides(color = "none", 
         alpha = "none") + 
  scale_y_continuous(breaks = c(-1.0, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00), 
                     limits = c(-1.0, 1.0)) + 
  scale_x_discrete(breaks = c("PICRUSt2,  Permuted", "PICRUSt2,  Unpermuted", 
                              "Tax4Fun2,  Permuted", "Tax4Fun2,  Unpermuted"), 
                   labels = c("PICRUSt2, \nPermuted", "PICRUSt2, \nUnpermuted", 
                              "Tax4Fun2, \nPermuted", "Tax4Fun2, \nUnpermuted")) + 
  theme_light(base_size = 9) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
plot_ra_cluster_nolegend



# generate plot of transformed p value (comparing preterm birth cases to term birth controls)
# correlations from cluster stratified analysis

plot_pcorr_cluster_legend <- ggplot() + 
  geom_boxplot(data = plot_pcorr_cluster[plot_pcorr_cluster$perm == "Permuted", ], 
               aes(x = as.factor(method), y = corr, color = cluster),
               outlier.shape = NA, fill = "white", lwd = 0.35, show.legend = F) + 
  geom_jitter(data = plot_pcorr_cluster[plot_pcorr_cluster$perm == "Permuted", ], 
              aes(x = as.factor(method), y = corr, color = cluster, alpha = 0.01, shape = perm), 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15, jitter.height = 0, seed = 1), 
              cex = 1, show.legend = F) + 
  geom_point(data = plot_pcorr_cluster[plot_pcorr_cluster$perm == "Unpermuted", ], 
             aes(x = as.factor(method), y = corr, color = as.factor(cluster), shape = perm), 
             cex = 3, position = position_dodge(width = 0.75)) + 
  xlab("") + 
  ylab("Inference Spearman rho") + 
  scale_color_manual(values = c("black", "#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c("All", 
                              expression(paste(italic("L. crispatus"), " dominated")), 
                              expression(paste(italic("L. iners"), " dominated")), 
                              "Mixed"), 
                     "Cluster") + 
  scale_shape_manual(values = c(20, 17), "Permutation") + 
  theme_light(base_size = 9) + 
  theme(legend.key.size = unit(0.7, "char"), 
        legend.key.height = unit(0.7, "char"), 
        legend.text.align = 0) + 
  guides(color = guide_legend(order = 1, override.aes = list(size = 2)), 
         shape = guide_legend(override.aes = list(fill = "black", size = c(3, 4)), order = 2, reverse = T))
plot_pcorr_cluster_legend



### figure 4 patchwork of relative abundance correlation plot and p correlation plot 

(plot_ra_cluster_nolegend | plot_pcorr_cluster_legend) +
  plot_annotation(tag_level = "A") +
  plot_layout(design = c(area(1, 1, 20, 10),
                         area(1, 11, 20, 17.5))) # TLBR






##########################################################################
### GENERATE FIGURE 5 WEIGHTED NSTI, PERCENT READS DISCARDED HISTOGRAM ###
##########################################################################



### prepare dataframe for plotting

# merge wnsti, discard, and cluster_f from data
hist <- left_join(wnsti, discard)
# hist should have 72 obs of 4 vars

# convert hist to long for two-panel plot with panels for wNSTI and % reads discarded
hist <- reshape2::melt(hist,
                       id.vars = c("sample",
                                   "cluster_f"),
                       measure.vars = c("wnsti",
                                        "reads"))
hist$variable <- ifelse(hist$variable == "wnsti", "PICRUSt2 weighted NSTI", 
                        "Tax4Fun2 % reads discarded")
# plot should have 144 obs of 4 vars



# generate plot
histogram <- ggplot(hist, aes(x = value, color = cluster_f, fill = cluster_f)) +
  geom_histogram(bins = 50, alpha = 0.65) +
  scale_color_manual(values = c("#ca0020","#ffd92f","#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus")," dominated")),
                                expression(paste(italic("L. iners")," dominated")),
                                "Mixed"),
                     "Cluster") +
  scale_fill_manual(values = c("#ca0020","#ffd92f","#0571b0"), 
                    labels = c(expression(paste(italic("L. crispatus")," dominated")),
                               expression(paste(italic("L. iners")," dominated")),
                               "Mixed"),
                    "Cluster") +
  facet_wrap(~variable, scales = "free_y") +
  ylab("Count") +
  xlab("Value") +
  theme_light(base_size = 10) +
  theme(legend.text.align = 0,
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"))
histogram






########################################################################
### GENERATE FIGURE 6 KO FUNCTIONAL CATEGORY STRATIFIED RESULTS PLOT ###
########################################################################



### generate plot of relative abundance correlations from KO functional category stratified analysis

plot_ra_ko_nolegend <- ggplot(data = plot_ra_ko, 
                              aes(x = as.factor(method_perm), y = corr, color = cat)) + 
  geom_boxplot(outlier.shape = NA, fill = "white", lwd = 0.35,
               position = position_dodge(width = 0.85)) + 
  stat_boxplot(geom = "errorbar", width = 0.5, lwd = 0.35, position = position_dodge(width = 0.85)) + 
  geom_jitter(data = plot_ra_ko, 
              aes(x = as.factor(method_perm), y = corr, color = cat, alpha = 0.0001), 
              shape = 20, cex = 1,
              position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.05, jitter.height = 0, seed = 1), 
              show.legend = F) + 
  xlab("") + 
  ylab("Spearman rho") + 
  scale_color_manual(values = c("black", 
                                "#CC0000", 
                                "#FF7733", 
                                "#FFD500", 
                                "#41ab5d",
                                "#00FFFF",
                                "#3377FF", 
                                "#FF66B3", 
                                "#8c510a"), 
                     "KO functional category", 
                     labels = c("All", 
                                "Metabolism", 
                                "Genetic information\nprocessing", 
                                "Environmental information\nprocessing", 
                                "Cellular processes", 
                                "Organismal systems", 
                                "Human diseases", 
                                "BRITE hierarchies", 
                                "Not included in\npathway or BRITE")) + 
  scale_x_discrete(breaks = c("PICRUSt2,  Permuted", "PICRUSt2,  Unpermuted", 
                              "Tax4Fun2,  Permuted", "Tax4Fun2,  Unpermuted"), 
                   labels = c("PICRUSt2, \nPermuted", "PICRUSt2, \nUnpermuted", 
                              "Tax4Fun2, \nPermuted", "Tax4Fun2, \nUnpermuted")) + 
  guides(color = "none") + 
  scale_y_continuous(breaks = c(-0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00), 
                     limits = c(-0.75, 1.00)) + 
  theme_light(base_size = 9) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1), 
        legend.key.size = unit(.8, "char"), legend.key.height = unit(0.8, "char"))
# plot_ra_ko_nolegend



# generate plot of transformed p value (comparing preterm birth cases to term birth controls)
# correlations from KO functional category stratified analysis

plot_pcorr_ko_legend <- ggplot() + 
  geom_boxplot(data = plot_pcorr_ko[plot_pcorr_ko$perm == "Permuted", ], 
               aes(x = as.factor(method), y = corr, color = cat), 
               outlier.shape = NA, fill = "white", lwd = 0.35, show.legend = F,
               position = position_dodge(width = 0.75)) + 
  geom_jitter(data = plot_pcorr_ko[plot_pcorr_ko$perm == "Permuted", ], 
              aes(x = as.factor(method), y = corr, color = cat, alpha = 0.01, shape = perm), 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1, jitter.height = 0, seed = 1), 
              show.legend = F, cex = 1) + 
  geom_point(data = plot_pcorr_ko[plot_pcorr_ko$perm == "Unpermuted", ], 
             aes(x = as.factor(method), y = corr, color = cat, shape = perm), cex = 3, 
             position = position_dodge(width = 0.75)) + 
  xlab("") + 
  ylab("Inference Spearman rho") + 
  scale_color_manual(values = c("black",
                                "#CC0000",
                                "#FF7733", 
                                "#FFD500",
                                "#41ab5d", 
                                "#00FFFF",
                                "#3377FF",
                                "#FF66B3", 
                                "#8c510a"), 
                     "KO functional category", 
                     labels = c("All",
                                "Metabolism", 
                                "Genetic information\nprocessing", 
                                "Environmental information\nprocessing", 
                                "Cellular processes",
                                "Organismal systems", 
                                "Human diseases", 
                                "BRITE hierarchies", 
                                "Not included in\npathway or BRITE")) + 
  scale_y_continuous(breaks = c(-0.5, -0.25, 0.0, 0.25, 0.50), 
                     limits = c(-0.5, 0.5)) + 
  scale_shape_manual(values = c(20, 17), "Permutation") + 
  theme_light(base_size = 9) + 
  theme(legend.key.size = unit(0.7, "char"), legend.key.height = unit(0.7, "char")) + 
  guides(color = guide_legend(order = 1, override.aes = list(size = 2)), 
         shape = guide_legend(override.aes = list(fill = "black", size = c(3, 4)), order = 2, reverse = T))
# plot_pcorr_ko_legend




### figure 6 patchwork of relative abundance correlation plot and p correlation plot 
(plot_ra_ko_nolegend / plot_pcorr_ko_legend) +
  plot_annotation(tag_level = "A") +
  plot_layout(design = c(area(1, 1, 10.5, 10),
                         area(11, 1, 20, 10)), # TLBR
              guides = "collect")






##########################################################
### GENERATE SUPPLEMENTAL FIGURE 4 P VS P SCATTER PLOT ###
##########################################################



### generate scatter plot of transformed p values (comparing preterm birth cases to term 
### birth controls) from cluster stratified analysis

plot_p_cluster <- ggplot(data = plot_pvp_cluster, 
                         aes(x = p_obs, y = p_pred, color = cluster)) +  
  geom_point(alpha = 0.17) +  
  geom_smooth(se = FALSE, method = lm, lwd = 1.5) +  
  facet_wrap(plot_pvp_cluster$method) +  
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     "Cluster", 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")), 
                                expression(paste(italic("L. iners"), " dominated")), 
                                "Mixed")) +  
  xlab("Transformed p value from observed metagenome") +  
  ylab("Transformed p value from predicted metagenome") +  
  theme_light(base_size = 9) +  
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        legend.text.align = 0)
# plot_p_cluster



### generate scatter plot of transformed p values (comparing preterm birth cases to term 
### birth controls) from KO category stratified analysis

### prepare dataframe

# when use factor order used in other plots for the KO functional category stratified
# analysis, BRITE hierarchies is plotted second to last and so on top most of the data
# and the plot looks overwhelmingly pink because BRITE hierarchies category is much
# larger than others
# to fix, will level the factor in reverse order, and then plot legend in reverse of the
# reversed order
plot_pvp_ko$cat <- factor(plot_pvp_ko$cat,
                          levels = c("Not included in pathway or BRITE",
                                     "BRITE hierarchies",
                                     "Human diseases",
                                     "Organismal systems",
                                     "Cellular processes",
                                     "Environmental information processing",
                                     "Genetic information processing",
                                     "Metabolism"))

# make long dataframe for plottingm, need to double length of p_obs, cat, gene
# need to stack p_pi and p_tax and add variable indicating pi or tax
p_obs <- rep(plot_pvp_ko$p_obs, 2)
gene <- rep(plot_pvp_ko$gene, 2)
cat <- rep(plot_pvp_ko$cat, 2)
p_pred <- c(plot_pvp_ko$p_pi, plot_pvp_ko$p_tax)
method <- c(rep("PICRUSt2", nrow(plot_pvp_ko)),
            rep("Tax4Fun2", nrow(plot_pvp_ko)))
plot_pvp_ko <- data.frame(p_obs)
plot_pvp_ko$gene <- gene
plot_pvp_ko$cat <- cat
plot_pvp_ko$p_pred <- p_pred
plot_pvp_ko$method <- factor(method)



### generate scatter plot of transformed p values (comparing preterm birth cases to term 
### birth controls) from KO functional category stratified analysis
plot_p_ko <- ggplot(data = plot_pvp_ko, 
                    aes(x = p_obs, y = p_pred, color = cat, alpha = cat)) + 
  geom_point() + 
  scale_alpha_manual(values = c(0.5, 0.35, rep(0.5, 6)), 
                     guide = "none") + 
  geom_smooth(se = FALSE, method = lm, lwd = 1.5) + 
  facet_wrap(plot_pvp_ko$method) + 
  scale_color_manual(values = c("#8c510a", 
                                "#FF66B3",
                                "#3377FF",
                                "#00FFFF", 
                                "#41ab5d", 
                                "#FFD500",
                                "#FF7733",
                                "#CC0000"), 
                     "KO functional category", 
                     labels = c("Not included in\npathway or BRITE", 
                                "BRITE hierarchies", 
                                "Human diseases", 
                                "Organismal systems", 
                                "Cellular processes",
                                "Environmental information\nprocessing", 
                                "Genetic information\nprocessing", 
                                "Metabolism")) + 
  xlab("Transformed p value from observed metagenome") + 
  ylab("Transformed p value from predicted metagenome") + 
  theme_light(base_size = 9) + 
  guides(color = guide_legend(override.aes = list(shape = 16), reverse = T), 
         axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)), 
         axis.title.x = element_text(margin = margin(t = 4, r = 0, b = 0, l = 0)))
# plot_p_ko



### patchworl final plot
(plot_p_cluster / plot_p_ko) +
  plot_annotation(tag_level = "A") +
  plot_layout(design = c(area(1, 1, 10.5, 10),
                         area(11, 1, 20, 10))) # TLBR






##################################################
### GENERATE FIGURE 7 CORRELATION SCATTER PLOT ###
##################################################



### generate dataframe with sample, cluster, L. crispatus relative abundance, L. iners relative abundance, 
### ratio of L. crispatus relative abundance:L. iners relative abundance, and
### level-1 KO functional category relative abundances

plot_scatter <- data[, colnames(data) %in% c("sample",
                                             "cluster",
                                             "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.s__Lactobacillus_crispatus_ra",
                                             "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.s__Lactobacillus_iners_ra",
                                             "met_ra",
                                             "gen_ra",
                                             "env_ra",
                                             "cel_ra",
                                             "org_ra",
                                             "hum_ra",
                                             "bri_ra",
                                             "not_ra")]
# plot_scatter should have 72 obs of 12 vars

# make variable for ratio of L. crispatus:L. iners relative abundance
plot_scatter <- rename(plot_scatter, 
                       crisp = d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.s__Lactobacillus_crispatus_ra)
plot_scatter <- rename(plot_scatter, 
                       iners = d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.s__Lactobacillus_iners_ra)
plot_scatter$crisp_iners <- plot_scatter$crisp / plot_scatter$iners
# plot_scatter should have 72 obs of 13 vars



### generate long dataframe for plotting

plot_scatter <- reshape2::melt(plot_scatter,
                       id.vars = c("sample", 
                                   "cluster",
                                   "crisp",
                                   "iners",
                                   "crisp_iners"),
                       measure.vars = c("met_ra",
                                        "gen_ra",
                                        "env_ra",
                                        "cel_ra",
                                        "org_ra",
                                        "hum_ra",
                                        "bri_ra",
                                        "not_ra"))
# plot_scatter should have 576 obs of 7 vars

# generate a factor KO functional category variable, factor cluster
plot_scatter <- plot_scatter %>% mutate(cat = case_when(variable == "met_ra" ~ "Metabolism",
                                                        variable == "gen_ra" ~ "Genetic information\nprocessing",
                                                        variable == "env_ra" ~ "Environmental information\nprocessing",
                                                        variable == "cel_ra" ~ "Cellular processes",
                                                        variable == "org_ra" ~ "Organismal systems",
                                                        variable == "hum_ra" ~ "Human diseases",
                                                        variable == "bri_ra" ~ "BRITE hierarchies",
                                                        variable == "not_ra" ~ "Not included in\npathway or BRITE"))
# plot_scatter should have 576 obs of 8 vars

plot_scatter$cat <- factor(plot_scatter$cat,levels = c("Metabolism",
                                                       "Genetic information\nprocessing",
                                                       "Environmental information\nprocessing",
                                                       "Cellular processes",
                                                       "Organismal systems",
                                                       "Human diseases",
                                                       "BRITE hierarchies",
                                                       "Not included in\npathway or BRITE"))

plot_scatter$cluster <- factor(plot_scatter$cluster,
                               levels = c(3, 2, 1))



### estimate r squared from correlations between:
### L. crispatus relative abundance and KO functional category relative abundance
### L. iners relative abundance and KO functional category relative abundance
### ratio of L. crispatus relative abundance:L. iners relative abundance and KO functional 
### category relative abundance
### put r squared in dataframe along with coordinates to plot the values, 
### merge these into plot_scatter

# for L. crispatus realtive abundances
r2_crisp <- c()
variable <- c()
for(each in unique(plot_scatter$variable)){
  model <- lm(value ~ crisp, data = plot_scatter[plot_scatter$variable == each, ])
  r2_crisp <- append(r2_crisp, scales::percent(summary(model)$r.squared, accuracy = 0.1))
  variable <- append(variable, each)
}
r2_crisp <- as.data.frame(cbind(r2_crisp, variable))
r2_crisp$x_crisp <- c(0.5)
r2_crisp$y_crisp <- c(0.18, 0.20, 0.04, 0.05, 0.030, 0.06, 0.41, 0.05)
plot_scatter <- left_join(plot_scatter, r2_crisp)
# plot_scatter should have 576 obs of 11 vars

# for L. iners relative abundances
r2_iners <- c()
variable <- c()
for(each in unique(plot_scatter$variable)){
  model <- lm(value ~ iners, data = plot_scatter[plot_scatter$variable == each, ])
  r2_iners <- append(r2_iners, scales::percent(summary(model)$r.squared, accuracy = 0.1))
  variable <- append(variable, each)
}
r2_iners <- as.data.frame(cbind(r2_iners, variable))
r2_iners$x_iners <- c(0.875, 0.125, 0.875, 0.625, 0.500, 0.50, 0.500, 0.875)
r2_iners$y_iners <- c(0.225, 0.225, 0.030, 0.050, 0.035, 0.06, 0.405, 0.05)
plot_scatter <- left_join(plot_scatter, r2_iners)
# plot_scatter should have 576 obs of 14 vars

# for ratio of L. crispatus relative abundances:L. iners relative abundances
r2_ratio <- c()
variable <- c()
for(each in unique(plot_scatter$variable)){
  model <- lm(value ~ crisp_iners, data = plot_scatter[plot_scatter$variable == each, ])
  r2_ratio <- append(r2_ratio, scales::percent(summary(model)$r.squared, accuracy = 0.1))
  variable <- append(variable, each)
}
r2_ratio <- as.data.frame(cbind(r2_ratio, variable))
r2_ratio$x_ratio <- c(10)
r2_ratio$y_ratio <- c(0.165, 0.225, 0.03, 0.05, 0.0325, 0.06, 0.41, 0.03)
plot_scatter <- left_join(plot_scatter, r2_ratio)
# plot_scatter should have 576 obs of 17 vars



### generate scatter plots of L. crispatus relative abundance v KO functional category 
### relative abundance, faceted by KO functional category, with legend for patchwork

plot_crisp_scat_legend <- ggplot(data = plot_scatter, aes(x = crisp, y = value, color = factor(cluster))) + 
  geom_point() + 
  facet_wrap(plot_scatter$cat, scales = "free", nrow = 2) + 
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")), 
                              expression(paste(italic("L. iners"), " dominated")), 
                              "Mixed"), 
                     "Cluster") + 
  xlab(expression(paste(italic("L. crispatus"), " relative abundance"))) + 
  ylab("KO functional category\nrelative abundance") + 
  geom_text(mapping = aes(x = x_crisp, y = y_crisp, label = r2_crisp), color = "black", size = 2.5, fontface = "italic", alpha = 0.03) + 
  theme_light(base_size = 8) + 
  theme(legend.text.align = 0, legend.text = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.direction = "horizontal", 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0))) + 
  guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5))
# plot_crisp_scat_legend



### generate scatter plots of L. iners relative abundance v KO functional category 
### relative abundance, faceted by KO functional category

plot_iners_scat <- ggplot(data = plot_scatter, aes(x = iners, y = value, color = factor(cluster))) + 
  geom_point() + 
  facet_wrap(plot_scatter$cat, scales = "free", nrow = 2) + 
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")), 
                              expression(paste(italic("L. iners"), " dominated")), 
                              "Mixed"), 
                     "Cluster") + 
  xlab(expression(paste(italic("L. iners"), " relative abundance"))) + 
  ylab("KO functional category\nrelative abundance") + 
  labs(tag  = "B") + 
  geom_text(mapping = aes(x = x_iners, y = y_iners, label = r2_iners), color = "black", size = 2.5, fontface = "italic", alpha = 0.03) + 
  guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5)) + 
  theme_light(base_size = 8) + 
  theme(legend.position = "none", 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)), 
        plot.tag = element_text(size = 12))
# plot_iners_scat



### generate scatter plots of ratio of L. crispatus relative abundance:L. iners relative abundance
### v KO functional category relative abundance, faceted by KO functional category

plot_ratio_scat <- ggplot(data = plot_scatter, aes(x = crisp_iners, y = value, color = factor(cluster))) + 
  geom_point() + 
  facet_wrap(plot_scatter$cat, scales = "free", nrow = 2) + 
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")), 
                              expression(paste(italic("L. iners"), " dominated")), 
                              "Mixed"), 
                     "Cluster") + 
  xlab(expression(paste(italic("L. crispatus"), " relative abundance : ", italic("L. iners"), " relative abundance"))) + 
  ylab("KO functional category\nrelative abundance") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x)10^x), labels = trans_format("log10", math_format(10^.x))) + 
  labs(tag  = "C") + 
  geom_text(mapping = aes(x = x_ratio, y = y_ratio, label = r2_ratio), color = "black", size = 2.5, fontface = "italic", alpha = 0.03) + 
  guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5)) + 
  theme_light(base_size = 8) + 
  theme(legend.position = "none", 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)), 
        plot.tag = element_text(size = 12))
# plot_ratio_scat



### generate figure 7, patchwork crisp, iners, ratio scatter plots

plot_crisp_scat_legend / plot_iners_scat / plot_ratio_scat +
  guide_area() +
  plot_annotation(tag_level = "A") +
  plot_layout(guides = "collect",
              design = c(area(1, 1, 10, 10), 
                         area(11, 1, 20, 10), 
                         area(21, 1, 30, 10),
                         area(31, 1, 32, 10))) + # TLBR
  theme(plot.tag = element_text(size = 11), legend.position = "bottom")
  
  




##################################################################################
### GENERATE SUPPLEMENTAL FIGURE 4 STACKED BAR PLOT ORDERED BY WMGS DENDROGRAM ###
##################################################################################



### will use plotdata_16s and plotdata_meta from figure 1 

### will order samples by label variable in leaf_meta

# check to make sure sample variable in plotdata_16S, plotdata_meta, and leaf_meta$label are formatted the same
unique(plotdata_16s$sample)
unique(plotdata_meta$sample)
leaf_meta$label
# leaf_meta$label has leading "X" in sample names, need to drop
leaf_meta$label <- substr(leaf_meta$label, 2, 6)



# generate 16S composition stacked bar plot ordered by WMGS dendrogram
# will pull legend for later patchworked final plot
(stacked_bar_16s_meta <- ggplot(data = plotdata_16s, 
                                aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
   geom_bar(stat = "identity",
            aes(fill = factor(order, levels = c("16",
                                                "15",
                                                "14",
                                                "13",
                                                "12",
                                                "11",
                                                "10",
                                                "9",
                                                "8",
                                                "7",
                                                "6",
                                                "5",
                                                "4",
                                                "3",
                                                "2",
                                                "1")))) +
   scale_fill_manual(values = c("gray48",
                                "#006837",
                                "#31a354",
                                "#c2e699",
                                "#ffffcc",
                                "#fbb4b9",
                                "#f768a1",
                                "#c51b8a",
                                "#7a0177",
                                
                                "#c6dbef",
                                "#9ecae1",
                                "#6baed6",
                                "#2171b5",
                                "#08519c",
                                "#08306b",
                                "#031737"),
                     labels = c("Other",
                                expression(paste(italic("Ureaplasma"), "")),
                                expression(paste(italic("Streptococcus agalactiae"), "")),
                                expression(paste(italic("Sneathia sanguinegens"), "")),
                                expression(paste(italic("Megasphaera"), " unidentified")),
                                expression(paste(italic("Gardnerella"), " uncultured")),
                                expression(paste(italic("Gardnerella"), "")),
                                expression(paste(italic("Fannyhessea vaginae"), "")),
                                expression(paste(italic("Bifidobacteriaceae"), "")),
                                expression(paste(italic("Lactobacillus"), " metagenome")),
                                expression(paste(italic("Lactobacillus"), "")),
                                expression(paste(italic("Lactobacillus paragasseri"), "")),
                                expression(paste(italic("Lactobacillus jensenii"), "")),
                                expression(paste(italic("Lactobacillus iners"), "")),
                                expression(paste(italic("Lactobacillus cripatus")," ")),
                                expression(paste(italic("Lactobacillus antri"), ""))),
                     "Bacterial taxa") +
   guides(fill = guide_legend(reverse = T,
                              override.aes = list(size = 4))) +
   scale_x_discrete(limits = rev, labels = NULL) +
   labs(x="", y = "Taxon relative abundance") +
   geom_rect(aes(xmin = 0.5, xmax = 72.5, ymin = 1.0275, ymax = 1.03), 
             fill = "black") + # long vertical bar in brackets
   geom_rect(aes(xmin = 0.5, xmax = 0.35, ymin = 1.01, ymax = 1.03), 
             fill = "black") + # short horizontal bar at bottom
   geom_rect(aes(xmin = 30.43, xmax = 30.58, ymin = 1.01, ymax = 1.03), 
             fill = "black") + # short horizontal bar at boundary of clusters
   geom_rect(aes(xmin = 72.35, xmax = 72.5, ymin = 1.01, ymax = 1.03), 
             fill = "black") + # short horizontal bar at top
   annotate(geom = "text", x = 51, y = 1.08, 
            label = "WMGS cluster 1\n(metabolism,  uncharacterized KO enriched)", 
            color = "black", size = 2, angle = 90) +
   annotate(geom = "text", x = 15, y = 1.08, 
            label = "WMGS cluster 2\n(genetic information processing KO enriched)", 
            color = "black", size = 2, angle = 90) +
   coord_flip() +
   theme_minimal(base_size = 9) +
   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
         legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text.align = 0,
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# pull legend for patchworking
legend_16_meta <- as_ggplot(get_legend(stacked_bar_16s_meta))



# generate 16S composition stacked bar plot ordered by WMGS dendrogram 
# with tiles for case status, no legend
# will use plot for later patchworked final plot
(stacked_bar_16s_case_nolegend_meta <- ggplot(data = plotdata_16s, 
                                              aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
    geom_bar(stat = "identity",
             aes(fill = factor(order, levels = c("16",
                                                 "15",
                                                 "14",
                                                 "13",
                                                 "12",
                                                 "11",
                                                 "10",
                                                 "9",
                                                 "8",
                                                 "7",
                                                 "6",
                                                 "5",
                                                 "4",
                                                 "3",
                                                 "2",
                                                 "1")))) +
    scale_fill_manual(values = c("gray48",
                                 "#006837",
                                 "#31a354",
                                 "#c2e699",
                                 "#ffffcc",
                                 "#fbb4b9",
                                 "#f768a1",
                                 "#c51b8a",
                                 "#7a0177",
                                 
                                 "#c6dbef",
                                 "#9ecae1",
                                 "#6baed6",
                                 "#2171b5",
                                 "#08519c",
                                 "#08306b",
                                 "#031737"),
                      labels = c("Other",
                                 expression(paste(italic("Ureaplasma"), "")),
                                 expression(paste(italic("Streptococcus agalactiae"), "")),
                                 expression(paste(italic("Sneathia sanguinegens"), "")),
                                 expression(paste(italic("Megasphaera"), " unidentified")),
                                 expression(paste(italic("Gardnerella"), " uncultured")),
                                 expression(paste(italic("Gardnerella"), "")),
                                 expression(paste(italic("Fannyhessea vaginae"), "")),
                                 expression(paste(italic("Bifidobacteriaceae"), "")),
                                 expression(paste(italic("Lactobacillus"), " metagenome")),
                                 expression(paste(italic("Lactobacillus"), "")),
                                 expression(paste(italic("Lactobacillus paragasseri"), "")),
                                 expression(paste(italic("Lactobacillus jensenii"), "")),
                                 expression(paste(italic("Lactobacillus iners"), "")),
                                 expression(paste(italic("Lactobacillus cripatus")," ")),
                                 expression(paste(italic("Lactobacillus antri"), ""))),
                      "Bacterial taxa") +
    guides(fill = guide_legend(reverse = T,
                               override.aes = list(size = 4),
                               order = 3)) +
    new_scale("fill") +
    geom_tile(data = plotdata_16s, aes(x = sample, y = case_y,
                                       fill = factor(C_CASE2, levels = c(1, 0)),
                                       height = 0.05, width = 0.875)) +
    scale_fill_manual(values = c("red3", "black"),
                      labels = c("Preterm case", "Term control"),
                      "Birth outcome") +
    guides(fill = guide_legend(reverse = F, 
                               override.aes = list(size = 4),
                               order = 2)) +
    scale_x_discrete(limits = rev, labels = NULL) +
    labs(x = "", y = "Taxon relative abundance") +
    geom_rect(aes(xmin = 0.5, xmax = 72.5, ymin = 1.0275, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 0.5, xmax = 0.35, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 30.43, xmax = 30.58, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 72.35, xmax = 72.5, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    annotate(geom = "text", x = 51, y = 1.08, 
             label = "WMGS cluster 1\n(metabolism,  uncharacterized KO enriched)", 
             color = "black", size = 2, angle = 90) +
    annotate(geom = "text", x = 15, y = 1.08, 
             label = "WMGS cluster 2\n(genetic information processing KO enriched)", 
             color = "black", size = 2, angle = 90) +coord_flip() +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



# generate plot of tiles for case status ordered by WMGS dendrogram
# will pull legend for later patchworked final plot
(tiles_meta <- ggplot(data = plotdata_16s, 
                      aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
    geom_tile(data = plotdata_16s, aes(x = sample, y = case_y,
                                       fill = factor(C_CASE2, levels = c(1, 0)),
                                       height = 0.05, width = 0.875)) +
    scale_fill_manual(values = c("red3", "black"),
                      labels = c("Preterm case", "Term control"),
                      "Birth outcome") +
    guides(fill = guide_legend(reverse = F, 
                               override.aes = list(size = 4),
                               order = 2)) +
    scale_x_discrete(limits = rev, labels = NULL) +
    labs(x = "", y = "") +
    coord_flip() +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text.align = 0,
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# pull legend for patchworking
legend_tile_meta <- as_ggplot(get_legend(tiles_meta))



# generate plot of WMGS dendrogram 
# to use in final patchworked plot
(dendrogram_meta <- ggplot(segment(dendro_data(dendro_meta, type = "rectangle"))) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lwd = .35) +
    coord_flip() + 
    scale_y_reverse() +
    theme_void() +
    scale_x_reverse(expand = expansion(mult = 0.01)))



# generate plot spacer 
# to use in final patchworked plot
plot_spacer <- plot_spacer()



# generate observed metagenome KO level 1 stacked bar plot ordered by WMGS dendrogram
# will pull legend for final patchworked plot
(stacked_bar_meta_meta <- ggplot(data = plotdata_meta, 
                                 aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
    geom_bar(stat = "identity", 
             aes(fill = variable)) +
    scale_fill_manual("Functional category",
                      breaks=c("met_ra", 
                               "gen_ra", 
                               "env_ra", 
                               "cel_ra", 
                               "org_ra", 
                               "hum_ra", 
                               "bri_ra", 
                               "not_ra"),
                      values=c("#CC0000",
                               "#FF7733",
                               "#FFD500",
                               "#41ab5d",
                               "#00FFFF",
                               "#3377FF",
                               "#FF66B3",
                               "#8c510a"),
                      labels=c("Metabolism",
                               "Genetic information processing",
                               "Environmental information processing",
                               "Cellular processes",
                               "Organismal systems",
                               "Human diseases",
                               "BRITE hierarchies",
                               "Not included in pathway or BRITE")) +
    guides(fill = guide_legend(reverse = F, override.aes = list(size = 4))) +
    scale_x_discrete(limits = rev, labels = NULL) +
    labs(x = "", y = "KO functional category relative abundance") +
    geom_rect(aes(xmin = 0.5, xmax = 72.5, ymin = 1.0275, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 0.5, xmax = 0.35, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 30.43, xmax = 30.58, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 72.35, xmax = 72.5, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    annotate(geom = "text", x = 51, y = 1.08, 
             label = "WMGS cluster 1\n(metabolism,  uncharacterized KO enriched)", 
             color = "black", size = 2, angle = 90) +
    annotate(geom = "text", x = 15, y = 1.08, 
             label = "WMGS cluster 2\n(genetic information processing KO enriched)", 
             color = "black", size = 2, angle = 90) +coord_flip() +
    coord_flip() +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# pull legend for patchworking
legend_meta_meta <- as_ggplot(get_legend(stacked_bar_meta_meta))



# generate observed metagenome KO level 1 stacked bar plot ordered by WMGS dendrogram
# without legend
# will use for final patchworked plot
(stacked_bar_meta_nolegend_meta <- ggplot(data = plotdata_meta, 
                                          aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
    geom_bar(stat = "identity", 
             aes(fill = variable)) +
    scale_fill_manual("Functional category",
                      breaks=c("met_ra", 
                               "gen_ra", 
                               "env_ra", 
                               "cel_ra", 
                               "org_ra", 
                               "hum_ra", 
                               "bri_ra", 
                               "not_ra"),
                      values=c("#CC0000",
                               "#FF7733",
                               "#FFD500",
                               "#41ab5d",
                               "#00FFFF",
                               "#3377FF",
                               "#FF66B3",
                               "#8c510a"),
                      labels=c("Metabolism",
                               "Genetic information processing",
                               "Environmental information processing",
                               "Cellular processes",
                               "Organismal systems",
                               "Human diseases",
                               "BRITE hierarchies",
                               "Not included in pathway or BRITE")) +
    guides(fill = guide_legend(reverse = F, override.aes = list(size = 4))) +
    scale_x_discrete(limits = rev, labels = NULL) +
    labs(x = "", y = "KO functional category relative abundance") +
    geom_rect(aes(xmin = 0.5, xmax = 72.5, ymin = 1.0275, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 0.5, xmax = 0.35, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 30.43, xmax = 30.58, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    geom_rect(aes(xmin = 72.35, xmax = 72.5, ymin = 1.01, ymax = 1.03), 
              fill = "black") + 
    annotate(geom = "text", x = 51, y = 1.08, 
             label = "WMGS cluster 1\n(metabolism,  uncharacterized KO enriched)", 
             color = "black", size = 2, angle = 90) +
    annotate(geom = "text", x = 15, y = 1.08, 
             label = "WMGS cluster 2\n(genetic information processing KO enriched)", 
             color = "black", size = 2, angle = 90) +coord_flip() +
    coord_flip() +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



# generate supplemental figure 5 patchworked plot
(stacked_bars_patch_meta <- (dendrogram_meta | 
                        stacked_bar_16s_case_nolegend_meta | 
                        stacked_bar_meta_nolegend_meta | 
                        legend_tile_meta | 
                        legend_meta_meta | 
                        legend_16_meta | 
                        plot_spacer)) +
  plot_layout(design = c(area(1, 1, 20, 10), #dendrogram area TLBR
                         area(1, 9, 20, 58), #16S stacked bar and case status tile area
                         area(1, 58, 20, 107), #metagenome stacked bar area
                         area(0, 114, 2), #case status legend area
                         area(4.5, 123, 8), #metagenome legend area
                         area(11, 119, 20), #16S legend area
                         area(1, 135, 20, 138))) #spacer area






#################################################################################################################
### GENERATE SUPPLEMENTAL FIGURE 5 OBSERVED METAGENOME KO LEVEL 2 STACKED BAR PLOT ORDERED BY WMGS DENDROGRAM ###
#################################################################################################################



### using plotdata_meta2 from supplemental figure 2

# plot ordered by WMGS hierarchical clustering dendrogram with tiles for birth
(stacked_bar_meta2_meta <- ggplot(data = plotdata_meta2, 
                                  aes(x = factor(sample, levels = leaf_meta$label), y = value)) +
   geom_bar(stat = "identity",
            aes(fill = variable)) +
   scale_fill_manual("Functional category",
                     values=c("#521414", 
                              "#660000", 
                              "#8A0F0F", 
                              "#CC0000", 
                              "#FF3333", 
                              "#EB4747", 
                              "#FF6666", 
                              "#E87D7D", 
                              "#FF9999", 
                              "#FFCCCC", 
                              "#F7D4D4",
                              
                              "#CC4400",
                              "#FF5500",
                              "#FF7733",
                              "#FF9966",
                              "#FFBB99",
                              
                              "#FFD500",
                              "#FFEE99",
                              
                              "#005a32",
                              "#238b45",
                              "#41ab5d",
                              "#74c476",
                              "#a1d99b",
                              "#c7e9c0",
                              
                              "#1F4747",
                              "#006666",
                              "#178282",
                              "#1FADAD",
                              "#00CCCC",
                              "#00FFFF",
                              "#66FFFF",
                              "#99FFFF",
                              "#CCFFFF",
                              
                              "#002266",
                              "#003399",
                              "#0055FF",
                              "#3377FF",
                              "#6699FF",
                              "#99BBFF",
                              "#CCDDFF",
                              
                              "#330066",
                              "#6600CC",
                              "#9933FF",
                              "#CC99FF",
                              "#E5CCFF",
                              
                              "#FF3399",
                              "#FF66B3",
                              "#FF99CC",
                              "#FFCCE6",
                              
                              "#543005",
                              "#8c510a",
                              "#bf812d",
                              "#dfc27d"),
                     labels=c("Carbohydrate metabolism",
                              "Energy metabolism",
                              "Lipid metabolism",
                              "Nucleotide metabolism",
                              "Amino acid metabolism",
                              "Metabolism of other amino acids",
                              "Glycan biosynthesis and metabolism",
                              "Metabolism of cofactors and vitamins",
                              "Metabolism of terpenoids and polyketides",
                              "Biosynthesis of other secondary metabolites",
                              "Xenobiotics degradation and metabolism",
                              
                              "Transcription",
                              "Translation",
                              "Folding, sorting and degradation",
                              "Replication and repair",
                              "Information processing in viruses",
                              
                              "Membrane transport",
                              "Signal transduction",
                              
                              "Transport and catabolism",
                              "Cell motility",
                              "Cell growth and death",
                              "Cellular community - eukaryotes",
                              "Cellular community - prokaryotes",
                              "Aging",
                              
                              "Immune system",
                              "Endocrine system",
                              "Circulatory system",
                              "Digestive system",
                              "Excretory system",
                              "Nervous system",
                              "Sensory system",
                              "Development and regulation",
                              "Environmental adaptation",
                              
                              "Cancer: overview",
                              "Cancer: specific types",
                              "Immune disease",
                              "Neurodegenerative disease",
                              "Substance dependence",
                              "Cardiovascular disease",
                              "Endocrine and metabolic disease",
                              
                              "Infectious disease: bacterial",
                              "Infectious disease: viral",
                              "Infectious disease: parasitic",
                              "Drug resistance: antimicrobial",
                              "Drug resistance: antineoplastic",
                              
                              "Protein families: metabolism",
                              "Protein families: genetic information processing",
                              "Protein families: signaling and cellular processes",
                              "Viral protein families",
                              
                              "Unclassified: metabolism",
                              "Unclassified: genetic information processing",
                              "Unclassified: signaling and cellular processes",
                              "Poorly characterized")) +
   guides(fill = guide_legend(reverse = F,override.aes = list(size = 1), 
                              order = 3, ncol = 2)) +
   new_scale("fill") +
   geom_tile(data = plotdata_meta2, aes(x = sample, y = case_y,
                                        fill = factor(C_CASE2, levels = c(1,0)),
                                        height = 0.05, width = 0.875)) +
   scale_fill_manual(values = c("red3", "black"),
                     labels = c("Preterm case", "Term control"),
                     "Birth outcome") +
   guides(fill = guide_legend(reverse = F, override.aes = list(size = 1), order = 2)) +
   scale_x_discrete(limits = rev, labels = NULL) +
   labs(x = "", y = "KO functional category relative abundance") +
   geom_rect(aes(xmin = 0.5, xmax = 72.5, ymin = 1.0275, ymax = 1.03), 
             fill = "black") + 
   geom_rect(aes(xmin = 0.5, xmax = 0.35, ymin = 1.01, ymax = 1.03), 
             fill = "black") + 
   geom_rect(aes(xmin = 30.43, xmax = 30.58, ymin = 1.01, ymax = 1.03), 
             fill = "black") + 
   geom_rect(aes(xmin = 72.35, xmax = 72.5, ymin = 1.01, ymax = 1.03), 
             fill = "black") + 
   annotate(geom = "text", x = 51, y = 1.08, 
            label = "WMGS cluster 1\n(metabolism,  uncharacterized KO enriched)", 
            color = "black", size = 2, angle = 90) +
   annotate(geom = "text", x = 15, y = 1.08, 
            label = "WMGS cluster 2\n(genetic information processing KO enriched)", 
            color = "black", size = 2, angle = 90) +coord_flip() +
   coord_flip() +
   theme_minimal(base_size = 9) +
   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
         legend.key.size = unit(0.9, "char"), legend.key.height = unit(0.9, "char"), legend.text = element_text(size = 5),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



# patchwork final plot
(dendrogram_meta | plot_spacer | stacked_bar_meta2_meta) +
  plot_layout(widths = c(1.5, -1.3, 9.5))






#############################################################################
### GENERATE SUPPLEMENTAL FIGURE 6 OBSERVED METAGENOME KO ALPHA DIVERSITY ###
#############################################################################



### make long dataset for plotting

# add sample variable to alpha_meta
alpha_meta$sample <- as.numeric(rownames(alpha_meta))
# alpha_meta should have 72 obs of 5 vars

# merge C_CASE2 and cluster_meta from data with alpha_meta
plot_alpha_meta <- data[, colnames(data) %in% c("sample",  "cluster_f", "cluster_meta")]
plot_alpha_meta <- left_join(plot_alpha_meta, alpha_meta)
# plot_alpha_meta should have 72 obs of 7 vars

# reshape to long
plot_alpha_meta <- reshape2::melt(plot_alpha_meta, 
                     id.vars = c("sample",  "cluster_f", "cluster_meta"), 
                     measure.vars = c("Observed",  "Shannon",  "Simpson",  "InvSimpson"))
colnames(plot_alpha_meta) <- c("sample",  "cluster_f", "cluster_meta", "metric",  "alpha")
# plot_alpha_meta should have 288 obs of 5 vars



### generate plot colored by WMGS cluster
plot_alpha_meta_meta <- ggplot(data = plot_alpha_meta, 
                               aes(x = factor(cluster_meta), y = alpha, color = factor(cluster_meta)))  + 
  geom_violin(lwd = 0.75) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  geom_jitter(shape = 1) +
  scale_shape(solid = F) + 
  facet_wrap(factor(plot_alpha_meta$metric, levels = c("Observed", "Shannon", "Simpson", "InvSimpson")), 
             scales = "free") + 
  xlab("") + 
  ylab("WMGS Alpha diversity") + 
  theme_light(base_size = 11) +
  theme(legend.text.align = 0, 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
plot_alpha_meta_meta



### generate plot colored by WMGS cluster
plot_alpha_meta_16s <- ggplot(data = plot_alpha_meta, 
                               aes(x = cluster_f, y = alpha, color = cluster_f))  + 
  geom_violin(lwd = 0.75) + 
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  geom_jitter(shape = 1) +
  scale_shape(solid = F) + 
  facet_wrap(factor(plot_alpha_meta$metric, levels = c("Observed", "Shannon", "Simpson", "InvSimpson")), 
             scales = "free") + 
  xlab("") + 
  ylab("WMGS Alpha diversity") + 
  theme_light(base_size = 11) +
  theme(legend.text.align = 0, 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
plot_alpha_meta_16s



### patchwork final plot
plot_alpha_meta_meta / plot_alpha_meta_16s +
  plot_annotation(tag_level = "A") +
  theme(plot.tag = element_text(size = 11))






##################################################################
### GENERATE SUPPLEMENTAL FIGURE 7 WMGS VS 16S ALPHA DIVERSITY ###
##################################################################



### make long dataset for plotting

# add sample variable to alpha
alpha$sample <- as.numeric(rownames(alpha))
# alpha should have 72 obs of 5 vars

# add "_meta" to colnames for WMGS alpha diversity metrics
colnames(alpha_meta) <- c("Observed_meta", "Shannon_meta", "Simpson_meta", "InvSimpson_meta", "sample")

# merge alpha, alpha_meta, and cluster variables
plot_alpha_vs <- left_join(alpha, alpha_meta)
plot_alpha_vs <- left_join(plot_alpha_vs, data[, colnames(data) %in% c("sample", "cluster_f", "cluster_meta_f")])
# plot_alpha_vs should have 72 obs of 11 vars

# lengthen over 16S alpha diversity
plot_alpha_vs <- reshape2::melt(plot_alpha_vs,
                                id.vars = c("sample", "cluster_f", "cluster_meta_f", 
                                            "Observed_meta", "Shannon_meta",
                                            "Simpson_meta", "InvSimpson_meta"),
                                measure.vars = colnames(plot_alpha_vs)[colnames(plot_alpha_vs) %!in% c("sample", "cluster_f", "cluster_meta_f", 
                                                                                                       "Observed_meta", "Shannon_meta",
                                                                                                       "Simpson_meta", "InvSimpson_meta")])
plot_alpha_vs <- dplyr::rename(plot_alpha_vs, 
                               metric_16s = variable,
                               val_16s = value)
# plot_alpha_vs should have 288 obs of 9 vars

# lengthen over WMGS alpha diversity
plot_alpha_vs <- reshape2::melt(plot_alpha_vs,
                                id.vars = c("sample", "cluster_f", "cluster_meta_f", 
                                            "metric_16s", "val_16s"),
                                measure.vars = colnames(plot_alpha_vs)[colnames(plot_alpha_vs) %!in% c("sample", "cluster_f", "cluster_meta_f", 
                                                                                                       "metric_16s", "val_16s")])
plot_alpha_vs <- dplyr::rename(plot_alpha_vs,
                               metric_meta = variable,
                               val_meta = value)
# plot_alpha_vs should have 1152 obs of 7 vars

# only want rows with same metric for 16S and WMGS data 
# drop rows with mismatched 16s and WMGS metrics
plot_alpha_vs <- plot_alpha_vs[substr(plot_alpha_vs$metric_16s, 1, 3) == substr(plot_alpha_vs$metric_meta, 1, 3),]
# plot_alpha_vs should have 288 obs of 7 vars



### generate plot colored by 16S cluster status with cluster-specific loess
plot_alpha_vs_16s <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta,
                                                      color = cluster_f, fill = cluster_f)) + 
  geom_smooth() +
  geom_point(aes(shape = cluster_meta_f)) + 
  scale_color_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  scale_fill_manual(values = c("#ca0020", "#ffd92f", "#0571b0"), 
                    labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                               expression(paste(italic("L. iners"), " dominated")),
                               "Mixed"),
                    "16S Cluster") +
  scale_shape_manual(values = c(16, 8),
                     labels = c(1, 2),
                     "WMGS Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 4, r = 0, b = 0, l = 0)),
        legend.text.align = 0) 
# plot_alpha_vs_16s



### generate plot colored by WMGS cluster status with cluster-specific loess
plot_alpha_vs_meta <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta,
                                                       color = cluster_meta_f, fill = cluster_meta_f)) + 
  geom_smooth() +
  geom_point(aes(shape = cluster_f)) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  scale_fill_manual(values = c("#e7298a",  "#e6ab02"),  
                    labels = c("1", "2"), 
                    "WMGS cluster") +
  scale_shape_manual(values = c(16, 8, 18),
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 4, r = 0, b = 0, l = 0)),
        legend.text.align = 0) 
# plot_alpha_vs_meta



### generate supplemental figure 7, patchwork
plot_alpha_vs_16s / plot_alpha_vs_meta +
  plot_annotation(tag_level = "A") +
  theme(plot.tag = element_text(size = 11))






##################################################################
### GENERATE SUPPLEMENTAL FIGURE 7 WMGS VS 16S ALPHA DIVERSITY ###
##################################################################



### using plot_alpha_vs from supplemental figure 8

### generate plot colored by WMGS cluster status with no loess
### only shape legend to pull for patchwork
plot_alpha_vs_meta <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta)) + 
  geom_point(aes(color = cluster_meta_f, fill = cluster_meta_f, shape = cluster_f)) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  scale_fill_manual(values = c("#e7298a",  "#e6ab02"),  
                    labels = c("1", "2"), 
                    "WMGS cluster") +
  scale_shape_manual(values = c(16, 8, 18),
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.direction = "horizontal", legend.text.align = 0) +
  guides(color = "none", fill = "none",
         shape = guide_legend(title.position = "top", title.hjust = 0.5)) 
plot_alpha_vs_meta

# pull legend for patchwork
legend_alpha_vs_meta_shape <- as_ggplot(get_legend(plot_alpha_vs_meta))



### generate plot colored by 16S cluster status with no loess
### no legend for patchwork
plot_alpha_vs_meta_nolegend <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta)) + 
  geom_point(aes(color = cluster_meta_f, fill = cluster_meta_f, shape = cluster_f)) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  scale_fill_manual(values = c("#e7298a",  "#e6ab02"),  
                    labels = c("1", "2"), 
                    "WMGS cluster") +
  scale_shape_manual(values = c(16, 8, 18),
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  labs(tag  = "A") +
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none")
plot_alpha_vs_meta_nolegend



### generate plot colored by 16S cluster status with overall loess
### no legend for patchwork
plot_alpha_vs_meta_loess <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta)) + 
  geom_smooth(color = "grey30") +
  geom_point(aes(color = cluster_meta_f, fill = cluster_meta_f, shape = cluster_f)) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  scale_fill_manual(values = c("#e7298a",  "#e6ab02"),  
                    labels = c("1", "2"), 
                    "WMGS cluster") +
  scale_shape_manual(values = c(16, 8, 18),
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  labs(tag  = "B") +
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none")
plot_alpha_vs_meta_loess





# pull legend for patchwork
legend_alpha_vs_meta_color <- as_ggplot(get_legend(plot_alpha_vs_meta_loess_cluster))



### generate plot colored by 16S cluster status with cluster-specific loess
### only color legend to pull for patchwork
plot_alpha_vs_meta_loess_cluster_nolegend <- ggplot(data = plot_alpha_vs, aes(x = val_16s, y = val_meta,
                                                                              color = cluster_meta_f, fill = cluster_meta_f)) + 
  geom_smooth() +
  geom_point(aes(shape = cluster_f)) + 
  scale_color_manual(values = c("#e7298a",  "#e6ab02"),  
                     labels = c("1", "2"), 
                     "WMGS cluster") +
  scale_fill_manual(values = c("#e7298a",  "#e6ab02"),  
                    labels = c("1", "2"), 
                    "WMGS cluster") +
  scale_shape_manual(values = c(16, 8, 18),
                     labels = c(expression(paste(italic("L. crispatus"), " dominated")),
                                expression(paste(italic("L. iners"), " dominated")),
                                "Mixed"),
                     "16S Cluster") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  facet_wrap(~plot_alpha_vs$metric_16s, scales = "free") +
  xlab("16S Alpha diversity") + 
  ylab("WMGS Alpha diversity") + 
  labs(tag  = "C") +
  theme_light(base_size = 9) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none")
plot_alpha_vs_meta_loess_cluster_nolegend



### generate supplemental figure 7, patchwork
plot_alpha_vs_meta_nolegend / plot_alpha_vs_meta_loess / plot_alpha_vs_meta_loess_cluster_nolegend / legend_alpha_vs_meta_shape / legend_alpha_vs_meta_color +
  guide_area() +
  plot_layout(design = c(area(1, 1, 10, 10), 
                         area(11, 1, 20, 10), 
                         area(21, 1, 30, 10),
                         area(31, 3, 33),
                         area(31, 8, 33))) + # TLBR
  theme(plot.tag = element_text(size = 11))





