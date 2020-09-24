####read data
library(readxl)
#url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
metadata_filename <- "PBMC8_metadata.xlsx"
#download.file(paste0(url, "/", metadata_filename), destfile = metadata_filename,
              #mode = "wb")
md <- read_excel(metadata_filename)

## Make sure condition variables are factors with the right levels
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
head(md)

## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00")
names(color_conditions) <- levels(md$condition)

###unzip data
#fcs_filename <- "PBMC8_fcs_files.zip"
#download.file(paste0(url, "/", fcs_filename), destfile = fcs_filename, 
              #mode = "wb")
#unzip(fcs_filename)

###load the content of the .fcs files into R

#if (!requireNamespace("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
#BiocManager::install("flowCore")
library(flowCore)
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, 
                        truncate_max_range = FALSE)
fcs_raw

pData(fcs_raw)

#####
panel_filename <- "PBMC8_panel.xlsx"
#download.file(paste0(url, "/", panel_filename), destfile = panel_filename, 
              #mode = "wb")
panel <- read_excel(panel_filename)
head(data.frame(panel))

# Replace problematic characters 
panel$Antigen <- gsub("-", "_", panel$Antigen)

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)

# Replace problematic characters 
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)

# Lineage markers
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Functional markers
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Spot checks
all(lineage_markers %in% panel_fcs$desc)

all(functional_markers %in% panel_fcs$desc)

##Data transformation
## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(lineage_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

##Diagnostic plots
## Generate sample IDs corresponding to each cell in the `expr` matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", 
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = expression, color = condition, 
                 group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions)


# Get the median marker expression per sample
library(dplyr)

expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% 
  summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions) 

###基于FlowSOM和ConensusClusterPlus的细胞群体识别

library(FlowSOM)

fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

## Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

#####必须加载这2个包
library(pheatmap)
library(RColorBrewer)
###########
plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = FALSE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend)
  
}

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers], 
                                expr01 = expr01[, lineage_markers], 
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)


##### Visual representation with t-SNE

## Find and skip duplicates
dups <- which(!duplicated(expr[, lineage_markers]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids) 

## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 1000)  

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]

## Run t-SNE
library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)

## Plot t-SNE colored by CD4 expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], 
                 expr[tsne_inds, lineage_markers])

ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD4)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn("CD4", 
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

#############
dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

## Plot t-SNE colored by clusters
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

########## Cluster merging and annotation
##Manual cluster merging and annotating based on heatmaps
cluster_merging1_filename <- "PBMC8_cluster_merging1.xlsx"
url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
download.file(paste0(url, "/", cluster_merging1_filename), 
              destfile = cluster_merging1_filename, mode = "wb")
cluster_merging1 <- read_excel(cluster_merging1_filename)
data.frame(cluster_merging1)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

mm <- match(code_clustering1, cluster_merging1$original_cluster)
code_clustering1m <- cluster_merging1$new_cluster[mm]
####
dr$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])
gg_tsne_cl1m <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
gg_tsne_cl1m

## Facet per sample
gg_tsne_cl1m + facet_wrap(~ sample_id) 

## Facet per condition
gg_tsne_cl1m + facet_wrap(~ condition)

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers],
                                expr01 = expr01[, lineage_markers], cell_clustering = cell_clustering1,
                                color_clusters = color_clusters, cluster_merging = cluster_merging1)

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers],
                                expr01 = expr01[, lineage_markers], cell_clustering = cell_clustering1m,
                                color_clusters = color_clusters)

###Differential analysis
library(lme4)
library(multcomp)
## Model formula without random effects
model.matrix( ~ condition, data = md)

## Create contrasts
contrast_names <- c("BCRXLvsRef")
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
K
FDR_cutoff <- 0.05

##### Differential cell population abundance
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100

counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props), 
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

## Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])

ggplot(ggdf, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_clusters) 

##########
ggdf$patient_id <- factor(md$patient_id[mm])

ggplot(ggdf) +
  geom_boxplot(aes(x = condition, y = proportion, color = condition, 
                   fill = condition),  position = position_dodge(), alpha = 0.5, 
               outlier.color = NA) +
  geom_point(aes(x = condition, y = proportion, color = condition, 
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2))

formula_glmer_binomial1 <- y/total ~ condition + (1|sample_id)
formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)

differential_abundance_wrapper <- function(counts, md, formula, K){
  ## Fit the GLMM for each cluster separately
  ntot <- colSums(counts)
  
  fit_binomial <- lapply(1:nrow(counts), function(i){
    
    data_tmp <- data.frame(y = as.numeric(counts[i, md$sample_id]), 
                           total = ntot[md$sample_id], md)
    
    fit_tmp <- glmer(formula, weights = total, family = binomial, 
                     data = data_tmp)
    
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
    })
    
    return(t(out))
  })
  
  ### Extract fitted contrast coefficients
  coeffs <- lapply(fit_binomial, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  colnames(coeffs) <- paste0("coeff_", contrast_names)
  rownames(coeffs) <- rownames(counts)
  
  ### Extract p-values
  pvals <- lapply(fit_binomial, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(counts)
  
  ### Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  
  return(list(coeffs = coeffs, pvals = pvals, adjp = adjp))
}

da_out1 <- differential_abundance_wrapper(counts, md = md, 
                                          formula = formula_glmer_binomial1, K = K)
apply(da_out1$adjp < FDR_cutoff, 2, table)

da_out2 <- differential_abundance_wrapper(counts, md = md, 
                                          formula = formula_glmer_binomial2, K = K)
apply(da_out2$adjp < FDR_cutoff, 2, table)

da_output2 <- data.frame(cluster = rownames(props),  props, 
                         da_out2$coeffs, da_out2$pvals, da_out2$adjp, row.names = NULL)
print(head(da_output2), digits = 2)
##################
normalization_wrapper <- function(expr, th = 2.5){
  expr_norm <- apply(expr, 1, function(x){ 
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0){
      x <- (x - mean(x, na.rm = TRUE))
    }else{ 
      x <- (x - mean(x, na.rm = TRUE)) / sdx
    }
    x[x > th] <- th
    x[x < -th] <- -th
    return(x)
  })
  expr_norm <- t(expr_norm)
}

plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp, 
                                              condition, color_conditions, th = 2.5){
  ## Order samples by condition
  oo <- order(condition)
  condition <- condition[oo]
  expr_norm <- expr_norm[, oo, drop = FALSE]
  
  ## Create the row labels with adj p-values and other objects for pheatmap
  labels_row <- paste0(rownames(expr_norm), " (", sprintf( "%.02e", sign_adjp), ")")
  labels_col <- colnames(expr_norm)
  annotation_col <- data.frame(condition = factor(condition))
  rownames(annotation_col) <- colnames(expr_norm)
  annotation_colors <- list(condition = color_conditions)
  color <- colorRampPalette(c("#87CEFA", "#56B4E9", "#56B4E9", "#0072B2", 
                              "#000000", "#D55E00", "#E69F00", "#E69F00", "#FFD700"))(100)
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  gaps_col <- as.numeric(table(annotation_col$condition))
  
  pheatmap(expr_norm, color = color, breaks = breaks, 
           legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, 
           labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, 
           annotation_col = annotation_col, annotation_colors = annotation_colors, 
           fontsize = 8)
}
## Apply the arcsine-square-root transformation
asin_table <- asin(sqrt((t(t(counts_table) / colSums(counts_table)))))
asin <- as.data.frame.matrix(asin_table)
## Keep significant clusters and sort them by significance
sign_clusters <- names(which(sort(da_out2$adjp[, "adjp_BCRXLvsRef"]) < FDR_cutoff))
## Get the adjusted p-values
sign_adjp <- da_out2$adjp[sign_clusters , "adjp_BCRXLvsRef", drop=FALSE]
## Keep samples for condition A and normalize to mean = 0 and sd = 1
asin_norm <- normalization_wrapper(asin[sign_clusters, ])
mm <- match(colnames(asin_norm), md$sample_id)
plot_differential_heatmap_wrapper(expr_norm = asin_norm, sign_adjp = sign_adjp, 
                                  condition = md$condition[mm], color_conditions = color_conditions)

########## Differential analysis of marker expression stratified by cell population

## Get median marker expression per sample and cluster
expr_median_sample_cluster_tbl <- data.frame(expr[, functional_markers], 
                                             sample_id = sample_ids, cluster = cell_clustering1m) %>%
  group_by(sample_id, cluster) %>% 
  summarize_all(funs(median))
## Melt
expr_median_sample_cluster_melt <- melt(expr_median_sample_cluster_tbl, 
                                        id.vars = c("sample_id", "cluster"), value.name = "median_expression", 
                                        variable.name = "antigen")
## Rearange so the rows represent clusters and markers
expr_median_sample_cluster <- dcast(expr_median_sample_cluster_melt, 
                                    cluster + antigen ~ sample_id,  value.var = "median_expression")
rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$cluster, 
                                               "_", expr_median_sample_cluster$antigen)
## Eliminate clusters with low frequency
clusters_keep <- names(which((rowSums(counts < 5) == 0)))
keepLF <- expr_median_sample_cluster$cluster %in% clusters_keep
expr_median_sample_cluster <- expr_median_sample_cluster[keepLF, ]
## Eliminate cases with zero expression in all samples
keep0 <- rowSums(expr_median_sample_cluster[, md$sample_id]) > 0
expr_median_sample_cluster <- expr_median_sample_cluster[keep0, ]
#####
ggdf <- expr_median_sample_cluster_melt[expr_median_sample_cluster_melt$cluster 
                                        %in% clusters_keep, ]
## Add info about samples
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$patient_id <- factor(md$patient_id[mm])
ggplot(ggdf) +
  geom_boxplot(aes(x = antigen, y = median_expression, 
                   color = condition, fill = condition), 
               position = position_dodge(), alpha = 0.5, outlier.color = NA) +
  geom_point(aes(x = antigen, y = median_expression, color = condition, 
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge(), 
             size = 0.7) +
  facet_wrap(~ cluster, scales = "free_y", ncol=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2))


differential_expression_wrapper <- function(expr_median, md, model = "lmer", formula, K){
  
  ## Fit LMM or LM for each marker separately
  fit_gaussian <- lapply(1:nrow(expr_median), function(i){
    
    data_tmp <- data.frame(y = as.numeric(expr_median[i, md$sample_id]), md)
    
    switch(model, 
           lmer = {
             fit_tmp <- lmer(formula, data = data_tmp)
           },
           lm = {
             fit_tmp <- lm(formula, data = data_tmp)
           })
    
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
    })
    
    return(t(out))
  })
  
  ### Extract fitted contrast coefficients
  coeffs <- lapply(fit_gaussian, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  colnames(coeffs) <- paste0("coeff_", contrast_names)
  rownames(coeffs) <- rownames(expr_median)
  
  ### Extract p-values
  pvals <- lapply(fit_gaussian, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(expr_median)
  
  ### Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  
  return(list(coeffs = coeffs, pvals = pvals, adjp = adjp))
}

formula_lm <- y ~ condition
formula_lmer <- y ~ condition + (1|patient_id)

de_out1 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster, 
                                           md = md, model = "lm", formula = formula_lm, K = K)
apply(de_out1$adjp < FDR_cutoff, 2, table)

de_out2 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster, 
                                           md = md, model = "lmer", formula = formula_lmer, K = K)
apply(de_out2$adjp < FDR_cutoff, 2, table)

de_output2 <- data.frame(expr_median_sample_cluster, 
                         de_out2$coeffs, de_out2$pvals, de_out2$adjp, row.names = NULL)
print(head(de_output2), digits = 2)


## Keep the significant markers, sort them by significance and group by cluster
sign_clusters_markers <- names(which(de_out2$adjp[, "adjp_BCRXLvsRef"] < FDR_cutoff))
oo <- order(expr_median_sample_cluster[sign_clusters_markers, "cluster"], 
            de_out2$adjp[sign_clusters_markers, "adjp_BCRXLvsRef"])
sign_clusters_markers <- sign_clusters_markers[oo]

## Get the significant adjusted p-values
sign_adjp <- de_out2$adjp[sign_clusters_markers , "adjp_BCRXLvsRef"]

## Normalize expression to mean = 0 and sd = 1
expr_s <- expr_median_sample_cluster[sign_clusters_markers,md$sample_id]
expr_median_sample_cluster_norm <- normalization_wrapper(expr_s)

mm <- match(colnames(expr_median_sample_cluster_norm), md$sample_id)
plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_cluster_norm, 
                                  sign_adjp = sign_adjp, condition = md$condition[mm],
                                  color_conditions = color_conditions)
