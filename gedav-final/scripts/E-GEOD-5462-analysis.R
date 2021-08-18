# Transcription profiling of human sequential breast cancer biopsies during letrozole treatment ----
# https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-5462
# In the present investigation, we have exploited the opportunity provided by
# neoadjuvant treatment of a group of postmenopausal women with large operable
# or locally advanced breast cancer (in which therapy is given with the primary
# tumour remaining within the breast) to take sequential biopsies of the same
# cancers before and after 10-14 days treatment with letrozole. RNA extracted
# from the biopsies has been subjected to Affymetrix microarray analysis and the
# data from paired biopsies interrogated to discover genes whose expression is
# most influenced by oestrogen deprivation. Experiment Overall Design: biopsies
# were taken from the same subjects both pretreatment and after 10-14 days
# Letrozol, 2.5 mg/day, oral

# LOAD LIBRARIES ----------------------------------------------------------
library(glue)
library(stringr)
library(dendextend)
library(pheatmap)
library(marray)
library(limma)
library(genefilter)
library(oligo)
library(pd.hg.u133a)
library(hgu133a.db)
library(gplots)
library(affycoretools)
library(RSQLite)
library(DBI)
library(RColorBrewer)
library(multtest)
library(MASS)
source("scripts/clust_hm.R")


# SET CONSTANTS -----------------------------------------------------------
set_id <- "E-GEOD-5462"
data_dir <- file.path("data", set_id)
raw_dir <- file.path(data_dir, "raw")
processed_dir <- file.path(data_dir, "RData_objs")
img_dir <- file.path("img")
ann_url <-
  glue("https://www.ebi.ac.uk/arrayexpress/files/{set_id}/{set_id}.sdrf.txt")
ann_file <- file.path(data_dir, glue("{set_id}.sdrf.txt"))

# SET UP DIRECTORY STRUCTURE ----------------------------------------------
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir)
}
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir)
}
if (!dir.exists("img")) {
  dir.create("img")
}

# DOWNLOAD & FORMAT ANNOTATION --------------------------------------------
# Download annotation file
if (!file.exists(ann_file)) {
  download.file(ann_url, ann_file)
}
ann <- read.delim(ann_file)
ann.filt <- ann[, c(
  "Source.Name",
  "Hybridization.Name",
  "Array.Data.File",
  "Comment..ArrayExpress.FTP.file."
)]

splitFactors <- as.data.frame(t(sapply(ann.filt$Hybridization.Name, 
                                       function(x) {sort(str_split(x, ";", simplify = TRUE))})))
colnames(splitFactors) <- c("PatientID", "Tissue", "Gender", "Treatment")
splitFactors$PatientID <- as.numeric(str_replace(splitFactors$PatientID, "A|B", ""))

ann.filt.final <- data.frame(ann.filt[, -2], splitFactors)
row.names(ann.filt.final) <- ann.filt.final$Array.Data.File

# Reorder so samples are in the same order (PatientID, Tmt)
ann.filt.final <- ann.filt.final[order(ann.filt.final$PatientID, ann.filt.final$Treatment),]
ann.filt.final$Treatment.Factors <- factor(ann.filt.final$Treatment, 
                                           levels = unique(ann.filt.final$Treatment),
                                           labels = c("Letrozole", "Pretreatment"))
ann.filt.final$Sample.Name <- paste(ann.filt.final$Treatment.Factors,
                                    ann.filt.final$PatientID,
                                    sep = "_")


# READ IN RAW DATA --------------------------------------------------------
# Create annotated data frame of the annotation file
ann.filt.ADF <- AnnotatedDataFrame(ann.filt.final)

# Download/extract raw data celfiles
if (!all(ann$Array.Data.File %in% list.files(raw_dir))) {
  raw_zips <- unique(ann$Comment..ArrayExpress.FTP.file.)
  for (i in 1:length(raw_zips)) {
    destzip <- file.path(data_dir, glue("raw_{i}.zip"))
    if (!file.exists(destzip)) {
      options(timeout = 180)
      download.file(raw_zips[i], destfile = destzip)
    }
    unzip(destzip, exdir = raw_dir, overwrite = TRUE)
  }
}

# Read in the raw CEL files
if (!file.exists(file.path(processed_dir, "raw_data.RData"))) {
  message("Importing CEL files. Saving an RData object for quicker loading.")
  raw_data <-
    oligo::read.celfiles(
      filenames = file.path(raw_dir, rownames(ann.filt.final)),
      phenoData = ann.filt.ADF,
      verbose = TRUE
    )
  saveRDS(object = raw_data,
          file = file.path(processed_dir, "raw_data.RData"))
  message(glue("Saved as: '{file.path(processed_dir, 'raw_data.RData')}'"))
} else {
  message(glue("Reading from file: '{file.path(processed_dir, \"raw_data.RData\")}'"))
  raw_data <-
    readRDS(file = file.path(processed_dir, "raw_data.RData"))
}

# QUALITY ANALYSIS --------------------------------------------------------
# How many of each sample type are there?
table(pData(raw_data)$Treatment.Factors)
# They are also paired samples
table(table(pData(raw_data)$PatientID) == 2)

# log2 data boxplot
set.seed(1234)
tmt <- pData(raw_data)$Treatment.Factors
labelCols <- c("green4", "yellow3")[tmt]
png(file = "img/1-raw-boxplot.png", width = 1115, height = 796)
oligo::boxplot(
  raw_data,
  main = "Boxplot of raw data (log2 transformed)",
  ylab = "log2(Intensity)",
  col = labelCols,
  las = 2,
  cex.axis = 0.5,
  names = pData(raw_data)$Sample.Name
)
legend(
  "topleft",
  legend = levels(pData(raw_data)$Treatment.Factors),
  pch = 15,
  col = unique(labelCols),
  horiz = T,
  box.lty = 0,
  bg = "transparent"
)
dev.off()

# RLE boxplot
if (!file.exists(file.path(processed_dir, "rle_data.RData"))) {
  message("Performing RMA *without* quantile normalization.")
  rle_data <- rma(raw_data, normalize = FALSE)
  saveRDS(object = rle_data,
          file = file.path(processed_dir, "rle_data.RData"))
  message(glue("Saved as: '{file.path(processed_dir, 'rle_data.RData')}'"))
} else {
  message(glue("Reading from file: '{file.path(processed_dir, 'rle_data.RData')}'"))
  rle_data <-
    readRDS(file = file.path(processed_dir, "rle_data.RData"))
}
eset_row_medians <- apply(exprs(rle_data), 1, median)
RLE_prerma <-
  apply(exprs(rle_data), 2, function(x)
    x - eset_row_medians)

png(file = "img/2-rle-boxplot.png", width = 1115, height = 796)
oligo::boxplot(
  RLE_prerma,
  transfo = identity,
  main = "RLE Boxplot (non-normalized data)",
  ylab = "Relative log expression",
  col = labelCols,
  las = 2,
  cex.axis = 0.5,
  names = rle_data$Sample.Name,
  outline = FALSE
)
legend(
  "topleft",
  legend = levels(tmt),
  pch = 15,
  col = unique(labelCols),
  horiz = T,
  box.lty = 0,
  bg = "transparent"
)
dev.off()

# OUTLIER DETECTION -------------------------------------------------------
exprs_rle <- exprs(rle_data)
colnames(exprs_rle) <- rle_data$Sample.Name
tmt <- rle_data$Treatment.Factors
labelCols <- c("green4", "yellow3")[tmt]

# Check how the samples separate
pcs <- prcomp(t(exprs_rle))$x
pc_var <- round(apply(pcs, 2, var) / sum(apply(pcs, 2, var)) * 100, 2)
png(file = "img/3-sample-pca.png", width = 1115, height = 796)
plot(
  pcs[, 2] ~ pcs[, 1],
  type = "n",
  ylab = glue("PC2 ({pc_var[2]}%)"),
  xlab = glue("PC1 ({pc_var[1]}%)"),
  main = paste("PC2 vs PC1 (all samples)", "Sample type colored", sep = "\n")
)
text(
  pcs[, 2] ~ pcs[, 1],
  labels = colnames(exprs_rle),
  col = labelCols,
  font = 2
)
legend(
  "topleft",
  legend = levels(tmt),
  text.col = unique(labelCols),
  text.font = 2,
  inset = 0.02
)
dev.off()
# It doesn't appear that there is much within patient sample bias

# Correlation
corr <- cor(exprs_rle, use = "pairwise.complete.obs")

## Correlation heatmaps
# without clustering
png(file = "img/4-noclust-heatmap.png", width = 1115, height = 796)
hm_noclust <- clust_hm(corr, factor = tmt, label_colors = labelCols)
dev.off()

# with clustering
png(file = "img/5-clust-heatmap.png", width = 1115, height = 796)
hm_clust <- clust_hm(corr, clust = TRUE, factor = tmt, label_colors = labelCols)
dev.off()

# get outliers based on cluster dendrogram
out_clust <- cutree(hclust(dist(corr)), k = 3) >= 2
(hm_outliers <- names(out_clust[out_clust]))

## Average correlation plots
avg_corr <- apply(corr, 1, mean)
png(file = "img/6-avg-corr.png", width = 1115, height = 796)
plot(
  range(avg_corr) ~ range(1, length(avg_corr)),
  type = "n",
  axes = F,
  ylab = "Average Sample-sample Correlation",
  xlab = "",
  main = glue("Average correlation plot\n{set_id}"),
  ylim = c(0.75, 0.95)
)
box()
points(
  avg_corr,
  pch = 21,
  bg = ifelse(avg_corr < 0.9, "red", "lightblue")
)
abline(h = 0.9, col = rgb(1, 0, 0, 0.5))
axis(1, at = seq(1, 116), labels = colnames(exprs_rle), las = 2, cex.axis = 0.5)
axis(2, las = 2, cex.axis = 0.6, at = seq(0.75, 0.95, by = 0.025))
legend("bottomleft", cex = 0.8, 
       legend = c("Above Threshold", "Potential Outliers", "Threshold (90% correlation)"),
       pt.bg = c("lightblue", "red", NA), pch = c(21, 21, NA), 
       lty = c(0, 0, 1), col = c("black", "black", "red"), 
       inset = 0.02, bty = "n")
dev.off()

# Grab outliers from this analysis
ac_outliers <- names(avg_corr[avg_corr < 0.9])

# Find samples who cluster differently and whose average correlation falls below
# the threshold
outliers <- intersect(hm_outliers, ac_outliers)

# Outliers which show up in both tests are considered true outliers
# extract the patient IDs
(outlierPatientID <- unique(as.numeric(str_extract(outliers, "\\d+"))))
filter <- rle_data$PatientID %in% outlierPatientID
table(filter)

# NORMALIZATION -------------------------------------------------------
if (!file.exists(file.path(processed_dir, "rma_data.RData"))) {
  message("Normalizing data with RMA, then saving for future use.")
  # RMA normalization after outlier removal
  rma_data <- rma(raw_data[, !filter])
  saveRDS(object = rma_data,
          file = file.path(processed_dir, "rma_data.RData"))
  message(glue("Saved as: \"{file.path(processed_dir, 'rma_data.RData')}\""))
} else {
  message(glue("Reading from file: '{file.path(processed_dir, 'rma_data.RData')}'"))
  rma_data <-
    readRDS(file = file.path(processed_dir, "rma_data.RData"))
}
# RMA boxplot
set.seed(1234)
tmt <- rma_data$Treatment.Factors
labelCols <- c("green4", "yellow3")[tmt]
png(file = "img/6-rma-boxplot.png", width = 1115, height = 796)
oligo::boxplot(
  rma_data,
  transfo = identity,
  main = "RMA Boxplot (normalized data)",
  ylab = "log2(Normalized Intensity)",
  las = 2, col = cols_rma,
  cex.axis = 0.5,
  names = rma_data$Sample.Name
)
legend(
  "topleft",
  legend = levels(tmt),
  pch = 15,
  col = unique(labelCols),
  horiz = T,
  box.lty = 0,
  bg = "transparent"
)
dev.off()

# See how normalization affected RLE
rma_row_medians <- apply(exprs(rma_data), 1, median)
RLE_postrma <-
  apply(exprs(rma_data), 2, function(x)
    x - rma_row_medians)

png(file = "img/7-rle-postrma-boxplot.png", width = 1115, height = 796)
oligo::boxplot(
  RLE_postrma,
  transfo = identity,
  main = "RLE Boxplot (normalized data)",
  ylab = "Relative log expression",
  las = 2, col = labelCols,
  cex.axis = 0.5,
  names = rma_data$Sample.Name,
  outline = FALSE
)
legend(
  "topleft",
  legend = levels(tmt),
  pch = 15,
  col = unique(labelCols),
  horiz = T,
  box.lty = 0,
  bg = "transparent"
)
dev.off()

exprs_rma <- exprs(rma_data)
colnames(exprs_rma) <- rma_data$Sample.Name
tmt <- rma_data$Treatment.Factors
labelCols <- c("green4", "yellow3")[tmt]

# Dendrogram
# NOT USED
png(file = "img/7-dendrogram.png", width = 1115, height = 796)
dend <- dist(t(exprs_rma)) %>%
  hclust() %>%
  as.dendrogram() %>%
  set("labels_cex", 0.6) %>%
  set("labels_colors", labelCols[order.dendrogram(.)]) %>%
  set("branches_lwd", 1.5)
plot(dend,
     axes = FALSE,
     ylab = "Height",
     font = 2)
axis(2, lwd = 2, las = 2)
title(paste("Cluster Dendrogram", glue("{set_id} data set"), sep = "\n"))
dev.off()

# Sample PCA
# NOT USED
pcs <- prcomp(t(exprs_rma))$x
pc_var <- round(apply(pcs, 2, var) / sum(apply(pcs, 2, var)) * 100, 2)
png(file = "img/8-rma-pca.png", width = 1115, height = 796)
plot(
  pcs[, 2] ~ pcs[, 1],
  type = "n",
  ylab = glue("PC2 ({pc_var[2]}%)"),
  xlab = glue("PC1 ({pc_var[1]}%)"),
  main = paste("PC2 vs PC1 (normalized data)", "Sample type colored", sep = "\n")
)
text(pcs[, 2] ~ pcs[, 1],
     labels = colnames(exprs_rma),
     col = labelCols, font = 2)
legend(
  "topleft",
  legend = levels(tmt),
  text.col = c("green4", "yellow3"),
  text.font = 2,
  inset = 0.02
)
dev.off()

# FILTERING GENES -----------------------------------------
# Remove lowly expressed genes
# 1) Calculate row cv of each gene for each group
row_cv <- apply(rma_data, 1, function(thisRow) {sd(thisRow)/mean(thisRow)})

# 2) Choose a cutoff value to remove enrichment of invariant transcripts
cv.cutoff <- shorth(row_cv)
png(file = "img/9-hist_row_cv.png", width = 1115, height = 796)
hist(
  x = row_cv,
  breaks = 100,
  col = "coral4",
  xlab = "Coefficient of Variation",
  main = glue("Histogram of CVs\n(normalized data, all genes)\n{set_id}")
)
abline(v = cv.cutoff, lwd = 2.5, col = "red")
text(x = cv.cutoff, y = -150, srt = -45, adj = c(0, 0),
     label = glue("cutoff: {round(cv.cutoff,2)}"), 
     xpd = NA, font = 2, col = "red")
dev.off()

# 3) Find which rows contain transcript levels greater than cutoff in at least the
# minimum number of samples
keep <- row_cv >= cv.cutoff

# 4) Filter out low expressing genes from eset
rma_dat_keep <- rma_data[keep,]

# Histograms pre-/post-filtering
png(file = "img/10-hist_mean_exp.png", width = 1115, height = 796)
par(mfrow = c(2, 1))
hist(
  x = apply(rma_data, 1, mean),
  breaks = 100,
  freq = FALSE,
  col = "coral4",
  xlab = "Row mean expression",
  xlim = c(5, 14),
  main = glue("Mean expression, all genes (normalized)\n{set_id}")
)
hist(
  x = apply(exprs(rma_dat_keep), 1, mean),
  breaks = 100,
  freq = FALSE,
  col = "coral4",
  xlim = c(5, 14),
  xlab = "Row mean expression",
  main = glue("Mean expression, filtered genes (normalized)\n{set_id}")
)
dev.off()

# Annotate filtered data
rma_dat_keep <- annotateEset(rma_dat_keep, hgu133a.db, columns = c("PROBEID", "SYMBOL", "GENENAME"))

# FEATURE SELECTION ------------------------------------------------------------
(tmt <- rma_dat_keep$Treatment.Factors)
(tmt <- relevel(tmt, ref = "Pretreatment"))
design <- model.matrix(~tmt)

# Fit linear model
fit <- lmFit(rma_dat_keep, design)

# Define contrasts and perform empirical Bayes
fit2 <- eBayes(fit)
tt <- topTable(fit2, number = Inf, coef = "tmtLetrozole", sort.by = "none")

# Verify that the p-values extracted are in the same order
all(fit2$p.value[,2] == tt$P.Value)

# MULTIPLICITY ADJUSTMENT -------------------------------------------------
# Extract raw p-values
raw.p <- tt$P.Value
names(raw.p) <- rownames(tt)
  
# Extract log fold change
lfc <- tt$logFC
names(lfc) <- rownames(tt)
abs.fc <- 2^abs(lfc)

# Adjust for multiple comparisons
adj.p <- p.adjust(raw.p, method = "fdr")
pv.cutoff <- 0.05

# Retained genes associated with p-value cutoff
table(adj.p < pv.cutoff)
p.trans <- -log10(adj.p)

png(file = "img/11-histogram-fdrs.png", width = 1115, height = 796)
# Histogram of all adjusted p-values (FDRs)
par(mfrow = c(2,1))
hist(p.trans, breaks = 100, col = "lightblue", axes = FALSE, 
     main = "Histogram of all FDRs", xlab = "-log10(FDR)")
axis(1, cex.axis = 0.8)
axis(2, cex.axis = 0.8, las = 2)
abline(v = -log10(pv.cutoff), lwd = 1.5, col = "red")
# Histogram of significant FDRs
hist(p.trans[adj.p < pv.cutoff], 100, col = "lightblue", axes = FALSE, 
     main = glue("Histogram of significant FDRs < {pv.cutoff}"), xlab="-log10(FDR)")
axis(1, cex.axis = 0.8)
axis(2, cex.axis = 0.8, las = 2)
dev.off()

# Fold change cutoff
fc.cutoff <- 1.25
table(abs.fc > fc.cutoff)
table(adj.p < pv.cutoff)

# Volcano plot of retained genes
png(file = "img/12_1-volcano-plot.png", width = 1115, height = 796)
plot(p.trans ~ lfc, pch = 16, cex = 0.6,
     xlab = "Log2(FC)", ylab = "-log10(p-value)",
     main = "Volcano Plot of Selected Features")
points(p.trans[p.trans > -log10(pv.cutoff) & lfc > log2(fc.cutoff)] ~ lfc[p.trans > -log10(pv.cutoff) & lfc > log2(fc.cutoff)], pch = 21, bg = "green")
points(p.trans[p.trans > -log10(pv.cutoff) & lfc < -log2(fc.cutoff)] ~ lfc[p.trans > -log10(pv.cutoff) & lfc < -log2(fc.cutoff)], pch = 21, bg = "red")
abline(h = -log10(pv.cutoff), v = c(log2(fc.cutoff), -log2(fc.cutoff)))
dev.off()

# Subset data
cutoff.pv <- adj.p[adj.p < pv.cutoff]
cutoff.fc <- abs.fc[abs.fc > fc.cutoff]

# Find genes which meet the cutoff 
final_genes <- intersect(names(cutoff.pv), names(cutoff.fc))
final_gene_rows <- which(rownames(exprs(rma_dat_keep)) %in% final_genes)
filtered_rma_dat_keep <- rma_dat_keep[final_gene_rows, ]

# How many features are left?
dim(filtered_rma_dat_keep) # 125 features remain

# They're all ones with significant pvalues and absolute fold changes > 1.5
all(fData(filtered_rma_dat_keep)$PROBEID == final_genes)

# HCA CLUSTERING POST-FILTER --------------------------------------------------
exprs_filt_rma <- exprs(filtered_rma_dat_keep)
colnames(exprs_filt_rma) <- filtered_rma_dat_keep$Sample.Name
ann_col <- data.frame(Tmt.Status = tmt)
rownames(ann_col) <- colnames(exprs_filt_rma)
png(file = "img/12_2-clust-post-filt.png", width = 1115, height = 796)
pheatmap(exprs_filt_rma,
         main = paste("Clustered Heatmap", glue("(Genes with FDR < {pv.cutoff} and |FC| > {fc.cutoff})")),
         border_color = NA,
         color = colorRampPalette(c("red", "black", "green"))(100),
         annotation_col = ann_col,
         annotation_colors = list(Tmt.Status=c(Pretreatment = "yellow3", Letrozole = "green4")),
         scale = "row",
         cutree_rows = 3,
         cutree_cols = 2,
         angle_col = 90,
         fontsize_col = 6,
         fontsize_row = 6
         )
dev.off()

# CLASSIFICATION WITH FILTERED GENES --------------------------------------
datx <- as.data.frame(t(exprs(filtered_rma_dat_keep)))
clas <- as.character(filtered_rma_dat_keep$Treatment.Factors)
datx <- data.frame(clas, datx)

# Truncate the classification labels
datx.rownames <- filtered_rma_dat_keep$Sample.Name
datx.rownames <- gsub("^Letrozole", "let", datx.rownames)
datx.rownames <- gsub("^Pretreatment", "pre", datx.rownames)

colnames(datx) <- c("tmt", fData(filtered_rma_dat_keep)$PROBEID)
rownames(datx) <- datx.rownames
datx <- datx[order(rownames(datx)), ]

# create train/test sets on 60%/40% of data
set.seed(11)
pre <- grep("^Pre", datx$tmt)
post <- grep("^Let", datx$tmt)
train <- c(sample(pre, 0.6*length(pre)),
           sample(post, 0.6*length(post)))

# subset based on their classification
train_set <- datx[train, ]
test_set <- datx[-train, ]

# pull out their classes and remove from data frame
classes <- list(
  train = train_set[,1],
  test = test_set[,1]
)
train_set <- train_set[,-1]
test_set <- test_set[,-1]

# train and test model
datx.train <- lda(classes$train ~ ., train_set)
datx.pred <- predict(datx.train, test_set)

# confusion matrix
(tab <- table(datx.pred$class, classes$test))

# correct classification
sum(diag(tab))
# error rate
(1-sum(diag(tab))/sum(tab)) * 100

# discriminant function plot
actual.class <- c(21, 22)[as.factor(classes$test)]
predicted.class <- c("red2", "green3")[as.factor(datx.pred$class)]
png(file = "img/13-classification.png", width = 1115, height = 796)
plot(datx.pred$x,
     pch = actual.class,
     bg = predicted.class, 
     axes = FALSE, 
     ylab = "Discriminant Function",
     xlab = "",
     cex = 1.5,
     main = glue("Discriminant Function\n({set_id} data set)"))
axis(side=1, las = 2,
     labels = rownames(datx.pred$x),
     cex = 0.6,
     at = 1:length(rownames(datx.pred$x)))
axis(side=2)
abline(v = 0.5*(nrow(datx.pred$x)+1), col = "blue", lwd = 1.5)
box()  
legend(
  "topleft",
  legend = c(
    "Point shape" ,
    "Actual: pretreatment",
    "Actual: letrozole",
    "",
    "Point color",
    "Prediction: pretreatment",
    "Prediction: letrozole",
    "",
    "Classification divider"
  ),
  pch = c(NA, unique(actual.class), NA,  NA, 18, 18, NA, NA),
  pt.bg = c(NA, "black", "black", NA, NA, NA, NA, NA, NA),
  col = c(NA, NA, NA, NA, NA, unique(predicted.class), NA, "blue"),
  lty = c(0, 0, 0, 0, 0, 0, 0, 0, 1),
  lwd = c(0, 0, 0, 0, 0, 0, 0, 0, 1.5),
  bty = "n",
  pt.cex = 1.5,
  inset = 0.02
)
dev.off()

# FUNCTIONAL ANALYSIS -----------------------------------------------------
# get top 5 upregulated genes
lfc_up <- sort(lfc, decreasing = TRUE)[1:5]
2^lfc_up
# get top 5 downregulated genes
lfc_down <- sort(lfc, decreasing = FALSE)[1:5]
2^abs(lfc_down)

# concatenate them into a vector
top_10_up_down <- c(lfc_up, lfc_down)

# write top 10 genes vector to a file for functional analysis
write.table(
  x = names(top_10_up_down),
  file = file.path(data_dir, "top10_up_down.txt"),
  quote = FALSE,
  sep = "\n",
  row.names = FALSE,
  col.names = FALSE
)
