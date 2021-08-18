library(MASS)
source("scripts/k.speClust2.R")

# Look at the directory where our data is stored
dir("data")

# Load data and annotation files with header/rownames set
dat <- read.table("data/sotiriou_data.txt",
                  header = TRUE,
                  row.names = 1)
ann <- read.table("data/sotiriou_annotations.txt",
                  header = TRUE,
                  row.names = 1)

# Make factor of "site" column
ann$site <- as.factor(ann$site)
str(ann)

# Principal component calculation
dat.pca <-prcomp(t(dat))
pcs <- dat.pca$x

# Plot PCA PC2 vs PC1
plot(
    pcs[, 2] ~ pcs[, 1],
    bg = c("red", "blue")[ann$site],
    pch = 21, cex = 0.8, xlab = "PC1", ylab = "PC2",
    main = "PCA of Sotiriou data\nPC2 vs PC1"
)

# Add legend to display factor and color
legend(
    "bottomleft", legend = levels(ann$site),
    pt.bg = c("red", "blue"),
    pch = 21, cex = 0.8, inset = 0.02
)

# Get variation % of all eigenvalues
dat.pca.scree <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2) * 100, 2)

# Filter out those explaining less than 1% variance
dat.pca.scree <- dat.pca.scree[dat.pca.scree > 1]

# Scree plot
plot(
    x = 1:length(dat.pca.scree),
    y = dat.pca.scree,
    type = "b", pch = 21,
    bg = "purple",
    ylab = "Percent variance",
    xlab = "# of Components"
)
title("Scree plot of Sotiriou data
      (components accounting for >=1% of variance)")

# Add up variance percent of first two eigenvalues
first2eigenvalues <- sum(dat.pca.scree[1:2])

# Graphical parameters
par(mfrow = c(1, 2), cex = 0.64, las = 2)

# Distance matrix
dat.dist <- dist(t(dat))

# Classic MDS
dat.loc <- cmdscale(dat.dist)
plot(dat.loc[, 2] ~ dat.loc[, 1], xlab = "PC1", ylab = "PC2",
     pch = 21, bg = c("red", "blue")[ann$site],
     main = "Classic MDS plot\nSotiriou data")
legend("bottomleft", legend = levels(ann$site),
       inset = 0.02, pch = 21, pt.bg = c("red", "blue"))

# Non-metric MDS
dat.mds <- isoMDS(dat.dist, trace = FALSE)
plot(dat.mds$points[,2] ~ dat.mds$points[,1], xlab = "PC1", ylab = "PC2",
     pch = 21, bg = c("red", "blue")[ann$site],
     main = paste0("Non-metric MDS plot\nSotiriou data (Stress = ",
                   round(dat.mds$stress, 2), ")"))
legend("bottomleft", legend = levels(ann$site),
       inset = 0.02, pch = 21, pt.bg = c("red", "blue"))

# Scale and center the rows of our matrix
temp <- scale(t(dat), center = TRUE, scale = TRUE)

# Calculate/plot 2D embedding of weighted graph Laplacian
phi <- k.speClust2(X = t(temp))
plot(
    phi[, 2] ~ phi[, 1], pch = 21, cex = 0.8,
    xlab = "phi1", ylab = "phi2",
    bg = c("red", "blue")[ann$site],
    main="Weighted Graph Laplacian\nSotiriou data"
)
legend(
    "bottomleft", legend = levels(ann$site),
    inset = 0.02, pch = 21, pt.bg = c("red", "blue")
)

