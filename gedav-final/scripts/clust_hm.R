clust_hm <- function(corr_matrix, clust = FALSE, label_colors, factor, labels = colnames(corr_matrix)) {
    cut_cols <- "black"
    legend <- levels(factor)
    text.col <- rep("black", 3)
    pch <- rep(15, 3)
    col <- unique(label_colors)
    main <- "Correlation Heatmap"
    
    if (clust) {
        cut_cols <- ifelse(cutree(hclust(dist(corr_matrix)), k = 3) >= 2, "red", "black")
        legend <- c(levels(factor), "Potential Outlier")
        text.col <- c(rep("black", 3), "red")
        pch <- c(15, 15, 15, NULL)
        col <- c(unique(label_colors), NA)
        main <- "Correlation Heatmap\n(Clustered)"
    }
    
    hm <- heatmap.2(
        corr_matrix,
        dendrogram = ifelse(clust, "row", "none"),
        Rowv = ifelse(clust, clust, FALSE),
        Colv = clust,
        symm = TRUE,
        trace = "none",
        col = brewer.pal(11, "PiYG"),
        density.info = "none",
        margins = c(0.5, 5.5),
        labRow = labels,
        labCol = NA,
        cexRow = 0.6,
        lhei = c(1, 4),
        main = main,
        RowSideColors = label_colors,
        colRow = cut_cols,
        offsetRow = 0.02
    )
    
    legend(
        x = 0.55,
        y = 0.95,
        legend = legend,
        text.col = text.col,
        # cex = 0.6,
        pch = pch,
        col = col,
        xpd = NA,
        horiz = TRUE,
        xjust = 0.5,
        title = "Cancer Type"
    )
    return(hm)
} 
