norm_and_plot2 <- function(array, method=c("none", "median", "loess", "printTipLoess"), plot = TRUE) {
    # named list of normalized data
    norms <- list(
        NotNormalized = maNorm(array, norm = method[1]),
        MedianNormalized = maNorm(array, norm = method[2]),
        LoessNormalized = maNorm(array, norm = method[3]),
        PrintTipLoessNormalized = maNorm(array, norm = method[4])
    )
    if (plot) {
        # graph normalized data
        par(mfrow = c(2, 2), lwd = 2)
        for (i in seq_along(norms)) {
            d.method <- names(norms)[i]
            d.name <- deparse(substitute(array))
            maPlot(
                m = norms[[i]],
                lines.func = NULL,
                legend.func = NULL,
                main = paste(d.name, d.method, sep = ": ")
            )
        }
    }
    else {
        return(norms)
    }
}
