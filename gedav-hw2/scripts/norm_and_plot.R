# define function to reduce redundancy
norm_and_plot <- function(array, plot = TRUE) {
    # named list of normalized data
    norms <- list(
        notNormalized = maNorm(array, norm = "n"),
        medianNormalized = maNorm(array, norm = "m"),
        loessNormalized = maNorm(array, norm = "l"),
        printTipLoess = maNorm(array, norm = "p")
    )
    if (plot) {
        # graph normalized data
        par(mfrow = c(2, 2))
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