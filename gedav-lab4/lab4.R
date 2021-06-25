# Load the marray library and the swirl data set
library("marray")
data("swirl")

# Plot an MvA plot of array 3 without any stratified lines.
maPlot(swirl[, 3],
       lines.func = NULL,
       legend.func = NULL,
       main = "Print-tip Loess pre-normalization")
dev.off()

# Normalize array 3 by global median location normalization.
median.norm.swirl <- maNorm(swirl[, 3], norm = "median")

# Plot an MvA plot of the normalized array without the stratified lines or
# legend.
maPlot(
    median.norm.swirl,
    lines.func = NULL,
    legend.func = NULL,
    main = "Print-tip Loess post-normalization (median)"
)
dev.off()

# Repeat #3 and #4 applying loess global intensity normalization.
loess.norm.swirl <- maNorm(swirl[, 3], norm = "loess")
maPlot(
    loess.norm.swirl,
    lines.func = NULL,
    legend.func = NULL,
    main = "Print-tip Loess post-normalization (loess)"
)
dev.off()

# Go to the course website and retrieve the compressed file called ‘GenePix
# files’. Open it up and put the contents in a directory.
dir(path = "./data")
system(command = "unzip -o data/GenePix_files.zip -d data/gp")
dir(path = "./data/gp")

a.cdna <-
    read.GenePix(
        path = "data/gp",
        name.Gf = "F532 Median",
        name.Gb = "B532 Median",
        name.Rf = "F635 Median",
        name.Rb = "B635 Median",
        name.W = "Flags"
    )

norm_and_plot <- function(patient) {
    
    # normalize
    n.norm <- maNorm(patient, norm = "n")
    p.norm <- maNorm(patient, norm = "p")
    s.norm <- maNorm(patient, norm = "s")
    
    # 3 plots on one page
    par(mfrow=c(3,1))
    
    # plots
    maPlot(m = patient, lines.func = NULL, legend.func = NULL, main = "none normalization")
    maPlot(m = p.norm, lines.func = NULL, legend.func = NULL, main = "printTipLoess normalization")
    maPlot(m = s.norm, lines.func = NULL, legend.func = NULL, main = "scalePrintTipMAD normalization")
    mtext(text = paste(deparse(substitute(patient)), "plots"), outer = TRUE, line = -1.5, font = 2)
}

norm_and_plot(a.cdna[,1])
norm_and_plot(a.cdna[,2])

p.norm <- maNorm(a.cdna, norm = "p")
s.norm <- maNorm(a.cdna, norm = "s")

p.norm.df <- data.frame(maW(p.norm))
rownames(p.norm.df) <- make.names(probes, unique = TRUE)
colnames(p.norm.df) <- c("patient1", "patient2")

s.norm.df <- data.frame(maW(s.norm))
rownames(s.norm.df) <- make.names(probes, unique = TRUE)
colnames(s.norm.df) <- c("patient1", "patient2")

pkgs <- c("affy", "limma", "affyPLM", "fpc")
sa <- "simpleaffy"
(not_installed <- pkgs[!(pkgs %in% installed.packages())])
if (length(not_installed)) {
    install.packages(not_installed, dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages(sapply(pkgs, library, character.only = TRUE, quietly = TRUE))


    
    