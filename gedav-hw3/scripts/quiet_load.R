quiet_load <- function(x) {
    if (x %in% installed.packages()) {
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
}