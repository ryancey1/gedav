px <- function(x, y, col = par("col"), pch = par("pch"), lm.col = "red", ...) {
    points(x, y, pch = pch, col = col, ...)
    abline(0, 1, col = lm.col)
}