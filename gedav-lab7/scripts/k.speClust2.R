k.speClust2 <- function (X, qnt = NULL) {
    dist2full <- function(dis) {
        n <- attr(dis, "Size")
        full <- matrix(0, n, n)
        full[lower.tri(full)] <- dis
        full + t(full)
    }
    dat.dis <- dist(t(X), "euc") ^ 2
    if (!is.null(qnt)) {
        eps <- as.numeric(quantile(dat.dis, qnt))
    }
    if (is.null(qnt)) {
        eps <- min(dat.dis[dat.dis != 0])
    }
    kernal <- exp(-1 * dat.dis / (eps))
    K1 <- dist2full(kernal)
    diag(K1) <- 0
    D = matrix(0, ncol = ncol(K1), nrow = ncol(K1))
    tmpe <- apply(K1, 1, sum)
    tmpe[tmpe > 0] <- 1 / sqrt(tmpe[tmpe > 0])
    tmpe[tmpe < 0] <- 0
    diag(D) <- tmpe
    L <- D %*% K1 %*% D
    X <- svd(L)$u
    Y <- X / sqrt(apply(X ^ 2, 1, sum))
}