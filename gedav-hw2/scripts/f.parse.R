# qRT-PCR file formatting and calculation of fold changes (continues on next slide)

f.parse <- function(path = pa,
                    file = fi,
                    out = out.fi) {
    d <- read.table(
        paste(path, file, sep = ""),
        skip = 11,
        sep = ",",
        header = T
    )
    u <- as.character(unique(d$Name))
    u <- u[u != ""]
    u <- u[!is.na(u)]
    
    ref <- unique(as.character(d$Name[d$Type == "Reference"]))
    u <- unique(c(ref, u))
    hg <- c("B-ACTIN", "GAPDH", "18S")
    hg <- toupper(hg)
    p <- unique(toupper(as.character(d$Name.1)))
    p <- sort(setdiff(p, c("", hg)))
    
    mat <- matrix(0, nrow = length(u), ncol = length(p))
    dimnames(mat) <- list(u, p)
    for (i in 1:length(u)) {
        print(paste(i, ": ", u[i], sep = ""))
        tmp <- d[d$Name %in% u[i], c(1:3, 6, 9)]
        g <- toupper(unique(as.character(tmp$Name.1)))
        g <- sort(setdiff(g, c("", hg)))
        
        for (j in 1:length(g)) {
            v <- tmp[toupper(as.character(tmp$Name.1)) %in% g[j], 5]
            v <- v[v != 999]
            v <-
                v[((v / mean(v)) < 1.5) & ((v / mean(v)) > 0.67)]	#gene j vector
            
            hv3 <- NULL
            for (k in 1:length(hg)) {
                #housekeeping gene vector (each filtered by reps)
                hv <-
                    tmp[toupper(as.character(tmp$Name.1)) %in% hg[k], 5]
                hv <- hv[hv != 999]
                hv3 <-
                    c(hv3, hv[((hv / mean(hv)) < 1.5) & ((hv / mean(hv)) > 0.67)])
            }
            
            
            # qRT-PCR file formatting and calculation of fold changes (cont)
            
            
            sv <-
                mean(as.numeric(v)) - mean(as.numeric(hv3))	#scaled value for gene j
            
            if (i == 1) {
                #reference sample only
                mat[u[i], g[j]] <- sv
                next
            }
            
            mat[u[i], g[j]] <- sv - mat[u[1], g[j]]
        }
    }
    
    mat[1, ][!is.na(mat[1, ])] <- 0
    fc <- 2 ^ (-1 * mat)
    write.table(
        t(c("Subject", dimnames(mat)[[2]])),
        paste(path, out, sep = ""),
        quote = F,
        sep = "\t",
        col.names = F,
        row.names = F
    )
    write.table(
        round(fc, 3),
        paste(path, out, sep = ""),
        quote = F,
        sep = "\t",
        append = T,
        col.names = F
    )
}

pa <- "data/"
fi <- "qRT-PCR.csv"
out.fi <- "out.fi"

f.parse(pa, fi, out.fi)
