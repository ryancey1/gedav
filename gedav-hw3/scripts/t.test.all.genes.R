t.test.all.genes <- function(x,s1,s2) {
      x1 <- x[s1]
      x2 <- x[s2]
      x1 <- as.numeric(x1)
      x2 <- as.numeric(x2)
      t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
      out <- as.numeric(t.out$p.value)
      return(out)
  }