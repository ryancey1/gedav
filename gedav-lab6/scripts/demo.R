# Number1 -----------------------------------------------------------------
# Download the GEO Brain Aging study from the class website. Also obtain the
# annotation file for this data frame.
dir("./data")

# Number2 -----------------------------------------------------------------
# Load into R, using read.table() function and the header=T/row.names=1
# arguments for each data file.
brain.dat <- read.table(file = "data/agingStudy11FCortexAffy.txt",
                        header = TRUE,
                        row.names = 1)
brain.ann <- read.table(file = "data/agingStudy1FCortexAffyAnn.txt",
                        header = TRUE,
                        row.names = 1)
# Number3 -----------------------------------------------------------------
# Prepare 2 separate vectors for comparison. The first is a comparison between
# male and female patients. The current data frame can be left alone for this,
# since the males and females are all grouped together. The second vector is
# comparison between patients >= 50 years of age and those < 50 years of age.
# 
# To do this, you must use the annotation file and logical operators to isolate
# the correct arrays/samples.
male <- brain.ann$Gender == "M" # males
over50 <- brain.ann$Age >= 50 # equal to or older than 50

# Number4 -----------------------------------------------------------------
# Run the t.test function from the notes using the first gene vector below for
# the gender comparison. Then use the second gene vector below for the age
# comparison. Using these p-values, use either p.adjust in the base library or
# mt.rawp2adjp in the multtest library to adjust the values for multiple
# corrections with the Holm's method.
library("multtest")
source("scripts/t.test.all.genes.R")

# gene of interest vectors (gender and age, respectively)
g.g <- c(1394, 1474, 1917, 2099, 2367, 2428, 2625, 3168, 3181, 3641, 
         3832, 4526, 4731, 4863, 6062, 6356, 6684, 6787, 6900, 7223, 
         7244, 7299, 8086, 8652, 8959, 9073, 9145, 9389, 10219, 11238, 
         11669, 11674, 11793)
g.a <- c(25, 302, 1847, 2324, 246, 2757, 3222, 3675, 4429, 4430, 
         4912, 5640, 5835, 5856, 6803, 7229, 7833, 8133, 8579, 8822, 
         8994, 10101, 11433, 12039, 12353, 12404, 12442, 67, 88, 100)

# raw p-values
rawp_gend <- apply(brain.dat[g.g, ], 1, t.test.all.genes, male, !male)
rawp_age <- apply(brain.dat[g.a, ], 1, t.test.all.genes, over50, !over50)

# adjusted p-values
mt.gend <- mt.rawp2adjp(rawp_gend, proc = "Holm")
mt.age <- mt.rawp2adjp(rawp_age, proc = "Holm")

# Number5 ----------------------------------------------------------------- 
# Sort the adjusted p-values and non-adjusted p-values and plot them vs. the
# x-axis of numbers for each comparison data set. Make sure that the two lines
# are different colors. Also make sure that the p-values are sorted before
# plotting.
### GENDER ###
gend.holm <- data.frame(mt.gend$adjp[order(mt.gend$index), ])

mt.plot(
    gend.holm,
    plottype = "pvsr",
    proc = colnames(gend.holm),
    leg = c(0, 0.099),
    lty = c(1:2),
    lwd = 1.5,
    col = c("red", "blue"),
    main = "Type I error rate vs. # of rejected hypotheses:\nHolm adjustment (GENDER)"
)

### AGE ###
age.holm <- mt.age$adjp[order(mt.age$index), ]
mt.plot(
    age.holm,
    plottype = "pvsr",
    proc = colnames(age.holm),
    leg = c(0.05, 1.03),
    lty = c(1:2),
    lwd = 1.5,
    col = c("red", "blue"),
    main = "Type I error rate vs. # of rejected hypotheses:\nHolm adjustment (AGE)"
)


# Number6 ----------------------------------------------------------------- 
# Repeat #4 and #5 with the Bonferroni method.
mt.gend <- mt.rawp2adjp(rawp_gend, proc = "Bonferroni")
mt.age <- mt.rawp2adjp(rawp_age, proc = "Bonferroni")

gend.bonf <- mt.gend$adjp[order(mt.gend$index), ]
mt.plot(
    gend.bonf,
    plottype = "pvsr",
    proc = colnames(gend.bonf),
    leg = c(0, 0.34),
    lty = c(1,3),
    lwd = 1.5,
    col = c("red", "blue"),
    main = "Type I error rate vs. # of rejected hypotheses:\nBonferroni adjustment"
)

age.bonf <- mt.age$adjp[order(mt.age$index), ]
mt.plot(
    age.bonf,
    plottype = "pvsr",
    proc = colnames(age.bonf),
    leg = c(0.05, 1.03),
    lty = c(1,3),
    lwd = 1.5,
    col = c("red", "blue"),
    main = "Type I error rate vs. # of rejected hypotheses:\nBonferroni & Holm adjustments"
)

# Clean out environment
rm(list = ls())

# Number7 ----------------------------------------------------------------- 
# Read in the log2 normalized fragments per kb per million mapped reads (FPKM)
# data matrix and annotation files. This is RNA-sequencing data that has
# normalized read counts on a similar scale to microarray intensities.
tcga.dat <- read.table(file = "data/tcga_brca_fpkm.txt",
                       header = TRUE,
                       row.names = 1,
                       sep = "\t")

tcga.ann <- read.table(file = "data/tcga_brca_fpkm_sam.txt",
                       header = TRUE,
                       row.names = 1,
                       sep = "\t")

# Number8 ----------------------------------------------------------------- 
# Use grep to subset the data matrix only by gene ‘GATA3’ and make sure to cast
# this vector to numeric.
gata3 <- grep("GATA3", rownames(tcga.dat))
tcga.dat.gata3 <- as.numeric(tcga.dat[gata3, ])

# Number9 ----------------------------------------------------------------- 
# Create a binary (1/0) vector for the patients where the upper 25% expression
# of GATA3 is coded as 1 and all other patients are coded as 0. Call this new
# variable ‘group’.
# We're only interested in the 75th percentile (upper 25%)
upper25th <- quantile(tcga.dat.gata3, prob = 0.75)

# Vectorize necessary variables to isolate them
group <- as.numeric(tcga.dat.gata3 >= upper25th)
status <- as.numeric(tcga.ann$vital_status == "DECEASED")
time <- as.numeric(tcga.ann$months_to_event)

# Number10 ---------------------------------------------------------------- 
# Create a data matrix with the ‘group’ variable you created in #9 and the
# remaining variables in the annotation file.
ann.dm <- data.frame(time, status, group)


# Number11 ---------------------------------------------------------------- 
# Run a Kaplan-Meier (KM) analysis to determine if a difference in survival
# experience exists between the two GATA3 expression groups using the survdiff
# function. Extract the p-value from the chi squared test output.
library(survival)
surv <- with(ann.dm, Surv(time, status))
(sdf <- survdiff(surv ~ group, data = ann.dm))
pval_km <- signif(pchisq(sdf$chisq, length(sdf$n) - 1, lower.tail = FALSE), 4)

# Number12 ---------------------------------------------------------------- 
# Now run a Cox proportion hazard (PH) regression model on just the grouping
# variable (i.e. no other covariates) and extract both the p-value and hazard
# ratio from the output.
cph <- coxph(surv ~ group, data = ann.dm)
summary(cph)
pval_cox <- signif(summary(cph)$coefficients[5], 4)
HR_cox <- signif(summary(cph)$coefficients[2], 4)


# Number13 ---------------------------------------------------------------- 
# Run the survfit() function only on the grouping variable (i.e. no other
# covariates) and plot the KM curves, being sure to label the two groups with a
# legend, two different colors for each line, and provide the KM p-value, Cox PH
# p-value, Cox PH hazard ratio, and sample sizes all in each of the two groups
# all on the plot.
sf <- survfit(surv ~ group, data = ann.dm)

pl <- ggsurvplot(sf, risk.table = TRUE)

addin <- c(
            paste0("KM p-value = ", pval_km),
            paste0("\nCoxPH p-value = ", pval_cox),
            paste0("\n\nCoxPH HR = ", HR_cox)
)

pl$plot <- pl$plot + 
    ggplot2::annotate(
        "text",
        x = 60,
        y = 0.95,
        label = addin,
        hjust = 0,
        vjust = 1,
        size = 4
    )

pl

# plot(
#     sf,
#     lty = 1:2,
#     lwd = 1.5,
#     col = c("blue", "red"),
#     xlab = "Months",
#     ylab = "Survival Probability",
#     main = "Kaplan-Meier Curve\nTCGA BrCa data grouped by GATA3"
# )
# legend(
#     x = "topright",
#     legend = c(
#         paste0("GATA3 low (", sf$n[[1]], ")"),
#         paste0("GATA3 high (", sf$n[[2]], ")"),
#         "",
#         paste0("KM p-value = ", pval_km),
#         paste0("CoxPH p-value = ", pval_cox),
#         paste0("CoxPH HR = ", HR_cox)
#     ),
#     col = c("blue", "red"),
#     lwd = 1.5,
#     lty = c(1, 2, 0, 0, 0, 0),
#     inset = 0.02,
#     cex = 0.8
# )

# Number14 ---------------------------------------------------------------- 
# Does this result agree with the Mehra et al, study result?

