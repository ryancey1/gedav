library(ggplot2)

# Number 1 ---------------------------------------------------------------- 
#Go to class website under Course Documents > Data Sets and download the SLE B
#cell data set (from Garaud et al).
system(command = "7z -y x sle_b_cell.7z")

# Number 2 ----------------------------------------------------------------
# Unzip the text file, and read into R (Hint: using the read.table()
# function with a “header=T” argument and “row.names=1” argument is one method
# to do this).
data = read.table(file = "sle_b_cell.txt", header = TRUE, row.names = 1)

# Number 3 ----------------------------------------------------------------
# Look at the dimensions of the data. There should be 26 samples. If you
# have 27 samples, you still have the row names in the first data column, so
# retry 2 to set the row names to these.
dim(data)

# Number 4 ----------------------------------------------------------------
# Print the sample names to screen.
colnames(data)

# Number 5 ----------------------------------------------------------------
# Plot the second SLE patient sample versus the first normal control samples in
# an xy scatter plot. Remember that the first argument is the x vector. Label
# the x and y-axes as 'Normal' and 'SLE', respectively. Title the plot, 'SLE B
# cell sample vs. Normal B cell sample – all probesets'. Add grey grid lines
# with the function grid().

# vectorize the desired columns
df_5 = data.frame(SLE = data$sle.2, Control = data$control.1, row.names = rownames(data))
# open device
png(filename = "plots/5_allprobesets.png")
# plot data
plot(x = df_5[,"Control"],
     y = df_5[,"SLE"],
     xlab = "Normal",
     ylab = "SLE",
     main = "SLE B cell sample vs. Normal B cell sample – all probesets")
grid()
# close device
dev.off()


# Number 6 ----------------------------------------------------------------
# Now do the same plot but pick only the first 20 probesets. Use the pch=15
# argument to change the shape and color the points blue with the col argument.
df_first_20 = data.frame(SLE = data[1:20, "sle.2"], Control = data[1:20, "control.1"], 
                  row.names = rownames(data)[1:20])
png(filename = "plots/6_20probesets.png")
plot(x = df_6[,"Control"], 
     y = df_6[,"SLE"],
     xlab = 'Control',
     ylab = 'SLE',
     main = 'SLE B cell sample vs. Normal B cell sample – first 20 probesets',
     pch = 15,
     col = 'blue')
grid()
dev.off()


# Number 7 ----------------------------------------------------------------
# Now plot the following gene in a gene profile plot, IGLJ3 (immunoglobulin
# lambda joining 3), which is probeset ID 211881_x_at. This type of plot has the
# sample indices across the x-axis and the intensities on the y-axis, so you can
# see a profile of the gene across experiments or arrays. First plot the ranges
# using the type=”n” argument and the plot() function, then add the genes with
# the lines() function call. Add grid lines.
IGLJ3_x = as.numeric(data["211881_x_at", ])
png(filename = "plots/7_geneprofile.png")
plot(range(1:26),
     range(IGLJ3_x),
     type = "n",
     xlab = "Sample Number",
     ylab="Intensity",
     main = "Gene profile plot of IGLJ3 (211881_x_at)")
lines(IGLJ3_x, type = "l", col = "red", lwd = 2)
points(IGLJ3_x, type = "p", pch = 16, col = c(rep("blue", 17), rep("green", 9)))
grid()
legend(x = "topright", legend = c("SLE", "Control"), pch = 16, col = c("blue", "green"), bg = "white")
dev.off()


# 8.) Finally, another way to visualize a gene profile across conditions is to
# graph a boxplot with a single distribution box per condition. To do this, we
# need to create a factor vector that indicates the disease or normal condition,
# then use this vector with the expression vector for IGLJ3 in the boxplot
# function to create the graph.

# 8. ggplot ------------------------------------------------------------------
df <- data.frame(Sample.Type = c(rep("SLE", 17), rep("Control", 9)), IGLJ3.Intensity = IGLJ3_x)
ggplot(data = df, aes(x = Sample.Type, y = IGLJ3.Intensity, fill = Sample.Type)) +
    ylab("IGLJ3 Intensity") + 
    xlab("Sample Type") +
    ggtitle("Intensity distribution per condition") +
    geom_boxplot()
ggsave(filename = "plots/8_ggplot_boxplot.png", width = 5, height = 5)

# 8. Base plotting -----------------------------------------------------------
png(filename = "plots/8_boxplot_base.png")
boxplot(df$IGLJ3.Intensity ~ df$Sample.Type, 
        col=c("red", "blue"),
        ylab = "IGLJ3 Intensity", 
        xlab = "Sample Type",
        main = "Intensity Distribution by Sample Type")
legend(x = "topleft", 
       legend = c("Control", "SLE"), 
       col = c("red", "blue"),
       pch = 15)
dev.off()
