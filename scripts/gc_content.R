# 1. Install and load necessary package (only run install.packages once)
# install.packages("fitdistrplus")
library(fitdistrplus)

# 2. Read the data
# Assuming your file is named 'human_GC_SCFR_length.tsv' and is space-separated.
# The column names will be '13_seq_len' and '6_pct_gc' as per the header.
data <- read.table("human_GC_SCFR_length.tsv", 
                   header = TRUE, 
                   sep = " ")

# Rename columns for simpler access in the script (optional but good practice)
colnames(data) <- c("SCFR_Length", "GC_Content")

# --- Part 1: Plotting the Distribution of GC Content vs. Length ---

# 3. Create a scatter plot of GC Content vs. SCFR Length
# This shows the overall relationship between the two variables.
# Save plot to a file
png("GC_vs_Length_Scatter.png", width = 800, height = 600)
plot(data$GC_Content, data$SCFR_Length,
     xlab = "SCFR Length (13_seq_len)",
     ylab = "Percentage GC Content (6_pct_gc)",
     main = "Percentage GC Content vs. SCFR Length",
     pch = 19, # Solid dots
     col = "blue")
dev.off()
cat("Scatter plot saved as 'GC_vs_Length_Scatter.png'\n")

# 4. Create a histogram and density plot of GC Content (for visual distribution assessment)
# GC content is a proportion (0 to 1), so it's a continuous variable.
png("GC_Content_Distribution.png", width = 800, height = 600)
hist(data$GC_Content, 
     prob = TRUE, # Plot probabilities instead of frequencies
     main = "Distribution of Percentage GC Content",
     xlab = "Percentage GC Content (6_pct_gc)",
     col = "lightblue",
     border = "black")
lines(density(data$GC_Content), 
      col = "red", 
      lwd = 2)
legend("topright", 
       legend = c("Histogram", "Kernel Density Estimate"), 
       col = c("lightblue", "red"), 
       lty = c(NA, 1), 
       lwd = c(NA, 2), 
       pch = c(22, NA))
dev.off()
cat("Distribution plot saved as 'GC_Content_Distribution.png'\n")

# --- Part 2: Distribution Fitting using fitdistrplus ---

# 5. Fit common distributions to the GC_Content data
# GC content is a proportion (0-1), making the Beta distribution a biologically relevant candidate.
# We will also check Normal and Gamma distributions for comparison.

gc_data <- data$GC_Content

# a) Fit Normal distribution
fit_normal <- fitdist(gc_data, "norm")
cat("\n--- Normal Distribution Fit ---\n")
print(summary(fit_normal))

# b) Fit Gamma distribution (requires all data > 0, which is typically true for GC content)
fit_gamma <- fitdist(gc_data, "gamma")
cat("\n--- Gamma Distribution Fit ---\n")
print(summary(fit_gamma))

# c) Fit Beta distribution (most suitable for proportions/percentages)
# Need to supply starting parameters for the optimization.
# A common starting point is a Beta(1, 1) - i.e., Uniform distribution.
fit_beta <- fitdist(gc_data, "beta", start = list(shape1 = 1, shape2 = 1))
cat("\n--- Beta Distribution Fit ---\n")
print(summary(fit_beta))

# 6. Compare the goodness-of-fit statistics (AIC, BIC, Kolmogorov-Smirnov, etc.)
gof_stats <- gofstat(list(fit_normal, fit_gamma, fit_beta), 
                     fitnames = c("Normal", "Gamma", "Beta"))
cat("\n--- Goodness-of-Fit (GOF) Statistics ---\n")
print(gof_stats)

# 7. Visualize the fits (optional but recommended)
png("GC_Content_Fit_Pdfs.png", width = 800, height = 600)
par(mfrow = c(2, 2)) # Set up a 2x2 plotting area
plot.legend <- c("Normal", "Gamma", "Beta")
denscomp(list(fit_normal, fit_gamma, fit_beta), legendtext = plot.legend)
qqcomp(list(fit_normal, fit_gamma, fit_beta), legendtext = plot.legend)
cdfcomp(list(fit_normal, fit_gamma, fit_beta), legendtext = plot.legend)
ppcomp(list(fit_normal, fit_gamma, fit_beta), legendtext = plot.legend)
dev.off()
cat("Distribution comparison plots saved as 'GC_Content_Fit_Pdfs.png'\n")
