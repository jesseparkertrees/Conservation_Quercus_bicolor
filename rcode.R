
# Main Article Title: “Genomically informed seed orchard design for trailing-edge tree populations: a perspective from Quercus bicolor conservation”

# Jesse B. Parker (parkerjesseb@gmail.com)1 
# Sean Hoban (shoban@mortonarb.org)2,3
# Laura M. Thompson (lthompson@usgs.gov)1,4
# Scott E Schlarbaum (tenntip@utk.edu)1
# 1	School of Natural Resources, University of Tennessee, Knoxville TN 37996, USA
# 2	The Center for Tree Science, The Morton Arboretum, Lisle IL 60532, USA
# 3	Committee on Evolutionary Biology, University of Chicago, Chicago IL 60637, USA
# 4	U.S. Geological Survey, National Climate Adaptation Science Center, Reston VA 20192, USA
# Author for correspondence: parkerjesseb@gmail.com
# Keywords: Ex situ, germplasm, central-marginal hypothesis, introgression, genetic diversity 

# Disclaimer: Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.


#################################################################################################################
##################### I. LOAD PACKAGES #############################################################################
#################################################################################################################
library(dartR)
library(vcfR)
library(dplyr)
library(SNPfiltR)
library(ggplot2)
# Load additional packages for parallel processing
library(adegenet)
library(parallel)
library(doParallel)
library(abind)


#################################################################################################################
##################### II. SET UP DATA ###############################################################################
#################################################################################################################
# Set working directory to repository
setwd("~/Conservation_Quercus_bicolor")
# Load filtered .vcf (this has been filtered to only include biallelic SNPs with a call rate of 1)
vc<-read.vcfR("./qubi.vcf")
# Convert to genlight 
glv<-vcfR2genlight(vc)
glv<-gl.compliance.check(glv)
glv<-glv[2:135] # Filter out sample_1 (Q. bicolor x Q. stellata hybrid) 
# Load population data, admixture proportions, and other data aquired from Parker et al. (2025)
pop<-read.csv("./Qbicolor_data.csv")
pop<-pop[2:135,] # Filter out sample_1 (Q. bicolor x Q. stellata hybrid) and reference samples
names<-glv@ind.names
pop$ID == names # Check that names are the same
# Set population levels
glv@pop<-pop$pop %>% as.factor()
pop2<-pop$pop2 %>% as.data.frame()
colnames(pop2)<-"pop2"
glv@strata<-pop2

# Set up another genlight with less stringent filter (0.95)
# Load the vcf that has been filtered to a call rate of 0.95
vc2<-read.vcfR("qubi_95.vcf")
glv2<-vcfR2genlight(vc2)
glv2<-glv2[2:135] # Filter out sample_1 (Q. bicolor x Q. stellata hybrid) 
glv2<-gl.compliance.check(glv2)
names<-glv2@ind.names
pop$ID == names # Check that names are the same
glv2@pop<-pop$pop %>% as.factor()
glv2@strata<-pop2


#################################################################################################################
##################### III. ASSESS HYBRIDIZATION AND ESTABLISH THRESHOLD ##############################################
#################################################################################################################
# Filter to only the Tennessee individuals
pop_TN<-pop[1:69,]

# Plot hybrid count curve with both the sNMF and STRUCTURE ancestry proportions
thresh<-seq(0,1,0.01)
prop<-pop_TN$sNMF_V4 # sNMF ancestry component
counts <- sapply(thresh, function(t) sum(prop > t))
df<-data.frame(threshold = thresh, count = counts)
prop2<-pop_TN$X2spec_struc_swo_prop # STRUCTURE Q. bicolor ancestry component
counts2 <- sapply(thresh, function(t) sum(prop2 > t))
df2<-data.frame(threshold = thresh, count = counts2)
# Combine into one dataframe
df$method <- "sNMF ancestry proportion"
df2$method <- "STRUCTURE ancestry proportion"
combined_df <- rbind(df, df2)
my_colors <- c(
  "sNMF ancestry proportion" = "lightblue",
  "STRUCTURE ancestry proportion" = "olivedrab"
)
plot <- ggplot(combined_df, aes(x = threshold, y = count, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 46, linetype = "dashed", color = "black", size = 1) +
  annotate("text",
           x = max(combined_df$threshold) * 0.35,
           y = 46,
           label = "46 individual subset", vjust = -0.5, hjust = -0.1, color = "black") +
  scale_color_manual(values = my_colors, name = "") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
  labs(
    title = "",
    x = "Hybrid Threshold",
    y = "Genotype Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
plot

# Filter out hybrids below 0.84 snmf (equivalent to 0.9 STRUCTURE)
filt1<-pop[pop$X2spec_snmf_swo_prop>0.84,]
glv_filt1<-gl.keep.ind(glv, filt1$ID)
glv_filt1<-gl.filter.monomorphs(glv_filt1)

glv_filt95<-gl.keep.ind(glv2, filt1$ID)
glv_filt95<-gl.filter.monomorphs(glv_filt95)
edge_filt95<-glv_filt95[glv_filt95@strata$pop2=="EDGE"]
edge_filt95<-gl.filter.monomorphs(edge_filt95)
core_filt95<-glv_filt95[glv_filt95@strata$pop2=="CORE"]
core_filt95<-gl.filter.monomorphs(core_filt95)


edge_filt1<-glv_filt1[glv_filt1@strata$pop2=="EDGE"]
core_filt1<-glv_filt1[glv_filt1@strata$pop2=="CORE"]
edge_filt1<-gl.filter.monomorphs(edge_filt1)
core_filt1<-gl.filter.monomorphs(core_filt1)


#################################################################################################################
##################### IV. ASSESS HETEROZYGOSITY #####################################################################
#################################################################################################################
het<-gl.report.heterozygosity(glv2, method="ind") %>% as.data.frame() #use glv2 (less stringent call rate, more SNPs included)
het$rangepos<-0
het$rangepos[0:70]<-"TN"
het$rangepos[71:134]<-"CORE"
my_colors <- c(
  "CORE" = "lightblue",
  "TN" = "olivedrab"
)
# Histogram plot
plot <- ggplot(het, aes(x = Ho, fill = rangepos)) +
  geom_histogram(
    alpha = 0.9,
    position = "dodge",
    bins = 15,
    color = "black"
  ) +
  scale_fill_manual(values = my_colors, name = "Range position") +
  labs(
    title = "",
    x = "Observed Heterozygosity (Ho)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

plot # Long tail on right in TN due to hybridization

# Plot without hybrids
het<-gl.report.heterozygosity(glv_filt95, method="ind") %>% as.data.frame()
het$rangepos<-0
het$rangepos[0:46]<-"TN"
het$rangepos[47:111]<-"CORE"
# Histogram plot
plot <- ggplot(het, aes(x = Ho, fill = rangepos)) +
  geom_histogram(
    alpha = 0.9,
    position = "dodge",
    bins = 11,
    color = "black"
  ) +
  scale_fill_manual(values = my_colors, name = "Range position") +
  labs(
    title = "",
    x = "Observed Heterozygosity (Ho)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

plot

library("coin")
het$rangepos <- as.factor(het$rangepos)
wilcox_test(het$Ho ~ het$rangepos, distribution = "approximate")
library(car)
leveneTest(Ho ~ rangepos, data = het)

#################################################################################################################
##################### V. ALLELE RAREFACTION without hybrid filters #################################################
#################################################################################################################
# Divide into TN and core samples
glv_edge<-glv[glv@strata$pop2=="EDGE"]
glv_edge<-gl.filter.monomorphs(glv_edge)
glv_core<-glv[glv@strata$pop2=="CORE"]
glv_core<-gl.filter.monomorphs(glv_core)

################################# 1. Run for all pooled samples #################################
# Set parameters
num_reps <- 1000
num_cores <- detectCores() - 2  # Use available cores minus 2 for other tasks
# Set sample size
sample_sizes <- seq(10, 134, by = 1)
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(glv), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(glv)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)

summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df1 <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 2. Run for EDGE #################################
sample_sizes <- seq(10, 69, by = 1)
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv_edge), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(glv_edge), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(glv_edge)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)
summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df_edge <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 3. Run for CORE #################################
# Set parameters
sample_sizes <- seq(10, 65, by = 1)  # Sample sizes ranging from 5 to 110
allele_freqs <- colMeans(as.matrix(glv_core), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present

# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(glv_core), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(glv_core)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)  # Count alleles present in sample
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)

summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df_core <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 4. Plot all results ##############################
summary_df1$Population <- "All Samples"
summary_df_edge$Population <- "Tennessee Samples"
summary_df_core$Population <- "Range Core Samples"
# Combine both results into one data frame
combined_summary_df <- rbind(summary_df1, summary_df_edge, summary_df_core)
# Calculate thresholds for 70% and 90% of total alleles
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
threshold_70 <- total_alleles * 0.7
threshold_90 <- total_alleles * 0.9

my_colors <- c(
  "Tennessee Samples" = "#35B779FF",
  "Range Core Samples" = "#31688EFF",
  "All Samples" = "#FDE725FF"
)
plot <- ggplot(combined_summary_df, aes(x = Sample_Size, y = Mean_Alleles, color = Population, group = Population)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = threshold_70, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = threshold_90, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = max(sample_sizes) - 5, y = threshold_70, label = "70% alleles", vjust = -0.5, color = "black", hjust=-1) +
  annotate("text", x = max(sample_sizes) - 5, y = threshold_90, label = "90% alleles", vjust = -0.5, color = "black", hjust=-1) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(breaks = seq(5, 135, by = 10), labels = seq(5, 135, by = 10)) +
  labs(title = "", x = "Sample Size", y = "Mean Number of Alleles") +
  theme_minimal()
plot


#################################################################################################################
##################### VI. ALLELE RAREFACTION with hybrid filters ################################################
#################################################################################################################
# Set parameters
num_reps <- 1000
num_cores <- detectCores() - 2  # Use available cores minus 2 for other tasks

################################# 1. Run for all pooled samples #################################
# Set sample size
sample_sizes <- seq(10, length(glv_filt1@ind.names), by = 1)  
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv_filt1), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(glv_filt1), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(glv_filt1)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)  # Count alleles present in sample
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)
summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df1_filt1 <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 2. Run for EDGE #################################
# Set parameters
sample_sizes <- seq(10, length(edge_filt1@ind.names), by = 1)
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(edge_filt1), na.rm = TRUE)
total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(edge_filt1), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(edge_filt1)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)  # Count alleles present in sample
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)

summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df_edge_filt1 <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 3. Run for CORE #################################
# Set parameters
sample_sizes <- seq(10, length(core_filt1@ind.names), by = 1)  # Sample sizes ranging from 5 to 110

# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(core_filt1), na.rm = TRUE)

total_alleles <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present

# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Store results
results <- array(dim = c(length(sample_sizes), num_reps))

temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(core_filt1), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(core_filt1)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)  # Count alleles present in sample
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)

summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df_core_filt1 <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 4. Run for CoreHunter subset (see section VIII on how to generate) #################################
ch<-read.csv("./ch33.csv")
glv_ch<-gl.keep.ind(glv_filt1, ch$ID)
glv_ch<-gl.filter.monomorphs(glv_ch)
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv_ch), na.rm = TRUE)
total_alleles_ch <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
# Set parameters
sample_sizes <- seq(10, length(glv_ch@ind.names), by = 1)  # Sample sizes ranging from 5 to 110
# Parallel processing setup
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Store results
results <- array(dim = c(length(sample_sizes), num_reps))
temp <- foreach(nrep = 1:num_reps, .packages = 'adegenet') %dopar% {
  res <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    t <- sample_sizes[i]
    sampled_inds <- sample(1:nInd(glv_ch), t, replace = FALSE)
    sampled_alleles <- colSums(as.matrix(glv_ch)[sampled_inds, ], na.rm = TRUE)
    res[i] <- sum(sampled_alleles > 0, na.rm = TRUE)  # Count alleles present in sample
  }
  res
}
results[,] <- do.call(cbind, temp)  # Combine results into an array
stopCluster(cl)

summary_results <- apply(results, 1, mean, na.rm = TRUE)
summary_df_ch <- data.frame(Sample_Size = sample_sizes, Mean_Alleles = summary_results, Total_Alleles = total_alleles)

################################# 5. Determine alleles for grafted individuals #################################
# Extract grafted trees
graf<-pop[pop$Grafted.=="grafted",]
dim(graf)
glv_graf<-gl.keep.ind(glv_filt1, graf$ID)

glv_graf<-gl.keep.ind(glv_filt1, graf$ID)
glv_graf<-gl.filter.monomorphs(glv_graf)
# Extract allele frequencies and population information
allele_freqs <- colMeans(as.matrix(glv_graf), na.rm = TRUE)
total_alleles_graf <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present

################################# 6. Combine data, apply fitted curves, plot ###################################
summary_df1_filt1$Population <- "All Samples"
summary_df_edge_filt1$Population <- "Tennessee Samples"
summary_df_core_filt1$Population <- "Range Core Samples"
summary_df_ch$Population <- "CoreHunter Subset"

# Combine both results into one data frame
combined_summary_df_filt2 <- rbind(summary_df1_filt1, summary_df_core_filt1,summary_df_edge_filt1, summary_df_ch)
# Make sure factor levels are consistent
combined_summary_df_filt2$Population <- factor(
  combined_summary_df_filt2$Population,
  levels = c("All Samples", "Range Core Samples","Tennessee Samples","CoreHunter Subset")
)

# Add fitted curves (Michaelis-Menten model)
library(minpack.lm)
fit_mm_model <- function(df) {
  tryCatch({
    nlsLM(
      Mean_Alleles ~ (a * Sample_Size) / (b + Sample_Size),
      data = df,
      start = list(a = max(df$Mean_Alleles), b = 10)
    )
  }, error = function(e) {
    message("Model failed for population: ", unique(df$Population), "\nError: ", e$message)
    return(NA)
  })
}

# Fit to all populations
fits <- combined_summary_df_filt2 %>%
  group_by(Population) %>%
  group_split() %>%
  setNames(unique(combined_summary_df_filt2$Population)) %>%
  lapply(fit_mm_model)

prediction_df <- combined_summary_df_filt2 %>%
  group_by(Population) %>%
  group_split() %>%
  lapply(function(df) {
    pop <- unique(df$Population)
    model <- fits[[pop]]
    newdata <- data.frame(Sample_Size = 10:111)
    if (inherits(model, "nls") || inherits(model, "nlsLM")) {
      preds <- predict(model, newdata = newdata)
      data.frame(Population = pop, Sample_Size = 10:111, Predicted = preds)
    } else {
      data.frame(Population = pop, Sample_Size = 10:111, Predicted = NA)
    }
  }) %>%
  bind_rows()

prediction_df$Population <- factor(
  prediction_df$Population,
  levels = c("All Samples","Range Core Samples","Tennessee Samples","CoreHunter Subset")
)

# Extract parameter 'a' (asymptote) for each population fit
asymptotes <- sapply(fits, function(model) {
  if (inherits(model, "nls") || inherits(model, "nlsLM")) {
    coef(model)["a"]
  } else {
    NA
  }
})

# Combine and plot
summary_df1_filt1$Group <- "All Samples"
summary_df_edge_filt1$Group <- "Tennessee Samples"
summary_df_core_filt1$Group <- "Range Core Samples"
summary_df_ch$Group<- "Genomically Informed Selection"
# Run this for Population 2 (after running for Pop1)

combined_summary_df_filt2 <- rbind(summary_df1_filt1, summary_df_edge_filt1, summary_df_core_filt1, summary_df_ch)
# Make sure factor levels are consistent
combined_summary_df_filt2$Group <- factor(
  combined_summary_df_filt2$Group,
  levels = c("All Samples", "Range Core Samples","Tennessee Samples","Genomically Informed Selection")
)

# Calculate thresholds for 70% and 90% of total alleles
allele_freqs <- colMeans(as.matrix(glv_filt1), na.rm = TRUE)
total_alleles_filt1 <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
threshold_70_filt1 <- total_alleles_filt1 * 0.7
threshold_90_filt1 <- total_alleles_filt1 * 0.9
allele_freqs <- colMeans(as.matrix(glv_ch), na.rm = TRUE)
total_alleles_ch <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present
allele_freqs <- colMeans(as.matrix(edge_filt1), na.rm = TRUE)
edge_alleles_filt1 <- sum(allele_freqs > 0, na.rm = TRUE)  # Total unique alleles present

total_alleles_graf/thr_allsamp
5818/thr_allsamp
5818/edge_alleles_filt1

my_colors <- c(
  "Tennessee Samples" = "#35B779FF",
  "Range Core Samples" = "#31688EFF",
  "Genomically Informed Selection" = "#440154FF",
  "All Samples" = "#FDE725FF"
)
plot <- ggplot(combined_summary_df_filt2, aes(x = Sample_Size, y = Mean_Alleles, color = Group, group = Group)) +
  geom_line() +
  geom_point() +
  geom_line(data = prediction_df %>% filter(Population != "All Samples" & Population != "CoreHunter Subset"),
            aes(x = Sample_Size, y = Predicted, color = Population, group = Population),
            size = 0.8, alpha = 0.3) +
  geom_hline(yintercept = threshold_70_filt1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = threshold_90_filt1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = total_alleles_graf, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = total_alleles_filt1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = total_alleles_ch, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 75, y = threshold_70_filt1, label = "70% Alleles", vjust = -0.5, color = "black") +
  annotate("text", x = 75, y = threshold_90_filt1, label = "90% Alleles", vjust = -0.5, color = "black") +
  annotate("text", x = 75, y = total_alleles_graf, label = "Opportunistic Sample (72.6% Alleles)", vjust = -0.5, color = "black") +
  annotate("text", x = 75, y = total_alleles_ch, label = "Genomically Informed Selection (76.3% Alleles)", vjust = -0.5, color = "black") +
  annotate("text", x = 75, y = total_alleles_filt1, label = "100% Alleles", vjust = -0.5, color = "black") +
  scale_x_continuous(breaks = seq(5, 110, by = 10), labels = seq(5, 110, by = 10)) +
  scale_color_manual(values = my_colors) +
  labs(title = "", x = "Sample Size", y = "Mean Number of Alleles") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"))
plot
ggsave("allele_sampling_plot.png", plot, width = 10, height = 8, bg= "white")

#################################################################################################################
##################### VII. ASSESS RELATEDNESS AND USE COREHUNTER ################################################
#################################################################################################################
# Create network map of kinship estimates
glv_graf<-gl.keep.ind(glv_filt95, graf$ID)
glv_graf<-gl.filter.monomorphs(glv_graf)
grm_graf<-gl.grm(glv_graf)
grmnet_graf<-gl.grm.network(grm_graf, glv_graf, relatedness_factor=0.092) # Only pairs of individuals with relatedness ≥ 0.038 are connected.

grm_core<-gl.grm(core_filt95)
grmnet_core<-gl.grm.network(grm_core, core_filt1, relatedness_factor=0.092)

grm_edge<-gl.grm(edge_filt95)
grmnet_edge<-gl.grm.network(grm_edge, edge_filt1, relatedness_factor=0.092)

# Use CoreHunter on edge samples to select individuals
library(corehunter)
#drop individual with outlier heterozygosity (from section IV)
filt<-gl.drop.ind(edge_filt95, c("sample_6"))
filt<-gl.filter.monomorphs(filt)
grm_filt<-gl.grm(filt)
grmnet_filt<-gl.grm.network(grm_filt, filt, relatedness_factor=0.092)
# Strategically remove a few of the highly interrelated individuals
filt<-gl.drop.ind(filt, c("sample_7", "sample_50","sample_45","sample_28")) 
filt<-gl.filter.monomorphs(filt)
grafmat<-as.matrix(filt)
chg<-genotypes(grafmat, format='biparental')
sample_sizes <- 20:40
# Lists to store results
corehunt_results <- vector("list", length(sample_sizes))
names(corehunt_results) <- paste0("n", sample_sizes)
d_fer_results <- vector("list", length(sample_sizes))
names(d_fer_results) <- paste0("n", sample_sizes)
set.seed(67)
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  # Run Core Hunter
  ch_results <- sampleCore(chg, size = n, time = 45)
  corehunt_results[[i]] <- ch_results
  
  # Subset population using selected IDs
  d_fer_results[[i]] <- pop[pop$ID %in% ch_results$sel, ]
}

grmnet_results <- vector("list", length(d_fer_results))
names(grmnet_results) <- names(d_fer_results)
for (i in seq_along(d_fer_results)) {
  
  d_fer <- d_fer_results[[i]]
  n <- nrow(d_fer)
  glv_dfer <- gl.keep.ind(edge_filt95, d_fer$ID)
  glv_dfer <- gl.filter.monomorphs(glv_dfer)
  grm_dfer <- gl.grm(glv_dfer)
  
  gl.grm.network(grm_dfer, glv_dfer, relatedness_factor = 0.092)
}

d_fer_results$n33$ID %in% ch$ID
gl_ch33 <- gl.keep.ind(edge_filt95, d_fer_results$n33$ID)
gl_ch33 <- gl.filter.monomorphs(gl_ch33)
grm_dfer <- gl.grm(gl_ch33)
gl.grm.network(grm_dfer,gl_ch33, relatedness_factor = 0.092)
ch33_pop<-pop[pop$ID %in% d_fer_results$n33$ID,]
length(ch33_pop$ID)
mean(ch33_pop$sNMF_V4)
par(mar=c(2,2,2,2))
hist(chp_33$sNMF_V4, breaks=20)
ch33_het<-het[het$ind.name %in% d_fer_results$n33$ID,]
hist(ch33_het$Ho, breaks=20)
ch33<-cbind(ch33_pop,ch33_het$Ho)
#write.csv(ch33,"ch33.csv")

#png("grm_network_ch33.png", width = 3000, height = 3000, res = 300)
#grmnet_ch33<-gl.grm.network(grm_dfer, glv_dfer, relatedness_factor=0.092, title=NULL, node.label.size = 4)
#dev.off()
