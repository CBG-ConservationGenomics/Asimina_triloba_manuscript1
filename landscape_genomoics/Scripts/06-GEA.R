#######################################################################
# Genotype-Environment Association (GEA) with RDA and LFMM
# https://thewanglab.github.io/algatr/articles/RDA_vignette.html
#######################################################################

library(algatr)
rda_packages()
library(dplyr)
library(vegan)

set.seed(1234)

## ------------------------------------------------
## 0. Checks and paths
## ------------------------------------------------

data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
in_dir <- file.path(data_dir, "outputs", "00-preprocessing")
out_dir <- file.path(data_dir, "outputs", "06-gea")
plot_dir <- file.path(out_dir, "plots")
results_dir <- file.path(out_dir, "results")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

load_if_missing <- function(obj, rds_name) {
  if (!exists(obj, inherits = TRUE)) {
    p <- file.path(in_dir, rds_name)
    if (!file.exists(p)) stop("Missing required input: ", obj, " (expected ", p, ")")
    assign(obj, readRDS(p), envir = .GlobalEnv)
  }
}

load_if_missing("coord", "coord.rds")
load_if_missing("env_scaled", "env_scaled.rds")
load_if_missing("str_dos", "str_dos.rds")

#######################################################################
## 01. Redundancy analysis (RDA)
#######################################################################

# Scale the genomic data using Hellinger transformation
gen_h <- decostand(str_dos, "hellinger")

# Run a simple RDA with no variable selection
mod_full <- rda_run(gen_h, env_scaled, model = "full")
writeLines(capture.output(mod_full$call, summary(mod_full), RsquareAdj(mod_full)),
           file.path(results_dir, "rda_full_summary.txt"))

# Run RDA with variable selection
mod_best <- rda_run(gen_h, env_scaled,
                    model = "best",
                    Pin = 0.05,
                    R2permutations = 1000,
                    R2scope = TRUE
)
writeLines(capture.output(mod_best$call, mod_best$anova, RsquareAdj(mod_best)),
           file.path(results_dir, "rda_best_summary.txt"))

# Partial RDA with geography as covariable
mod_pRDA_geo <- rda_run(gen_h, env_scaled, coord,
                        model = "full",
                        correctGEO = TRUE,
                        correctPC = FALSE
)
writeLines(capture.output(anova(mod_pRDA_geo), RsquareAdj(mod_pRDA_geo), summary(mod_pRDA_geo)),
           file.path(results_dir, "rda_partial_geo_summary.txt"))

# Partial RDA with genetic structure and geography as covariables
mod_pRDA_gs <- rda_run(gen_h, env_scaled, coord,
                       model = "full",
                       correctGEO = TRUE,
                       correctPC = TRUE,
                       nPC = 2
)
writeLines(capture.output(summary(mod_pRDA_gs)), file.path(results_dir, "rda_partial_gs_summary.txt"))

png(file.path(plot_dir, "rda_partial_gs_plot.png"), width = 1200, height = 1200, res = 200)
plot(mod_pRDA_gs)
dev.off()

# Variance partitioning
varpart <- rda_varpart(gen_h, env_scaled, coord,
                       Pin = 0.05, R2permutations = 1000,
                       R2scope = TRUE, nPC = 2
)
writeLines(capture.output(str(varpart)),
           file.path(results_dir, "rda_varpart_table.txt"))

# Partial RDA with variable selection (for outlier detection)
mod_pRDA <- rda_run(gen_h, env_scaled, model = "best", correctPC = TRUE, nPC = 2)
writeLines(capture.output(mod_pRDA$anova), file.path(results_dir, "rda_partial_best_anova.txt"))

# Identifying outliers using the Z-scores method
png(file.path(plot_dir, "rda_plot_axes.png"), width = 1600, height = 1200, res = 200)
rda_plot(mod_pRDA, axes = "all", binwidth = 20)
dev.off()

rda_sig_z <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)
writeLines(paste0("Number of RDA outlier SNPs (z): ", length(rda_sig_z$rda_snps)),
           file.path(results_dir, "rda_outlier_count.txt"))

# Identifying outliers using the p-value method
rda_sig_p <- rda_getoutliers(mod_best, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, plot = FALSE)
    
# How many outlier SNPs were detected?
length(rda_sig_p$rda_snps)  

# Extract SNP names; choices is number of axes
snp_names <- rownames(scores(mod_best, choices = 2, display = "species"))

# Identify outliers that have q-values < 0.1
q_sig <-
  rda_sig_p$rdadapt %>%
  mutate(snp_names = snp_names) %>%
  filter(q.values <= 0.1)

# How many outlier SNPs were detected?
nrow(q_sig) 

# How many were in the intersect of all 3 methods
Reduce(intersect, list(
  q_sig$snp_names,
  rda_sig_p$rda_snps,
  rda_sig_z$rda_snps
))                      

# Visualizing RDA results with rda_plot()
rda_plot(mod_best, rda_sig_p$rda_snps, biplot_axes = c(1, 2), rdaplot = TRUE, manhattan = FALSE)

# RDA plot with outlier SNPs highlighted
png(file.path(plot_dir, "rda_plot_outliers.png"), width = 1600, height = 1200, res = 200)
rda_plot(mod_pRDA, rda_sig_z$rda_snps, rdaplot = TRUE, manhattan = FALSE, binwidth = 0.01)
dev.off()

# Manhattan plot with rda_sig_p results
rda_plot(mod_best, rda_sig_p$rda_snps, rda_sig_p$pvalues, rdaplot = FALSE, manhattan = TRUE)

# Interpreting RDA results
# Extract genotypes for outlier SNPs
rda_snps <- rda_sig_p$rda_snps
rda_gen <- gen_h[, rda_snps]
rda_sig_z 
rda_snps_z <- rda_sig_z$rda_snps
rda_gen_z <- gen_h[, rda_snps_z]

# Run correlation test
cor_df <- rda_cor(rda_gen, env_scaled)
cor_df_z <- rda_cor(rda_gen_z, env_scaled)

# Make a table from these results (displaying only the first 5 rows):
rda_table(cor_df, nrow = 5)
rda_table(cor_df_z, nrow = 5)

# Order by the strength of the correlation
rda_table(cor_df, order = TRUE, nrow = 5)
rda_table(cor_df_z, order = TRUE, nrow = 5)

# Only retain the top variable for each SNP based on the strength of the correlation
rda_table(cor_df, top = TRUE, nrow = 5)
rda_table(cor_df_z, top = TRUE, nrow = 5)

# Display results for only one environmental variable
rda_table(cor_df, var = "bio15", nrow = 5)
rda_table(cor_df_z, var = "bio15", nrow = 5)

# Save RDA correlation table (p < 0.05) to file
rda_sig_table <- rda_table(cor_df_z, sig = 0.05, order = TRUE)
write.csv(rda_sig_table, file.path(results_dir, "rda_p0.05.csv"), row.names = TRUE)

#######################################################################
## 02. Latent factor mixed models (LFMM)
#######################################################################

# K selection (run for reference; use chosen K in lfmm_run)
# select_K(gen_h, K_selection = "tracy_widom", criticalpoint = 2.0234)
# select_K(gen_h, K_selection = "quick_elbow", low = 0.08, max.pc = 0.90)
# select_K(gen_h, K_selection = "tess", coords = coord, Kvals = 1:10)
# select_K(gen_h, K_selection = "find_clusters", perc.pca = 90, max.n.clust = 10)

ridge_results <- lfmm_run(gen_h, env_scaled, K = 2, lfmm_method = "ridge")
lasso_results <- lfmm_run(gen_h, env_scaled, K = 2, lfmm_method = "lasso")

# Save LFMM tables to files
lfmm_ridge_tab <- lfmm_table(ridge_results$df)
write.table(lfmm_ridge_tab,
            file = file.path(results_dir, "lfmm_lasso_table.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

lfmm_lasso_tab <- lfmm_table(lasso_results$df)
write.table(lfmm_lasso_tab,
            file = file.path(results_dir, "lfmm_ridge_table.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# Save significant SNP tables as CSV for downstream use
write.csv(lasso_results$df, file.path(results_dir, "lfmm_lasso_results.csv"), row.names = TRUE)
write.csv(ridge_results$df, file.path(results_dir, "lfmm_ridge_results.csv"), row.names = TRUE)

# QQ plot
png(file.path(plot_dir, "lfmm_qqplot_lasso.png"), width = 1200, height = 1200, res = 200)
lfmm_qqplot(lasso_results$df)
dev.off()

# Manhattan plots
png(file.path(plot_dir, "lfmm_manhattan_lasso.png"), width = 1600, height = 800, res = 200)
lfmm_manhattanplot(lasso_results$df, sig = 0.05)
dev.off()

png(file.path(plot_dir, "lfmm_manhattan_ridge.png"), width = 1600, height = 800, res = 200)
lfmm_manhattanplot(ridge_results$df, sig = 0.05)
dev.off()



