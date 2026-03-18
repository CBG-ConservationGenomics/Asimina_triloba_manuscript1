# Hierarchical clustering from genetic distance (pipeline artifact)

data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
in_dir <- file.path(data_dir, "outputs", "00-preprocessing")
out_dir <- file.path(data_dir, "outputs", "03-hclust")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("ggplot2", "ggtree", "ape")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install once and re-run.\n",
    "For ggtree (Bioconductor):\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('ggtree')"
  )
}

library(ggplot2)
library(ggtree)
library(ape)

as_dist <- function(x) {
  if (inherits(x, "dist")) return(x)
  if (is.matrix(x)) return(stats::as.dist(x))
  if (is.data.frame(x)) return(stats::as.dist(as.matrix(x)))
  if (is.list(x)) {
    for (nm in c("dist", "dist_mat", "distance", "D")) {
      if (!is.null(x[[nm]])) return(as_dist(x[[nm]]))
    }
  }
  stop("Don't know how to convert `euc_gendist` to a `dist` object.")
}

dist_matrix <- as_dist(euc_gendist)
hc <- hclust(dist_matrix, method = "average")

newick_path <- file.path(out_dir, "hclust_average.newick")
ape::write.tree(ape::as.phylo(hc), file = newick_path)

png(file.path(out_dir, "hclust_dendrogram.png"), width = 1800, height = 1200, res = 200)
plot(hc, main = "Hierarchical clustering (Euclidean genetic distance)", xlab = "", sub = "")
dev.off()

p_circ <- ggtree(hc, layout = "circular") +
  geom_tiplab(size = 2) +
  ggtitle("Circular hierarchical clustering (Euclidean genetic distance)")

ggsave(
  filename = file.path(out_dir, "hclust_circular.png"),
  plot = p_circ,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)

