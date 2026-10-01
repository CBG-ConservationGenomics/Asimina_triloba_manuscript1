# ============================================================
# Make an UpSet plot showing orthogroup overlap among species
# from the current OrthoFinder analysis.
#
# Input:
#   orthofinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv
#
# Output:
#   upset_plot_all_orthogroups.pdf
#
# Description:
#   - Automatically detects the newest OrthoFinder Results_*
#     directory.
#   - Uses ALL orthogroups identified by OrthoFinder.
#   - Converts gene counts to presence/absence for each species.
#   - Generates an UpSet plot showing the number of orthogroups
#     shared among different combinations of species.
# ============================================================


# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(UpSetR)
  library(dplyr)
  library(readr)
})


# ------------------------------------------------------------
# 1. Auto-detect newest OrthoFinder results directory
# ------------------------------------------------------------

results_dir <- list.dirs(
  "orthofinder",
  recursive = FALSE,
  full.names = TRUE
)

results_dir <- results_dir[
  grepl("Results_", basename(results_dir))
]

if (length(results_dir) == 0) {
  stop("No orthofinder/Results_* directory found")
}

results_dir <- results_dir[
  order(file.info(results_dir)$mtime, decreasing = TRUE)
]

OF_DIR <- results_dir[1]

cat("Using OrthoFinder results directory:\n")
cat(OF_DIR, "\n\n")


# ------------------------------------------------------------
# 2. Define input/output files
# ------------------------------------------------------------

orthogroups_file <- file.path(
  OF_DIR,
  "Orthogroups",
  "Orthogroups.GeneCount.tsv"
)

output_png <- "upset_plot_all_orthogroups.png"

if (!file.exists(orthogroups_file)) {
  stop("Cannot find file: ", orthogroups_file)
}


# ------------------------------------------------------------
# 3. Read complete OrthoFinder orthogroup count table
# ------------------------------------------------------------

cat("Reading Orthogroups.GeneCount.tsv...\n")

orthogroups_df <- read_tsv(
  orthogroups_file,
  show_col_types = FALSE
)

# Remove trailing "_annotated" from species names, if present
orthogroups_df <- orthogroups_df %>%
  rename_with(~ sub("_annotated$", "", .x))


# Identify species columns
# Excludes Orthogroup identifier and Total column

species_cols <- setdiff(
  names(orthogroups_df),
  c("Orthogroup", "Total")
)

cat("Species included:\n")
print(species_cols)

cat("\nTotal orthogroups:", nrow(orthogroups_df), "\n")


# ------------------------------------------------------------
# 4. Convert gene counts to presence/absence
# ------------------------------------------------------------

# 1 = orthogroup contains at least one protein from species
# 0 = orthogroup absent from species

upset_input <- orthogroups_df %>%
  select(all_of(species_cols)) %>%
  mutate(
    across(
      everything(),
      ~ as.integer(.x > 0)
    )
  ) %>%
  as.data.frame()


# ------------------------------------------------------------
# 5. Calculate useful summary statistics
# ------------------------------------------------------------

# Orthogroups represented in every species

shared_all <- sum(rowSums(upset_input) == length(species_cols))

cat(
  "Orthogroups present in all",
  length(species_cols),
  "species:",
  shared_all,
  "\n"
)


# Species-specific orthogroups

cat("\nSpecies-specific orthogroups:\n")

for (sp in species_cols) {

  species_specific <- sum(
    upset_input[[sp]] == 1 &
    rowSums(upset_input) == 1
  )

  cat(sp, ":", species_specific, "\n")
}
# ------------------------------------------------------------
# Annonaceae-specific orthogroups
# ------------------------------------------------------------

annonaceae <- c(
  "Annona_cherimola",
  "Annona_montana",
  "Asimina_triloba"
)

non_annonaceae <- c(
  "Lindera_megaphylla",
  "Magnolia_kwangsiensis",
  "Persea_americana"
)

# Present in ALL three Annonaceae AND absent from all other species
annonaceae_shared_only <- rowSums(upset_input[, annonaceae]) == length(annonaceae) &
                          rowSums(upset_input[, non_annonaceae]) == 0

n_annonaceae_shared_only <- sum(annonaceae_shared_only)

cat(
  "\nOrthogroups shared by all three Annonaceae and absent from other species:",
  n_annonaceae_shared_only,
  "\n"
)


# ------------------------------------------------------------
# 6. Prepare publication-style species labels
# ------------------------------------------------------------

# Replace underscores with spaces for plotting
pretty_names <- gsub("_", " ", species_cols)

# Rename columns only in the plotting object
colnames(upset_input) <- pretty_names


# ------------------------------------------------------------
# 7. Set species order
# ------------------------------------------------------------

# IMPORTANT:
# UpSetR displays the first species in this vector at the TOP
# when keep.order = TRUE.
#
# Put Asimina triloba first so that it appears at the top.

species_order <- c(
  "Persea americana",
  "Lindera megaphylla",
  "Magnolia kwangsiensis",
  "Annona montana",
  "Annona cherimola",
  "Asimina triloba"
)

# ------------------------------------------------------------
# 8. Generate UpSet plot
# ------------------------------------------------------------

cat("\nGenerating UpSet plot...\n")

png(
  output_png,
  width = 3300,
  height = 1950,
  res = 300
)

upset(
  upset_input,

  # Six species
  nsets = length(species_order),

  # Show the largest intersections
  nintersects = 20,

  # Species order
  # Asimina triloba will appear at the top
  sets = species_order,
  keep.order = TRUE,

  # Largest intersections first
  order.by = "freq",

  # Horizontal numbers above bars
  number.angles = 0,

  # Do not display empty intersections
  empty.intersections = "off",

  # Give more room to intersection plot
  mb.ratio = c(0.60, 0.40),

  # Text sizing
  text.scale = c(
    1.3,   # Intersection Size title
    1.0,   # Intersection Size axis labels
    0.9,   # Numbers above bars
    1.2,   # Set Size title
    1.0,   # Set Size axis labels
    1.05   # Species names
  ),

  mainbar.y.label = "Orthogroups",
  sets.x.label = "Orthogroups per species"
)


# ------------------------------------------------------------
# 9. Italicize species names
# ------------------------------------------------------------

# UpSetR does not directly provide an argument for changing the
# font face of individual species labels, so modify the text
# grobs after the plot has been drawn.

library(grid)

# Force creation of all graphical objects
grid.force()

# Get names of graphical objects
grobs <- grid.ls(print = FALSE)$name

# Loop through species labels
for (sp in species_order) {

  # Find text grobs containing this species name
  matching_grobs <- grobs[
    sapply(
      grobs,
      function(g) {
        obj <- tryCatch(
          grid.get(g),
          error = function(e) NULL
        )

        if (is.null(obj) || is.null(obj$label))
          return(FALSE)

        any(as.character(obj$label) == sp)
      }
    )
  ]

  for (g in matching_grobs) {

    if (sp == "Asimina triloba") {

      # Highlight focal species:
      # bold + italic
      grid.edit(
        g,
        gp = gpar(
          fontface = "bold.italic"
        )
      )

    } else {

      # All other species:
      # italic only
      grid.edit(
        g,
        gp = gpar(
          fontface = "italic"
        )
      )
    }
  }
}


dev.off()


# ------------------------------------------------------------
# 10. Finish
# ------------------------------------------------------------

cat("\nDone!\n")
cat("UpSet plot written to:", output_png, "\n")
