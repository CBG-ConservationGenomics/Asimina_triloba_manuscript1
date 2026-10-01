#!/usr/bin/env Rscript

# ============================================================
# Plot gene family expansions and contractions across species
#
# Purpose:
#   Generate a diverging bar plot summarizing gene family
#   expansions and contractions inferred by CAFE5 across the
#   six species included in the comparative analysis.
#
#   - Expansions are plotted to the right.
#   - Contractions are plotted to the left.
#   - Light bars represent all inferred gene family changes.
#   - Dark overlays represent the subset classified as rapidly
#     changing based on CAFE5 branch probabilities < 0.05.
#
# Input:
#   Gene family counts summarized in Supplementary Table S17.
#
# Outputs:
#   gene_family_changes.png
#   gene_family_changes.pdf
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})


# ------------------------------------------------------------
# 2. Enter gene family counts
# ------------------------------------------------------------

df <- data.frame(
  Species = c(
    "Lindera megaphylla",
    "Persea americana",
    "Asimina triloba",
    "Annona cherimola",
    "Annona montana",
    "Magnolia kwangsiensis"
  ),

  Total_Expanded = c(
    1397, 2195, 359, 5660, 1241, 3570
  ),

  Total_Contracted = c(
    2464, 1170, 4276, 2071, 2602, 2859
  ),

  Rapid_Expanded = c(
    117, 86, 15, 738, 222, 97
  ),

  Rapid_Contracted = c(
    36, 35, 276, 66, 164, 25
  )
)


# ------------------------------------------------------------
# 3. Prepare data for plotting
# ------------------------------------------------------------

# Contractions are converted to negative values so that they
# extend to the left of zero. Expansions remain positive.

plot_df <- df %>%
  transmute(
    Species,
    Total_Expansion    = Total_Expanded,
    Rapid_Expansion    = Rapid_Expanded,
    Total_Contraction  = -Total_Contracted,
    Rapid_Contraction  = -Rapid_Contracted
  )


# ------------------------------------------------------------
# 4. Set species order
# ------------------------------------------------------------

# ggplot2 draws the first factor level at the bottom,
# so define the desired top-to-bottom order and reverse it.

species_order <- c(
  "Asimina triloba",
  "Annona cherimola",
  "Annona montana",
  "Lindera megaphylla",
  "Magnolia kwangsiensis",
  "Persea americana"
)

plot_df$Species <- factor(
  plot_df$Species,
  levels = rev(species_order)
)
# ------------------------------------------------------------
# 5. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(y = Species)
) +

  # Total contractions
  geom_col(
    aes(
      x = Total_Contraction,
      fill = "Total contraction"
    ),
    width = 0.70
  ) +

  # Total expansions
  geom_col(
    aes(
      x = Total_Expansion,
      fill = "Total expansion"
    ),
    width = 0.70
  ) +

  # Rapidly changing contractions
  # Drawn over the total contraction bars
  geom_col(
    aes(
      x = Rapid_Contraction,
      fill = "Rapid contraction"
    ),
    width = 0.70
  ) +

  # Rapidly changing expansions
  # Drawn over the total expansion bars
  geom_col(
    aes(
      x = Rapid_Expansion,
      fill = "Rapid expansion"
    ),
    width = 0.70
  ) +

  # Zero line separating contractions and expansions
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +

  # Colorblind-friendly colors
  scale_fill_manual(
    values = c(
      "Total contraction" = "#E69F73",
      "Rapid contraction" = "#A63603",
      "Total expansion"   = "#70C1A6",
      "Rapid expansion"   = "#006B4F"
    ),
    breaks = c(
      "Total contraction",
      "Rapid contraction",
      "Total expansion",
      "Rapid expansion"
    )
  ) +

  # Display contraction values as positive axis labels
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.12, 0.12))
  ) +

  # Axis labels
  labs(
    x = "Number of gene families",
    y = NULL,
    fill = NULL
  ) +

  # Publication-style theme
  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.text.y = element_text(
      face = "italic",
      size = 11
    ),

    axis.title.x = element_text(
      size = 12
    ),

    legend.position = "bottom",

    legend.text = element_text(
      size = 10
    )
  ) +

  # Arrange legend in two rows
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )


# ------------------------------------------------------------
# 6. Display figure
# ------------------------------------------------------------

print(p)


# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "gene_family_changes.png",
  plot = p,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = "gene_family_changes.pdf",
  plot = p,
  width = 9,
  height = 5.5
)


# ------------------------------------------------------------
# 8. Finish
# ------------------------------------------------------------

cat("\nDone!\n")
cat("Created:\n")
cat("  gene_family_changes.png\n")
cat("  gene_family_changes.pdf\n")

~
~
