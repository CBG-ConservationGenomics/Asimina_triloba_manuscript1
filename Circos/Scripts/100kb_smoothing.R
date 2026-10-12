library(dplyr)

setwd("~/Desktop/Pawpaw_circos/data/")
dat <- read.table("SNPS_LDpruned_r0.9_n10_100kb.hist",
                  header = FALSE)
colnames(dat) <- c("chr","start","end","value")

k <- 5  # window size (odd)
half <- (k - 1) / 2

dat <- dat %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(value_smooth = zoo::rollapply(
    value, width = k, FUN = mean,
    fill = "extend", align = "center"
  )) %>%
  ungroup()

write.table(dat[,c("chr","start","end","value_smooth")],
            "SNPS_LDpruned_r0.9_n10_100kb.smooth.hist",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

###########
gene <- read.table("Asimina_triloba_gene_100kb.hist",
                  header = FALSE)
colnames(gene) <- c("chr","start","end","value")

k <- 5  # window size (odd)
half <- (k - 1) / 2

gene <- gene %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(value_smooth = zoo::rollapply(
    value, width = k, FUN = mean,
    fill = "extend", align = "center"
  )) %>%
  ungroup()

write.table(gene[,c("chr","start","end","value_smooth")],
            "Asimina_triloba_gene_100kb_smooth.hist",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##############
repeats <- read.table("Atri_repeats_100kb.bed",
                   header = FALSE)
colnames(repeats) <- c("chr","start","end","value")

k <- 5  # window size (odd)
half <- (k - 1) / 2

repeats <- repeats %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(value_smooth = zoo::rollapply(
    value, width = k, FUN = mean,
    fill = "extend", align = "center"
  )) %>%
  ungroup()

write.table(repeats[,c("chr","start","end","value_smooth")],
            "Atri_repeats_100kb_smooth.hist",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

