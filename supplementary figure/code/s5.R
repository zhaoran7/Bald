library(data.table)
library(ggplot2)

root <- "/mnt/d/bald"
out <- file.path(root, "sup")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

lead <- fread(file.path(root, "gu/lead/lead.assoc"))
reg <- fread(file.path(root, "gu/report/region_summary.tsv"))
sel <- fread(file.path(root, "gu/report/selected_region.tsv"))
match <- fread(file.path(root, "gu/report/core_archaic_match.tsv"))

steps <- data.table(
  y = 6:1,
  title = c("COJO lead variants", "LD-defined core blocks", "Risk haplotype alleles",
            "Archaic matching", "YRI outgroup filter", "Phylogenetic analysis"),
  n = c(nrow(lead), nrow(reg), uniqueN(paste(reg$trait, reg$id)),
        uniqueN(paste(match[n_match_risk > 0]$trait, match[n_match_risk > 0]$id)),
        nrow(sel), length(list.files(file.path(root, "gu/phy"), pattern = "full\\.phy_phyml_tree\\.txt$", recursive = TRUE))),
  note = c("Conditional lead variants from four GWAS traits",
           "High-LD SNPs define local risk-core haplotypes",
           "Alleles carried on lead-risk haplotypes",
           "Matched against Altai, Chagyrskaya, Vindija, Denisova and Denisova 25",
           "Regions retained when the risk-core allele is absent in YRI",
           "Maximum-likelihood trees for retained archaic-like regions")
)
wrap <- function(x, n = 62) vapply(strwrap(x, n, simplify = FALSE), paste, "", collapse = "\n")
steps[, lab := paste0(title, "\n", "n = ", format(n, big.mark = ","), "\n", wrap(note))]

edges <- data.table(x = 0, xend = 0, y = steps$y[-nrow(steps)] - 0.33, yend = steps$y[-1] + 0.33)

p <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.09, "inches")), linewidth = 0.35, color = "#4D4D4D") +
  geom_label(data = steps, aes(0, y, label = lab), label.size = 0.25, fill = "white",
             color = "black", size = 2.7, family = "Arial", lineheight = 0.95,
             label.padding = unit(0.28, "lines")) +
  scale_x_continuous(limits = c(-3.15, 3.15), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, 6.55), expand = c(0, 0)) +
  theme_void(base_size = 8, base_family = "Arial") +
  theme(plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(out, "S5.pdf"), p, width = 6.6, height = 7.2, device = cairo_pdf)
