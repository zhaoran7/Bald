suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

root <- "/mnt/f/bald/analysis/postgwas"
final_dir <- file.path(root, "magma_enrichment_final")
fig_dir <- file.path(final_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_res_f <- file.path(root, "result/magma_x/summary/magma_gene_results_chr1_22_X.tsv")
out_f <- file.path(fig_dir, "MAGMA_gene_manhattan_4traits.pdf")

trait_lab <- c(
  "bald.qt" = "Continuous",
  "bald12.bt" = "M shape",
  "bald13.bt" = "O shape",
  "bald14.bt" = "U shape"
)
traits <- names(trait_lab)
pt8 <- 8 / .pt

all <- fread(gene_res_f)
all <- all[trait %chin% traits & !is.na(CHR_INT) & is.finite(P)]
all[, trait_label := factor(unname(trait_lab[trait]), levels = unname(trait_lab))]
all[, logp := -log10(pmax(P, 1e-300))]
all[, LABEL := fifelse(!is.na(SYMBOL) & nzchar(SYMBOL), SYMBOL, LABEL)]

make_offsets <- function(d) {
  chr_info <- d[!is.na(CHR_INT), .(chr_len = max(MID, na.rm = TRUE)), by = CHR_INT]
  setorder(chr_info, CHR_INT)
  chr_info[, offset := shift(cumsum(chr_len), fill = 0)]
  chr_info[, center := offset + chr_len / 2]
  chr_info
}

chr_info <- make_offsets(all)
d <- merge(all, chr_info[, .(CHR_INT, offset)], by = "CHR_INT")
d[, BPcum := MID + offset]
d[, color_group := factor(CHR_INT %% 2)]
d[, bonf := -log10(0.05 / uniqueN(GENE)), by = trait]
d[, bonf_sig := bonferroni_pvalue < 0.05]

sig_wide <- dcast(
  unique(d[, .(GENE, SYMBOL, LABEL, CHR_INT, trait, bonf_sig)]),
  GENE + SYMBOL + LABEL + CHR_INT ~ trait,
  value.var = "bonf_sig",
  fill = FALSE
)
sig_wide[, n_bald_bonf := rowSums(.SD), .SDcols = traits]
shared_genes <- sig_wide[n_bald_bonf == length(traits), GENE]
d[, shared_all4_bonf := GENE %chin% shared_genes]

select_top_labels <- function(x, n_labels = 40L) {
  x <- x[!grepl("^ENSG", LABEL)]
  setorder(x, P)
  x <- head(x, min(nrow(x), n_labels))
  x[, label_rank_trait := seq_len(.N)]
  x[]
}

base_top <- rbindlist(lapply(traits, function(tr) {
  x <- select_top_labels(d[trait == tr], n_labels = 40L)
  x[, is_base_top40 := TRUE]
  x
}), fill = TRUE)

shared_to_sync <- unique(base_top[shared_all4_bonf == TRUE, GENE])
synced_shared <- d[trait %chin% traits & GENE %chin% shared_to_sync & !grepl("^ENSG", LABEL)]
synced_shared[, is_base_top40 := FALSE]
lab <- unique(rbindlist(list(base_top, synced_shared), fill = TRUE), by = c("trait", "GENE"))
lab[, is_synced_shared := shared_all4_bonf & GENE %chin% shared_to_sync]

layout_labels <- function(x, tr) {
  x <- copy(x)
  setorder(x, BPcum, P)
  n <- nrow(x)
  if (!n) return(x)
  x_range <- range(d$BPcum, na.rm = TRUE)
  x_pad <- diff(x_range) * 0.035
  x[, label_slot := seq_len(.N)]
  x[, label_x := seq(x_range[1] + x_pad, x_range[2] - x_pad, length.out = .N)]

  y_max <- max(d[trait == tr, logp], na.rm = TRUE)
  y_min <- min(d[trait == tr, logp], na.rm = TRUE)
  y_span <- max(6, y_max - y_min)
  x[, elbow_y := y_max + y_span * 0.085]
  x[, label_y := y_max + y_span * 0.240]
  x[]
}

lab <- rbindlist(lapply(traits, function(tr) layout_labels(lab[trait == tr], tr)), fill = TRUE)
lab[, label_color_group := fifelse(is_synced_shared, "shared", "label")]

p <- ggplot(d, aes(BPcum, logp)) +
  geom_point(aes(color = color_group), size = 0.38, alpha = 0.75) +
  geom_point(
    data = d[shared_all4_bonf == TRUE],
    aes(color = "shared"),
    size = 0.62,
    alpha = 0.92
  ) +
  geom_hline(aes(yintercept = bonf), linetype = "dashed",
             linewidth = 0.25, color = "#B24A4A") +
  geom_segment(
    data = lab,
    aes(x = BPcum, y = logp, xend = BPcum, yend = elbow_y,
        group = interaction(trait, GENE)),
    inherit.aes = FALSE,
    linewidth = 0.16,
    color = "#222222",
    alpha = 0.26
  ) +
  geom_segment(
    data = lab,
    aes(x = BPcum, y = elbow_y, xend = label_x, yend = label_y,
        group = interaction(trait, GENE)),
    inherit.aes = FALSE,
    linewidth = 0.16,
    color = "#222222",
    alpha = 0.26
  ) +
  geom_text(
    data = lab,
    aes(x = label_x, y = label_y, label = LABEL, color = label_color_group),
    inherit.aes = FALSE,
    size = pt8,
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    fontface = "plain"
  ) +
  facet_wrap(~trait_label, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_color_manual(
    values = c("0" = "#8FAAC8", "1" = "#D9D9D9",
               "shared" = "#B24A4A", "label" = "#222222"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = chr_info$center,
    labels = sub("23", "X", as.character(chr_info$CHR_INT)),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  coord_cartesian(clip = "on") +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "plain"),
    axis.text = element_text(size = 8, face = "plain"),
    axis.title = element_text(size = 9, face = "plain"),
    panel.spacing.y = unit(0.34, "in"),
    plot.margin = margin(14, 18, 8, 8)
  )

ggsave(out_f, p, width = 14, height = 14.8, device = cairo_pdf)
cat("Wrote", out_f, "\n")
