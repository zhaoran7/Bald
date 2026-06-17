suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

root <- "/mnt/f/bald/analysis/postgwas"
final_dir <- file.path(root, "magma_enrichment_final")
fig_dir <- file.path(final_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

enrich_f <- file.path(root, "result/magma_x/enrichment/magma_x_ora_msigdb_core.tsv")
out_f <- file.path(fig_dir, "MAGMA_enrichment_4traits_2x2.pdf")

trait_lab <- c(
  "bald.qt" = "Continuous",
  "bald12.bt" = "M shape",
  "bald13.bt" = "O shape",
  "bald14.bt" = "U shape"
)
traits <- names(trait_lab)
pt8 <- 8 / .pt

ora <- fread(enrich_f)
ora <- ora[trait %chin% traits & signal_set == "bonferroni" & n_overlap > 0]

category_map <- c(
  "GOBP" = "BP",
  "GOCC" = "CC",
  "GOMF" = "MF",
  "KEGG" = "KEGG",
  "REACTOME" = "REACTOME",
  "WP" = "WP",
  "HALLMARK" = "HALLMARK"
)
category_order <- c("BP", "CC", "MF", "KEGG", "REACTOME", "WP", "HALLMARK")
category_pal <- c(
  "BP" = "#8175A9",
  "CC" = "#D55E3F",
  "MF" = "#2F8DB6",
  "KEGG" = "#5AAE61",
  "REACTOME" = "#D6A33A",
  "WP" = "#7A6BA6",
  "HALLMARK" = "#666666"
)
category_short <- c(
  "BP" = "BP",
  "CC" = "CC",
  "MF" = "MF",
  "KEGG" = "KEGG",
  "REACTOME" = "REAC",
  "WP" = "WP",
  "HALLMARK" = "HM"
)

pretty_term <- function(x, width = 42L) {
  x <- sub("^(GOBP|GOCC|GOMF|KEGG|REACTOME|WP|HALLMARK)_", "", x)
  x <- tolower(gsub("_", " ", x))
  x <- gsub("\\bdna\\b", "DNA", x)
  x <- gsub("\\brna\\b", "RNA", x)
  x <- gsub("\\bwnt\\b", "WNT", x)
  x <- gsub("\\btcf\\b", "TCF", x)
  x <- gsub("\\bstat\\b", "STAT", x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

mix_with_white <- function(hex, amount) {
  amount <- pmin(1, pmax(0, amount))
  rgb <- grDevices::col2rgb(hex) / 255
  mixed <- 1 - amount * (1 - rgb)
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

prepare_category_plot <- function(tr, top_per_category = 3L) {
  z <- copy(ora[trait == tr])
  if (!nrow(z)) return(NULL)
  z[, category := unname(category_map[collection])]
  z <- z[category %chin% category_order]
  if (!nrow(z)) return(NULL)

  z <- z[order(fdr, p), head(.SD, top_per_category), by = category]
  z[, category := factor(category, levels = category_order)]

  pieces <- list()
  bars <- list()
  y_cursor <- 0
  row_step <- 1.42
  for (cat in rev(category_order)) {
    zz <- z[as.character(category) == cat]
    if (!nrow(zz)) next
    setorder(zz, -fdr, -p)
    zz[, y := y_cursor + row_step * seq_len(.N)]
    pieces[[cat]] <- zz
    bars[[cat]] <- data.table(
      category = cat,
      ymin = y_cursor + 0.45,
      ymax = y_cursor + row_step * nrow(zz) + 0.55,
      y_mid = y_cursor + row_step * (nrow(zz) + 1) / 2
    )
    y_cursor <- y_cursor + row_step * nrow(zz) + 1.55
  }

  dplot <- rbindlist(pieces, fill = TRUE)
  bplot <- rbindlist(bars, fill = TRUE)
  if (!nrow(dplot)) return(NULL)

  dplot[, score := -log10(pmax(fdr, 1e-300))]
  dplot[, score_scaled := {
    rng <- range(score, na.rm = TRUE)
    if (diff(rng) == 0) rep(0.65, .N) else (score - rng[1]) / diff(rng)
  }, by = category]
  dplot[, point_fill := mapply(
    mix_with_white,
    category_pal[as.character(category)],
    0.30 + 0.70 * score_scaled
  )]
  dplot[, label := pretty_term(set, 42L)]

  max_count <- max(dplot$n_overlap, na.rm = TRUE)
  x_left <- -max(2.5, max_count * 0.33)
  x_bar0 <- x_left * 0.62
  x_bar1 <- x_left * 0.50
  x_cat <- x_left * 0.78
  x_text <- max_count + max(0.5, max_count * 0.035)
  x_right <- max_count + max(12, max_count * 0.75)

  bplot[, `:=`(
    fill_col = category_pal[category],
    category_label = category_short[category],
    x_bar0 = x_bar0,
    x_bar1 = x_bar1,
    x_cat = x_cat
  )]
  dplot[, `:=`(
    x_text = x_text,
    x_right = x_right,
    x_left = x_left
  )]

  list(data = dplot, bars = bplot, x_left = x_left, x_right = x_right,
       max_count = max_count)
}

make_category_plot <- function(tr) {
  obj <- prepare_category_plot(tr)
  if (is.null(obj)) return(NULL)
  dplot <- obj$data
  bplot <- obj$bars
  breaks <- pretty(c(0, obj$max_count), n = 4)
  breaks <- breaks[breaks >= 0 & breaks <= obj$max_count]

  ggplot(dplot, aes(n_overlap, y)) +
    geom_rect(
      data = bplot,
      aes(xmin = x_bar0, xmax = x_bar1, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.95
    ) +
    geom_text(
      data = bplot,
      aes(x = x_cat, y = y_mid, label = category_label, color = category),
      inherit.aes = FALSE,
      size = pt8,
      hjust = 1
    ) +
    geom_point(
      aes(size = n_overlap, fill = point_fill, alpha = 0.55 + 0.45 * score_scaled),
      shape = 21,
      color = "#333333",
      stroke = 0.22
    ) +
    geom_text(
      aes(x = x_text, y = y, label = label, color = category),
      hjust = 0,
      size = pt8,
      lineheight = 0.92
    ) +
    scale_fill_identity() +
    scale_alpha_identity() +
    scale_color_manual(values = category_pal, guide = "none") +
    scale_size_continuous(range = c(1.2, 4.6), guide = "none") +
    scale_x_continuous(
      limits = c(obj$x_left, obj$x_right),
      breaks = breaks,
      expand = c(0, 0)
    ) +
    labs(
      title = unname(trait_lab[tr]),
      x = "Count",
      y = "Description"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "plain", size = 10, hjust = 0.02),
      axis.text.x = element_text(size = 8, face = "plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size = 8, face = "plain"),
      panel.grid.major.y = element_line(linewidth = 0.20, color = "grey88"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 18, 8, 5.5)
    )
}

plots <- lapply(traits, make_category_plot)
names(plots) <- traits
plots <- plots[!vapply(plots, is.null, logical(1))]
if (length(plots) != 4L) {
  stop("Expected four enrichment panels, got ", length(plots))
}

pdf(out_f, width = 16, height = 16, onefile = TRUE)
grid.newpage()
lay <- grid.layout(2, 2)
pushViewport(viewport(layout = lay))
for (i in seq_along(plots)) {
  row <- ifelse(i <= 2, 1, 2)
  col <- ifelse(i %% 2 == 1, 1, 2)
  print(plots[[i]], vp = viewport(layout.pos.row = row, layout.pos.col = col))
}
popViewport()
dev.off()

cat("Wrote", out_f, "\n")
