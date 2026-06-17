library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

root <- "/mnt/f/bald"
res <- file.path(root, "analysis/img/result/coloc_mri")
out_dir <- file.path(root, "analysis/img/plot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plink <- Sys.which("plink")
if (!nzchar(plink)) stop("plink not found")

ref <- file.path(res, "ld/coloc_locuszoom_ref")
genes_f <- "/mnt/f/Height2CVD/data/magma/ref/NCBI38.gene.loc"
ld_dir <- file.path(res, "ld")
icon_dir <- file.path(res, "brain_icons")

traits <- c("Continuous", "M shape", "O shape", "U shape")
bald_file <- c(
  "Continuous" = file.path(root, "data/gwas/bald/clean/bald/bald.gz"),
  "M shape"    = file.path(root, "data/gwas/bald/clean/bald12/bald12.gz"),
  "O shape"    = file.path(root, "data/gwas/bald/clean/bald13/bald13.gz"),
  "U shape"    = file.path(root, "data/gwas/bald/clean/bald14/bald14.gz")
)

cols <- c("#1D79B2", "#EBEA98", "#EB434A")
purple <- "#6A3D9A"
fs <- 8
gene_fs <- 6
track_h <- unit(2.5, "cm")
track_w <- unit(5.0, "cm")
gene_h <- unit(1.05, "cm")
head_h <- unit(0.55, "cm")
panel_gap_cm <- 0.32
col_gap_cm <- 0.65
panel_w_cm <- 6.1
fig_margin_cm <- 0.55

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
ok_rsid <- function(x) grepl("^rs[0-9]+$", x, ignore.case = TRUE)
fix_img_file <- function(f) sub("^/mnt/f/bald/img/data/gwas", file.path(root, "data/img/gwas"), f)
pretty_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Cortical volume ", "Cortical volume: ", x)
  x <- gsub("Cortical area ", "Cortical area: ", x)
  x <- gsub("Cortical thickness ", "Cortical thickness: ", x)
  x <- gsub("White matter ", "WM ", x)
  x <- gsub("Retrolenticular part of internal capsule", "RLIC", x)
  x <- gsub("Anterior limb of internal capsule", "ALIC", x)
  x
}

read_bald <- function(label, win) {
  cmd <- sprintf(
    "gzip -dc %s | awk -v c=%d -v s=%d -v e=%d 'NR==1 || ($2==c && $3>=s && $3<=e)'",
    shQuote(bald_file[[label]]), win$chr, win$start, win$end
  )
  x <- fread(cmd = cmd, showProgress = FALSE)
  if (!nrow(x)) return(data.table())
  setnames(x, c("SNP", "CHR", "POS", "EA", "NEA", "EAF", "N", "BETA", "SE", "P"))
  unique(x[ok_rsid(SNP), .(SNP, chr = as.integer(CHR), pos = as.integer(POS),
                           p = as.numeric(P))], by = "SNP")
}

read_bald_windows <- function(label, wins) {
  wf <- tempfile()
  fwrite(unique(wins[, .(chr, start, end)]), wf, sep = "\t", col.names = FALSE)
  cmd <- sprintf(
    "gzip -dc %s | awk 'NR==FNR{n++; c[n]=$1; s[n]=$2; e[n]=$3; next} FNR==1{print; next} {for(i=1;i<=n;i++) if($2==c[i] && $3>=s[i] && $3<=e[i]) {print; next}}' %s -",
    shQuote(bald_file[[label]]), shQuote(wf)
  )
  x <- fread(cmd = cmd, showProgress = FALSE)
  unlink(wf)
  if (!nrow(x)) return(data.table())
  setnames(x, c("SNP", "CHR", "POS", "EA", "NEA", "EAF", "N", "BETA", "SE", "P"))
  unique(x[ok_rsid(SNP), .(SNP, chr = as.integer(CHR), pos = as.integer(POS),
                           p = as.numeric(P))], by = "SNP")
}

lead_pos_from_bald <- function(snps) {
  tf <- tempfile(); writeLines(snps, tf)
  z <- rbindlist(lapply(bald_file, function(f) {
    out <- tempfile()
    system(sprintf("zgrep -w -F -f %s %s > %s || true", shQuote(tf), shQuote(f), shQuote(out)))
    if (!file.exists(out) || file.info(out)$size == 0) return(NULL)
    x <- fread(out, header = FALSE, showProgress = FALSE)
    unlink(out)
    setnames(x, c("SNP", "CHR", "POS", "EA", "NEA", "EAF", "N", "BETA", "SE", "P"))
    x[, .(lead_snp = SNP, chr38 = as.integer(CHR), lead_pos38 = as.integer(POS))]
  }), fill = TRUE)
  unlink(tf)
  unique(z, by = "lead_snp")
}

read_img <- function(f, snps, map) {
  f <- fix_img_file(f)
  sf <- tempfile(); writeLines(snps, sf)
  tf <- tempfile()
  system(sprintf("zgrep -w -F -f %s %s > %s || true", shQuote(sf), shQuote(f), shQuote(tf)))
  on.exit(unlink(c(sf, tf)), add = TRUE)
  if (!file.exists(tf) || file.info(tf)$size == 0) return(data.table())
  x <- fread(tf, header = FALSE, showProgress = FALSE)
  setnames(x, c("chr37", "SNP", "pos37", "a1", "a2", "beta", "se", "logp"))
  x <- merge(x[ok_rsid(SNP)], map, by = "SNP", all.x = FALSE)
  unique(x[, .(SNP, chr, pos, p = pmax(10^(-as.numeric(logp)), 1e-300))], by = "SNP")
}

ld_read <- function(f) {
  z <- fread(f)
  sc <- intersect(c("SNP_B", "ID_B"), names(z))[1]
  rc <- intersect(c("R2", "UNPHASED_R2"), names(z))[1]
  unique(z[, .(SNP = get(sc), r2 = as.numeric(get(rc)))], by = "SNP")
}

get_ld <- function(chr, lead, snps) {
  bim <- fread(paste0(ref, ".bim"), select = 2, col.names = "SNP")
  lead_ref <- if (lead %chin% bim$SNP) lead else snps[snps %chin% bim$SNP][1]
  if (is.na(lead_ref)) return(list(ref = lead, ld = data.table(SNP = character(), r2 = numeric())))
  old <- list.files(ld_dir, pattern = paste0("^.*", lead_ref, ".*\\.ld$"), full.names = TRUE)
  if (length(old)) {
    old <- old[order(file.info(old)$mtime, decreasing = TRUE)][1]
    z <- ld_read(old)
    return(list(ref = lead_ref, ld = unique(rbind(data.table(SNP = c(lead_ref, lead), r2 = 1), z), by = "SNP")))
  }
  key <- file.path(ld_dir, paste0("locus_", chr, "_", lead_ref, "_", safe(format(Sys.time(), "%H%M%S"))))
  xf <- paste0(key, ".snps")
  writeLines(unique(c(lead_ref, snps[snps %chin% bim$SNP])), xf)
  cmd <- sprintf(
    "%s --bfile %s --chr %s --extract %s --r2 --ld-snp %s --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0 --out %s >/dev/null 2>&1",
    shQuote(plink), shQuote(ref), chr, shQuote(xf), shQuote(lead_ref), shQuote(key)
  )
  system(cmd)
  if (!file.exists(paste0(key, ".ld"))) return(list(ref = lead_ref, ld = data.table(SNP = c(lead_ref, lead), r2 = 1)))
  list(ref = lead_ref, ld = unique(rbind(data.table(SNP = c(lead_ref, lead), r2 = 1), ld_read(paste0(key, ".ld"))), by = "SNP"))
}

crop_icon <- function(icon, hemi = NA_character_) {
  if (!file.exists(icon) || !requireNamespace("png", quietly = TRUE)) return(NULL)
  im <- png::readPNG(icon)
  w <- dim(im)[2]
  q <- floor(w / 4)
  if (identical(tolower(hemi), "right")) {
    im[, (3 * q + 1):w, , drop = FALSE]
  } else {
    im[, 1:q, , drop = FALSE]
  }
}

track_plot <- function(d, xlim, lead_pos, lead_snp, label, icon = NULL, hemi = NA_character_) {
  ymax <- max(d$y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  p <- ggplot(d, aes(pos, y))
  ic <- crop_icon(icon %||% "", hemi)
  if (!is.null(ic)) {
    p <- p + annotation_custom(
      rasterGrob(ic, interpolate = TRUE),
      xmin = xlim[2] - diff(xlim) * 0.245, xmax = xlim[2] - diff(xlim) * 0.005,
      ymin = ymax * 0.50, ymax = ymax * 0.96
    )
  }
  p <- p +
    geom_vline(xintercept = lead_pos, color = "#EB434A", linewidth = 0.24) +
    geom_point(aes(color = r2), size = 0.75, alpha = 0.88) +
    annotate("text", x = xlim[1] + diff(xlim) * 0.015, y = ymax * 1.01,
             label = label, hjust = 0, vjust = 1, family = "sans", size = fs / ggplot2::.pt) +
    scale_color_gradientn(colours = cols, limits = c(0, 1), guide = "none") +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, ymax * 1.05), clip = "off") +
    labs(x = NULL, y = expression(-log[10](italic(p)))) +
    theme_classic(base_size = fs, base_family = "sans") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.32),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = fs), axis.title.y = element_text(size = fs, margin = margin(r = 2)),
          axis.line = element_blank(), plot.margin = margin(0.5, 1, 0.5, 1))
  p + geom_point(data = d[SNP == lead_snp], aes(pos, y), inherit.aes = FALSE,
                 shape = 23, fill = purple, color = purple, size = 1.05, stroke = 0.25)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

gene_plot <- function(z, xlim) {
  if (nrow(z)) {
    z <- z[order(start, end)]
    z[, `:=`(x1 = pmax(start, xlim[1]), x2 = pmin(end, xlim[2]))]
    z <- z[x2 > x1]
    z[, y := 1L]
    if (nrow(z) > 1) {
      for (i in seq_len(nrow(z))) {
        while (any(z[seq_len(i - 1)]$end >= z$start[i] - diff(xlim) * .035 &
                   z[seq_len(i - 1)]$y == z$y[i])) z$y[i] <- z$y[i] + 1L
      }
    }
    z[, width := x2 - x1]
    lab <- z[order(-width)][seq_len(min(.N, 4))]
  } else lab <- data.table()
  ggplot() +
    {if (nrow(z)) geom_segment(data = z, aes(x = x1, xend = x2, y = y, yend = y),
                               color = "#2B2678", linewidth = 0.34,
                               arrow = arrow(length = unit(0.025, "inches"), type = "closed"))} +
    {if (nrow(lab)) geom_text(data = lab, aes(x = (x1 + x2) / 2, y = y + 0.16, label = gene),
                              color = "#2B2678", family = "sans", size = gene_fs / ggplot2::.pt,
                              check_overlap = TRUE)} +
    scale_x_continuous(limits = xlim, breaks = round(c(xlim[1], mean(xlim), xlim[2])),
                       labels = scales::comma, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.55, max(1.7, if (nrow(z)) max(z$y) + 0.55 else 1.7)), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = fs, base_family = "sans") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.32),
          axis.text.x = element_text(size = gene_fs), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line = element_blank(),
          plot.margin = margin(0.5, 1, 0.5, 1))
}

legend_grob <- function() {
  pal <- grDevices::colorRampPalette(cols)(256)
  gTree(children = gList(
    textGrob(expression(r^2~"with lead SNP"), x = 0.5, y = 0.93,
             gp = gpar(fontfamily = "sans", fontsize = gene_fs)),
    rasterGrob(matrix(pal, nrow = 1), x = 0.5, y = 0.55,
               width = 0.96, height = 0.28, interpolate = TRUE),
    segmentsGrob(x0 = 0.02, x1 = 0.98, y0 = 0.34, y1 = 0.34,
                 gp = gpar(linewidth = 0.35)),
    segmentsGrob(x0 = c(0.02, 0.50, 0.98), x1 = c(0.02, 0.50, 0.98),
                 y0 = 0.34, y1 = 0.40, gp = gpar(linewidth = 0.35)),
    textGrob(c("0", ".5", "1"), x = c(0.02, 0.50, 0.98), y = 0.16,
             gp = gpar(fontfamily = "sans", fontsize = gene_fs))
  ))
}

heatmap_grob <- function() {
  hm <- best[, .(PP.H4 = max(PP.H4)), by = .(trait, group_id)]
  hm <- merge(CJ(trait = traits, group_id = groups$group_id, unique = TRUE),
              hm, by = c("trait", "group_id"), all.x = TRUE)
  hm <- merge(hm, groups[, .(group_id, chr, lead_snp, lead_pos38, standard_name)], by = "group_id")
  hm[, trait := factor(trait, levels = rev(traits))]
  hm[, locus := factor(group_id, levels = groups$group_id)]
  locus_lab <- setNames(
    paste0(pretty_name(groups$standard_name), "\nchr", groups$chr, " | ", groups$lead_snp),
    groups$group_id
  )
  ggplot(hm, aes(locus, trait, fill = PP.H4)) +
    geom_tile(color = "#9A9A9A", linewidth = 0.18) +
    scale_fill_gradient(low = "#F6EFF7", high = "#7A0177", limits = c(0.8, 1),
                        na.value = "white", name = "PP.H4") +
    scale_x_discrete(labels = locus_lab, position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = fs, base_family = "sans") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = gene_fs, angle = 45, hjust = 0, vjust = 0, color = "black"),
          axis.text.y = element_text(size = fs, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = gene_fs, color = "black"),
          legend.text = element_text(size = gene_fs, color = "black"),
          plot.margin = margin(1, 1, 1, 1))
}

align <- function(gs) {
  mw <- do.call(unit.pmax, lapply(gs, `[[`, "widths"))
  lapply(gs, function(g) { g$widths <- mw; g })
}

best <- fread(file.path(res, "coloc_best_by_pair.tsv"))[PP.H4 > 0.8]
best[, file := fix_img_file(file)]
best[, group_id := paste(img, chr, start, end, sep = "|")]
lead_pos <- lead_pos_from_bald(unique(best$lead_snp))
groups <- best[, .SD[which.max(PP.H4)], by = group_id]
groups <- merge(groups, lead_pos, by = "lead_snp", all.x = TRUE, sort = FALSE)
groups[is.na(lead_pos38), `:=`(chr38 = chr, lead_pos38 = lead_pos)]
groups[, `:=`(chr = as.integer(chr38), start38 = pmax(1L, as.integer(lead_pos38 - 500000L)),
              end38 = as.integer(lead_pos38 + 500000L))]
cnt <- best[, .(n_bald = uniqueN(trait), traits = paste(unique(trait), collapse = ";")), by = group_id]
groups <- merge(groups, cnt, by = "group_id")
groups <- groups[order(-n_bald, chr, lead_pos38, standard_name)]
groups[, panel_h_cm := 0.55 + (n_bald + 1) * 2.5 + 1.05]
wins <- groups[, .(chr, start = start38, end = end38)]
message("Caching bald GWAS windows")
bald_cache <- lapply(traits, read_bald_windows, wins = wins)
names(bald_cache) <- traits
map_all <- unique(rbindlist(lapply(bald_cache, \(x) x[, .(SNP, chr, pos)]), fill = TRUE), by = "SNP")
message("Caching imaging GWAS windows")
img_cache <- setNames(lapply(unique(groups$file), \(f) read_img(f, map_all$SNP, map_all)),
                      unique(groups$file))

genes <- fread(genes_f, header = FALSE)
setnames(genes, c("gene_id", "chr", "start", "end", "strand", "gene"))
genes[, chr := as.integer(chr)]

make_panel <- function(i) {
  g <- groups[i]
  win <- list(chr = g$chr, start = g$start38, end = g$end38)
  xlim <- c(win$start, win$end)
  bald <- lapply(traits, \(tr) bald_cache[[tr]][chr == g$chr & pos >= xlim[1] & pos <= xlim[2]])
  names(bald) <- traits
  img <- img_cache[[g$file]][chr == g$chr & pos >= xlim[1] & pos <= xlim[2]]
  mem <- best[group_id == g$group_id]
  pdat <- rbindlist(c(list(data.table(track = "MRI", img[, .(SNP, pos, p)])),
                      lapply(unique(mem$trait), \(tr) data.table(track = tr, bald[[tr]][, .(SNP, pos, p)]))),
                    fill = TRUE)
  pdat <- unique(pdat[ok_rsid(SNP) & is.finite(pos) & is.finite(p)], by = c("track", "SNP"))
  ld <- get_ld(g$chr, g$lead_snp, unique(pdat$SNP))
  pdat <- merge(pdat, ld$ld, by = "SNP", all.x = TRUE)
  pdat[SNP == g$lead_snp, r2 := 1]
  pdat <- pdat[pos >= xlim[1] & pos <= xlim[2] & is.finite(r2)]
  pdat[, y := -log10(pmax(p, 1e-300))]
  if (!nrow(pdat)) return(nullGrob())

  order <- c("MRI", traits[traits %chin% unique(mem$trait)])
  icon <- file.path(icon_dir, paste0(g$img, ".png"))
  plots <- lapply(order, function(tr) {
    lab <- if (tr == "MRI") pretty_name(g$standard_name) else tr
    track_plot(pdat[track == tr], xlim, g$lead_pos38, g$lead_snp, lab,
               icon = if (tr == "MRI") icon else NULL, hemi = g$hemi)
  })
  gp <- gene_plot(genes[chr == g$chr & start <= xlim[2] & end >= xlim[1]], xlim)
  grobs <- align(lapply(c(plots, list(gp)), ggplotGrob))
  ttl <- textGrob(paste0("chr", g$chr, ":", scales::comma(g$lead_pos38), " | ", g$lead_snp),
                  gp = gpar(fontfamily = "sans", fontsize = fs))
  arrangeGrob(grobs = c(list(ttl), grobs), ncol = 1,
              heights = unit.c(head_h, rep(track_h, length(order)), gene_h))
}

panels <- lapply(seq_len(nrow(groups)), make_panel)
groups[, panel_i := .I]

p2 <- groups[n_bald == 1, panel_i]
p3 <- groups[n_bald == 2, panel_i]
p4 <- groups[n_bald == 3, panel_i]
p5 <- groups[n_bald == 4, panel_i]
layout_ids <- list(
  c(p2[1:3], p3[1]),
  c(p5, p4),
  c(p2[4:5], p3[2]),
  c(p2[6:7], p3[3]),
  c(p2[8:9], p3[4])
)
layout_ids <- lapply(layout_ids, \(x) x[!is.na(x)])

col_height <- function(ids) sum(groups$panel_h_cm[ids]) + max(length(ids) - 1, 0) * panel_gap_cm
max_col_cm <- max(vapply(layout_ids, col_height, numeric(1)))
col_grob <- lapply(layout_ids, function(ids) {
  gs <- list()
  hs <- list()
  for (j in seq_along(ids)) {
    gs[[length(gs) + 1]] <- panels[[ids[j]]]
    hs[[length(hs) + 1]] <- unit(groups$panel_h_cm[ids[j]], "cm")
    if (j < length(ids)) {
      gs[[length(gs) + 1]] <- nullGrob()
      hs[[length(hs) + 1]] <- unit(panel_gap_cm, "cm")
    }
  }
  rest <- max_col_cm - col_height(ids)
  if (rest > 0) {
    gs[[length(gs) + 1]] <- nullGrob()
    hs[[length(hs) + 1]] <- unit(rest, "cm")
  }
  arrangeGrob(grobs = gs, ncol = 1, heights = do.call(unit.c, hs))
})
body_grobs <- list()
body_widths <- list()
for (j in seq_along(col_grob)) {
  body_grobs[[length(body_grobs) + 1]] <- col_grob[[j]]
  body_widths[[length(body_widths) + 1]] <- unit(panel_w_cm, "cm")
  if (j < length(col_grob)) {
    body_grobs[[length(body_grobs) + 1]] <- nullGrob()
    body_widths[[length(body_widths) + 1]] <- unit(col_gap_cm, "cm")
  }
}
body <- arrangeGrob(grobs = body_grobs, ncol = length(body_grobs), widths = do.call(unit.c, body_widths))

pdf_file <- file.path(out_dir, "s13.pdf")
body_w_cm <- 5 * panel_w_cm + 4 * col_gap_cm
fig_w_cm <- body_w_cm + 2 * fig_margin_cm
fig_h_cm <- max_col_cm + 2 * fig_margin_cm
col_x0 <- (0:4) * (panel_w_cm + col_gap_cm)
rest_cm <- max_col_cm - vapply(layout_ids, col_height, numeric(1))
hm_w_cm <- 3 * panel_w_cm + 2 * col_gap_cm
hm_h_cm <- max(2.2, min(rest_cm[3:5]) - 0.45)
hm_h_cm <- min(hm_h_cm, 5.7)
leg_w_cm <- 2.65
leg_h_cm <- 0.95

pdf(pdf_file, width = fig_w_cm / 2.54, height = fig_h_cm / 2.54, useDingbats = FALSE, bg = "white")
grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))
grid.text("A", x = unit(fig_margin_cm, "cm"), y = unit(fig_h_cm - 0.18, "cm"),
          hjust = 0, vjust = 1, gp = gpar(fontfamily = "sans", fontsize = 10, fontface = "bold"))
pushViewport(viewport(x = unit(fig_margin_cm + body_w_cm / 2, "cm"),
                      y = unit(fig_margin_cm + max_col_cm / 2, "cm"),
                      width = unit(body_w_cm, "cm"),
                      height = unit(max_col_cm, "cm")))
grid.draw(body)
popViewport()
pushViewport(viewport(x = unit(fig_margin_cm + col_x0[2] + panel_w_cm / 2, "cm"),
                      y = unit(fig_margin_cm + leg_h_cm / 2, "cm"),
                      width = unit(leg_w_cm, "cm"),
                      height = unit(leg_h_cm, "cm")))
grid.draw(legend_grob())
popViewport()
grid.text("B", x = unit(fig_margin_cm + col_x0[3], "cm"),
          y = unit(fig_margin_cm + hm_h_cm + 0.24, "cm"),
          hjust = 0, vjust = 0, gp = gpar(fontfamily = "sans", fontsize = 10, fontface = "bold"))
pushViewport(viewport(x = unit(fig_margin_cm + col_x0[3] + hm_w_cm / 2, "cm"),
                      y = unit(fig_margin_cm + hm_h_cm / 2, "cm"),
                      width = unit(hm_w_cm, "cm"),
                      height = unit(hm_h_cm, "cm")))
grid.draw(ggplotGrob(heatmap_grob()))
popViewport()
dev.off()

manifest <- groups[, .(panel = panel_i, chr, lead_snp, lead_pos_grch38 = lead_pos38,
                       imaging_trait = standard_name, bald_traits = traits,
                       n_colocalized_bald_traits = n_bald, PP.H4_max = PP.H4)]
fwrite(manifest, file.path(out_dir, "s13.tsv"), sep = "\t")
message(pdf_file)
