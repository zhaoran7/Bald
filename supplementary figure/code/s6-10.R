library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

base <- "/mnt/f/bald"
out <- file.path(base, "manuscript/sup")
cache <- file.path(tempdir(), "locuszoom_block_cache")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

lift <- "/mnt/d/bald/tool/liftOver"
chain <- "/mnt/d/bald/tool/hg38ToHg19.over.chain.gz"
refgene_file <- "/mnt/f/bald/analysis/archaic/plot/causal/resource/refGene.hg19.txt.gz"
if (!file.exists(refgene_file)) refgene_file <- "/mnt/d/bald/gu/ref_tracks/refGene.txt.gz"
vcf_dir <- "/mnt/f/refGen/1kg_phase3/vcf"

traits <- data.table(
  trait = c("bald", "bald12", "bald13", "bald14"),
  label = c("Continuous", "M shape", "O shape", "U shape"),
  gwas = file.path(base, "data/gwas/bald/clean",
                   c("bald", "bald12", "bald13", "bald14"),
                   c("bald.gz", "bald12.gz", "bald13.gz", "bald14.gz"))
)

regions <- data.table(
  region = paste0("region", 1:5),
  trait = c("bald13", "bald12", "bald", "bald", "bald13"),
  id = c("1.rs12405323.170418964", "6.rs9349320.45269814",
         "8.rs1041791.109112070", "12.rs417915.52842074",
         "23.rs6525167.66423003"),
  chr = c("1", "6", "8", "12", "X"),
  chr_num = c(1L, 6L, 8L, 12L, 23L),
  lead = c("rs12405323", "rs9349320", "rs1041791", "rs417915", "rs6525167"),
  lead_bp = c(170418964L, 45269814L, 109112070L, 52842074L, 66423003L)
)

sel <- fread(file.path(base, "analysis/archaic/result/locus/report/selected_region.tsv"))
regions <- merge(
  regions,
  sel[, .(trait, id, core_start, core_end, n_ld_snp, core_size_bp, best_lineage)],
  by = c("trait", "id"), all.x = TRUE, sort = FALSE
)
regions[, `:=`(plot_start = lead_bp - 500000L, plot_end = lead_bp + 500000L)]
regions[core_start < plot_start, plot_start := core_start - 50000L]
regions[core_end > plot_end, plot_end := core_end + 50000L]
regions[, `:=`(pull_start = plot_start - 2000000L, pull_end = plot_end + 2000000L,
               risk_label = paste0("Risk haplotype ", seq_len(.N), " (",
                                   sprintf("%.1f", core_size_bp / 1000), " kb)"))]

fmt_bp <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
run <- function(cmd) {
  s <- system(cmd)
  if (!identical(s, 0L)) stop("Command failed: ", cmd)
}

lw <- .18
pal_r2 <- colorRampPalette(c("#FFFFFF", "#FEE5E5", "#FCAE91", "#FB6A4A", "#CB181D"))(101)
pal_dp <- colorRampPalette(c("#FFFFFF", "#E5ECF6", "#9ECAE1", "#3182BD", "#08519C"))(101)
ld_col <- function(x, pal) pal[pmax(1L, pmin(101L, round(pmax(0, pmin(1, x)) * 100) + 1L))]
scale_var <- function(x) {
  z <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(z)
  rg <- range(x[ok])
  z[ok] <- if (diff(rg) == 0) 1 else (x[ok] - rg[1]) / diff(rg)
  z
}
scale_r2 <- function(x) pmax(0, pmin(1, (x - .95) / .05))
scale_dp <- function(x) pmax(0, pmin(1, (x - .95) / .05))

read_ld <- function(r, tr) {
  d <- file.path(base, "analysis/archaic/result/locus/ld", tr,
                 paste0("chr", r$chr))
  f <- file.path(d, paste0(r$lead, ".vcor"))
  if (!file.exists(f)) return(data.table(pos37 = integer(), r2 = numeric()))
  x <- fread(f)
  setnames(x, c("#CHROM_A", "POS_A", "ID_A", "CHROM_B", "POS_B", "ID_B", "UNPHASED_R2"),
           c("CHR_A", "POS_A", "ID_A", "CHR_B", "pos37", "SNP", "r2"), skip_absent = TRUE)
  unique(rbind(
    x[, .(SNP, pos37 = as.integer(pos37), r2 = pmax(0, pmin(1, as.numeric(r2))))],
    data.table(SNP = r$lead, pos37 = r$lead_bp, r2 = 1)
  ), by = "pos37")
}

lift_gwas <- function(r, tr) {
  cf <- file.path(cache, sprintf("%s_%s.tsv", r$region, tr$trait))
  if (file.exists(cf) && file.info(cf)$size > 0) return(fread(cf))
  raw <- tempfile(tmpdir = cache, fileext = ".raw.tsv")
  bed <- tempfile(tmpdir = cache, fileext = ".bed")
  mapped <- tempfile(tmpdir = cache, fileext = ".mapped.bed")
  unmapped <- tempfile(tmpdir = cache, fileext = ".unmapped.bed")
  run(sprintf(
    "zcat -f %s | awk -v c=%s -v s=%d -v e=%d 'NR==1 || ($2==c && $3>=s && $3<=e)' > %s",
    shQuote(tr$gwas), ifelse(r$chr == "X", 23, as.integer(r$chr)),
    r$pull_start, r$pull_end, shQuote(raw)
  ))
  x <- fread(raw)
  if (nrow(x) == 0) return(data.table())
  x[, POS := as.integer(POS)]
  b <- x[is.finite(POS), .(
    chrom = ifelse(CHR %in% c("23", 23), "chrX", paste0("chr", CHR)),
    start = pmax(POS - 1L, 0L),
    end = POS,
    SNP
  )]
  fwrite(b, bed, sep = "\t", col.names = FALSE)
  run(sprintf("%s %s %s %s %s >/dev/null 2>&1",
              shQuote(lift), shQuote(bed), shQuote(chain), shQuote(mapped), shQuote(unmapped)))
  m <- fread(mapped, col.names = c("chrom37", "start37", "end37", "SNP"))
  m[, pos37 := as.integer(end37)]
  x <- merge(x, m[, .(SNP, pos37)], by = "SNP")
  x <- x[pos37 >= r$plot_start & pos37 <= r$plot_end]
  if (nrow(x) == 0) return(data.table())
  x[, `:=`(trait = tr$trait, trait_label = tr$label,
           logp = -log10(pmax(as.numeric(P), .Machine$double.xmin)))]
  fwrite(x, cf, sep = "\t")
  x
}

read_refgene <- function() {
  g <- fread(refgene_file, header = FALSE)
  setnames(g, 1:13, c("bin", "tx", "chrom", "strand", "txStart", "txEnd",
                      "cdsStart", "cdsEnd", "exonCount", "exonStarts",
                      "exonEnds", "score", "gene"))
  g[, len := txEnd - txStart]
  g[order(gene, -len), .SD[1], by = gene][order(chrom, txStart)]
}

place_genes <- function(g, r) {
  span <- r$plot_end - r$plot_start
  g[, w := pmax(x1 - x0, nchar(gene) * span * .0065)]
  g[, xlab := pmin(r$plot_end - w / 2, pmax(r$plot_start + w / 2, (x0 + x1) / 2))]
  g[, `:=`(lab0 = pmax(r$plot_start, xlab - w / 2),
           lab1 = pmin(r$plot_end, xlab + w / 2))]
  last <- numeric()
  rr <- integer(nrow(g))
  for (i in seq_len(nrow(g))) {
    k <- which(last + span * .010 < g$lab0[i])[1]
    if (is.na(k)) { last <- c(last, -Inf); k <- length(last) }
    rr[i] <- k
    last[k] <- max(last[k], g$lab1[i])
  }
  g[, row := rr]
  g
}

risk_band <- function(r) {
  dummy <- "Continuous"
  ggplot(data.table(dummy = dummy, xmin = r$core_start / 1e6, xmax = r$core_end / 1e6,
                    x = (r$core_start + r$core_end) / 2e6, y0 = .34, y1 = .34,
                    ylab = .76, lab = r$risk_label), aes(x, y0)) +
    facet_grid(dummy ~ ., switch = "y") +
    geom_segment(aes(x = xmin, xend = xmax, y = y0, yend = y1),
                 color = "#8B5A2B", linewidth = .42, lineend = "butt", inherit.aes = FALSE) +
    geom_text(aes(x = x, y = ylab, label = lab), size = 2.2, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(r$plot_start, r$plot_end) / 1e6,
                       breaks = seq(r$plot_start, r$plot_end, length.out = 3) / 1e6,
                       labels = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(y = expression(-log[10](P)), x = NULL) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 7, colour = "transparent"),
          axis.title.y = element_text(size = 7, colour = "transparent"),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.margin = margin(2, 2, 2, 2))
}

gene_track <- function(r, refgene) {
  dummy <- "Continuous"
  basep <- ggplot(data.table(x = c(r$plot_start, r$plot_end) / 1e6, y = 0, dummy = dummy), aes(x, y)) +
    geom_blank() +
    facet_grid(dummy ~ ., switch = "y") +
    scale_x_continuous(limits = c(r$plot_start, r$plot_end) / 1e6,
                       breaks = seq(r$plot_start, r$plot_end, length.out = 3) / 1e6,
                       labels = fmt_bp(seq(r$plot_start, r$plot_end, length.out = 3))) +
    labs(y = expression(-log[10](P)), x = NULL)
  g <- refgene[chrom == paste0("chr", r$chr) & txEnd >= r$plot_start & txStart <= r$plot_end]
  if (!nrow(g)) {
    p <- basep + scale_y_continuous(limits = c(-1, 0), expand = c(0, 0)) +
      theme_classic(base_size = 7) +
      theme(strip.background = element_blank(), strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, size = 7, colour = "transparent"),
            axis.title.y = element_text(size = 7, colour = "transparent"),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank(), axis.text.x = element_text(size = 6),
            plot.margin = margin(2, 2, 2, 2))
    attr(p, "n_gene_rows") <- 1L
    return(p)
  }
  g[, `:=`(x0 = pmax(txStart, r$plot_start), x1 = pmin(txEnd, r$plot_end))]
  g <- place_genes(g[x1 > x0][order(x0, x1)], r)
  g[, dummy := dummy]
  ex <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
    st <- as.integer(strsplit(g$exonStarts[i], ",", fixed = TRUE)[[1]])
    en <- as.integer(strsplit(g$exonEnds[i], ",", fixed = TRUE)[[1]])
    data.table(dummy = dummy, gene = g$gene[i], row = g$row[i],
               x0 = pmax(st, r$plot_start), x1 = pmin(en, r$plot_end))
  }), fill = TRUE)
  ex <- ex[x1 > x0]
  nr <- max(g$row)
  ts <- if (nr > 12) 1.25 else if (nr > 8) 1.45 else if (nr > 5) 1.65 else 1.9
  p <- basep +
    geom_segment(data = g, aes(x = x0 / 1e6, xend = x1 / 1e6, y = -row, yend = -row),
                 linewidth = lw, color = "#2B246D", inherit.aes = FALSE) +
    geom_rect(data = ex, aes(xmin = x0 / 1e6, xmax = x1 / 1e6,
                             ymin = -row - .09, ymax = -row + .09),
              fill = "#2B246D", color = NA, inherit.aes = FALSE) +
    geom_text(data = g, aes(x = xlab / 1e6, y = -row + .27, label = gene),
              size = ts, color = "#2B246D", fontface = "italic", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(-nr - .55, -.35), expand = c(0, 0)) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 7, colour = "transparent"),
          axis.title.y = element_text(size = 7, colour = "transparent"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.text.x = element_text(size = 6),
          plot.margin = margin(2, 2, 2, 2))
  attr(p, "n_gene_rows") <- nr
  p
}

stack_left <- function(p0, p1, p2, nr) {
  g0 <- ggplotGrob(p0)
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  w <- Reduce(unit.pmax, list(g0$widths, g1$widths, g2$widths))
  g0$widths <- w
  g1$widths <- w
  g2$widths <- w
  arrangeGrob(g0, g1, g2, ncol = 1, heights = c(.55, 5.1, max(1.05, .50 + .20 * nr)))
}

locus_plot <- function(r, refgene) {
  ld <- read_ld(r, r$trait)
  z <- rbindlist(lapply(seq_len(nrow(traits)), function(i) {
    x <- lift_gwas(r, traits[i])
    if (!nrow(x)) return(NULL)
    x <- merge(x, ld[, .(pos37, r2)], by = "pos37", all.x = TRUE)
    x[is.na(r2), r2 := 0]
    x[, r2 := pmax(0, pmin(1, r2))]
    x
  }), fill = TRUE)
  z[, lead_hit := pos37 == r$lead_bp | SNP == r$lead]
  p1 <- ggplot() +
    geom_vline(xintercept = r$lead_bp / 1e6, color = "#D62728", linewidth = lw) +
    geom_point(data = z[lead_hit != TRUE & r2 < .2], aes(pos37 / 1e6, logp, color = r2),
               size = .48, alpha = .72) +
    geom_point(data = z[lead_hit != TRUE & r2 >= .2 & r2 < .5], aes(pos37 / 1e6, logp, color = r2),
               size = .55, alpha = .82) +
    geom_point(data = z[lead_hit != TRUE & r2 >= .5 & r2 < .8], aes(pos37 / 1e6, logp, color = r2),
               size = .62, alpha = .92) +
    geom_point(data = z[lead_hit != TRUE & r2 >= .8], aes(pos37 / 1e6, logp, color = r2),
               size = .72, alpha = 1) +
    geom_point(data = z[lead_hit == TRUE], aes(pos37 / 1e6, logp),
               shape = 23, size = 1.35, stroke = lw,
               color = "#7B3294", fill = "#7B3294") +
    facet_grid(trait_label ~ ., scales = "free_y", switch = "y") +
    scale_color_gradientn(colors = c("#2B7BBA", "#9ECAE1", "#F4D35E", "#F89C4B", "#D73027"),
                          limits = c(0, 1), breaks = c(0, .5, 1),
                          name = expression(r^2)) +
    scale_x_continuous(limits = c(r$plot_start, r$plot_end) / 1e6,
                       breaks = seq(r$plot_start, r$plot_end, length.out = 3) / 1e6,
                       labels = NULL) +
    labs(y = expression(-log[10](P)), x = NULL) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 7),
          legend.position = c(.88, .95), legend.direction = "horizontal",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5),
          legend.key.width = unit(.42, "cm"), legend.key.height = unit(.10, "cm"),
          plot.margin = margin(2, 2, 2, 2))
  p0 <- risk_band(r)
  p2 <- gene_track(r, refgene)
  stack_left(p0, p1, p2, attr(p2, "n_gene_rows"))
}

vcf_gt <- function(r, pos) {
  pos <- sort(unique(as.integer(pos)))
  cf <- file.path(cache, sprintf("%s_%d_%d_%d_gt.tsv", r$region, length(pos), min(pos), max(pos)))
  if (!file.exists(cf) || file.info(cf)$size == 0) {
    bed <- tempfile(tmpdir = cache, fileext = ".bed")
    fwrite(data.table(chr = r$chr, start = pos - 1L, end = pos),
           bed, sep = "\t", col.names = FALSE)
    vcf0 <- file.path(base, "analysis/archaic/result/locus/coreVcf", r$trait, r$id, "kg.vcf.gz")
    vcf <- if (file.exists(vcf0)) vcf0 else file.path(vcf_dir, paste0("chr", r$chr, ".vcf.gz"))
    run(sprintf("bcftools query -R %s -f '%%CHROM\t%%POS\t%%ID\t%%REF\t%%ALT[\t%%GT]\n' %s > %s",
                shQuote(bed), shQuote(vcf), shQuote(cf)))
  }
  fread(cf, header = FALSE, fill = TRUE)
}

gt_diag <- function(gt, ref, alt, allele) {
  gt <- sub(":.*$", "", gt)
  gt <- gsub("\\|", "/", gt)
  unlist(lapply(strsplit(gt, "/", fixed = TRUE), function(v) {
    v <- c(v, rep(".", 2L))[1:2]
    a <- suppressWarnings(as.integer(v))
    b <- rep(NA_character_, 2L)
    b[which(a == 0)] <- ref
    b[which(a == 1)] <- alt
    as.integer(b == allele)
  }), use.names = FALSE)
}

pair_ld <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]
  b <- b[ok]
  if (length(a) < 20 || length(unique(a)) < 2 || length(unique(b)) < 2) return(c(r2 = NA, dp = NA))
  pa <- mean(a)
  pb <- mean(b)
  pab <- mean(a == 1 & b == 1)
  D <- pab - pa * pb
  den <- pa * (1 - pa) * pb * (1 - pb)
  r2 <- ifelse(den > 0, D^2 / den, NA_real_)
  dmax <- if (D >= 0) min(pa * (1 - pb), (1 - pa) * pb) else min(pa * pb, (1 - pa) * (1 - pb))
  dp <- ifelse(dmax > 0, abs(D / dmax), NA_real_)
  c(r2 = pmin(1, r2), dp = pmin(1, dp))
}

ld_block_plot <- function(r) {
  core <- fread(file.path(base, "analysis/archaic/result/locus/report/core_risk.tsv"))
  d <- core[trait == r$trait & id == r$id & is_diagnostic_archaic == TRUE &
              toupper(diagnostic_archaic_allele) %chin% c("A", "C", "G", "T")]
  d <- unique(d[order(pos), .(pos = as.integer(pos), diagnostic_archaic_allele = toupper(diagnostic_archaic_allele))], by = "pos")
  if (nrow(d) < 2) return(list(plot = ggplot() + theme_void(), r2_range = c(.95, 1), dp_range = c(.95, 1)))
  z <- vcf_gt(r, d$pos)
  if (!nrow(z)) return(list(plot = ggplot() + theme_void(), r2_range = c(.95, 1), dp_range = c(.95, 1)))
  setnames(z, 1:5, c("chr", "pos", "snp", "ref", "alt"))
  z[, pos := as.integer(pos)]
  z <- merge(z, d, by = "pos")
  z <- z[nchar(ref) == 1 & nchar(alt) == 1 & diagnostic_archaic_allele %chin% c(ref, alt)]
  z <- unique(z[order(pos)], by = "pos")
  n <- nrow(z)
  if (n < 2) return(list(plot = ggplot() + theme_void(), r2_range = c(.95, 1), dp_range = c(.95, 1)))
  gtm <- lapply(seq_len(n), function(i) {
    gt_diag(as.character(unlist(z[i, 6:(ncol(z) - 1), with = FALSE], use.names = FALSE)),
            z$ref[i], z$alt[i], z$diagnostic_archaic_allele[i])
  })
  mat <- rbindlist(lapply(seq_len(n), function(i) {
    rbindlist(lapply(seq_len(n), function(j) {
      ld <- if (i == j) c(r2 = 1, dp = 1) else pair_ld(gtm[[i]], gtm[[j]])
      type <- ifelse(j >= i, "r2", "dp")
      data.table(i = i, j = j, type = type, val = ifelse(type == "r2", ld["r2"], ld["dp"]))
    }))
  }))
  dp_range <- c(.95, 1)
  mat[, sc := NA_real_]
  mat[type == "r2", sc := scale_r2(val)]
  mat[type == "dp", sc := scale_dp(val)]
  mat[, col := fifelse(is.na(sc), "#FFFFFF",
                       fifelse(type == "r2", ld_col(sc, pal_r2), ld_col(sc, pal_dp)))]
  axis_pos <- 1 + (z$pos - min(z$pos)) / max(1, diff(range(z$pos))) * (n - 1)
  br <- pretty(range(z$pos) / 1e6, n = 3)
  br <- br[br >= min(z$pos) / 1e6 & br <= max(z$pos) / 1e6]
  con <- data.table(x = seq_len(n), xend = axis_pos)
  tick <- data.table(x = axis_pos, y0 = -1.18, y1 = -.82)
  p <- ggplot() +
    geom_tile(data = mat, aes(i, n - j + 1, fill = col), width = 1, height = 1,
              color = "grey82", linewidth = lw) +
    scale_fill_identity() +
    geom_rect(aes(xmin = .5, xmax = n + .5, ymin = .5, ymax = n + .5),
              fill = NA, color = "black", linewidth = lw) +
    geom_segment(data = con, aes(x = x, xend = xend, y = .46, yend = -.82),
                 linewidth = lw, color = "black") +
    geom_segment(data = tick, aes(x = x, xend = x, y = y0, yend = y1),
                 linewidth = lw, color = "#F04B4B") +
    geom_segment(aes(x = 1, xend = n, y = -1.35, yend = -1.35), linewidth = lw) +
    annotate("text", x = 1 + (br - min(z$pos) / 1e6) / max(1e-9, diff(range(z$pos) / 1e6)) * (n - 1),
             y = -1.80, label = sprintf("%.1f", br), size = 1.9) +
    annotate("text", x = n / 2, y = -2.28,
             label = paste0("Chromosome ", r$chr, " coordinate (Mb)"), size = 2.4) +
    coord_fixed(ratio = 1, xlim = c(.5, n + .5), ylim = c(-4.80, n + .5), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 2, 6, 2))
  list(plot = p, r2_range = c(.95, 1), dp_range = dp_range)
}

add_ld_legend <- function(r2_range, dp_range) {
  grid.text(expression(Correlation~(R^2)), x = .30, y = .092, gp = gpar(fontsize = 7))
  for (k in 0:9) grid.rect(x = .19 + k * .025, y = .060, width = .024, height = .018,
                           gp = gpar(fill = ld_col(k / 9, pal_r2), col = "grey70", lwd = lw))
  grid.text(sprintf("%.2f", r2_range[1]), x = .20, y = .028, gp = gpar(fontsize = 6))
  grid.text(sprintf("%.2f", r2_range[2]), x = .42, y = .028, gp = gpar(fontsize = 6))
  grid.text("Correlation (D')", x = .70, y = .092, gp = gpar(fontsize = 7))
  for (k in 0:9) grid.rect(x = .59 + k * .025, y = .060, width = .024, height = .018,
                           gp = gpar(fill = ld_col(k / 9, pal_dp), col = "grey70", lwd = lw))
  grid.text(sprintf("%.2f", dp_range[1]), x = .60, y = .028, gp = gpar(fontsize = 6))
  grid.text(sprintf("%.2f", dp_range[2]), x = .82, y = .028, gp = gpar(fontsize = 6))
}

refgene <- read_refgene()

for (i in seq_len(nrow(regions))) {
  r <- regions[i]
  block_obj <- ld_block_plot(r)
  pdf(file.path(out, paste0(r$region, ".pdf")), width = 11.2, height = 6.8, useDingbats = FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    3, 4,
    heights = unit.c(unit(.5, "cm"), unit(1, "null"), unit(.5, "cm")),
    widths = unit.c(unit(.5, "cm"), unit(1.10, "null"), unit(1, "null"), unit(.5, "cm"))
  )))
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.draw(locus_plot(r, refgene))
  grid.text("A", x = unit(.02, "npc"), y = unit(.98, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 12, fontface = "bold"))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  pushViewport(viewport(width = .96, height = .98))
  print(block_obj$plot, newpage = FALSE)
  add_ld_legend(block_obj$r2_range, block_obj$dp_range)
  grid.text("B", x = unit(.02, "npc"), y = unit(.98, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 12, fontface = "bold"))
  popViewport(3)
  dev.off()
}
