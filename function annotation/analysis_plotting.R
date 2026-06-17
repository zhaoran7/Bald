library(data.table)
library(ggplot2)
for (p in c("officer", "rvg", "openxlsx")) if (!requireNamespace(p, quietly = TRUE)) stop("Please install: ", p)

base <- if (.Platform$OS.type == "windows") "F:/bald/analysis/archaic" else "/mnt/f/bald/analysis/archaic"
root <- file.path(base, "result/locus")
ref <- file.path(base, "ref_tracks")
if (!file.exists(file.path(ref, "refGene.txt.gz"))) {
  ref <- if (.Platform$OS.type == "windows") "D:/bald/gu/ref_tracks" else "/mnt/d/bald/gu/ref_tracks"
}
if (!file.exists(file.path(ref, "refGene.txt.gz"))) stop("missing ref_tracks: ", ref)
out <- file.path(base, "plot")
cache <- file.path(root, "plot/.Fig3_cache")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

sel <- fread(file.path(root, "report/selected_region.tsv"))
main_regions <- data.table(
  trait = c("bald13", "bald12", "bald", "bald", "bald13"),
  id = c("1.rs12405323.170418964",
         "6.rs9349320.45269814",
         "8.rs1041791.109112070",
         "12.rs417915.52842074",
         "23.rs6525167.66423003"),
  title = c("O shape", "M shape", "Continuous; U shape", "Continuous", "Continuous; M shape; O shape"),
  lead_label = c("rs12405323", "rs9349320", "rs1041791", "rs417915", "rs143054933 / rs5919284 / rs6525167")
)
regions <- merge(main_regions, sel, by = c("trait", "id"), all.x = TRUE, sort = FALSE)
regions[, archaics := strsplit(matched_archaics, ";", fixed = TRUE)]
setorder(regions, lead_chr, lead_bp)

cell_key <- c(GM12878 = "Gm12878", H1hESC = "H1hesc", HSMM = "Hsmm",
              HUVEC = "Huvec", K562 = "K562", NHEK = "Nhek", NHLF = "Nhlf")
mark_key <- c(H3K27ac = "H3k27ac", H3K4me1 = "H3k4me1", H3K4me3 = "H3k4me3")
cell_col <- c(GM12878 = "#ef777d", H1hESC = "#facd7c", HSMM = "#87cdbd",
              HUVEC = "#83c5e0", K562 = "#6b6bb5", NHEK = "#9b65b2", NHLF = "#df79ad")
fs6 <- 6 / ggplot2::.pt
fs7 <- 7 / ggplot2::.pt
variant_col <- "#2f78b7"

nice_top <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(1)
  m <- max(v)
  s <- if (m < 50) 5 else 10
  z <- ceiling(m / s) * s
  if (z <= m) z <- z + s
  z
}

refs <- c("refGene", "wgEncodeRegDnaseClusteredV3", "wgEncodeRegTfbsClusteredV3", "wgEncodeOpenChromDnaseNhekPk")
for (cl in cell_key) for (mk in mark_key) refs <- c(refs, paste0("wgEncodeBroadHistone", cl, mk, "StdPk"))
for (x in unique(refs)) {
  f <- file.path(ref, paste0(x, ".txt.gz"))
  if (!file.exists(f) || file.info(f)$size == 0)
    download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/", x, ".txt.gz"),
                  f, mode = "wb", quiet = TRUE)
}

read_region <- function(tr, rid) fread(file.path(root, "report/region_summary.tsv"))[trait == tr & id == rid][1]
read_core <- function(tr, rid) fread(file.path(root, "report/core_risk.tsv"))[trait == tr & id == rid]
read_ld <- function(tr, snp, bp) {
  f <- file.path(root, "ld", tr, "ld.tsv")
  if (!file.exists(f)) return(data.table(pos = integer(), lead_r2 = numeric()))
  x <- fread(f)[lead_snp == snp & lead_bp == bp, .(pos, lead_r2 = R2)]
  unique(x, by = "pos")
}
read_gene <- function(chr, st, en) {
  g <- fread(file.path(ref, "refGene.txt.gz"), header = FALSE)
  setnames(g, 1:13, c("bin", "tx", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                      "exonCount", "exonStarts", "exonEnds", "score", "gene"))
  g <- g[chrom == paste0("chr", chr) & txEnd >= st & txStart <= en]
  g[, len := txEnd - txStart]
  g[order(gene, -len), .SD[1], by = gene][order(txStart)]
}

read_peak <- function(file, chr, st, en) {
  f <- file.path(ref, paste0(file, ".txt.gz"))
  if (!file.exists(f) || file.info(f)$size == 0) return(data.table())
  x <- fread(f, header = FALSE, sep = "\t", quote = "", fill = Inf)
  if (ncol(x) < 4) return(data.table())
  setnames(x, 1:4, c("bin", "chrom", "start", "end"))
  x[, signal := 1]
  if (ncol(x) >= 8) x[, signal := suppressWarnings(as.numeric(V8))]
  if (ncol(x) >= 6) x[!is.finite(signal), signal := suppressWarnings(as.numeric(V6))]
  x <- x[chrom == paste0("chr", chr) & end >= st & start <= en, .(start, end, signal)]
  x[, `:=`(start = pmax(start, st), end = pmin(end, en))]
  x[is.finite(signal) & end > start]
}

read_signal <- function(file, class, cell, mark, chr, st, en) {
  cf <- file.path(cache, sprintf("signal_%s_chr%s_%d_%d.tsv", file, chr, st, en))
  if (file.exists(cf)) return(fread(cf))
  if (.Platform$OS.type == "windows") return(data.table())
  bw <- if (class == "histone mark") {
    file.path(ref, paste0(sub("StdPk$", "StdSig", file), ".bigWig"))
  } else if (class == "skin open chromatin" && cell == "NHEK") {
    file.path(ref, "wgEncodeUwDnaseNhekRawRep1.bigWig")
  } else NA_character_
  exe <- file.path(ref, "bigWigToBedGraph")
  if (is.na(bw) || !file.exists(bw) || !file.exists(exe)) return(data.table())
  tmp <- tempfile(fileext = ".bedGraph")
  ok <- system2(exe, c(bw, tmp, paste0("-chrom=chr", chr), paste0("-start=", st), paste0("-end=", en)),
                stdout = FALSE, stderr = FALSE)
  if (!identical(ok, 0L) || !file.exists(tmp)) return(data.table())
  x <- fread(tmp, col.names = c("chrom", "start", "end", "signal"))
  x[, `:=`(start = pmax(start, st), end = pmin(end, en))]
  x <- x[end > start & is.finite(signal), .(start, end, signal)]
  fwrite(x, cf, sep = "\t")
  x
}

empty_shape <- function() data.table(x = numeric(), y0 = numeric(), y = numeric(), grp = character(),
                                     col = character(), axis_top = numeric(), h = numeric())

ribbon <- function(x, y0, h, id, col, top = NA_real_) {
  if (!nrow(x)) return(empty_shape())
  x <- rbind(x[, .(x = start, signal)], x[, .(x = end, signal)])[order(x)]
  x <- x[, .(signal = mean(signal, na.rm = TRUE)), by = x]
  raw <- x$signal
  k <- min(61L, max(9L, 2L * floor(nrow(x) / 70L) + 1L))
  if (nrow(x) >= k) {
    sm <- frollmean(raw, n = k, align = "center", fill = NA_real_)
    x[, signal := fifelse(is.na(sm), raw, sm)]
  }
  top <- if (is.finite(top) && top > 0) top else nice_top(x$signal)
  x[, `:=`(y0 = y0, y = y0 + h * pmin(pmax(signal, 0), top) / top,
           grp = id, col = col, axis_top = top, h = h)]
  x
}

peak_ribbon <- function(x, y0, h, id, col, top = NA_real_) {
  if (!nrow(x)) return(empty_shape())
  top <- if (is.finite(top) && top > 0) top else nice_top(x$signal)
  rbindlist(lapply(seq_len(nrow(x)), function(i) {
    z <- x[i]
    data.table(x = seq(z$start, z$end, length.out = 11),
               y0 = y0,
               y = y0 + h * pmin(pmax(z$signal, 0), top) / top * c(0, .12, .35, .75, .9, 1, .85, .55, .26, .08, 0),
               grp = paste(id, i, sep = "_"), col = col, axis_top = top, h = h)
  }))
}

axis_table <- function(x, xat, len) {
  empty_axis <- data.table(x = numeric(), x0 = numeric(), xlab = numeric(),
                           y0 = numeric(), y_mid = numeric(), y_top = numeric(),
                           lab_top = character())
  if (!nrow(x)) return(empty_axis)
  a <- x[, .(h = max(h, na.rm = TRUE), axis_top = max(axis_top, na.rm = TRUE)), by = y0]
  a <- a[is.finite(axis_top) & axis_top > 0]
  if (!nrow(a)) return(empty_axis)
  a[, `:=`(x = xat, x0 = xat - len, xlab = xat - len * 1.45,
           y_mid = y0 + h / 2, y_top = y0 + h,
           lab_top = format(round(axis_top), big.mark = ",", scientific = FALSE, trim = TRUE))]
  a
}

overlap_pos <- function(pos, peak) {
  if (!nrow(peak)) return(rep(FALSE, length(pos)))
  vapply(pos, function(p) any(peak$start <= p & peak$end >= p), logical(1))
}

add_track <- function(meta, plist, chr, st, en, core, file, label, class, cell, mark, col, weight) {
  x <- read_peak(file, chr, st, en)
  id <- paste0("T", nrow(meta) + 1L)
  if (nrow(x)) x[, `:=`(cell = cell, mark = mark)]
  plist[[id]] <- x
  meta <- rbind(meta, data.table(id, file, label, class, cell, mark, col,
                                 n_peak = nrow(x),
                                 n_core_overlap = if (nrow(x)) sum(overlap_pos(core$pos, x)) else 0L,
                                 max_signal = if (nrow(x)) max(x$signal, na.rm = TRUE) else 0),
                fill = TRUE)
  meta[, score := n_core_overlap * 100 + log1p(max_signal) * 4 + n_peak * .02 + weight]
  list(meta = meta, plist = plist)
}

catalog <- function(chr, st, en, core) {
  meta <- data.table()
  plist <- list()
  z <- add_track(meta, plist, chr, st, en, core, "wgEncodeRegDnaseClusteredV3", "ENCODE\nDNase I", "open chromatin", "ENCODE", "DNase", "#2b75b8", 30); meta <- z$meta; plist <- z$plist
  z <- add_track(meta, plist, chr, st, en, core, "wgEncodeRegTfbsClusteredV3", "ENCODE\nTF binding", "TF binding", "ENCODE", "TFBS", "#222222", 20); meta <- z$meta; plist <- z$plist
  z <- add_track(meta, plist, chr, st, en, core, "wgEncodeOpenChromDnaseNhekPk", "NHEK\nDNase I", "skin open chromatin", "NHEK", "DNase", "#8FAAC8", 80); meta <- z$meta; plist <- z$plist
  for (cl in names(cell_key)) for (mk in names(mark_key)) {
    wt <- c(NHEK = 45, NHLF = 25, HSMM = 10, HUVEC = 10, GM12878 = 8, H1hESC = 8, K562 = 8)[cl]
    f <- paste0("wgEncodeBroadHistone", cell_key[cl], mark_key[mk], "StdPk")
    z <- add_track(meta, plist, chr, st, en, core, f, paste(cl, mk, sep = "\n"),
                   "histone mark", cl, mk, cell_col[cl], wt)
    meta <- z$meta
    plist <- z$plist
  }
  list(meta = meta[order(-score)], plist = plist)
}

pick_main <- function(meta) {
  rbind(meta[class == "open chromatin" & cell == "ENCODE" & mark == "DNase"][1],
        meta[class == "TF binding" & cell == "ENCODE" & mark == "TFBS"][1],
        meta[class == "skin open chromatin" & cell == "NHEK" & mark == "DNase"][1],
        meta[class == "histone mark" & cell == "NHEK" & mark == "H3K27ac"][1],
        fill = TRUE)
}

read_arch <- function(tr, rid, core, archaics) {
  cf <- file.path(cache, paste0("arch_", tr, "_", gsub("[^A-Za-z0-9]+", "_", rid), ".tsv"))
  if (file.exists(cf)) return(fread(cf))
  if (.Platform$OS.type == "windows") stop("Run once in WSL to create archaic cache: ", cf)
  nm <- c(Altai = "altai", Chagyr = "chagyr", Vindija = "vindija", Denisova = "denisova", Denisova25 = "denisova25")
  d <- file.path(root, "coreVcf", tr, rid)
  x <- rbindlist(lapply(archaics, function(a) {
    v <- fread(cmd = sprintf("bcftools query -f '%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s",
                             shQuote(file.path(d, paste0(nm[a], ".vcf.gz")))),
               col.names = c("pos", "ref", "alt", "gt"), fill = TRUE)
    v <- merge(core[, .(pos, risk_core_allele)], v, by = "pos", all.x = TRUE, sort = FALSE)
    v[, present := !is.na(ref) & !is.na(alt) & !is.na(gt)]
    v[, match := mapply(function(ref, alt, gt, risk) {
      if (anyNA(c(ref, alt, gt, risk))) return(FALSE)
      alle <- c(ref, strsplit(alt, ",", fixed = TRUE)[[1]])
      idx <- unique(suppressWarnings(as.integer(strsplit(gsub("\\|", "/", sub(":.*", "", gt)), "/", fixed = TRUE)[[1]])))
      any(alle[idx + 1L] == risk, na.rm = TRUE)
    }, ref, alt, gt, risk_core_allele)]
    v[, .(pos, archaic = a, present, match)]
  }), fill = TRUE)
  fwrite(x, cf, sep = "\t")
  x
}

exon_table <- function(g, x0, x1) {
  if (!nrow(g)) return(data.table(gene = character(), y = numeric(), start = numeric(), end = numeric()))
  rbindlist(lapply(seq_len(nrow(g)), function(i) {
    a <- as.integer(strsplit(g$exonStarts[i], ",", fixed = TRUE)[[1]])
    b <- as.integer(strsplit(g$exonEnds[i], ",", fixed = TRUE)[[1]])
    data.table(gene = g$gene[i], y = g$y[i],
               start = pmax(a[!is.na(a)], x0), end = pmin(b[!is.na(b)], x1))
  }), fill = TRUE)[end > start]
}

arrow_table <- function(g, span) {
  if (!nrow(g)) return(data.table(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
  rbindlist(lapply(seq_len(nrow(g)), function(i) {
    z <- g[i]
    a <- z$tx0 + span * .02
    b <- z$tx1 - span * .02
    if (!is.finite(a) || !is.finite(b) || b <= a) return(NULL)
    xs <- seq(a, b, by = max(span * .015, 600))
    if (!length(xs)) return(NULL)
    d <- ifelse(z$strand == "+", 1, -1)
    rbind(data.table(x = xs, xend = xs + d * span * .004, y = z$y, yend = z$y + .035),
          data.table(x = xs, xend = xs + d * span * .004, y = z$y, yend = z$y - .035))
  }), fill = TRUE)
}

track_shape <- function(z, plist, chr, st, en, y0, h) {
  sig <- read_signal(z$file, z$class, z$cell, z$mark, chr, st, en)
  if (nrow(sig)) ribbon(sig, y0, h, z$id, z$col) else peak_ribbon(plist[[z$id]], y0, h, z$id, z$col)
}

hist_summary <- function(sel, meta, plist, chr, st, en, h, tops) {
  drop <- sel[class == "histone mark", paste(cell, mark)]
  rbindlist(lapply(names(mark_key), function(mk) {
    y0 <- c(H3K27ac = 2.85, H3K4me1 = 1.60, H3K4me3 = 0.35)[mk]
    top <- tops[[mk]]
    rbindlist(lapply(names(cell_key), function(cl) {
      if (paste(cl, mk) %in% drop) return(data.table())
      f <- paste0("wgEncodeBroadHistone", cell_key[cl], mark_key[mk], "StdPk")
      sig <- read_signal(f, "histone mark", cl, mk, chr, st, en)
      pk <- plist[[meta[file == f, id][1]]]
      if (is.null(pk)) pk <- data.table()
      id <- paste(mk, cl)
      if (nrow(sig)) ribbon(sig, y0, h, id, cell_col[cl], top = top) else peak_ribbon(pk, y0, h, id, cell_col[cl], top = top)
    }), fill = TRUE)
  }), fill = TRUE)
}

plot_region <- function(r, full = FALSE) {
  reg <- read_region(r$trait, r$id)
  core <- read_core(r$trait, r$id)
  ld_all <- read_ld(r$trait, reg$lead_snp, reg$lead_bp)
  if (!reg$lead_bp %in% ld_all$pos) ld_all <- rbind(ld_all, data.table(pos = reg$lead_bp, lead_r2 = 1), fill = TRUE)
  ld_all[pos == reg$lead_bp & !is.finite(lead_r2), lead_r2 := 1]
  core <- merge(core, ld_all, by = "pos", all.x = TRUE, sort = FALSE)
  core[pos == reg$lead_bp & !is.finite(lead_r2), lead_r2 := 1]
  core[!is.finite(lead_r2), lead_r2 := 0.98]
  gwas_dot <- unique(rbind(ld_all[, .(pos, lead_r2)], core[, .(pos, lead_r2)], fill = TRUE), by = "pos")
  gwas_dot <- gwas_dot[is.finite(lead_r2) & (lead_r2 >= 0.8 | pos == reg$lead_bp)]
  gwas_dot[, lead := pos == reg$lead_bp]
  chr <- ifelse(reg$lead_chr == 23, "X", as.character(reg$lead_chr))
  st <- reg$core_start
  en <- reg$core_end
  cmn <- min(c(core$pos, gwas_dot$pos), na.rm = TRUE)
  cmx <- max(c(core$pos, gwas_dot$pos), na.rm = TRUE)
  span <- max(cmx - cmn, (en - st) * .08)
  x0 <- cmn - span * .03
  x1 <- cmx + span * .07
  span <- x1 - x0

  arch <- merge(read_arch(r$trait, r$id, core, unlist(r$archaics)), core[, .(pos)], by = "pos")
  arch[, y := 10.15 - .16 * (match(archaic, unlist(r$archaics)) - 1)]

  g <- read_gene(chr, st, en)
  if (nrow(g)) g[, `:=`(y = 9.42 - .16 * (.I - 1), tx0 = pmax(txStart, x0), tx1 = pmin(txEnd, x1))]
  if (!nrow(g)) g <- data.table(gene = character(), y = numeric(), tx0 = numeric(), tx1 = numeric(),
                                 strand = character(), exonStarts = character(), exonEnds = character())
  ex <- exon_table(g, x0, x1)
  ar <- arrow_table(g, span)

  cc <- catalog(chr, st, en, core)
  sel <- if (full) cc$meta[order(class, mark, cell)] else pick_main(cc$meta)[order(-score)]
  sel[, y := if (full) 8.05 - (.I - 1) * .95 else c(8.00, 6.70, 5.40, 4.10)[.I]]

  track_h <- if (full) .62 else .56 * 1.2 * 1.5
  hist_h <- .46 * 1.2 * 1.5

  trk <- rbindlist(lapply(seq_len(nrow(sel)), function(i) track_shape(sel[i], cc$plist, chr, st, en, sel$y[i], track_h)), fill = TRUE)
  hist_tops <- if (r$title == "O shape") c(H3K27ac = 60, H3K4me1 = 40, H3K4me3 = 10) else c(H3K27ac = 120, H3K4me1 = 30, H3K4me3 = 100)
  hs <- if (full) empty_shape() else hist_summary(sel, cc$meta, cc$plist, chr, st, en, hist_h, hist_tops)
  ax <- rbind(axis_table(trk, x0, span * .005), axis_table(hs, x0, span * .005), fill = TRUE)

  label_x <- x0 - span * .020
  lbl <- rbind(data.table(x = label_x, y = c(10.65, 10.15), label = c("GWAS", reg$best_lineage)),
               sel[, .(x = label_x, y = y + track_h * .52, label)],
               fill = TRUE)
  if (!full) lbl <- rbind(lbl, data.table(x = label_x,
                                          y = c(2.85, 1.60, 0.35) + hist_h * .52,
                                          label = c("ENCODE\nH3K27ac", "ENCODE\nH3K4me1", "ENCODE\nH3K4me3")),
                          fill = TRUE)

  lead_lab <- data.table(pos = reg$lead_bp, label = r$lead_label)
  lead_lab[, side := ifelse(pos > x0 + span * .72, -1, 1)]
  lead_lab[, `:=`(x2 = pos + side * span * .026,
                  y2 = 11.08,
                  angle = side * 34)]
  ar_lab <- data.table(x = numeric(), xend = numeric(), y = numeric(), yend = numeric(), label = character())
  if (chr == "X") {
    ar_dist <- max(0, 66763806 - reg$core_end)
    ar_lab <- data.table(x = x1 + span * .012, xend = x1 + span * .070,
                         y = 9.28, yend = 9.28,
                         label = sprintf("AR, %.0f kb downstream", ar_dist / 1000))
  }
  hleg <- data.table(cell = names(cell_col), x = seq(x0 + span * .20, x0 + span * .82, length.out = length(cell_col)), y = 3.90)
  fill_cols <- unique(na.omit(c(trk$col, hs$col)))

  ggplot() +
    geom_segment(data = gwas_dot, aes(pos, 10.66, xend = pos, yend = 10.94), linewidth = .20) +
    geom_point(data = gwas_dot[lead == FALSE], aes(pos, 10.80), color = variant_col, size = 1.85) +
    geom_point(data = gwas_dot[lead == TRUE], aes(pos, 10.80), shape = 23, fill = "#6a3d9a", color = "#6a3d9a", size = 1.95) +
    geom_segment(data = lead_lab, aes(pos, 10.88, xend = x2, yend = y2), color = "#6a3d9a", linewidth = .22) +
    geom_text(data = lead_lab, aes(x2, y2, label = label), color = "#6a3d9a",
              hjust = ifelse(lead_lab$side > 0, 0, 1), vjust = -0.2, angle = lead_lab$angle, size = fs6) +
    geom_segment(data = arch[present == TRUE], aes(pos, y - .14, xend = pos, yend = y + .14), linewidth = .20) +
    geom_point(data = arch[match == TRUE], aes(pos, y), color = variant_col, size = 1.85) +
    geom_text(data = unique(arch[, .(archaic, y)]), aes(x1, y, label = archaic), hjust = 0, size = fs6) +
    geom_segment(data = g, aes(tx0, y, xend = tx1, yend = y), color = "#332b77", linewidth = .55) +
    geom_segment(data = ar, aes(x, y, xend = xend, yend = yend), color = "#332b77", linewidth = .22) +
    geom_rect(data = ex, aes(xmin = start, xmax = end, ymin = y - .06, ymax = y + .06), fill = "#332b77") +
    geom_text(data = g, aes(tx1 + span * .015, y, label = gene), color = "#332b77", hjust = 0, size = fs6, fontface = "italic") +
    geom_segment(data = ar_lab, aes(x, y, xend = xend, yend = yend), color = "#332b77", linewidth = .25,
                 arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
    geom_text(data = ar_lab, aes(xend + span * .010, y, label = label), color = "#332b77",
              hjust = 0, size = fs6, fontface = "italic") +
    geom_ribbon(data = trk, aes(x, ymin = y0, ymax = y, group = grp, fill = col), color = NA, alpha = .82) +
    geom_ribbon(data = hs, aes(x, ymin = y0, ymax = y, group = grp, fill = col), color = NA, alpha = .70) +
    geom_hline(yintercept = c(sel$y, if (!full) c(2.85, 1.60, 0.35)), color = "grey30", linewidth = .22) +
    geom_segment(data = ax, aes(x, y0, xend = x, yend = y_top), linewidth = .20) +
    geom_segment(data = ax, aes(x0, y0, xend = x, yend = y0), linewidth = .20) +
    geom_segment(data = ax, aes(x0, y_mid, xend = x, yend = y_mid), linewidth = .20) +
    geom_segment(data = ax, aes(x0, y_top, xend = x, yend = y_top), linewidth = .20) +
    geom_text(data = ax, aes(xlab, y0, label = "0"), hjust = 1, size = fs6) +
    geom_text(data = ax, aes(xlab, y_top, label = lab_top), hjust = 1, size = fs6) +
    {if (!full) geom_rect(data = hleg, aes(xmin = x - span * .01, xmax = x + span * .01, ymin = y - .08, ymax = y + .08, fill = cell), color = NA)} +
    {if (!full) geom_text(data = hleg, aes(x + span * .017, y, label = cell), hjust = 0, size = fs6)} +
    {if (!full) annotate("text", x = x0 + span * .08, y = 3.90, label = "H3 marks:", hjust = 0, size = fs6)} +
    geom_text(data = lbl, aes(x, y, label = label), hjust = 1, size = fs7) +
    scale_fill_manual(values = c(cell_col, setNames(fill_cols, fill_cols)), guide = "none") +
    scale_x_continuous(labels = function(x) format(round(x), big.mark = ","), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = sprintf("%s | chromosome %s: %s-%s", r$title, chr, format(st, big.mark = ","), format(en, big.mark = ",")),
         x = paste0("chr", chr, ":"), y = NULL) +
    coord_cartesian(xlim = c(x0, x1), ylim = c(if (full) min(sel$y) - .45 else 0, 11.05), clip = "off") +
    theme_classic(base_size = 6) +
    theme(plot.title = element_text(size = 7, hjust = .5),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 6, hjust = 0),
          axis.text.x = element_text(size = 6),
          plot.margin = margin(8, 75, 8, 90))
}

write_ppt <- function(ps, file, full = FALSE) {
  ppt <- officer::read_pptx()
  if (full || length(ps) != 2L) {
    for (p in ps) {
      ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
      ppt <- officer::ph_with(ppt, rvg::dml(ggobj = p), officer::ph_location(left = .15, top = .10, width = 12.9, height = 7.2))
    }
  } else {
    ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
    ppt <- officer::ph_with(ppt, rvg::dml(ggobj = ps[[1]]), officer::ph_location(left = .15, top = .05, width = 12.9 * .8, height = 3.60))
    ppt <- officer::ph_with(ppt, rvg::dml(ggobj = ps[[2]]), officer::ph_location(left = .15, top = 3.85, width = 12.9 * .8, height = 3.60))
  }
  target <- file.path(out, file)
  if (file.exists(target) && !file.remove(target)) target <- sub("\\.pptx$", "_new.pptx", target)
  print(ppt, target = target)
  message("Saved: ", target)
}

source_data <- function(file, selected_only = TRUE) {
  region_tab <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    z <- read_region(r$trait, r$id)
    z[, `:=`(figure_panel = r$title, trait = r$trait, region = r$id, archaics = paste(unlist(r$archaics), collapse = ", "))]
    z
  }), fill = TRUE)

  core_tab <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    reg <- read_region(r$trait, r$id)
    z <- read_core(r$trait, r$id)
    z[, `:=`(figure_panel = r$title, trait = r$trait, region = r$id, lead_snp = reg$lead_snp, lead_chr = reg$lead_chr, lead_bp = reg$lead_bp)]
    z
  }), fill = TRUE)

  arch_tab <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    reg <- read_region(r$trait, r$id)
    core <- read_core(r$trait, r$id)
    z <- read_arch(r$trait, r$id, core, unlist(r$archaics))
    z[, `:=`(figure_panel = r$title, trait = r$trait, region = r$id, lead_snp = reg$lead_snp, lead_chr = reg$lead_chr, lead_bp = reg$lead_bp)]
    z
  }), fill = TRUE)

  gene_tab <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    reg <- read_region(r$trait, r$id)
    chr <- ifelse(reg$lead_chr == 23, "X", as.character(reg$lead_chr))
    z <- read_gene(chr, reg$core_start, reg$core_end)
    if (nrow(z)) z[, `:=`(figure_panel = r$title, trait = r$trait, region = r$id)]
    z
  }), fill = TRUE)

  track_tab <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    reg <- read_region(r$trait, r$id)
    core <- read_core(r$trait, r$id)
    chr <- ifelse(reg$lead_chr == 23, "X", as.character(reg$lead_chr))
    cc <- catalog(chr, reg$core_start, reg$core_end, core)
    cc$meta[, `:=`(figure_panel = r$title, trait = r$trait, region = r$id, lead_snp = reg$lead_snp,
                   selected = id %in% pick_main(cc$meta)$id)]
    cc$meta
  }), fill = TRUE)

  openxlsx::write.xlsx(
    list(region_summary = region_tab,
         core_variants = core_tab,
         archaic_calls = arch_tab,
         gene_models = gene_tab,
         tracks = if (selected_only) track_tab[selected == TRUE] else track_tab),
    file.path(out, file),
    overwrite = TRUE
  )
}

main <- lapply(seq_len(nrow(regions)), function(i) plot_region(regions[i], FALSE))
full <- lapply(seq_len(nrow(regions)), function(i) plot_region(regions[i], TRUE))
write_ppt(main, "Fig3.pptx", FALSE)
write_ppt(full, "Fig3.full.pptx", TRUE)
message(out)
