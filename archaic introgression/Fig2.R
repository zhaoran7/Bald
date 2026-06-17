library(data.table)
library(ape)
library(ggplot2)
library(officer)
library(rvg)
suppressPackageStartupMessages(library(circlize))

root <- "/mnt/f/bald/analysis/archaic/result/locus"
panel_f <- "/mnt/f/refGen/1kg_phase3/integrated_call_samples_v3.20130502.ALL.panel"
out_f <- "/mnt/f/bald/analysis/archaic/plot/Fig2.pptx"
dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)

trait_lab <- c(bald = "Continuous", bald12 = "M shape", bald13 = "O shape", bald14 = "U shape")
sp_col <- c(AFR = "#9D76B1", AMR = "#73B28E", EAS = "#BA4E4F", EUR = "#55728E", SAS = "#CA997A")
lin_col <- c(Neanderthal = "#31859C", Denisovan = "#984807")
trait_col <- c("Continuous" = "#079662", "M shape" = "#cc0404",
               "O shape" = "#1d79b2", "U shape" = "#5817d1")
anc_col <- "#2F6B3F"; ring_col <- "#FCB4A4"; hi_col <- "#31859C"
txt <- 9.4 / 12; count_txt <- 9 / 12; title_txt <- 12 / 12; ideo_txt <- 10 / 12
pt8 <- 8 / 12

read_panel <- function(f) {
  x <- fread(f, fill = TRUE)
  if (all(c("sample", "pop", "super_pop") %in% names(x))) x <- x[, .(sample, pop, super_pop)] else {
    x <- fread(f, header = FALSE, fill = TRUE, select = 1:3)
    setnames(x, c("sample", "pop", "super_pop"))
    x <- x[sample != "sample"]
  }
  x[, sample_i := .I][]
}

read_tr <- function(f) read.tree(text = readLines(f, warn = FALSE)[grep("^\\(", readLines(f, warn = FALSE))[1]])
tips_below <- function(tr, node) {
  k <- tr$edge[tr$edge[, 1] == node, 2]
  unlist(lapply(k, function(i) if (i <= Ntip(tr)) i else tips_below(tr, i)), use.names = FALSE)
}
parent <- function(tr, node) {
  x <- tr$edge[tr$edge[, 2] == node, 1]
  if (length(x)) x[1] else NA_integer_
}
tip_type <- function(meta, labs) {
  x <- rep("modern", length(labs)); names(x) <- labs
  i <- match(labs, meta$label); x[!is.na(i)] <- meta$type[i[!is.na(i)]]
  x[grepl("vindija|altai|chagyr|denisova|neand|archaic", labs, TRUE)] <- "archaic"
  x[labs == "Ancestral"] <- "ancestral"; x
}
span <- function(a) {
  a <- sort((a + 2*pi) %% (2*pi))
  g <- c(diff(a), a[1] + 2*pi - a[length(a)])
  j <- which.max(g); s <- a[(j %% length(a)) + 1]; e <- a[j]
  if (e < s) e <- e + 2*pi
  c(s, e)
}
best_clade <- function(tr, meta, archs) {
  nt <- Ntip(tr); labs <- tr$tip.label; tp <- tip_type(meta, labs)
  archs <- intersect(archs, labs)
  if (!length(archs)) archs <- labs[tp == "archaic"]
  z <- rbindlist(lapply((nt + 1):(nt + tr$Nnode), function(n) {
    t <- labs[tips_below(tr, n)]
    b <- suppressWarnings(as.numeric(tr$node.label[n - nt]))
    if (!all(archs %chin% t) || !any(tp[t] == "modern")) return(NULL)
    data.table(node = n, n = length(t), boot = b, modern = sum(tp[t] == "modern"))
  }), fill = TRUE)
  if (!nrow(z)) return(NA_integer_)
  z[, boot_rank := fifelse(is.finite(boot), boot, -Inf)]
  z[order(n, -boot_rank, -modern)]$node[1]
}
wedge <- function(a0, a1, r0, r1, col) {
  th <- seq(a0, a1, length.out = 18)
  polygon(c(r0*cos(th), rev(r1*cos(th))), c(r0*sin(th), rev(r1*sin(th))), col = col, border = NA)
}

ring_dat <- function(trait0, id0, panel) {
  h <- fread(file.path(root, "report", "hap_match.tsv"))[trait == trait0 & id == id0]
  h[, label := hap_id]
  d <- rbindlist(lapply(seq_len(nrow(h)), function(i) {
    idx <- as.integer(sub("^s([0-9]+)_[12]$", "\\1", unlist(strsplit(h$copies[i], ";", TRUE))))
    panel[sample_i %in% idx, .N, by = super_pop][, .(label = h$label[i], super_pop, prop = N / sum(N))]
  }), fill = TRUE)
  d <- dcast(d, label ~ super_pop, value.var = "prop", fill = 0)
  for (sp in setdiff(names(sp_col), names(d))) d[, (sp) := 0]
  merge(h[, .(label, n)], d, by = "label", all.x = TRUE)
}

tree_rot <- function(trait, id) {
  tr <- read_tr(file.path(root, "phy", trait, paste0(id, ".main.phy_phyml_tree.txt")))
  tr$tip.label <- trimws(tr$tip.label)
  if ("Ancestral" %in% tr$tip.label) tr <- tryCatch(root(tr, "Ancestral", resolve.root = TRUE), error = function(e) tr)
  tr <- ladderize(tr)
  tf <- tempfile(fileext = ".png"); png(tf, 2, 2, units = "in", res = 72)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  plot.phylo(tr, "fan", use.edge.length = TRUE, show.tip.label = FALSE, no.margin = TRUE)
  p <- get("last_plot.phylo", .PlotPhyloEnv); i <- which(tr$tip.label == "Ancestral")
  90 - atan2(p$yy[i], p$xx[i]) * 180 / pi
}

draw_tree <- function(region, panel) {
  trait0 <- region$trait[1]; id0 <- region$id[1]
  tr <- read_tr(file.path(root, "phy", trait0, paste0(id0, ".main.phy_phyml_tree.txt")))
  tr$tip.label <- trimws(tr$tip.label)
  if ("Ancestral" %in% tr$tip.label) tr <- tryCatch(root(tr, "Ancestral", resolve.root = TRUE), error = function(e) tr)
  if ("Ancestral" %in% tr$tip.label && !is.null(tr$edge.length)) {
    ei <- which(tr$edge[, 2] == match("Ancestral", tr$tip.label))
    if (length(ei)) tr$edge.length[ei] <- 0
  }
  tr <- ladderize(tr)
  meta <- fread(file.path(root, "phy", trait0, paste0(id0, ".main.meta.tsv")))
  archs <- unlist(strsplit(region$matched_archaics[1], ";", TRUE)); archs <- archs[nzchar(archs)]
  ring <- ring_dat(trait0, id0, panel)
  arch_col <- lin_col[region$best_lineage[1]]
  if (is.na(arch_col)) arch_col <- lin_col["Denisovan"]

  nt <- Ntip(tr); labs <- tr$tip.label; tp <- tip_type(meta, labs)
  dep <- node.depth.edgelength(tr)[1:nt]; ia <- which(labs == "Ancestral")
  r <- max(dep[setdiff(seq_len(nt), ia)], na.rm = TRUE); if (!is.finite(r)) r <- max(dep, na.rm = TRUE)
  rr <- r * c(lab = 1.06, text = 1.15, pop0 = 1.54, pop1 = 1.68,
              bar0 = 1.83, bar1 = 2.00, lim = 2.18)

  par(mar = c(.08, .08, .50, .08), xpd = NA, pty = "s")
  plot.phylo(tr, "fan", use.edge.length = TRUE, show.tip.label = FALSE, no.margin = TRUE,
             edge.width = .30, rotate.tree = tree_rot(trait0, id0),
             x.lim = c(-rr["lim"], rr["lim"]), y.lim = c(-rr["lim"], rr["lim"]))
  p <- get("last_plot.phylo", .PlotPhyloEnv); a <- atan2(p$yy[1:nt], p$xx[1:nt])
  h <- best_clade(tr, meta, archs)
  clade_labs <- if (!is.na(h)) labs[tips_below(tr, h)] else character()
  clade_modern_n <- sum(tp == "modern" & labs %chin% clade_labs)
  if (!is.na(h)) {
    s <- span(atan2(p$yy[tips_below(tr, h)], p$xx[tips_below(tr, h)]))
    th <- seq(s[1], s[2], length.out = 500); r0 <- sqrt(p$xx[h]^2 + p$yy[h]^2)
    polygon(c(r0*cos(th), rev(rr["lab"]*cos(th))), c(r0*sin(th), rev(rr["lab"]*sin(th))),
            col = adjustcolor(arch_col, .20), border = NA)
  }
  keep <- tp != "ancestral"
  segments(p$xx[seq_len(nt)][keep], p$yy[seq_len(nt)][keep],
           rr["lab"]*cos(a[keep]), rr["lab"]*sin(a[keep]), lty = 3, col = "grey55", lwd = .10)
  symbols(p$xx[seq_len(nt)][keep], p$yy[seq_len(nt)][keep], circles = rep(r*.0023, sum(keep)),
          inches = FALSE, add = TRUE, bg = "black", fg = "black", lwd = .10)
  symbols(0, 0, circles = r*.010, inches = FALSE, add = TRUE, bg = anc_col, fg = anc_col, lwd = .10)

  deg <- a * 180/pi; flip <- deg < -90 | deg > 90
  label_keep <- keep
  for (i in which(label_keep)) {
    col <- if (tp[i] == "archaic") arch_col else if (tp[i] == "ancestral") anc_col else "black"
    label_r <- rr["text"] - if (tp[i] %chin% c("archaic", "ancestral")) .035*r else 0
    text(label_r*cos(a[i]), label_r*sin(a[i]), labs[i],
         srt = ifelse(flip[i], deg[i] + 180, deg[i]), adj = if (flip[i]) c(1, .5) else c(0, .5),
         cex = txt, col = col)
  }
  anc_a <- pi/2
  segments(0, 0, rr["lab"] * cos(anc_a), rr["lab"] * sin(anc_a), lty = 3, col = "grey55", lwd = .16)
  text((rr["text"] - .035*r) * cos(anc_a), (rr["text"] - .035*r) * sin(anc_a), "Ancestral",
       srt = 90, adj = c(0, .5), cex = txt, col = anc_col)
  if (!is.na(h)) {
    b <- suppressWarnings(as.numeric(tr$node.label)); nd <- nt + seq_len(tr$Nnode)
    kids <- tr$edge[tr$edge[, 1] == h, 2]
    arch_kids <- kids[kids > nt & vapply(kids, function(k) any(labs[tips_below(tr, k)] %chin% archs), logical(1))]
    j <- match(unique(c(h, arch_kids)), nd); j <- j[!is.na(j) & !is.na(b[j])]
    if (length(j)) nodelabels(b[j], node = nd[j], frame = "n", cex = txt)
  }

  for (x in rr[c("pop0","pop1","bar0","bar1")]) {
    th <- seq(0, 2*pi, length.out = 720); lines(x*cos(th), x*sin(th), col = "#D2D2D2", lwd = .18)
  }
  max_n <- max(ring$n, na.rm = TRUE); da <- 2*pi/nt*.36
  for (i in seq_len(nt)) {
    z <- ring[label == labs[i]][1]
    if (tp[i] == "modern" && nrow(z)) {
      r0 <- rr["pop0"]
      for (sp in names(sp_col)) {
        psp <- as.numeric(z[[sp]]); if (is.na(psp)) psp <- 0
        if (psp > 0) wedge(a[i]-da, a[i]+da, r0, r0 + psp*(rr["pop1"]-rr["pop0"]), sp_col[sp])
        r0 <- r0 + psp*(rr["pop1"]-rr["pop0"])
      }
      h1 <- rr["bar0"] + z$n/max_n*(rr["bar1"]-rr["bar0"])
      wedge(a[i]-da, a[i]+da, rr["bar0"], h1, ring_col)
      if (labs[i] %chin% clade_labs) {
        text((h1 + .035*r)*cos(a[i]), (h1 + .035*r)*sin(a[i]), z$n,
             srt = ifelse(flip[i], deg[i] + 180, deg[i]),
             adj = if (flip[i]) c(1, .5) else c(0, .5), cex = count_txt)
      }
    } else if (tp[i] == "archaic") {
      cx <- mean(rr[c("pop0","pop1")]) * cos(a[i]); cy <- mean(rr[c("pop0","pop1")]) * sin(a[i])
      symbols(cx, cy, circles = r*.014, inches = FALSE, add = TRUE, bg = arch_col, fg = "#666666", lwd = .10)
    }
  }
  title(sprintf("%s | chr%s:%s-%s", region$trait_label[1], region$lead_chr[1],
                format(region$core_start[1], big.mark = ","), format(region$core_end[1], big.mark = ",")),
        cex.main = title_txt, line = -.78, font.main = 2)
}

draw_ideogram <- function(regions) {
  cyto <- as.data.table(read.cytoband(species = "hg19")[[1]])
  setnames(cyto, 1:5, c("chr", "start", "end", "name", "gieStain"))
  chr_order <- paste0("chr", c(1:22, "X"))
  cyto <- cyto[chr %chin% chr_order]
  cyto[, chr := factor(chr, levels = chr_order)]
  setorder(cyto, chr, start)
  chr_df <- cyto[, .(chr_len = max(end)), by = chr]
  chr_df[, chr := as.character(chr)]
  stain_col <- c(gneg = "#FFFFFF", gpos25 = "#D9D9D9", gpos50 = "#A6A6A6",
                 gpos75 = "#6E6E6E", gpos100 = "#000000", gvar = "#CFCFCF",
                 stalk = "#BDBDBD", acen = "#D95F02")
  cyto[, fill := stain_col[gieStain]]
  cyto[is.na(fill), fill := "#FFFFFF"]
  reg <- copy(regions)
  reg[, chr := paste0("chr", fifelse(lead_chr == 23, "X", as.character(lead_chr)))]
  reg[, lineage_col := fifelse(best_lineage == "Neanderthal", lin_col["Neanderthal"], lin_col["Denisovan"])]

  par(mai = c(0, 0, 0, 0), xpd = NA)
  circos.clear()
  on.exit(circos.clear(), add = TRUE)
  circos.par(start.degree = 88, gap.degree = c(rep(0.7, 22), 10),
             cell.padding = c(0, 0, 0, 0), track.margin = c(0.001, 0.001),
             points.overflow.warning = FALSE, canvas.xlim = c(-1.10, 1.10),
             canvas.ylim = c(-1.10, 1.10))
  circos.initialize(chr_df$chr, xlim = cbind(0, chr_df$chr_len))
  circos.trackPlotRegion(ylim = c(0, 1.55), track.height = 0.105, bg.border = NA,
                         panel.fun = function(x, y) {
    s <- CELL_META$sector.index
    d <- cyto[as.character(chr) == s]
    for (i in seq_len(nrow(d))) {
      if (d$gieStain[i] == "acen") {
        k <- sum(d$gieStain[seq_len(i)] == "acen")
        xx <- if (k %% 2) c(d$start[i], d$end[i], d$start[i]) else c(d$end[i], d$start[i], d$end[i])
        circos.polygon(xx, c(0.10, 0.53, 0.96), col = d$fill[i], border = NA)
      } else {
        circos.rect(d$start[i], 0.10, d$end[i], 0.96, col = d$fill[i], border = NA)
      }
    }
    circos.rect(CELL_META$xlim[1], 0.10, CELL_META$xlim[2], 0.96, col = NA, border = "#555555", lwd = 0.25)
    circos.text(CELL_META$xcenter, 1.32, sub("chr", "", CELL_META$sector.index),
                facing = "bending.inside", niceFacing = TRUE, cex = ideo_txt, col = "black")
  })
  circos.trackPlotRegion(ylim = c(-0.45, 1), track.height = 0.285, bg.border = NA,
                         panel.fun = function(x, y) {
    z <- reg[chr == CELL_META$sector.index]
    if (nrow(z)) {
      for (i in seq_len(nrow(z))) {
        y1 <- 0.08 - 0.28*((i - 1) %% 3)
        circos.segments(z$lead_bp[i], 1.00, z$lead_bp[i], y1 + 0.08, col = z$lineage_col[i], lwd = 0.8)
        circos.points(z$lead_bp[i], y1 + 0.04, pch = 16, cex = 0.35, col = z$lineage_col[i])
        circos.text(z$lead_bp[i], y1, z$lead_snp[i], facing = "clockwise",
                    niceFacing = TRUE, cex = ideo_txt, col = z$lineage_col[i])
      }
    }
  })
  legend("center", c("Neanderthal-like", "Denisovan-like"), lty = 1, lwd = 1.4,
         col = lin_col[c("Neanderthal", "Denisovan")], bty = "n", cex = pt8, y.intersp = 1.25)
  title("Inherited candidate regions", cex.main = title_txt, line = -0.6, font.main = 2)
}

freq_table <- function(regions) {
  core <- fread(file.path(root, "report", "core_risk.tsv"),
                select = c("trait", "id", "pos", "diagnostic_archaic_allele", "is_diagnostic_archaic"))
  one_region <- function(r) {
    cr <- core[trait == r$trait & id == r$id & is_diagnostic_archaic == TRUE &
                 !is.na(diagnostic_archaic_allele) & diagnostic_archaic_allele %chin% c("A","C","G","T"),
               .(pos, diagnostic_archaic_allele)]
    if (!nrow(cr)) return(NULL)
    kg <- fread(file.path(root, "mat", r$trait, r$id, "kg.tsv"), header = FALSE)
    sample_cols <- paste0("s", seq_len(ncol(kg) - 5))
    setnames(kg, c("chr", "pos", "ref", "alt", "aa", sample_cols))
    kg <- merge(kg, cr, by = "pos")
    if (!nrow(kg)) return(NULL)
    samp <- copy(panel[seq_along(sample_cols)])
    samp[, col := sample_cols]
    pop_den <- samp[, .(den = 2L * .N), by = .(pop, super_pop)]
    gt_to_count <- function(gt, ref, alt, diag) {
      sp <- strsplit(gsub("\\|", "/", gt), "/", fixed = TRUE)
      a1 <- vapply(sp, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
      a2 <- vapply(sp, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
      a1 <- fifelse(a1 == "0", ref, fifelse(a1 == "1", alt, NA_character_))
      a2 <- fifelse(a2 == "0", ref, fifelse(a2 == "1", alt, NA_character_))
      as.integer(a1 == diag) + as.integer(a2 == diag)
    }
    site_freq <- rbindlist(lapply(seq_len(nrow(kg)), function(i) {
      cnt <- vapply(samp$col, function(cc) gt_to_count(kg[[cc]][i], kg$ref[i], kg$alt[i], kg$diagnostic_archaic_allele[i]),
                    integer(1))
      z <- copy(samp)
      z[, cnt := cnt]
      z[, .(num = sum(cnt)), by = .(pop, super_pop)][pop_den, on = c("pop","super_pop")][,
        .(region = r$freq_region, pop, super_pop, freq = fifelse(den > 0, num / den, NA_real_))]
    }))
    site_freq[, .(freq = max(freq, na.rm = TRUE)), by = .(region, pop, super_pop)]
  }
  d <- rbindlist(lapply(seq_len(nrow(regions)), function(i) one_region(regions[i])), fill = TRUE)
  d[, set := "All individuals"]
  d[]
}

draw_freq_heatmap <- function(regions, panel) {
  heat_txt <- 10 / ggplot2::.pt
  d <- freq_table(regions)
  fwrite(d, file.path(root, "plot", "Fig2.frequency.tsv"), sep = "\t")
  pop_order <- c("ACB","ASW","ESN","GWD","LWK","MSL","YRI",
                 "CLM","MXL","PEL","PUR",
                 "CDX","CHB","CHS","JPT","KHV",
                 "CEU","FIN","GBR","IBS","TSI",
                 "BEB","GIH","ITU","PJL","STU")
  pop_dt <- unique(d[pop %chin% pop_order, .(pop, super_pop)])
  pop_dt[, super_pop := factor(super_pop, levels = names(sp_col))]
  pop_dt[, pop := factor(pop, levels = pop_order)]
  setorder(pop_dt, super_pop, pop)
  pop_dt[, x := .I]
  y_dt <- data.table(region = rev(regions$freq_region), y = seq_len(nrow(regions)))
  d <- merge(d, pop_dt, by = c("pop", "super_pop"))
  d <- merge(d, y_dt, by = "region")
  grp <- pop_dt[, .(xmin = min(x) - .5, xmax = max(x) + .5), by = super_pop]
  y_top0 <- nrow(regions) + .62
  y_top1 <- nrow(regions) + .94
  ggplot() +
    geom_rect(data = grp, aes(xmin = xmin, xmax = xmax, ymin = y_top0, ymax = y_top1),
              inherit.aes = FALSE, fill = sp_col[as.character(grp$super_pop)], color = "white", linewidth = .25) +
    geom_text(data = grp, aes(x = (xmin + xmax)/2, y = (y_top0 + y_top1)/2, label = super_pop),
              inherit.aes = FALSE, family = "Arial", size = heat_txt, color = "white") +
    geom_tile(data = d, aes(x = x, y = y, fill = freq), color = "white", linewidth = .30, width = .98, height = .88) +
    geom_text(data = d, aes(x = x, y = y, label = ifelse(is.na(freq), "NA", sprintf("%.1f", 100*freq))),
              family = "Arial", size = heat_txt, color = "#202020") +
    scale_fill_gradientn(colours = c("#F8FBFF", "#E6F2F8", "#C7E1EE", "#9EC7DF", "#6BAED6"),
                         limits = c(0, max(d$freq, na.rm = TRUE)),
                         labels = function(x) sprintf("%.0f", 100*x),
                         name = "Diagnostic\nhaplotype (%)") +
    scale_x_continuous(breaks = pop_dt$x, labels = as.character(pop_dt$pop),
                       expand = expansion(mult = c(.005, .005))) +
    scale_y_continuous(breaks = y_dt$y, labels = y_dt$region,
                       limits = c(.5, y_top1 + .02), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial", base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
          axis.text.y = element_text(size = 10), axis.title = element_blank(),
          legend.title = element_text(size = 8), legend.text = element_text(size = 8),
          legend.position = "right", plot.margin = margin(2, 2, 2, 2))
}

draw_region_schematic <- function(regions) {
  lab_size <- 10 / .pt
  fmt_bp <- function(x) format(round(x), big.mark = ",", scientific = FALSE)
  fmt_p <- function(x) ifelse(x == 0, "0", formatC(x, format = "e", digits = 2))
  r <- copy(regions)
  r[, chr := fifelse(lead_chr == 23, "X", as.character(lead_chr))]
  r[, hap := paste0("Haplotype ", .I)]
  r[, display_snp := lead_snp]
  r[chr == "X", display_snp := "rs143054933 / rs5919284 / rs6525167"]
  r[, coord_lab := paste0("chr", chr, ": ", fmt_bp(core_start), "-", fmt_bp(core_end))]
  r[, info_lab := paste0(n_ld_snp, " SNPs\n", sprintf("%.1f kb", (core_end-core_start+1)/1000), "\n", "p-ILS=", fmt_p(p_ils))]
  r[, lineage_col := lin_col[best_lineage]]
  r[, x := seq(9, 91, length.out = .N)]
  r[, y := c(36, 57, 38, 53, 51)[seq_len(.N)]]
  r[, zoom_x := x + 6.0]
  r[, zoom_y := y + c(-1.0, 7.5, -3.2, -3.2, 6.6)[seq_len(.N)]]
  r[, zoom_h := 10 * ((core_end-core_start+1)/max(core_end-core_start+1))^.65 + 3.8]
  r[, rs_y := y + 7.2]
  r[, coord_y := y + 2.5]

  chr_order <- unique(r$chr)
  cyto <- as.data.table(read.cytoband(species = "hg19")[[1]])
  setnames(cyto, c("chrom", "start", "end", "band", "stain"))
  cyto[, chr := sub("^chr", "", chrom)]
  cyto <- cyto[chr %chin% chr_order]
  chr_len <- cyto[, .(chr_len = max(end)), by = chr]
  chr_len[, h := 51 * (chr_len / max(chr_len))^.90]
  r <- merge(r, chr_len[, .(chr, chr_len, h)], by = "chr", all.x = TRUE)
  cyto <- merge(cyto, chr_len, by = "chr")
  cyto <- merge(cyto, r[, .(chr, x, y)], by = "chr")
  cyto[, ymin := y - h/2 + start/chr_len*h]
  cyto[, ymax := y - h/2 + end/chr_len*h]
  band_col <- c(gneg="white", gpos25="#D9D9D9", gpos50="#9A9A9A", gpos75="#4F4F4F",
                gpos100="#111111", acen="#314EAE", gvar="#EFEFEF", stalk="#C7C7C7")
  cyto[, fill := band_col[stain]]
  cyto[is.na(fill), fill := "#EEEEEE"]

  dots <- r[, .(trait=unlist(strsplit(trait_label, "; ", fixed=TRUE))), by=.(chr, x, y, display_snp)]
  dots[, dot_x := x - 2.7]
  dots[, dot_y := y + seq(-1.8, 1.8, length.out=.N), by=.(chr, display_snp)]

  leg_l <- data.table(x=c(42,66), y=7.3, lab=c("Neanderthal-like", "Denisovan-like"), lineage=c("Neanderthal","Denisovan"))
  leg_t <- data.table(x=c(38,53,66,79), y=4.3, trait=c("Continuous", "M shape", "O shape", "U shape"))

  ggplot() +
    geom_segment(data=cyto, aes(x=x, xend=x, y=ymin, yend=ymax), linewidth=5.8, color="#9A9A9A", lineend="butt") +
    geom_segment(data=cyto, aes(x=x, xend=x, y=ymin, yend=ymax), linewidth=4.9, color=cyto$fill, lineend="butt") +
    geom_segment(data=r, aes(x=x, xend=x, y=y-h/2, yend=y+h/2), linewidth=.28, color="#333333", lineend="round") +
    geom_text(data=r, aes(x=x, y=y+h/2+3.3, label=paste0("chr", chr)), family="Arial", size=lab_size) +
    geom_segment(data=r, aes(x=x-2.0, xend=x+2.0, y=y, yend=y, color=best_lineage), linewidth=.55, lineend="round") +
    geom_point(data=r, aes(x=x, y=y, fill=best_lineage), shape=21, size=1.7, color="white", stroke=.16) +
    geom_polygon(data=r[, .(x=c(x+2.3, zoom_x-.85, zoom_x-.85), y=c(y, zoom_y-zoom_h/2, zoom_y+zoom_h/2), grp=.GRP, best_lineage=best_lineage), by=seq_len(nrow(r))],
                 aes(x=x, y=y, group=grp, fill=best_lineage), alpha=.13, color=NA) +
    geom_segment(data=r, aes(x=zoom_x, xend=zoom_x, y=zoom_y-zoom_h/2, yend=zoom_y+zoom_h/2, color=best_lineage), linewidth=4.8, lineend="round") +
    geom_segment(data=r, aes(x=zoom_x, xend=zoom_x, y=zoom_y-zoom_h/2, yend=zoom_y+zoom_h/2), linewidth=.8, color="white", alpha=.35, lineend="round") +
    geom_text(data=r, aes(x=x-3.0, y=rs_y, label=display_snp), hjust=1, family="Arial", size=lab_size) +
    geom_text(data=r, aes(x=x-3.0, y=coord_y, label=coord_lab), hjust=1, family="Arial", size=lab_size, lineheight=.86) +
    geom_text(data=r, aes(x=zoom_x+1.1, y=zoom_y+zoom_h/2+2.0, label=hap), hjust=0, family="Arial", fontface="bold", size=lab_size) +
    geom_text(data=r, aes(x=zoom_x+1.1, y=zoom_y, label=info_lab, color=best_lineage), hjust=0, family="Arial", size=lab_size, lineheight=.86, parse=FALSE) +
    geom_point(data=dots, aes(dot_x, dot_y, fill=trait), shape=21, size=1.75, color="white", stroke=.16) +
    geom_segment(data=leg_l, aes(x=x, xend=x+6.5, y=y, yend=y, color=lineage), linewidth=2.5, lineend="round") +
    geom_text(data=leg_l, aes(x=x+7.2, y=y, label=lab), hjust=0, family="Arial", size=lab_size) +
    geom_point(data=leg_t, aes(x=x, y=y, fill=trait), shape=21, size=1.8, color="white", stroke=.16) +
    geom_text(data=leg_t, aes(x=x+1.1, y=y, label=trait), hjust=0, family="Arial", size=lab_size) +
    scale_color_manual(values=lin_col, guide="none") +
    scale_fill_manual(values=c(lin_col, trait_col), guide="none") +
    coord_cartesian(xlim=c(0, 103), ylim=c(0, 75), clip="off") +
    labs(title="Candidate archaic-like MPB haplotypes") +
    theme_void(base_family="Arial") +
    theme(plot.title=element_text(size=10, face="plain", hjust=.02), plot.margin=margin(0, 1, 0, 1))
}

set_size <- function(pptx, w, h) {
  td <- tempfile(); dir.create(td); unzip(pptx, exdir = td)
  f <- file.path(td, "ppt", "presentation.xml")
  x <- gsub('<p:sldSz cx="[0-9]+" cy="[0-9]+"[^>]*/>',
            sprintf('<p:sldSz cx="%d" cy="%d" type="custom"/>', round(w*914400), round(h*914400)),
            readLines(f, warn = FALSE))
  writeLines(x, f, useBytes = TRUE)
  old <- setwd(td); on.exit(setwd(old), add = TRUE)
  tmp <- tempfile(fileext = ".pptx")
  system2("zip", c("-qr9X", tmp, list.files(td, recursive = TRUE, all.files = TRUE, no.. = TRUE)))
  file.copy(tmp, pptx, TRUE)
  unlink(td, TRUE)
}
panel_label <- function(ppt, lab, x, y) {
  ph_with(ppt, fpar(ftext(lab, prop = fp_text(font.size = 12, bold = TRUE, font.family = "Arial")),
                    fp_p = fp_par(text.align = "right")),
          ph_location(x, y, .45, .35))
}

panel <- read_panel(panel_f)
sel <- fread(file.path(root, "report", "selected_region.tsv"))
sel[, trait_name := unname(trait_lab[trait])]
pick_region <- function(tr, id0, lab = NULL, merged_label = NULL) {
  x <- copy(sel[trait == tr & id == id0][1])
  if (!nrow(x)) stop("missing selected region: ", tr, " ", id0)
  x[, trait_label := if (is.null(lab)) trait_name else lab]
  x[, freq_region := if (is.null(merged_label)) {
    sprintf("%s | chr%s:%s-%s", trait_label, fifelse(lead_chr == 23, "X", as.character(lead_chr)),
            format(core_start, big.mark = ","), format(core_end, big.mark = ","))
  } else merged_label]
  x
}
regions <- rbindlist(list(
  pick_region("bald13", "1.rs12405323.170418964", "O shape"),
  pick_region("bald12", "6.rs9349320.45269814", "M shape"),
  pick_region("bald", "8.rs1041791.109112070", "Continuous; U shape"),
  pick_region("bald", "12.rs417915.52842074", "Continuous"),
  pick_region("bald13", "23.rs6525167.66423003", "Continuous; M shape; O shape",
              "Continuous; M shape; O shape | chrX:66,186,241-66,423,003")
), fill = TRUE)
regions[, region_label := sprintf("%s | chr%s:%s-%s | %s", trait_label,
                                  fifelse(lead_chr == 23, "X", as.character(lead_chr)),
                                  format(core_start, big.mark = ","),
                                  format(core_end, big.mark = ","),
                                  best_lineage)]

ppt <- read_pptx()
slide_w <- 28.0; slide_h <- 27.0
cell_w <- 7.85; cell_h <- 7.85
gap_x <- .34; gap_y <- .34
x0 <- (slide_w - 3 * cell_w - 2 * gap_x) / 2
top_y <- .32; top_h <- 4.05
y0 <- top_y + top_h + .48
panel_w <- 3 * cell_w + 2 * gap_x
bot_y <- y0 + 2 * cell_h + gap_y + .46
bot_h <- 5.35
ppt <- add_slide(ppt, "Blank", "Office Theme")
ppt <- panel_label(ppt, "A", x0 - .28, top_y - .08)
ppt <- ph_with(ppt, dml(ggobj = draw_region_schematic(regions)), ph_location(x0, top_y, panel_w, top_h))
ppt <- panel_label(ppt, "B", x0 - .28, y0 - .08)
locs <- data.table(
  i = 1:5,
  x = c(x0 + (0:2)*(cell_w + gap_x),
        x0 + .5*(cell_w + gap_x) + (0:1)*(cell_w + gap_x)),
  y = c(rep(y0, 3), rep(y0 + cell_h + gap_y, 2))
)
for (i in seq_len(nrow(regions))) {
  ppt <- ph_with(ppt, dml(code = draw_tree(regions[i], panel)), ph_location(locs$x[i], locs$y[i], cell_w, cell_h))
}
ppt <- panel_label(ppt, "C", x0 - .28, bot_y - .08)
ppt <- ph_with(ppt, dml(ggobj = draw_freq_heatmap(regions, panel)), ph_location(x0, bot_y, panel_w, bot_h))
print(ppt, target = out_f)
set_size(out_f, slide_w, slide_h)
message("Saved: ", out_f)
