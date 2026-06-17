pacman::p_load(data.table, circlize, GenomicRanges, AnnotationDbi,
  TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, VennDiagram, grid)

# path
out  <- "D:/bald/res/Fig. 1.pdf"
dir0 <- "D:/bald"
dirg <- file.path(dir0, "data/gwas")
dirc <- file.path(dir0, "data/cojo")

# const
traits <- c("bald.qt","bald12.bt","bald13.bt","bald14.bt")
track_lab <- c("bald.qt"="Continuous","bald12.bt"="M shape","bald13.bt"="O shape","bald14.bt"="U shape")

cut_bg <- 1e-3; sig_cut <- 5e-8; max_log10 <- 50

col_ns  <- "#f2b21a"; col_sig <- "#9e423c"
hi_col <- c(shared_all="black", specific_bald.qt="#079662", specific_bald12.bt="#cc0404", specific_bald13.bt="#1d79b2", specific_bald14.bt="#5817d1")
trait_col <- c( "bald.qt" = hi_col["specific_bald.qt"], "bald12.bt" = hi_col["specific_bald12.bt"], "bald13.bt" = hi_col["specific_bald13.bt"], "bald14.bt" = hi_col["specific_bald14.bt"])

chr_df <- data.table(
  CHR = 1:23,
  chr_len = c(
    248956422,242193529,198295559,190214555,181538259,170805979,
    159345973,145138636,138394717,133797422,135086622,133275309,
    114364328,107043718,101991189,90338345,83257441,80373285,
    58617616,64444167,46709983,50818468,156040895
  )
)
chr_df[, chr := paste0("chr", fifelse(CHR == 23L, "X", as.character(CHR)))]
chr_map <- setNames(chr_df$chr_len, chr_df$CHR)

# function
cap_p <- function(p) {
  p <- as.numeric(p)
  p[!is.na(p) & p <= 0] <- 10^(-max_log10)
  p[p < 10^(-max_log10)] <- 10^(-max_log10)
  p
}

parse_chr_pos <- function(s) {
  s <- as.character(s)
  chr0 <- sub(":.*", "", s)
  pos0 <- sub("^[^:]+:([0-9]+).*", "\\1", s)
  chr_num <- suppressWarnings(as.integer(chr0))
  chr_num[chr0 %chin% c("X","x","23")] <- 23L
  data.table(CHR = chr_num, POS = suppressWarnings(as.integer(pos0)))
}

getv <- function(x, nm, default = NA) if (nm %in% names(x)) x[[nm]] else rep(default, nrow(x))

# gene
gene_db <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
gene_db <- gene_db[as.character(seqnames(gene_db)) %in% paste0("chr", c(1:22, "X"))]

map_gene_label <- function(entrez) {
  m <- as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(as.character(entrez)),
    keytype = "ENTREZID",
    columns = c("SYMBOL","ALIAS","GENENAME")
  )))
  m[, label := fifelse(
    !is.na(SYMBOL) & SYMBOL != "", SYMBOL,
    fifelse(!is.na(ALIAS) & ALIAS != "", ALIAS,
            fifelse(!is.na(GENENAME) & GENENAME != "", GENENAME, paste0("ENTREZ:", ENTREZID)))
  )]
  m <- m[, .(label = label[1L]), by = ENTREZID]
  setNames(m$label, m$ENTREZID)
}

annot_gene <- function(dt) {
  qgr <- GRanges(
    seqnames = paste0("chr", fifelse(dt$CHR == 23L, "X", as.character(dt$CHR))),
    ranges = IRanges(dt$POS, dt$POS)
  )
  hit <- distanceToNearest(qgr, gene_db, ignore.strand = TRUE)
  out <- copy(dt[queryHits(hit)])
  out[, ENTREZID := as.character(mcols(gene_db)$gene_id[subjectHits(hit)])]
  lab <- map_gene_label(out$ENTREZID)
  out[, gene := unname(lab[ENTREZID])]
  out[is.na(gene) | gene == "", gene := paste0("ENTREZ:", ENTREZID)]
  out[, .(pos_id, CHR, POS, gene)]
}

# read
read_gwas <- function(tr) {
  x <- fread(file.path(dirg, paste0(tr, ".gz")), showProgress = FALSE)
  y <- parse_chr_pos(x$SNP)
  z <- data.table(
    trait = tr,
    SNP   = as.character(x$SNP),
    CHR   = fcoalesce(fifelse(as.character(x$CHR) %chin% c("X","x","23"), 23L, suppressWarnings(as.integer(as.character(x$CHR)))), y$CHR),
    POS   = fcoalesce(suppressWarnings(as.integer(x$POS)), y$POS),
    EAF   = as.numeric(x$EAF),
    P_raw = as.numeric(x$P)
  )
  z <- z[CHR %between% c(1L,23L) & !is.na(POS) & !is.na(P_raw)]
  z <- z[EAF > 0.005 & EAF < 0.995]
  z <- z[P_raw < cut_bg]
  z[, LOGP := -log10(cap_p(P_raw))]
  z[, chr := paste0("chr", fifelse(CHR == 23L, "X", as.character(CHR)))]
  z[, .(trait, SNP, CHR, POS, P_raw, LOGP, chr)]
}

read_cojo <- function(tr) {
  x <- fread(file.path(dirc, tr, "jma.cojo"), fill = TRUE, quote = "", showProgress = FALSE)
  y <- parse_chr_pos(x$SNP)
  chr_col <- names(x)[names(x) %chin% c("Chr","CHR")][1]
  pos_col <- names(x)[names(x) %chin% c("bp","BP")][1]
  chr_raw <- as.character(x[[chr_col]])
  chr_num <- suppressWarnings(as.integer(chr_raw))
  chr_num[chr_raw %chin% c("X","x","23")] <- 23L
  
  z <- data.table(
    trait = tr,
    SNP   = as.character(x$SNP),
    CHR   = fcoalesce(chr_num, y$CHR),
    POS   = fcoalesce(suppressWarnings(as.integer(x[[pos_col]])), y$POS),
    p     = as.numeric(x$p),
    pJ    = as.numeric(x$pJ),
    b     = as.numeric(x$b),
    se    = as.numeric(x$se),
    bJ    = as.numeric(x$bJ),
    bJ_se = as.numeric(x$bJ_se)
  )
  z <- z[CHR %between% c(1L,23L) & !is.na(POS) & !is.na(pJ)]
  z[, LOGPJ := -log10(cap_p(pJ))]
  z[, chr := paste0("chr", fifelse(CHR == 23L, "X", as.character(CHR)))]
  z[, pos_id := paste0("chr", fifelse(CHR == 23L, "X", as.character(CHR)), ":", POS)]
  z[, .(trait, SNP, CHR, POS, p, pJ, LOGPJ, b, se, bJ, bJ_se, chr, pos_id)]
}

# venn
make_venn_grob <- function(c0) {
  s1 <- unique(c0[trait == "bald.qt",   pos_id])
  s2 <- unique(c0[trait == "bald12.bt", pos_id])
  s3 <- unique(c0[trait == "bald13.bt", pos_id])
  s4 <- unique(c0[trait == "bald14.bt", pos_id])
  
  u <- unique(c(s1, s2, s3, s4))
  dt <- data.table(id = u, A = u %chin% s1, B = u %chin% s2, C = u %chin% s3, D = u %chin% s4)
  
  grobTree(VennDiagram::draw.quad.venn(
    area1 = sum(dt$A), area2 = sum(dt$B), area3 = sum(dt$C), area4 = sum(dt$D),
    n12 = sum(dt$A & dt$B), n13 = sum(dt$A & dt$C), n14 = sum(dt$A & dt$D),
    n23 = sum(dt$B & dt$C), n24 = sum(dt$B & dt$D), n34 = sum(dt$C & dt$D),
    n123 = sum(dt$A & dt$B & dt$C), n124 = sum(dt$A & dt$B & dt$D),
    n134 = sum(dt$A & dt$C & dt$D), n234 = sum(dt$B & dt$C & dt$D),
    n1234 = sum(dt$A & dt$B & dt$C & dt$D),
    category = unname(track_lab[traits]),
    fill = c(hi_col["specific_bald.qt"], hi_col["specific_bald12.bt"], hi_col["specific_bald13.bt"], hi_col["specific_bald14.bt"]),
    alpha = rep(0.2, 4),
    col = NA, lwd = 0, lty = "blank",
    cex = 0.50, fontfamily = "sans", fontface = "plain",
    cat.cex = 0.50,
    cat.fontfamily = "sans", cat.fontface = "plain",
    cat.col = c(hi_col["specific_bald.qt"], hi_col["specific_bald12.bt"], hi_col["specific_bald13.bt"], hi_col["specific_bald14.bt"]),
    margin = 0.01,
    ind = FALSE
  ))
}

# legend
draw_legend <- function() {
  pushViewport(viewport(x = 0.75, y = 0.05, width = 0.15, height = 0.1, just = c("left","bottom")))
  grid.text("Labels", x = 0.00, y = 0.92, just = "left", gp = gpar(fontsize = 6))
  yy <- c(0.80, 0.67, 0.54, 0.41, 0.28)
  nm <- c("shared_all","specific_bald.qt","specific_bald12.bt","specific_bald13.bt","specific_bald14.bt")
  tx <- c("Shared","Continuous","M shape","O shape","U shape")
  for (i in 1:5) {
    grid.points(0.05, yy[i], pch = 4, size = unit(1.8, "mm"), gp = gpar(col = hi_col[nm[i]], lwd = 1))
    grid.text(tx[i], x = 0.11, y = yy[i], just = "left", gp = gpar(fontsize = 5.5, col = hi_col[nm[i]]))
  }
  grid.text("GWAS", x = 0.56, y = 0.92, just = "left", gp = gpar(fontsize = 6))
  grid.points(0.6, 0.78, pch = 16, size = unit(1.8, "mm"), gp = gpar(col = col_ns, fill = NA, lwd = 1))
  grid.text(expression(italic(p) > 5 %*% 10^-8), x = 0.66, y = 0.8, just = "left", gp = gpar(fontsize = 5.5))
  grid.points(0.6, 0.63, pch = 16, size = unit(1.8, "mm"), gp = gpar(col = col_sig))
  grid.text(expression(italic(p) <= 5 %*% 10^-8), x = 0.66, y = 0.65, just = "left", gp = gpar(fontsize = 5.5))
  popViewport()
}

# gwas
# gwas_list <- setNames(lapply(traits, read_gwas), traits)
# gwas_all  <- rbindlist(gwas_list)

# cojo
cojo_list <- setNames(lapply(traits, read_cojo), traits)
cojo_all  <- rbindlist(cojo_list)

# main
g <- copy(gwas_all)[trait %chin% traits & CHR %in% 1:23,
                    .(trait = as.character(trait), CHR = as.integer(CHR), POS = as.integer(POS), P_raw = as.numeric(P_raw), LOGP = as.numeric(LOGP))]
g <- g[!is.na(POS) & !is.na(P_raw) & !is.na(LOGP)]
g[, L := chr_map[as.character(CHR)]][POS >= 1 & POS <= L]
g[, `:=`(chr = paste0("chr", fifelse(CHR == 23L, "X", as.character(CHR))), logp = pmin(LOGP, max_log10))][, L := NULL]

c0 <- copy(cojo_all)[trait %chin% traits & CHR %in% 1:23,
                     .(trait = as.character(trait), SNP = as.character(SNP), CHR = as.integer(CHR), POS = as.integer(POS), pJ = as.numeric(pJ), chr, pos_id)]
c0 <- c0[!is.na(POS) & !is.na(pJ) & pJ < sig_cut]
c0[, L := chr_map[as.character(CHR)]][POS >= 1 & POS <= L][, L := NULL]

ann <- annot_gene(unique(c0[, .(pos_id, CHR, POS)]))
x <- merge(c0[, .(trait, SNP, pos_id, CHR, POS, chr, pJ)], ann, by = c("pos_id","CHR","POS"), all.x = TRUE)

pres <- dcast(unique(x[, .(pos_id, trait)]), pos_id ~ trait, value.var = "trait", fun.aggregate = length, fill = 0)
for (v in traits) if (!v %in% names(pres)) pres[, (v) := 0L]

site0 <- x[, .(
  SNP = SNP[which.min(pJ)],
  CHR = CHR[1L],
  POS = POS[1L],
  chr = chr[1L],
  gene = gene[1L],
  locus_pJ = min(pJ, na.rm = TRUE)
), by = pos_id]

site <- merge(site0, pres, by = "pos_id", all.x = TRUE)

site[, n_trait := rowSums(.SD), .SDcols = traits]
site[, cls := fifelse(
  n_trait == 4, "shared_all",
  fifelse(n_trait == 1 & `bald.qt`   > 0, "specific_bald.qt",
          fifelse(n_trait == 1 & `bald12.bt` > 0, "specific_bald12.bt",
                  fifelse(n_trait == 1 & `bald13.bt` > 0, "specific_bald13.bt",
                          fifelse(n_trait == 1 & `bald14.bt` > 0, "specific_bald14.bt", NA_character_)))))
]
site <- site[!is.na(cls)]
site[, col := hi_col[cls]]

lead_map <- rbindlist(c(
  lapply(traits, function(tr) site[cls == "shared_all", .(trait = tr, pos_id, CHR, POS, chr, gene, col)]),
  list(
    site[cls == "specific_bald.qt",   .(trait = "bald.qt",   pos_id, CHR, POS, chr, gene, col)],
    site[cls == "specific_bald12.bt", .(trait = "bald12.bt", pos_id, CHR, POS, chr, gene, col)],
    site[cls == "specific_bald13.bt", .(trait = "bald13.bt", pos_id, CHR, POS, chr, gene, col)],
    site[cls == "specific_bald14.bt", .(trait = "bald14.bt", pos_id, CHR, POS, chr, gene, col)]
  )
), use.names = TRUE)

lead <- merge(lead_map, g[, .(trait, CHR, POS, P_raw, logp)], by = c("trait","CHR","POS"), all.x = TRUE)
lead <- merge(lead, c0[, .(trait, CHR, POS, pJ)], by = c("trait","CHR","POS"), all.x = TRUE)
lead[is.na(P_raw), P_raw := pJ]
lead[is.na(logp), logp := pmin(-log10(pJ), max_log10)]

gene_disp <- site[, .(
  n_show = .N,
  cls_show = if (.N == 1L) cls[1L] else "mixed"
), by = gene]

lab <- merge(
  site[, .SD[which.min(locus_pJ)], by = gene][, .(gene, CHR, POS, chr)],
  gene_disp,
  by = "gene",
  all.x = TRUE
)

lab[, label_col := "black"]
lab[n_show == 1L & grepl("^specific_", cls_show), label_col := unname(hi_col[cls_show])]

lab <- lab[!grepl("^LOC|^LINC", gene, ignore.case = TRUE)]

setorder(lab, CHR, POS)
lab[, L := chr_len_map[as.character(CHR)]]
lab[, h := pmax(1e5, round(L * 0.0012))]
lab_df <- lab[, .(chr, start = pmax(1L, POS - h), end = pmin(L, POS + h), gene)]

cyto <- as.data.table(read.cytoband(species = "hg38")[[1]])
setnames(cyto, 1:5, c("chr","start","end","name","gieStain"))
cyto <- cyto[chr %chin% chr_df$chr]
stain_col <- c(gneg="#FFFFFF", gpos25="#D9D9D9", gpos50="#A6A6A6", gpos75="#6E6E6E", gpos100="#000000", gvar="#CFCFCF", stalk="#BDBDBD", acen="#E41A1C")
cyto[, fill := stain_col[gieStain]]
cyto[is.na(fill), fill := "#FFFFFF"]

venn <- make_venn_grob(c0)

# plot
pdf(out, width = 8, height = 8, useDingbats = FALSE)
par(mai = c(0,0,0,0), xpd = NA)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(rep(1, 22), 6),
  cell.padding = c(0,0,0,0),
  track.margin = c(0.0002, 0.0002),
  points.overflow.warning = FALSE,
  canvas.xlim = c(-0.78, 0.78),
  canvas.ylim = c(-0.78, 0.78)
)
circos.initialize(chr_df$chr, xlim = cbind(0, chr_df$chr_len))

circos.genomicLabels(
  gene_df, labels.column = 4, side = "outside",
  col = gene2$col, line_col = gene2$col, line_lwd = 0.28, cex = 0.4,
  labels_height = mm_h(40), connection_height = mm_h(7), padding = 0.25
)

circos.trackPlotRegion(ylim = c(0,1), track.height = 0.012, bg.border = NA, panel.fun = function(x, y) {
  s <- CELL_META$sector.index
  d <- cyto[chr == s]
  for (i in 1:nrow(d)) {
    if (d$gieStain[i] == "acen") {
      k <- sum(d$gieStain[1:i] == "acen")
      if (k %% 2 == 1) circos.polygon(c(d$start[i], d$end[i], d$start[i]), c(0.06,0.5,0.94), col = d$fill[i], border = NA)
      if (k %% 2 == 0) circos.polygon(c(d$end[i], d$start[i], d$end[i]), c(0.06,0.5,0.94), col = d$fill[i], border = NA)
    } else {
      circos.rect(d$start[i], 0.06, d$end[i], 0.94, col = d$fill[i], border = NA)
    }
  }
  circos.rect(CELL_META$xlim[1], 0.06, CELL_META$xlim[2], 0.94, col = NA, border = "#666666", lwd = 0.25)
})

circos.trackPlotRegion(ylim = c(0,1), track.height = 0.013, bg.border = NA, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 0.5, sub("chr", "", CELL_META$sector.index),
              facing = "bending.inside", niceFacing = TRUE, cex = 0.28, col = "black")
})

for (tr in traits) {
  gd <- g[trait == tr]
  ld <- lead[trait == tr]
  circos.trackPlotRegion(ylim = c(0, max_log10), track.height = 0.086, bg.border = NA, panel.fun = function(x, y) {
    s <- CELL_META$sector.index
    a <- gd[chr == s & P_raw > sig_cut]
    b <- gd[chr == s & P_raw <= sig_cut]
    z <- ld[chr == s]
    if (nrow(a)) circos.points(a$POS, a$logp, pch = 16, cex = 0.21, col = col_ns)
    if (nrow(b)) circos.points(b$POS, b$logp, pch = 16, cex = 0.21, col = col_sig)
    if (nrow(z)) circos.points(z$POS, z$logp, pch = 4, cex = 0.30, col = z$col)
    if (s == "chr1") {
      xr <- CELL_META$xlim
      circos.text(xr[1] - 0.18 * diff(xr), max_log10 * 0.50, track_lab[tr],
                  facing = "clockwise", niceFacing = FALSE,
                  adj = c(0.5, 1.1), cex = 0.50, col = trait_col[tr])
    }
  })
}

pushViewport(viewport(x = 0.50, y = 0.49, width = 0.23, height = 0.23))
grid.draw(venn)
popViewport()
draw_legend()

circos.clear()
dev.off()
message(out)