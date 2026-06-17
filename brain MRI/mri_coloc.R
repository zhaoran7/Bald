library(data.table)

to_wsl <- function(x) if (.Platform$OS.type == "unix" && grepl("^[A-Za-z]:", x)) {
  paste0("/mnt/", tolower(substr(x, 1, 1)), sub("^[A-Za-z]:", "", x))
} else x

if (!requireNamespace("coloc", quietly = TRUE)) {
  stop("R package 'coloc' is not installed. Install it with:\n  conda install -y -c conda-forge r-coloc\n", call. = FALSE)
}

root <- to_wsl(Sys.getenv("IMG_ROOT", "D:/bald/img"))
img_dir <- to_wsl(Sys.getenv("IMG_GWAS_DIR", "F:/bald/img/data/gwas"))
bald_dir <- to_wsl(Sys.getenv("BALD_GWAS_DIR", "D:/bald/data/gwas"))
out <- file.path(root, "res", "coloc_mri")
cache <- file.path(out, "cache")
detail_dir <- file.path(out, "detail")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)

window_bp <- as.integer(Sys.getenv("COLOC_WINDOW", "500000"))
min_snps <- as.integer(Sys.getenv("COLOC_MIN_SNPS", "50"))
max_pairs <- as.integer(Sys.getenv("COLOC_MAX_PAIRS", "0"))
max_loci_per_pair <- as.integer(Sys.getenv("COLOC_MAX_LOCI_PER_PAIR", "0"))
img_n <- as.integer(Sys.getenv("IMG_N", "33000"))
sig_rule <- Sys.getenv("COLOC_SIG", "fdr")
p1 <- as.numeric(Sys.getenv("COLOC_P1", "1e-4"))
p2 <- as.numeric(Sys.getenv("COLOC_P2", "1e-4"))
p12 <- as.numeric(Sys.getenv("COLOC_P12", "1e-5"))

trait_map <- c("bald.qt" = "Continuous", "bald12.bt" = "M shape",
               "bald13.bt" = "O shape", "bald14.bt" = "U shape")
bald_files <- file.path(bald_dir, paste0(names(trait_map), ".gz"))
names(bald_files) <- names(trait_map)
bald_case_frac <- c(
  "bald12.bt" = 48940 / (48940 + 67636),
  "bald13.bt" = 56531 / (56531 + 67636),
  "bald14.bt" = 38836 / (38836 + 67636)
)

ok_rsid <- function(x) grepl("^rs[0-9]+$", x, ignore.case = TRUE)
uc <- function(x) toupper(as.character(x))
comp <- function(x) chartr("ACGT", "TGCA", x)
pal <- function(a, b) paste0(a, b) %chin% c("AT", "TA", "CG", "GC")
msg <- function(...) message(format(Sys.time(), "[%F %T] "), ...)

read_img <- function(f) {
  x <- fread(cmd = paste("gzip -dc", shQuote(f)), showProgress = FALSE)
  setnames(x, c("chr", "SNP", "pos", "a1", "a2", "beta", "se", "logp"))
  x[, `:=`(
    chr = as.integer(chr), pos = as.integer(pos),
    ea = uc(a2), oa = uc(a1),
    p = pmax(10^(-as.numeric(logp)), 1e-300)
  )]
  x[chr %in% 1:22 & ok_rsid(SNP), .(SNP, chr, pos, ea, oa, beta, se, p)]
}

read_bald <- function(f) {
  x <- fread(cmd = paste("gzip -dc", shQuote(f)), showProgress = FALSE)
  setnames(x, c("CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P"),
           c("chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"))
  x[, `:=`(chr = as.integer(chr), pos = as.integer(pos), ea = uc(ea), oa = uc(oa))]
  x[chr %in% 1:22 & ok_rsid(SNP), .(SNP, chr, pos, ea, oa, eaf, N, beta, se, p)]
}

keep_windows <- function(x, w) {
  if (!nrow(x) || !nrow(w)) return(x[0])
  xi <- x[, .(chr, start = pos, end = pos, SNP, pos, ea, oa, eaf = if ("eaf" %in% names(x)) eaf else NA_real_,
              N = if ("N" %in% names(x)) N else NA_real_, beta, se, p)]
  wi <- w[, .(chr, start, end)]
  setkey(xi, chr, start, end)
  setkey(wi, chr, start, end)
  unique(foverlaps(xi, wi, nomatch = 0L)[, names(x), with = FALSE], by = "SNP")
}

merge_windows <- function(inst, pair_id) {
  inst <- copy(inst)
  inst[, `:=`(chr = suppressWarnings(as.integer(chr)),
              pos = suppressWarnings(as.integer(pos)))]
  inst <- inst[chr %in% 1:22 & is.finite(pos)]
  w <- unique(inst[, .(chr = chr,
                     start = pmax(1L, pos - window_bp),
                     end = pos + window_bp,
                     lead_snp = SNP,
                     lead_pos = pos)])
  if (!nrow(w)) return(w)
  setorder(w, chr, start, end)
  outw <- list()
  for (cc in unique(w$chr)) {
    z <- w[chr == cc]
    cur <- z[1]
    leads <- cur$lead_snp
    if (nrow(z) >= 2) {
      for (i in 2:nrow(z)) {
        if (z$start[i] <= cur$end) {
          cur$end <- max(cur$end, z$end[i])
          leads <- c(leads, z$lead_snp[i])
        } else {
          cur$lead_snp <- paste(unique(leads), collapse = ";")
          outw[[length(outw) + 1L]] <- copy(cur)
          cur <- z[i]
          leads <- cur$lead_snp
        }
      }
    }
    cur$lead_snp <- paste(unique(leads), collapse = ";")
    outw[[length(outw) + 1L]] <- copy(cur)
  }
  ans <- rbindlist(outw, fill = TRUE)
  ans[, `:=`(pair_id = pair_id, locus_index = seq_len(.N))]
  if (max_loci_per_pair > 0 && nrow(ans) > max_loci_per_pair) ans <- ans[seq_len(max_loci_per_pair)]
  ans[]
}

read_bald_inst <- function(trait) {
  cl <- list.files(file.path(root, "res", "mr", "tmp_mri272"),
                   pattern = paste0("^bald_", gsub("\\.", "\\\\.", trait), "\\.chr[0-9]+\\.clumps$"),
                   full.names = TRUE)
  ids <- unique(unlist(lapply(cl, function(f) {
    z <- tryCatch(fread(f, fill = TRUE), error = function(e) NULL)
    if (is.null(z)) return(character())
    id <- intersect(c("ID", "SNP"), names(z))[1]
    if (length(id)) z[[id]] else character()
  })))
  if (!length(ids)) {
    msg("No cached clumps for ", trait, "; using genome-wide significant SNPs as coloc loci.")
    return(read_bald(bald_files[[trait]])[p <= 5e-8, .(SNP, chr, pos, ea, oa, eaf, beta, se, p)])
  }
  read_bald(bald_files[[trait]])[SNP %chin% ids, .(SNP, chr, pos, ea, oa, eaf, beta, se, p)]
}

harmonise_region <- function(img, bald) {
  x <- merge(img, bald, by = "SNP", suffixes = c(".img", ".bald"))
  if (!nrow(x)) return(x)
  if (!"eaf.bald" %in% names(x) && "eaf" %in% names(x)) setnames(x, "eaf", "eaf.bald")
  if (!"N.bald" %in% names(x) && "N" %in% names(x)) setnames(x, "N", "N.bald")
  x <- x[!pal(ea.img, oa.img)]
  same <- x$ea.img == x$ea.bald & x$oa.img == x$oa.bald
  swap <- x$ea.img == x$oa.bald & x$oa.img == x$ea.bald
  csame <- x$ea.img == comp(x$ea.bald) & x$oa.img == comp(x$oa.bald)
  cswap <- x$ea.img == comp(x$oa.bald) & x$oa.img == comp(x$ea.bald)
  keep <- same | swap | csame | cswap
  x <- x[keep]
  if (!nrow(x)) return(x)
  flip <- (swap | cswap)[keep]
  x[, beta.bald := fifelse(flip, -beta.bald, beta.bald)]
  x[, maf := pmin(eaf.bald, 1 - eaf.bald)]
  x <- x[is.finite(beta.img) & is.finite(se.img) & is.finite(beta.bald) & is.finite(se.bald) &
      is.finite(maf) & maf > 0 & maf < 0.5]
  setorder(x, SNP, pos.img, pos.bald)
  unique(x, by = "SNP")
}

run_coloc <- function(d, bald) {
  if (nrow(d) < min_snps) return(NULL)
  ds1 <- list(beta = d$beta.img, varbeta = d$se.img^2, snp = d$SNP,
              MAF = d$maf, N = img_n, type = "quant")
  ds2 <- list(beta = d$beta.bald, varbeta = d$se.bald^2, snp = d$SNP,
              MAF = d$maf, N = round(median(d$N.bald, na.rm = TRUE)))
  if (bald == "bald.qt") {
    ds2$type <- "quant"
  } else {
    ds2$type <- "cc"
    ds2$s <- unname(bald_case_frac[[bald]])
  }
  suppressWarnings(coloc::coloc.abf(ds1, ds2, p1 = p1, p2 = p2, p12 = p12))
}

mr <- fread(file.path(root, "res", "mri_mr.tsv"))[
  status == "ok" & method != "" & !is.na(p)
]
if (sig_rule == "bonferroni") {
  mr <- mr[bonferroni < 0.05]
} else {
  mr <- mr[fdr_mri272 < 0.05]
}
mr[, `:=`(
  img = selected_img,
  field = selected_field,
  phenotype_name = selected_phenotype_name,
  bald = names(trait_map)[match(trait, trait_map)],
  mr_beta = beta,
  mr_se = se,
  mr_p = p,
  mr_fdr = fdr_mri272,
  mr_bonferroni = bonferroni
)]
manifest <- fread(file.path(root, "res", "mri272_component_manifest.tsv"))[, .(img, field, file, pheno_id)]
mr <- merge(mr, manifest, by = c("img", "field"), all.x = TRUE)
mr <- mr[!is.na(bald) & file.exists(file)]
setorder(mr, mr_fdr, mr_p)
if (max_pairs > 0 && nrow(mr) > max_pairs) mr <- mr[seq_len(max_pairs)]
mr[, pair_id := sprintf("pair_%04d", .I)]
fwrite(mr, file.path(out, "coloc_pairs.tsv"), sep = "\t")

img_inst <- readRDS(file.path(root, "res", "mr", "mri272_img_instruments.rds"))
bald_inst <- lapply(names(bald_files), read_bald_inst)

loci <- rbindlist(lapply(seq_len(nrow(mr)), function(i) {
  z <- mr[i]
  inst <- if (z$direction == "img_to_bald") img_inst[[z$img]] else bald_inst[[z$bald]]
  if (is.null(inst) || !nrow(inst)) return(NULL)
  w <- merge_windows(inst, z$pair_id)
  if (!nrow(w)) return(NULL)
  cbind(z[, .(direction, bald, trait, img, field, file, pheno_id, phenotype_name,
              mri_id, mri_class, atlas, standard_name, measure, region, hemi,
              mr_beta, mr_se, mr_p, mr_fdr, mr_bonferroni)], w)
}), fill = TRUE)
fwrite(loci, file.path(out, "coloc_loci_input.tsv"), sep = "\t")

if (!nrow(loci)) stop("No coloc loci were generated.")

msg("cache bald GWAS subsets")
bald_cache <- lapply(names(bald_files), function(bn) {
  cf <- file.path(cache, paste0(bn, ".tsv.gz"))
  if (file.exists(cf)) return(fread(cf))
  x <- keep_windows(read_bald(bald_files[[bn]]), unique(loci[bald == bn, .(chr, start, end)]))
  fwrite(x, cf, sep = "\t")
  x
})
names(bald_cache) <- names(bald_files)

msg("cache imaging GWAS subsets")
img_cache <- list()
for (im in unique(loci$img)) {
  cf <- file.path(cache, paste0(im, ".tsv.gz"))
  if (file.exists(cf)) {
    img_cache[[im]] <- fread(cf)
  } else {
    f <- unique(loci[img == im, file])[1]
    x <- keep_windows(read_img(f), unique(loci[img == im, .(chr, start, end)]))
    fwrite(x, cf, sep = "\t")
    img_cache[[im]] <- x
  }
}

msg("run coloc")
res <- list()
for (i in seq_len(nrow(loci))) {
  z <- loci[i]
  img <- img_cache[[z$img]][chr == z$chr & pos >= z$start & pos <= z$end]
  bald <- bald_cache[[z$bald]][chr == z$chr & pos >= z$start & pos <= z$end]
  d <- harmonise_region(img, bald)
  if (nrow(d) < min_snps) {
    rr <- data.table(status = "too_few_snps", n_snps = nrow(d))
  } else {
    co <- tryCatch(run_coloc(d, z$bald), error = function(e) e)
    if (inherits(co, "error") || is.null(co)) {
      rr <- data.table(status = "coloc_failed", n_snps = nrow(d), error = conditionMessage(co))
    } else {
      sm <- as.list(co$summary)
      cres <- as.data.table(co$results)
      pp_col <- intersect(c("SNP.PP.H4", "SNP.PP.H4.abf"), names(cres))[1]
      if (is.na(pp_col)) stop("Could not find SNP-level posterior probability column in coloc results.")
      top <- cres[which.max(cres[[pp_col]])]
      rr <- data.table(status = "ok", n_snps = nrow(d),
                       nsnps_coloc = sm$nsnps,
                       PP.H0 = sm$PP.H0.abf, PP.H1 = sm$PP.H1.abf,
                       PP.H2 = sm$PP.H2.abf, PP.H3 = sm$PP.H3.abf, PP.H4 = sm$PP.H4.abf,
                       top_snp = top$snp, top_snp_pp_h4 = top[[pp_col]],
                       top_pos = d[SNP == top$snp, pos.bald][1])
      det <- merge(cres[, .(SNP = snp, SNP.PP.H4 = .SD[[1]]), .SDcols = pp_col], d, by = "SNP")
      det[, `:=`(pair_id = z$pair_id, locus_index = z$locus_index, bald = z$bald, img = z$img)]
      fwrite(det[order(-SNP.PP.H4)], file.path(detail_dir, paste0(z$pair_id, "_locus", z$locus_index, ".tsv.gz")), sep = "\t")
    }
  }
  res[[i]] <- cbind(z, rr)
  if (i %% 50 == 0) msg("coloc ", i, "/", nrow(loci))
}

ans <- rbindlist(res, fill = TRUE)
setorder(ans, -PP.H4, status, mr_fdr)
fwrite(ans, file.path(out, "coloc_loci.tsv"), sep = "\t")

best <- ans[status == "ok"][order(-PP.H4), .SD[1], by = pair_id]
setorder(best, -PP.H4)
fwrite(best, file.path(out, "coloc_best_by_pair.tsv"), sep = "\t")

clean_result <- copy(best)
if (nrow(clean_result)) {
  clean_result[, evidence := fifelse(PP.H4 >= 0.8, "strong",
                              fifelse(PP.H4 >= 0.5, "moderate",
                              fifelse(PP.H4 >= 0.2, "weak", "none")))]
  clean_result[, locus := paste0("chr", chr, ":", start, "-", end)]
  clean_result <- clean_result[, .(
    direction, trait, mri_class, atlas, standard_name, measure, region, hemi,
    locus, chr, start, end, lead_snp, n_snps,
    PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, evidence,
    top_snp, top_snp_pp_h4, top_pos,
    mr_beta, mr_se, mr_p, mr_fdr, mr_bonferroni
  )]
  fwrite(clean_result, file.path(out, "coloc_results.tsv"), sep = "\t")
  fwrite(clean_result[evidence %chin% c("strong", "moderate")],
         file.path(out, "coloc_signals.tsv"), sep = "\t")
  fwrite(clean_result[, .N, by = .(trait, mri_class, evidence)][order(trait, mri_class, evidence)],
         file.path(out, "coloc_signal_counts.tsv"), sep = "\t")
}

if (nrow(best)) {
  p <- ggplot2::ggplot(best, ggplot2::aes(x = PP.H4, y = -log10(mr_fdr), color = trait)) +
    ggplot2::geom_point(size = 1.8, alpha = .85) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::labs(x = "coloc PP.H4", y = expression(-log[10]("MR FDR")), color = NULL)
  ggplot2::ggsave(file.path(out, "coloc_pph4_summary.pdf"), p, width = 5.2, height = 3.6)

  hm <- best[PP.H4 >= 0.5]
  if (nrow(hm)) {
    hm[, locus := paste0("chr", chr, ":", format(start, big.mark = ",", scientific = FALSE),
                         "-", format(end, big.mark = ",", scientific = FALSE))]
    hm[, img_label := paste0(standard_name, " (", mri_class, ")")]
    hm[, locus_label := paste0(trait, "\n", locus)]
    hm[, trait := factor(trait, levels = c("Continuous", "M shape", "O shape", "U shape"))]
    setorder(hm, trait, chr, start, -PP.H4)
    xlev <- unique(hm$locus_label)
    ylev <- unique(hm[order(PP.H4)]$img_label)
    hm[, `:=`(locus_label = factor(locus_label, levels = xlev),
              img_label = factor(img_label, levels = ylev))]
    ph <- ggplot2::ggplot(hm, ggplot2::aes(x = locus_label, y = img_label)) +
      ggplot2::geom_tile(ggplot2::aes(fill = PP.H4), color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", PP.H4)), size = 2.2) +
      ggplot2::scale_fill_gradient(low = "#EEF5FA", high = "#1D79B2", limits = c(0.5, 1), name = "PP.H4") +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                     axis.title = ggplot2::element_blank(),
                     legend.position = "right") +
      ggplot2::labs(title = "MRI-baldness colocalization signals")
    ggplot2::ggsave(file.path(out, "coloc_heatmap.pdf"), ph,
                    width = max(8, min(18, 0.42 * length(xlev) + 3)),
                    height = max(4, min(12, 0.22 * length(ylev) + 2)))
    if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
      doc <- officer::read_pptx()
      doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
      doc <- officer::ph_with(doc, rvg::dml(ggobj = ph),
                              location = officer::ph_location(left = 0.3, top = 0.3,
                                                              width = max(8, min(13, 0.42 * length(xlev) + 3)),
                                                              height = max(4, min(7, 0.22 * length(ylev) + 2))))
      print(doc, target = file.path(out, "coloc_heatmap.pptx"))
    }
  }
}

msg("done: ", out)
