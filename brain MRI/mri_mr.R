library(data.table)

to_wsl <- function(x) {
  if (.Platform$OS.type == "unix" && grepl("^[A-Za-z]:", x)) {
    paste0("/mnt/", tolower(substr(x, 1, 1)), sub("^[A-Za-z]:", "", x))
  } else x
}

check_dir <- function(x, label) {
  if (!dir.exists(x)) {
    stop(
      label, " is not readable: ", x, "\n",
      "If this is an F: drive path in WSL, remount it first:\n",
      "  sudo umount -l /mnt/f 2>/dev/null || true\n",
      "  sudo mkdir -p /mnt/f\n",
      "  sudo mount -t drvfs F: /mnt/f -o metadata,uid=$(id -u),gid=$(id -g)\n",
      call. = FALSE
    )
  }
}

root  <- to_wsl(Sys.getenv("IMG_ROOT", "D:/bald/img"))
name_f <- to_wsl(Sys.getenv("IMG_NAME", "D:/data/ukb/phe/common/img.lst"))
img0  <- to_wsl(Sys.getenv("IMG_GWAS_DIR", "F:/bald/img/data/gwas"))
gwas0 <- to_wsl(Sys.getenv("BALD_GWAS_DIR", "D:/bald/data/gwas"))
ld0   <- to_wsl(Sys.getenv("LD_PGEN_DIR", "F:/refGen/1kg_phase3"))
out0  <- file.path(root, "res")
mr0   <- file.path(out0, "mr")
tmp0  <- file.path(mr0, "tmp_mri272")
dir.create(mr0, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp0, recursive = TRUE, showWarnings = FALSE)
check_dir(img0, "IMG_GWAS_DIR")
check_dir(ld0, "LD_PGEN_DIR")

trait_map <- c("bald.qt" = "Continuous", "bald12.bt" = "M shape",
               "bald13.bt" = "O shape", "bald14.bt" = "U shape")
trait_levels <- c("Continuous", "M shape", "O shape", "U shape")
bald_files <- file.path(gwas0, paste0(names(trait_map), ".gz"))
names(bald_files) <- names(trait_map)

p_inst <- as.numeric(Sys.getenv("P_INST", "5e-8"))
clump_r2 <- Sys.getenv("CLUMP_R2", "0.001")
clump_kb <- Sys.getenv("CLUMP_KB", "10000")
download_missing <- Sys.getenv("DOWNLOAD_MISSING", "TRUE") == "TRUE"
drop_pal <- Sys.getenv("DROP_PALINDROMIC", "TRUE") == "TRUE"

std_region <- function(x) trimws(gsub("\\s+", " ", gsub("\\+", " and ", gsub("[-_]+", " ", x))))
ok_rsid <- function(x) grepl("^rs[0-9]+$", x, ignore.case = TRUE)
uc <- function(x) toupper(as.character(x))
comp <- function(x) chartr("ACGT", "TGCA", x)
pal <- function(a, b) paste0(a, b) %chin% c("AT", "TA", "CG", "GC")
pz <- function(z) 2 * pnorm(-abs(z))
msg <- function(...) message(format(Sys.time(), "[%F %T] "), ...)

mri272 <- function() {
  nm <- fread(name_f, header = FALSE)
  setnames(nm, c("field", "phenotype_name", "all"))
  nm[, `:=`(field = as.character(field), field_num = sub("^p", "", as.character(field)), img = paste0("img_", field))]

  cort <- nm[grepl("^aparc-Desikan_", phenotype_name)]
  cort[, hemi := fifelse(grepl("_lh_", phenotype_name), "left", "right")]
  cort[, measure := sub("^aparc-Desikan_[lr]h_([^_]+)_.*", "\\1", phenotype_name)]
  cort[, region := std_region(sub("^aparc-Desikan_[lr]h_[^_]+_", "", phenotype_name))]
  cort[, `:=`(
    mri_id = paste("cortex", measure, hemi, region, sep = "_"),
    mri_class = "cortex", atlas = "Desikan-Killiany",
    standard_name = paste("Cortical", measure, hemi, region),
    component_n = 1L, collapse_method = "single"
  )]

  sub_keep <- c("Lateral-Ventricle", "Thalamus-Proper", "Caudate", "Putamen",
                "Pallidum", "Hippocampus", "Amygdala", "Accumbens-area")
  sub <- nm[grepl("^aseg_[lr]h_volume_", phenotype_name)]
  sub <- sub[rowSums(sapply(sub_keep, function(k) grepl(k, sub$phenotype_name))) > 0]
  sub[, hemi := fifelse(grepl("^aseg_lh_", phenotype_name), "left", "right")]
  sub[, region := std_region(sub("^aseg_[lr]h_volume_", "", phenotype_name))]
  sub[, `:=`(
    mri_id = paste("subcortex", hemi, region, sep = "_"),
    mri_class = "subcortex", atlas = "FreeSurfer aseg", measure = "volume",
    standard_name = paste("Subcortical volume", hemi, region),
    component_n = 1L, collapse_method = "single"
  )]

  wm <- nm[grepl("^IDP_dMRI_TBSS_(FA|MD)_", phenotype_name)]
  wm[, measure := sub("^IDP_dMRI_TBSS_([^_]+)_.*", "\\1", phenotype_name)]
  wm[, tract0 := sub("^IDP_dMRI_TBSS_[^_]+_", "", phenotype_name)]
  wm[, tract := sub("_[LR]$", "", tract0)]
  wm[, hemi := fifelse(grepl("_L$", tract0), "left", fifelse(grepl("_R$", tract0), "right", "midline"))]
  wm[, `:=`(
    mri_id = paste("white_matter", measure, std_region(tract), sep = "_"),
    mri_class = "white_matter", atlas = "JHU TBSS",
    region = std_region(tract), standard_name = paste("White matter", measure, std_region(tract)),
    collapse_method = "best_component_by_p"
  )]
  wm[, component_n := .N, by = mri_id]

  x <- rbindlist(list(
    cort[, .(mri_id, mri_class, atlas, standard_name, measure, region, hemi, img, field, field_num, phenotype_name, component_n, collapse_method)],
    sub[, .(mri_id, mri_class, atlas, standard_name, measure, region, hemi, img, field, field_num, phenotype_name, component_n, collapse_method)],
    wm[, .(mri_id, mri_class, atlas, standard_name, measure, region, hemi, img, field, field_num, phenotype_name, component_n, collapse_method)]
  ), fill = TRUE)
  stopifnot(uniqueN(x$mri_id) == 272L)
  x[]
}

read_big40 <- function(f = file.path(img0, "BIG40_IDPs.html")) {
  if (!file.exists(f)) stop("BIG40_IDPs.html is missing: ", f, call. = FALSE)
  s <- paste(readLines(f, warn = FALSE), collapse = "\n")
  rows <- strsplit(s, "<tr", fixed = TRUE)[[1]][-1]
  parse_row <- function(r) {
    cells <- regmatches(r, gregexpr("<t[dh][^>]*>.*?</t[dh]>", r, perl = TRUE))[[1]]
    cells <- gsub("<[^>]+>", "", cells)
    cells <- gsub("&nbsp;", " ", cells, fixed = TRUE)
    cells <- gsub("&amp;", "&", cells, fixed = TRUE)
    trimws(cells)
  }
  rbindlist(lapply(rows, function(r) {
    z <- parse_row(r)
    if (length(z) >= 4 && grepl("^[0-9]{4}$", z[2]) && grepl("^[0-9]+$", z[3])) {
      data.table(pheno_id = z[2], field_num = z[3], big40_short = z[4])
    } else NULL
  }), fill = TRUE)
}

ensure_gwas <- function(manifest) {
  base_url <- "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"
  manifest[, file := file.path(img0, sprintf("%s_%s.txt.gz", field, pheno_id))]
  manifest[, url := sprintf("%s/%s.txt.gz", base_url, pheno_id)]
  manifest[, ok := file.exists(file) & file.exists(paste0(file, ".ok")) & file.info(file)$size > 1e6]
  if (download_missing) {
    todo <- manifest[ok != TRUE & !is.na(pheno_id)]
    if (nrow(todo)) msg("download missing imaging GWAS files: ", nrow(todo))
    for (i in seq_len(nrow(todo))) {
      z <- todo[i]
      msg("download ", i, "/", nrow(todo), " ", basename(z$file))
      dir.create(dirname(z$file), recursive = TRUE, showWarnings = FALSE)
      unlink(paste0(z$file, ".ok"))
      code <- system2("curl", c("-L", "-C", "-", "--retry", "10", "--retry-delay", "20",
                                "-o", z$file, z$url))
      if (code == 0 && file.exists(z$file) && file.info(z$file)$size > 1e6) {
        writeLines(as.character(Sys.time()), paste0(z$file, ".ok"))
      }
    }
    manifest[, ok := file.exists(file) & file.exists(paste0(file, ".ok")) & file.info(file)$size > 1e6]
  }
  fwrite(manifest, file.path(out0, "mri_gwas_manifest.tsv"), sep = "\t")
  manifest
}

read_img <- function(f) {
  x <- fread(cmd = paste("gzip -dc", shQuote(f)))
  setnames(x, c("chr", "SNP", "pos", "a1", "a2", "beta", "se", "logp"))
  x[, `:=`(
    chr = as.character(as.integer(chr)),
    ea = uc(a2), oa = uc(a1),
    p = pmax(10^(-as.numeric(logp)), 1e-300),
    eaf = NA_real_
  )]
  x[chr %chin% as.character(1:22) & ok_rsid(SNP), .(SNP, chr, pos, ea, oa, eaf, beta, se, p)]
}

read_bald <- function(f, snps = NULL) {
  if (is.null(snps)) {
    x <- fread(cmd = paste("gzip -dc", shQuote(f)))
  } else {
    sf <- tempfile(tmpdir = tmp0)
    writeLines(unique(snps), sf)
    on.exit(unlink(sf), add = TRUE)
    cmd <- sprintf("gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{a[$1]=1; next} FNR==1 || ($1 in a)' %s -",
                   shQuote(f), shQuote(sf))
    x <- fread(cmd = cmd)
  }
  if (!nrow(x)) return(x)
  setnames(x, c("CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P"),
           c("chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"))
  x[, `:=`(chr = as.character(chr), ea = uc(ea), oa = uc(oa))]
  x[chr %chin% as.character(1:22) & ok_rsid(SNP)]
}

clump_local <- function(x, label) {
  x <- unique(x[is.finite(p) & p <= p_inst & chr %chin% as.character(1:22)], by = "SNP")
  if (!nrow(x)) return(x)
  plink2 <- Sys.which("plink2")
  if (!nzchar(plink2)) stop("plink2 not found")
  keep <- character()
  for (cc in sort(unique(x$chr))) {
    pbase <- file.path(ld0, paste0("ALL.chr", cc))
    if (!file.exists(paste0(pbase, ".pgen"))) next
    inf <- file.path(tmp0, paste0(label, ".chr", cc, ".clump.in"))
    outf <- file.path(tmp0, paste0(label, ".chr", cc))
    fwrite(x[chr == cc, .(SNP, P = p)], inf, sep = "\t")
    system2(plink2, c("--pfile", pbase, "--allow-extra-chr", "--clump", inf,
      "--clump-id-field", "SNP", "--clump-p-field", "P", "--clump-p1", p_inst,
      "--clump-r2", clump_r2, "--clump-kb", clump_kb, "--out", outf),
      stdout = TRUE, stderr = TRUE)
    cf <- c(paste0(outf, ".clumps"), paste0(outf, ".clumped"))
    cf <- cf[file.exists(cf)]
    if (length(cf)) {
      z <- tryCatch(fread(cf[1], fill = TRUE), error = function(e) NULL)
      id <- intersect(c("ID", "SNP"), names(z))[1]
      if (!is.na(id)) keep <- c(keep, z[[id]])
    }
  }
  x[SNP %chin% keep]
}

harmonise <- function(exposure, outcome) {
  x <- merge(exposure, outcome, by = "SNP", suffixes = c(".exposure", ".outcome"))
  if (!nrow(x)) return(x)
  x <- x[!drop_pal | !pal(ea.exposure, oa.exposure)]
  same <- x$ea.exposure == x$ea.outcome & x$oa.exposure == x$oa.outcome
  swap <- x$ea.exposure == x$oa.outcome & x$oa.exposure == x$ea.outcome
  csame <- x$ea.exposure == comp(x$ea.outcome) & x$oa.exposure == comp(x$oa.outcome)
  cswap <- x$ea.exposure == comp(x$oa.outcome) & x$oa.exposure == comp(x$ea.outcome)
  keep <- same | swap | csame | cswap
  x <- x[keep]
  x[, beta.outcome := fifelse((swap | cswap)[keep], -beta.outcome, beta.outcome)]
  x[, `:=`(bx = beta.exposure, by = beta.outcome, sx = se.exposure, sy = se.outcome)]
  x[is.finite(bx) & is.finite(by) & is.finite(sx) & is.finite(sy) & bx != 0]
}

mr_primary <- function(d) {
  ns <- uniqueN(d$SNP)
  if (!ns) return(data.table(nsnp = 0L, method = NA_character_, beta = NA_real_, se = NA_real_, p = NA_real_, q = NA_real_, q_p = NA_real_))
  r <- d$by / d$bx
  sr <- abs(d$sy / d$bx)
  w <- 1 / sr^2
  if (ns == 1) {
    b <- r[1]
    se <- sqrt(d$sy[1]^2 / d$bx[1]^2 + d$by[1]^2 * d$sx[1]^2 / d$bx[1]^4)
    return(data.table(nsnp = ns, method = "Wald ratio", beta = b, se = se, p = pz(b / se), q = NA_real_, q_p = NA_real_))
  }
  b <- sum(w * r) / sum(w)
  se <- sqrt(1 / sum(w))
  q <- sum(w * (r - b)^2)
  data.table(nsnp = ns, method = "IVW fixed", beta = b, se = se, p = pz(b / se), q = q, q_p = pchisq(q, ns - 1, lower.tail = FALSE))
}

metric <- mri272()
big40 <- read_big40()
manifest <- merge(metric, big40, by = "field_num", all.x = TRUE)
manifest <- ensure_gwas(manifest)
fwrite(manifest, file.path(out0, "mri272_component_manifest.tsv"), sep = "\t")
if (nrow(manifest[ok != TRUE | is.na(ok)])) {
  stop("Some 272-MRI component GWAS files are still missing or incomplete. Re-run the script to resume download, or check mri272_component_manifest.tsv.")
}

work <- manifest[ok == TRUE & !is.na(pheno_id)]
if (!nrow(work)) stop("No complete imaging GWAS file for the 272 MRI components")

msg("prepare bald instruments")
bald_inst <- lapply(names(bald_files), function(nm) clump_local(read_bald(bald_files[[nm]]), paste0("bald_", nm)))
names(bald_inst) <- names(bald_files)
union_bald <- unique(unlist(lapply(bald_inst, `[[`, "SNP")))

inst_cache <- file.path(mr0, "mri272_img_instruments.rds")
if (file.exists(inst_cache)) {
  msg("read cached imaging instruments")
  img_inst <- readRDS(inst_cache)
} else {
  img_inst <- vector("list", nrow(work))
  names(img_inst) <- work$img
  for (i in seq_len(nrow(work))) {
    z <- work[i]
    msg("clump imaging instruments ", i, "/", nrow(work), " ", z$img)
    img_inst[[z$img]] <- clump_local(read_img(z$file), paste0(z$field, "_", z$pheno_id))
  }
  saveRDS(img_inst, inst_cache)
}
union_img <- unique(unlist(lapply(img_inst, `[[`, "SNP")))
bald_out <- lapply(bald_files, read_bald, snps = union_img)

component_res <- list()
for (i in seq_len(nrow(work))) {
  z <- work[i]
  img_full <- NULL
  if (length(union_bald)) img_full <- read_img(z$file)[SNP %chin% union_bald]
  for (bn in names(bald_files)) {
    if (nrow(img_inst[[z$img]])) {
      d <- harmonise(img_inst[[z$img]], bald_out[[bn]])
      r <- mr_primary(d)
      r[, `:=`(direction = "img_to_bald", bald = bn, trait = trait_map[bn], img = z$img)]
      component_res[[length(component_res) + 1L]] <- r
    }
    if (nrow(bald_inst[[bn]]) && !is.null(img_full)) {
      d <- harmonise(bald_inst[[bn]], img_full)
      r <- mr_primary(d)
      r[, `:=`(direction = "bald_to_img", bald = bn, trait = trait_map[bn], img = z$img)]
      component_res[[length(component_res) + 1L]] <- r
    }
  }
}
comp_res <- rbindlist(component_res, fill = TRUE)
comp_res <- merge(manifest, comp_res, by = "img", all.y = TRUE)
fwrite(comp_res, file.path(mr0, "mri272_component_mr.tsv"), sep = "\t")

grid <- CJ(direction = c("img_to_bald", "bald_to_img"), trait = trait_levels, mri_id = unique(metric$mri_id))
collapse <- copy(comp_res)
collapse[, status := fifelse(is.na(method), "no_mr_result", fifelse(nsnp == 0, "no_harmonised_instrument", "ok"))]
setorderv(collapse, c("direction", "trait", "mri_id", "p"), na.last = TRUE)
collapse[, component_rank := seq_len(.N), by = .(direction, trait, mri_id)]
best <- collapse[component_rank == 1L]
best <- merge(grid, best, by = c("direction", "trait", "mri_id"), all.x = TRUE)
best[is.na(status), status := "missing_or_failed_gwas"]
best[, `:=`(
  selected_img = img,
  selected_field = field,
  selected_phenotype_name = phenotype_name,
  bonferroni = pmin(p * uniqueN(mri_id), 1),
  fdr_mri272 = p.adjust(p, "fdr")
), by = .(direction, trait)]
setorder(best, direction, trait, p)

out <- best[, .(
  direction, trait, mri_id, mri_class, atlas, standard_name, measure, region, hemi,
  component_n, collapse_method, selected_img, selected_field, selected_phenotype_name,
  status, method, nsnp, beta, se, p, bonferroni, fdr_mri272, q, q_p
)]
fwrite(out, file.path(out0, "mri_mr.tsv"), sep = "\t")
msg("done: ", file.path(out0, "mri_mr.tsv"))
