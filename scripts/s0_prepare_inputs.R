library(data.table)

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/"),
  error = function(e) NA_character_
)
if (is.na(script_file)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- if (length(cmd_file)) {
    normalizePath(sub("^--file=", "", cmd_file[1]), winslash = "/")
  } else {
    normalizePath("scripts/s0_prepare_inputs.R", winslash = "/", mustWork = FALSE)
  }
}

proj <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = FALSE)
data0 <- file.path(proj, "data")
cojo0 <- file.path(data0, "cojo")
gwas0 <- file.path(data0, "gwas")
lead0 <- file.path(data0, "lead")
res0 <- file.path(proj, "res")

traits <- c("bald.qt", "bald12.bt", "bald13.bt", "bald14.bt")

dir.create(lead0, recursive = TRUE, showWarnings = FALSE)
dir.create(res0, recursive = TRUE, showWarnings = FALSE)

pick1 <- function(nm, cand) {
  x <- cand[cand %chin% nm]
  if (length(x)) x[1] else NA_character_
}

norm_chr <- function(x) {
  suppressWarnings(as.integer(gsub("^chr", "", as.character(x), ignore.case = TRUE)))
}

read_cojo_leads <- function(f) {
  x <- fread(f)
  chr_col <- pick1(names(x), c("Chr", "CHR", "chr", "lead_chr"))
  snp_col <- pick1(names(x), c("SNP", "snp", "lead_snp", "rsid", "ID"))
  bp_col <- pick1(names(x), c("bp", "BP", "pos", "POS", "lead_bp"))

  if (is.na(chr_col) || is.na(snp_col) || is.na(bp_col)) {
    stop("COJO file must contain chromosome, SNP and bp columns: ", f)
  }

  y <- x[, .(
    lead_chr = norm_chr(get(chr_col)),
    lead_snp = as.character(get(snp_col)),
    lead_bp = suppressWarnings(as.integer(get(bp_col)))
  )]

  y[
    !is.na(lead_chr) &
      !is.na(lead_bp) &
      !is.na(lead_snp) &
      lead_snp != "" &
      lead_snp != "lead_snp"
  ]
}

read_gwas <- function(f) {
  nm <- names(fread(f, nrows = 0))
  snp_col <- pick1(nm, c("SNP", "rsid", "ID"))
  chr_col <- pick1(nm, c("CHR", "Chr", "chr"))
  bp_col <- pick1(nm, c("POS", "BP", "bp", "pos"))
  ea_col <- pick1(nm, c("EA", "A1", "effect_allele", "EffectAllele"))
  oa_col <- pick1(nm, c("NEA", "A2", "other_allele", "OtherAllele"))
  bt_col <- pick1(nm, c("BETA", "beta", "b"))

  if (anyNA(c(snp_col, chr_col, bp_col, ea_col, oa_col, bt_col))) {
    stop("GWAS file must contain SNP/CHR/POS/EA/NEA/BETA columns: ", f)
  }

  x <- fread(
    f,
    select = c(snp_col, chr_col, bp_col, ea_col, oa_col, bt_col),
    showProgress = FALSE
  )
  setnames(
    x,
    c(snp_col, chr_col, bp_col, ea_col, oa_col, bt_col),
    c("lead_snp", "gwas_chr", "gwas_bp", "effect_allele", "other_allele", "beta")
  )

  x[, `:=`(
    lead_snp = as.character(lead_snp),
    gwas_chr = norm_chr(gwas_chr),
    gwas_bp = suppressWarnings(as.integer(gwas_bp)),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele)),
    beta = suppressWarnings(as.numeric(beta))
  )]

  x[
    !is.na(lead_snp) &
      lead_snp != "" &
      !is.na(effect_allele) &
      !is.na(other_allele) &
      !is.na(beta)
  ]
}

build_risk_one_trait <- function(tr) {
  cojo_f <- file.path(cojo0, tr, "jma.cojo")
  gwas_f <- file.path(gwas0, paste0(tr, ".gz"))

  if (!file.exists(cojo_f)) stop("Missing COJO file: ", cojo_f)
  if (!file.exists(gwas_f)) stop("Missing GWAS file: ", gwas_f)

  lead <- unique(read_cojo_leads(cojo_f))
  lead[, row_id := .I]
  fwrite(
    lead[, .(lead_chr, lead_snp, lead_bp)],
    file.path(lead0, paste0(tr, ".lead.tsv")),
    sep = "\t"
  )

  gwas <- read_gwas(gwas_f)

  by_snp <- merge(
    lead,
    gwas,
    by = "lead_snp",
    all.x = TRUE,
    allow.cartesian = TRUE,
    sort = FALSE
  )
  by_snp[, match_rank := fifelse(!is.na(gwas_chr) & gwas_chr == lead_chr, 0L, 1L)]
  setorder(by_snp, row_id, match_rank)
  by_snp <- by_snp[, .SD[1], by = row_id]
  by_snp[, match_method := fifelse(!is.na(effect_allele), "SNP", NA_character_)]

  miss <- is.na(by_snp$effect_allele) | is.na(by_snp$other_allele) | is.na(by_snp$beta)
  if (any(miss)) {
    fallback <- merge(
      by_snp[miss, .(row_id, lead_chr, lead_bp)],
      gwas[, .(lead_chr = gwas_chr, lead_bp = gwas_bp, effect_allele, other_allele, beta)],
      by = c("lead_chr", "lead_bp"),
      all.x = TRUE,
      allow.cartesian = TRUE,
      sort = FALSE
    )
    fallback <- fallback[, .SD[1], by = row_id]
    m <- match(by_snp$row_id, fallback$row_id)
    fill <- miss & !is.na(m) & !is.na(fallback$effect_allele[m])
    by_snp[fill, `:=`(
      gwas_chr = lead_chr,
      gwas_bp = lead_bp,
      effect_allele = fallback$effect_allele[m[fill]],
      other_allele = fallback$other_allele[m[fill]],
      beta = fallback$beta[m[fill]],
      match_method = "CHR_BP"
    )]
  }

  by_snp[, `:=`(
    trait = tr,
    risk_allele = fifelse(beta >= 0, effect_allele, other_allele)
  )]

  bad <- by_snp[is.na(risk_allele) | risk_allele == ""]
  if (nrow(bad)) {
    miss_f <- file.path(res0, paste0(tr, ".risk.missing.tsv"))
    fwrite(
      bad[, .(trait, lead_chr, lead_snp, lead_bp, gwas_chr, gwas_bp, match_method)],
      miss_f,
      sep = "\t"
    )
    stop("Some COJO lead SNPs lack GWAS EA/NEA/BETA for ", tr, "; wrote: ", miss_f)
  }

  out <- by_snp[, .(
    trait,
    lead_chr,
    lead_snp,
    lead_bp,
    gwas_chr,
    gwas_bp,
    effect_allele,
    other_allele,
    beta,
    risk_allele,
    match_method
  )]
  setorder(out, trait, lead_chr, lead_bp, lead_snp)
  fwrite(out, file.path(res0, paste0(tr, ".risk.tsv")), sep = "\t")
  out
}

risk <- rbindlist(lapply(traits, build_risk_one_trait), fill = TRUE)
setorder(risk, trait, lead_chr, lead_bp, lead_snp)
fwrite(risk, file.path(res0, "risk.tsv"), sep = "\t")

message("Done:")
message("  lead files: ", lead0)
message("  risk file : ", file.path(res0, "risk.tsv"))
