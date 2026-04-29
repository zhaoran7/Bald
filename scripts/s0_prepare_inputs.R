library(data.table)

script_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
if (is.na(script_file)) {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- if (length(cmd_file)) normalizePath(sub("^--file=", "", cmd_file[1]), winslash = "/") else normalizePath("scripts/s0_prepare_inputs.R", winslash = "/")
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

read_gwas_header <- function(f) {
  x <- fread(f, nrows = 0)
  names(x)
}

for (tr in traits) {
  cojo_f <- file.path(cojo0, tr, "jma.cojo")
  gwas_f <- file.path(gwas0, paste0(tr, ".gz"))

  if (!file.exists(cojo_f)) stop("Missing COJO file: ", cojo_f)
  if (!file.exists(gwas_f)) stop("Missing GWAS file: ", gwas_f)

  cj <- fread(cojo_f)
  chr_col <- pick1(names(cj), c("Chr", "CHR", "chr"))
  bp_col <- pick1(names(cj), c("bp", "BP", "POS", "pos"))
  if (is.na(chr_col) || is.na(bp_col) || !"SNP" %in% names(cj)) {
    stop("COJO file must contain Chr/CHR, SNP and bp/BP columns: ", cojo_f)
  }

  lead <- cj[, .(
    lead_chr = as.integer(get(chr_col)),
    lead_snp = as.character(SNP),
    lead_bp = as.integer(get(bp_col))
  )]
  lead <- lead[!is.na(lead_chr) & !is.na(lead_bp) & nzchar(lead_snp)]
  fwrite(lead, file.path(lead0, paste0(tr, ".lead.tsv")), sep = "\t")

  gh <- read_gwas_header(gwas_f)
  keep <- unique(na.omit(c(
    pick1(gh, c("SNP", "rsid", "ID")),
    pick1(gh, c("CHR", "Chr", "chr")),
    pick1(gh, c("POS", "BP", "bp", "pos")),
    pick1(gh, c("EA", "A1", "effect_allele")),
    pick1(gh, c("NEA", "A2", "other_allele")),
    pick1(gh, c("BETA", "beta", "b"))
  )))
  gw <- fread(gwas_f, select = keep, showProgress = FALSE)
  setnames(gw, pick1(names(gw), c("SNP", "rsid", "ID")), "SNP")
  setnames(gw, pick1(names(gw), c("CHR", "Chr", "chr")), "CHR")
  setnames(gw, pick1(names(gw), c("POS", "BP", "bp", "pos")), "POS")
  setnames(gw, pick1(names(gw), c("EA", "A1", "effect_allele")), "effect_allele")
  setnames(gw, pick1(names(gw), c("NEA", "A2", "other_allele")), "other_allele")
  setnames(gw, pick1(names(gw), c("BETA", "beta", "b")), "beta")
  gw <- gw[, .(
    SNP = as.character(SNP),
    CHR = as.integer(CHR),
    POS = as.integer(POS),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele = toupper(as.character(other_allele)),
    beta = as.numeric(beta)
  )]

  z <- merge(
    lead[, .(lead_chr, lead_snp, lead_bp)],
    gw,
    by.x = c("lead_snp", "lead_chr", "lead_bp"),
    by.y = c("SNP", "CHR", "POS"),
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(z$effect_allele) || anyNA(z$other_allele) || anyNA(z$beta)) {
    stop("Some COJO lead SNPs were not found in GWAS or lack EA/NEA/BETA: ", tr)
  }

  z[, `:=`(
    trait = tr,
    risk_allele = fifelse(beta >= 0, effect_allele, other_allele)
  )]
  fwrite(
    z[, .(trait, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, risk_allele)],
    file.path(res0, paste0(tr, ".risk.tsv")),
    sep = "\t"
  )
}

risk <- rbindlist(lapply(traits, function(tr) fread(file.path(res0, paste0(tr, ".risk.tsv")))), fill = TRUE)
fwrite(risk, file.path(res0, "risk.tsv"), sep = "\t")
message("Done: ", file.path(res0, "risk.tsv"), " and ", lead0)
