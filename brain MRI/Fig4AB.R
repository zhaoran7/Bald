suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(VennDiagram)
  library(officer)
  library(rvg)
})

root <- if (dir.exists("D:/bald/img")) "D:/bald/img" else "/mnt/d/bald/img"
out <- file.path(root, "plot")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

traits <- c("Continuous", "M shape", "O shape", "U shape")
cols <- c("Continuous" = "#079662", "M shape" = "#cc0404", "O shape" = "#1d79b2", "U shape" = "#5817d1")
font_family <- if (.Platform$OS.type == "windows") "Arial" else "sans"

set_size <- function(pptx, w, h) {
  td <- tempfile(); dir.create(td); unzip(pptx, exdir = td)
  f <- file.path(td, "ppt", "presentation.xml")
  x <- gsub('<p:sldSz cx="[0-9]+" cy="[0-9]+"[^>]*/>',
            sprintf('<p:sldSz cx="%d" cy="%d" type="custom"/>', round(w * 914400), round(h * 914400)),
            readLines(f, warn = FALSE))
  writeLines(x, f, useBytes = TRUE)
  tmp <- tempfile(fileext = ".pptx")
  old <- setwd(td); on.exit({setwd(old); unlink(td, TRUE)}, add = TRUE)
  system2("zip", c("-qr9X", tmp, list.files(td, recursive = TRUE, all.files = TRUE, no.. = TRUE)))
  file.copy(tmp, pptx, overwrite = TRUE)
}

save_ppt <- function(file, expr, w, h) {
  ppt <- read_pptx()
  ppt <- add_slide(ppt, "Blank", "Office Theme")
  ppt <- ph_with(ppt, dml(code = expr), ph_location(0, 0, w, h))
  print(ppt, target = file)
  set_size(file, w, h)
}

lab_region <- function(region, hemi) {
  z <- c(
    bankssts = "Banks STS", caudalanteriorcingulate = "Caudal anterior cingulate",
    caudalmiddlefrontal = "Caudal middle frontal", cuneus = "Cuneus", entorhinal = "Entorhinal",
    fusiform = "Fusiform", inferiorparietal = "Inferior parietal", inferiortemporal = "Inferior temporal",
    insula = "Insula", isthmuscingulate = "Isthmus cingulate", lateraloccipital = "Lateral occipital",
    lateralorbitofrontal = "Lateral orbitofrontal", lingual = "Lingual", medialorbitofrontal = "Medial orbitofrontal",
    middletemporal = "Middle temporal", paracentral = "Paracentral", parahippocampal = "Parahippocampal",
    parsopercularis = "Pars opercularis", parsorbitalis = "Pars orbitalis", parstriangularis = "Pars triangularis",
    pericalcarine = "Pericalcarine", postcentral = "Postcentral", posteriorcingulate = "Posterior cingulate",
    precentral = "Precentral", precuneus = "Precuneus", rostralanteriorcingulate = "Rostral anterior cingulate",
    rostralmiddlefrontal = "Rostral middle frontal", superiorfrontal = "Superior frontal",
    superiorparietal = "Superior parietal", superiortemporal = "Superior temporal", supramarginal = "Supramarginal",
    transversetemporal = "Transverse temporal"
  )
  r <- fifelse(region %chin% names(z), unname(z[region]), tools::toTitleCase(region))
  h <- fifelse(hemi == "left", "L", fifelse(hemi == "right", "R", ""))
  trimws(paste(h, r))
}

x <- fread(file.path(root, "res", "mri_mr.tsv"))[direction == "img_to_bald" & status == "ok"]
x[, sig := fdr_mri272 < 0.05]
x[, trait := factor(trait, traits)]

sig_id <- unique(x[sig == TRUE, mri_id])
d <- x[mri_id %chin% sig_id]
info <- d[, .(
  label = lab_region(region[which.min(fdr_mri272)], hemi[which.min(fdr_mri272)]),
  group = {i <- which.min(fdr_mri272); if (mri_class[i] == "cortex") paste("Cortex", measure[i]) else if (mri_class[i] == "subcortex") "Subcortical volume" else paste("White matter", measure[i])},
  mfdr = min(fdr_mri272)
), by = mri_id][order(group, mfdr)]

d <- merge(d[sig == TRUE], info, by = "mri_id")
d[, `:=`(mri_id = factor(mri_id, info$mri_id), lFDR = pmin(-log10(fdr_mri272), 60))]
lim <- max(0.15, quantile(abs(d$beta), .98, na.rm = TRUE))

p_heat <- ggplot(d, aes(mri_id, trait)) +
  geom_point(aes(size = lFDR, fill = pmax(pmin(beta, lim), -lim)), shape = 21, color = "grey30", stroke = .15) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(info$label, info$mri_id)) +
  scale_y_discrete(limits = rev(traits)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), name = "MR beta") +
  scale_size(range = c(1.0, 4.3), name = expression(-log[10](FDR))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 5.6, base_family = font_family) +
  theme(strip.text = element_text(size = 5.6, face = "plain"),
        axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 5.6),
        axis.text.y = element_text(size = 5.6),
        legend.title = element_text(size = 5.6), legend.text = element_text(size = 5.6),
        legend.position = "bottom", legend.box = "horizontal",
        panel.grid = element_line(color = "grey90", linewidth = .22), panel.spacing.x = unit(.45, "mm"))

venn_grob <- local({
  s <- lapply(traits, function(tr) unique(x[trait == tr & sig == TRUE, mri_id]))
  u <- unique(unlist(s))
  z <- data.table(id = u, A = u %chin% s[[1]], B = u %chin% s[[2]], C = u %chin% s[[3]], D = u %chin% s[[4]])
  grobTree(draw.quad.venn(
    area1 = sum(z$A), area2 = sum(z$B), area3 = sum(z$C), area4 = sum(z$D),
    n12 = sum(z$A & z$B), n13 = sum(z$A & z$C), n14 = sum(z$A & z$D),
    n23 = sum(z$B & z$C), n24 = sum(z$B & z$D), n34 = sum(z$C & z$D),
    n123 = sum(z$A & z$B & z$C), n124 = sum(z$A & z$B & z$D),
    n134 = sum(z$A & z$C & z$D), n234 = sum(z$B & z$C & z$D), n1234 = sum(z$A & z$B & z$C & z$D),
    category = traits, fill = unname(cols[traits]), alpha = rep(.24, 4), col = NA, lwd = 0, lty = "blank",
    cex = .75, cat.cex = .75, fontface = "plain", cat.fontface = "plain",
    fontfamily = font_family, cat.fontfamily = font_family, cat.col = unname(cols[traits]), margin = .02, ind = FALSE
  ))
})

save_ppt(file.path(out, "Fig3_heatmap.pptx"), print(p_heat), 12, 4)
save_ppt(file.path(out, "Fig3_venn.pptx"), grid.draw(venn_grob), 3, 3)
message(file.path(out, "Fig3_heatmap.pptx"))
message(file.path(out, "Fig3_venn.pptx"))
