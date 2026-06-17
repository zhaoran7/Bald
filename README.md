# Bald manuscript analysis files

This folder contains the manuscript-facing code and main result files for four analysis blocks.
Large raw inputs, genome-wide intermediate files, downloaded reference panels and local caches are not included.

## archaic introgression

- `code/archaic_locus_scripts`: locus-level archaic introgression workflow.
- `code/ibdmix`: full IBDmix scripts configured for local paths.
- `code/main_figure`: Fig. 2 plotting scripts.
- `results/locus_lead`: COJO lead-variant input tables.
- `results/locus_report`: core haplotype, archaic matching, YRI filter, selected-region and phylogeny summary tables.
- `results/phyml`: PhyML input, tree and metadata files.
- `results/causal_followup`: compact causal follow-up summaries and scripts; large LD/resource files were excluded.
- `figures`: Fig. 2 PDFs.

## function annotation

- `code/archaic_regulatory_annotation`: regulatory annotation scripts for archaic-like haplotypes.
- `code/main_figure`: Fig. 3 plotting scripts.
- `code/magma_enrichment_final` and `code/postgwas`: MAGMA/pathway enrichment code.
- `results/archaic_regulatory_annotation`: candidate regulatory annotation inputs/results.
- `results/Fig3_cache` and `results/Fig3_support`: Fig. 3 support tracks and cached annotation tracks.
- `results/magma_x_summary`, `results/magma_x_enrichment`, `results/pathway`: MAGMA and enrichment result summaries.
- `figures`: Fig. 3, function-region supplement PDFs and MAGMA/pathway figures.

## brain MRI

- `code/img_scripts`: MRI MR, coloc and brain-plot scripts.
- `code/main_figure`: Fig. 4 plotting scripts.
- `results/mr`: primary MR result tables.
- `results/coloc_mri`: coloc summary tables, PPH4 plot and brain icons.
- `results`: MRI manifest and Fig. 4 source tables.
- `figures`: Fig. 4 and related MRI/coloc PDFs.

## supplementary figure

- `code`: supplementary figure plotting scripts from the manuscript plus supporting GWAS, post-GWAS and MRI plotting scripts.
- `results`: supplementary tables.
- `figures`: generated supplementary figure PDFs.
