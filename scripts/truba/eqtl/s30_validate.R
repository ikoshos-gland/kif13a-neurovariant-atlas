# Reference-implementation check (CRAN coloc 5.2.3 + susieR 0.14.2) of the Python coloc.abf, and
# SuSiE-based multi-signal colocalisation for KIF13A. Compute node only.
#   1) coloc's bundled example (coloc_test_data D1,D2): R output + inputs exported for Python
#   2) the real KIF13A inputs written by s20_coloc.py export: coloc.abf in R on identical SNPs
#   3) runsusie on eQTL z (1000G EUR LD) and on BIG40 L3 z; coloc.susie where both have CSs
suppressMessages({
  .libPaths(c("/arf/scratch/avural/kif13a_eqtl/Rlib", .libPaths()))
  library(coloc); library(susieR)
})
out <- "/arf/scratch/avural/kif13a_eqtl/results/r_validation"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
inp <- "/arf/scratch/avural/kif13a_eqtl/results/r_input"
cat("coloc", as.character(packageVersion("coloc")), " susieR", as.character(packageVersion("susieR")), "\n")

## 1) bundled example
data(coloc_test_data)
D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
for (nm in c("D1", "D2")) {
  D <- coloc_test_data[[nm]]
  write.table(data.frame(snp = D$snp, beta = D$beta, varbeta = D$varbeta, MAF = D$MAF),
              file.path(out, paste0("example_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(nm, "type", D$type, "N", D$N, "sdY", ifelse(is.null(D$sdY), NA, D$sdY), "\n")
}
r <- coloc.abf(D1, D2)
write.table(t(r$summary), file.path(out, "example_D1D2_R_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# p-value route on the same example (forces the MAF/N branch)
D1p <- list(pvalues = pnorm(-abs(D1$beta / sqrt(D1$varbeta))) * 2, MAF = D1$MAF, N = D1$N, type = "quant", sdY = D1$sdY, snp = D1$snp)
D2p <- list(pvalues = pnorm(-abs(D2$beta / sqrt(D2$varbeta))) * 2, MAF = D2$MAF, N = D2$N, type = "quant", sdY = D2$sdY, snp = D2$snp)
rp <- coloc.abf(D1p, D2p)
write.table(t(rp$summary), file.path(out, "example_D1D2_pvalue_R_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("example N:", D1$N, D2$N, " sdY:", D1$sdY, D2$sdY, "\n")
print(r$summary); print(rp$summary)

## 2) real inputs
res <- list(); sus <- list()
for (f in list.files(inp, pattern = "\\.tsv$")) {
  tag <- sub("\\.tsv$", "", f)
  d <- read.delim(file.path(inp, f))
  ld <- as.matrix(read.table(file.path(inp, paste0(tag, ".ld"))))
  dimnames(ld) <- list(d$snp, d$snp)
  e <- list(beta = d$beta_e, varbeta = d$se_e^2, sdY = 1, type = "quant", snp = d$snp)
  if (all(is.na(d$beta))) {
    g <- list(pvalues = d$p, MAF = d$maf_ref, N = median(d$n), sdY = 1, type = "quant", snp = d$snp)
  } else {
    g <- list(beta = d$beta, varbeta = d$se^2, sdY = 1, type = "quant", snp = d$snp)
  }
  rr <- suppressMessages(coloc.abf(g, e))
  res[[tag]] <- data.frame(input = tag, t(rr$summary))
  ## 3) SuSiE (signed z only)
  n_e <- median(1 / (2 * d$maf_ref * (1 - d$maf_ref) * d$se_e^2))
  se <- tryCatch(susie_rss(z = d$z_e_alt, R = ld, n = round(n_e), L = 5, coverage = 0.95,
                           min_abs_corr = 0.5, max_iter = 500, check_prior = FALSE),
                 error = function(err) NULL)
  ncs_e <- if (is.null(se)) NA else length(se$sets$cs)
  cs_e <- if (is.null(se) || ncs_e == 0) "" else paste(sapply(se$sets$cs, function(ix) d$rsid[ix[which.max(se$pip[ix])]]), collapse = ",")
  ncs_g <- NA; cs_g <- ""; h4s <- ""
  if (!all(is.na(d$beta))) {
    sg <- tryCatch(susie_rss(z = d$z_g_alt, R = ld, n = round(median(d$n)), L = 5, coverage = 0.95,
                             min_abs_corr = 0.5, max_iter = 500, check_prior = FALSE),
                   error = function(err) NULL)
    ncs_g <- if (is.null(sg)) NA else length(sg$sets$cs)
    cs_g <- if (is.null(sg) || ncs_g == 0) "" else paste(sapply(sg$sets$cs, function(ix) d$rsid[ix[which.max(sg$pip[ix])]]), collapse = ",")
    if (!is.null(se) && !is.null(sg) && ncs_e > 0 && ncs_g > 0) {
      cs <- suppressMessages(coloc.susie(sg, se))
      h4s <- paste(sprintf("%s~%s:H3=%.3f,H4=%.3f", cs$summary$hit1, cs$summary$hit2, cs$summary$PP.H3.abf, cs$summary$PP.H4.abf), collapse = ";")
    }
  }
  sus[[tag]] <- data.frame(input = tag, n_eff_eqtl = round(n_e), eqtl_cs = ncs_e, eqtl_cs_leads = cs_e,
                           gwas_cs = ncs_g, gwas_cs_leads = cs_g, coloc_susie = h4s)
  cat(tag, "nsnps", rr$summary["nsnps"], "H3", signif(rr$summary["PP.H3.abf"], 4), "H4", signif(rr$summary["PP.H4.abf"], 4),
      "| susie eQTL CS", ncs_e, cs_e, "| GWAS CS", ncs_g, cs_g, h4s, "\n")
}
write.table(do.call(rbind, res), file.path(out, "real_coloc_abf_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, sus), file.path(out, "real_susie_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
