#------------------------------------------------------------
# 2021-07-14
# 
# Run this Rscript to reestimate betas using LDpred2
#   "X_LDpred2_QCed.Rscript"
# in the current working directory.
#   *LDpred2-inf: infinitesimal model
#   *LDpred2(-grid): grid of models
#
# Input: [1] trait name - prefix
#        [2] file name of the matched variants
#        [3] number of cores
# Output: a table with recaluclated betas
# 
# 
# Please perform QC before running LDpred2:
# (1) match GWAS summary stats with LD reference data 
#    (HapMap3 common variants)
# (2) match with GEL variants
# (3) additional QC on GWAS
#
# bsub -J ${trait} -q long -o bjob_LDpred2_${trait}.out -e bjob_LDpred2_${trait}.error -R"select[mem>10000] rusage[mem=10000]" -M10000 -n15 -R "span[hosts=1]" "/software/R-3.6.1/bin/Rscript X_LDpred2_QCed.Rscript ${trait} ${trait}_matched_variants.RDS 15"
#------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  stop("Please provide (1) trait, used as the prefix, (2) file name of the matched variants and (3) number of cores", call. = F)
}

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")
library(bigsnpr)
library(data.table)


trait = args[1]
matchedvariant = args[2]
NCORES = args[3]

write(paste0("*** GWAS Trait: ", trait), file = paste0(trait, "_LDpred2.log"))


## Load data
# LD panel provdied by LDpred2 software
map_ldref = readRDS("~/storage/Downloaded_data/LD_matrices/LDpred2/map.rds")

# matched variants
df_beta = readRDS(matchedvariant)
write(paste0("  Number of matched variants: ", nrow(df_beta)), 
      file = paste0(trait, "_LDpred2.log"), append = T)


## Run LDpred2
tmp = tempfile(tmpdir = paste0("tmp-data-", trait))

# calculate
for(chr in 1:22) {
  cat(chr, ".. ", sep = "")
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  
  corr_chr <- readRDS(paste0("~/storage/Downloaded_data/LD_matrices/LDpred2/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr_chr^2)   # this line was missing
    corr <- as_SFBM(corr_chr, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr_chr^2))   # this line was missing
    corr$add_columns(corr_chr, nrow(corr))
  }
}

saveRDS(ld, paste0(trait,"_interdata_ld.RDS"))
#ld = readRDS(paste0(trait,"_interdata_ld.RDS"))
saveRDS(corr, paste0(trait,"_interdata_correlation.RDS"))
#corr = readRDS(paste0(trait,"_interdata_correlation.RDS"))
write(paste0("  correlation has been calculated."), 
      file = paste0(trait, "_LDpred2.log"), append = T)


# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
#### is this correct?
#ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = length(ld),
                               chi2 = (beta / beta_se)^2,
                               sample_size = n_eff))
h2_est <- ldsc[["h2"]]
write(paste0("  LDSC heritability has been calculated: ", h2_est), 
      file = paste0(trait, "_LDpred2.log"), append = T)
write(ldsc, file = paste0(trait, "_LDpred2.log"), append = T)


## LDpred2-auto
#multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
#                               vec_p_init = seq_log(1e-4, 0.9, 30),
#                               ncores = NCORES)  # 13 hours
#beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)


## LDpred2-inf
beta_inf = snp_ldpred2_inf(corr, df_beta, h2 = h2_est) 
# beta output 
beta_inf_out = cbind(df_beta[, c("SNP_GEL","a0","a1")], beta_inf) 
write(paste0("  ...LDpred2-inf finished..."), 
      file = paste0(trait, "_LDpred2.log"), append = T)


## LDpred2-grid
h2_seq = round(h2_est * c(0.7, 1, 1.4), 4)
p_seq = signif(seq_log(1e-4, 1, length.out = 17), 2)
params = expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid = snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
beta_grid = as.data.frame(beta_grid) 
colnames(beta_grid) = paste(params$p,params$h2,params$sparse,sep='_')
write(paste0("  ...LDpred2-grid finished..."), 
      file = paste0(trait, "_LDpred2.log"), append = T)


# betas from the infinitesimal model and grid of models
betaall = cbind(beta_inf_out, beta_grid)
fwrite(betaall, file = "ldpred2_inf_grid_beta.txt", sep = "\t")
write(paste0("  ...Beta saved..."), 
      file = paste0(trait, "_LDpred2.log"), append = T)









