#------------------------------------------------------------
# 2021-07-16
# 
# Run this Rscript to reestimate betas using LDpred2
#   "X_LDpred2_QCed.Rscript"
# in the current working directory.
#   *LDpred2-inf: infinitesimal model
#   *LDpred2(-grid): grid of models
#
# Input: [1] trait name - prefix
#        [2] file name of the matched variants
# Output: a table with recaluclated betas
# 
# 
# Please perform QC before running LDpred2:
# (1) match GWAS summary stats with LD reference data 
#    (HapMap3 common variants)
# (2) match with GEL variants
# (3) additional QC on GWAS
#
#
# 2021-07-16: reinstalled bigsnpr
# 2021-07-19: solved the issue when using multiple cores 
#  in LDpred-auto, run this line before bsub
#$ export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/3.6:${R_LIBS_USER}"
#
# bsub -J ${trait} -q long -o bjob_LDpred2_${trait}.out -e bjob_LDpred2_${trait}.error -R"select[mem>60000] rusage[mem=60000]" -M60000 -n8 -R "span[hosts=1]" "/software/R-3.6.1/bin/Rscript X_LDpred2_QCed.Rscript ${trait} ${trait}_matched_variants.RDS"
#
# #run LDpred-auto only
# bsub -J ${trait} -q long -o bjob_LDpred2_${trait}_auto.out -e bjob_LDpred2_${trait}_auto.error -R"select[mem>60000] rusage[mem=60000]" -M60000 -n8 -R "span[hosts=1]" "/software/R-3.6.1/bin/Rscript X_LDpred2_auto_QCed.Rscript ${trait} ${trait}_matched_variants.RDS"
#
# using 7 cores, ~2-3 hours for LDpred-inf and -grid; ~7-8 hours for -auto
#------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Please provide (1) trait, used as the prefix and (2) file name of the matched variants", call. = F)
}

#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")
library(bigsnpr)
library(data.table)
#library(doParallel)

trait = args[1]
matchedvariant = args[2]

write(paste0("*** GWAS Trait: ", trait), file = paste0(trait, "_LDpred2.log"))



NCORES = nb_cores()
write(paste0("The max amount of cores detected is: ", NCORES), file = paste0(trait, "_LDpred2.log"), append = T)



## Load data
# LD panel provdied by LDpred2 software
map_ldref = readRDS("~/storage/Downloaded_data/LD_matrices/LDpred2/map.rds")

# matched variants
df_beta = readRDS(matchedvariant)
write(paste0("  Number of matched variants: ", nrow(df_beta)), 
      file = paste0(trait, "_LDpred2.log"), append = T)


## Run LDpred2
#tmp = tempfile(tmpdir = paste0("tmp-data-", trait))

# calculate correlation
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
    #corr <- as_SFBM(corr_chr, tmp)
    corr = as_SFBM(corr_chr, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

saveRDS(corr, paste0(trait,"_interdata_correlation.RDS"))
write(paste0("  correlation has been calculated."), 
      file = paste0(trait, "_LDpred2.log"), append = T)


## Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                               chi2 = (beta / beta_se)^2,
                               sample_size = n_eff))

h2_est <- ldsc[["h2"]]
write(paste0("  LDSC heritability has been calculated: ", h2_est), 
      file = paste0(trait, "_LDpred2.log"), append = T)
#        int     int_se         h2      h2_se
write(ldsc, file = paste0(trait, "_LDpred2.log"), append = T)



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
fwrite(betaall, file = paste0(trait, "_LDpred2_inf_grid_beta.txt"), sep = "\t")
write(paste0("  ...Beta saved..."), 
      file = paste0(trait, "_LDpred2.log"), append = T)



## LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)  

saveRDS(multi_auto, paste0(trait,"_LDpred2_auto.RDS"))
write(paste0("  ...LDpred2-auto finished..."), 
      file = paste0(trait, "_LDpred2.log"), append = T)

## beta output 
#beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
#saveRDS(beta_auto, paste0(trait,"_LDpred2_auto.RDS"))










