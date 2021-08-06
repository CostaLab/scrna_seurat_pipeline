mc_pK_identification <- function(scrna_sample){
  sweep_res_list <- paramSweep_v3(scrna_sample, PCs = 1:20, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- custom_findpk(sweep_stats)
  return(bcmvn)
}
