mc_doubletFinder_v3 <- function(sample, seurats, pN, pK_optimal, est_expected){
  result <- doubletFinder_v3(seurats[[sample]], PCs = 1:20, pN = pN, pK = pK_optimal[[sample]], nExp = est_expected[sample],
                             reuse.pANN = FALSE, sct = FALSE)
  result[["classifications"]] <- factor(result[[paste("DF.classifications", pN, pK_optimal[[sample]], est_expected[sample], sep = "_")]][,1],
                                        levels = c("Singlet", "Doublet"))
  return(result)
}
