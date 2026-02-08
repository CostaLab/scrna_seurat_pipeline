mc_doubletFinder <- function(sample, seurats, pN, pK_optimal, est_expected){
  pK_s <- pK_optimal[[sample]][1L]
  nExp_s <- est_expected[sample]
  result <- doubletFinder(seurats[[sample]], PCs = 1:20, pN = pN, pK = pK_s, nExp = nExp_s,
                          sct = FALSE)
  colname <- paste("DF.classifications", pN, pK_s, nExp_s, sep = "_")
  result[["classifications"]] <- factor(result[[colname]][,1], levels = c("Singlet", "Doublet"))
  return(result)
}
