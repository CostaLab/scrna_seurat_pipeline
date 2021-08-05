mc_est_expected <- function(x, doublet_rate, scrnas){
  # Assuming 7.5% doublet formation rate - tailor for your dataset
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
  nExp_poi <- round(doublet_rate[x]*nrow(scrnas[[x]]@meta.data))
  return(nExp_poi)
}
