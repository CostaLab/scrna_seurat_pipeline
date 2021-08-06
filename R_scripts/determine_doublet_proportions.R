determine_doublet_proportions <- function(scrna, doublet_lst){
  if(length(doublet_lst) == 1 & is.numeric(doublet_lst)){
      doublet_formation_rate <- rep(doublet_lst, length(unique(scrna$name)))
      names(doublet_formation_rate) <- names(doublet_lst)
    }
    if(length(doublet_lst) == 0 & is.null(doublet_lst)){
      number_of_cells <- table(scrna$name)
      doublet_formation_rate <- model_recovered(number_of_cells)
      names(doublet_formation_rate) <- names(number_of_cells)
    }
    if(length(doublet_lst) > 1){
      doublet_formation_rate <- c()
      for(i in 1:length(doublet_lst)){
        if(is.null(doublet_lst[[i]])){
          cells_sample <- table(scrna$name)[names(doublet_lst)[[i]]]
          doublet_formation_rate <- c(doublet_formation_rate, model_recovered(cells_sample))
        } else if(grepl(".csv$", doublet_lst[[i]])){
          cells_sample <- read.csv(doublet_lst[[i]])
          cells_sample <- as.numeric(gsub(",", "", cells_sample[1,1]))
          doublet_formation_rate <- c(doublet_formation_rate, model_recovered(cells_sample))
        } else if(is.numeric(doublet_lst[[i]]) & doublet_lst[[i]] < 1) {
            doublet_formation_rate <- c(doublet_formation_rate, doublet_lst[[i]])
        } else if(is.numeric(doublet_lst[[i]]) & doublet_lst[[i]] >= 1) doublet_formation_rate <- c(doublet_formation_rate, model_recovered(doublet_lst[[i]]))
      }
    }
  names(doublet_formation_rate) <- names(doublet_lst)
  return(doublet_formation_rate)
}
