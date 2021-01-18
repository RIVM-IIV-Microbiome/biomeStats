#' @title Create variable names
#' @param data.set dataset consisting of variables
#' @import graphics
#' @importFrom stats aggregate cor na.omit
#' @import utils
#' @import coin
#' @import dplyr
#' @import ggplot2
#' @import microbiome
#' @import phyloseq
#' @export

create.names <- function(data.set) {
  names.of.variables <- NULL
  for (j in (1:ncol(data.set))) {
    names.of.variables.j <- paste(unlist(strsplit(x = names(data.set)[j], split = "\\.")), collapse = " ")
    names.of.variables <- c(names.of.variables, paste(names.of.variables.j, sep = ""))
  }
  return(names.of.variables)
}

############################################################################################################

#' @title show.missing.by.variable
#' @param data.set dataset consisting of variables
#' @param plot Plot
#' @export

show.missing.by.variable <- function(data.set, plot) {
  missing.by.variable <- rep(0, ncol(data.set))
  names(missing.by.variable) <- names(data.set)
  for (j in (1:length(missing.by.variable))) {
    missing.by.variable[j] <- sum(is.na(data.set[, j]))
  }
  missing.by.variable <- missing.by.variable / nrow(data.set)
  if (plot) {
    par(mar = c(5, 12, 1, 2))
    barplot(missing.by.variable,
      horiz = TRUE, las = 1, xlab = "Proportion of missing observations",
      col = "lavender", cex.axis = 1.25, cex.lab = 1.25
    )
  }
  return(missing.by.variable)
}


############################################################################################################
#' @title check.strata
#' @param data.to.test dataset variable to specifically test
#' @param s Stratum that is created to check for in the analysis
#' @export
check.strata <- function(s, data.to.test) {
  stratum <- NULL
  aux.data.to.test <- subset(data.to.test, stratum == s)
  check <- min(length(unique(aux.data.to.test[, 1])),
               length(unique(aux.data.to.test[, 2])))
  return(check)
}
#

############################################################################################################

#' @title Wilson.interval
#' @param frequency frequency
#' @param n Number must be numeric
#' @param confidence confidence value
#' @importFrom stats as.formula qnorm sd
#' @importFrom grDevices dev.off pdf
#' @export
#
Wilson.interval <- function(frequency,n,confidence){
  kappa <- qnorm(1-(1-confidence)/2)
  estimate <- frequency/n
  aux1 <- (frequency+(kappa^2)/2+kappa*sqrt(n)*
             sqrt(estimate*(1-estimate)+(kappa^2)/(4*n)))/(n+(kappa^2))
  aux2 <- (frequency+(kappa^2)/2-kappa*sqrt(n)*sqrt(estimate*(1-estimate)+(kappa^2)/(4*n)))/
    (n+(kappa^2))
  return(c(max(0.0,aux2),min(1.0,aux1)))
}
#
############################################################################################################
