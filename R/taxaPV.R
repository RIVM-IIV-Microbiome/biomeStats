#' Calculate Proportional Variability For Taxa
#'
#' @name taxaPV
#'
#' @details Proportional Variability of taxa using relative abundance.
#'          The Proportional Variability is a non-parametric measure of variability.
#'          \url{https://doi.org/10.1371/journal.pone.0084074}. Here, I have
#'          incorporated subsampling approach to quantify Proportional Variability
#'          in microbiome studies.
#'
#' @param x A phyloseq object
#'
#' @param subsample Logical. Default=TRUE.
#'
#' @param lower_conf Lower confidence. Default=0.025
#' @param upper_conf Lower confidence. Default=0.975
#' @param sampling_n Number of subsampling. Default=9
#' @param round_pv_vals Return PV values by rounding to certain value. Default=3
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' library(microbiome)
#' data("FuentesIliGutData")
#' # Check for control samples how variable are core taxa
#' ps1 <- filterSampleData(FuentesIliGutData, ILI == "C") %>%
#'      microbiome::transform("compositional") %>%
#'      core(0.01, 75/100)
#' taxa_fd <- taxaPV(ps1)
#' taxa_fd
#'
#' @return A tibble with taxa ids, taxonomic information,
#'         two-group prevalence and fold change values.
#'
#' @author Original Author: Leo Lahti. Adapted by Sudarshan A. Shetty
#'
#' @references
#' Heath JP, Borowski P. (2013). Quantifying Proportional Variability.
#' \emph{PLoS ONE} 8(12): e84074.
#' \url{https://doi.org/10.1371/journal.pone.0084074}
#'
#' Heath J. (2006) Quantifying temporal variability in population abundances.
#' \emph{Oikos} 115, no. 3 (2006): 573-581.
#' \url{https://doi.org/10.1111/j.2006.0030-1299.15067.x}
#'
#' Shetty SA (2021). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#'
#' @export

taxaPV <- function(x,
                   subsample = TRUE,
                   lower_conf=0.025,
                   upper_conf=0.975,
                   sampling_n =9,
                   round_pv_vals = 3){

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if(subsample){
    all.pv <- .taxa_pv(x)
    sampled.pv <- .taxa_pv_subsampled(x,
                                      lower_conf = lower_conf,
                                      upper_conf = upper_conf,
                                      sampling_n = sampling_n)
    comb.pv <- cbind(all.pv, sampled.pv)
    comb.pv <- round(comb.pv, digits = round_pv_vals)
    comb.pv <- comb.pv %>%
      tibble::rownames_to_column("FeatureID")
    return(comb.pv)
  }

  all.pv <- .taxa_pv(x)
  all.pv <- round(all.pv, digits = round_pv_vals)
  all.pv <- all.pv %>%
    tibble::rownames_to_column("FeatureID")
  return(all.pv)
}


#' @importFrom phyloseq sample_names prune_samples nsamples
#' @importFrom stats quantile
.taxa_pv_subsampled <- function(x,
                                lower_conf = lower_conf,
                                upper_conf = upper_conf,
                                sampling_n = sampling_n){
  rand_sams <- ps.sub <- txvp <- sx <- cis_df <- sxi <- NULL
  s <- NULL
  for (i in seq_len(sampling_n)) {
    size_80 <- round(0.8*nsamples(x))
    rand_sams <- sample(sample_names(x), size= size_80, replace = TRUE)
    ps.sub <- prune_samples(sample_names(x) %in% rand_sams, x)
    #ps.sub <- prune_taxa(taxa_sums(ps.sub) > 0, ps.sub)
    txvp <- .taxa_pv(ps.sub)
    #rownames(txvp) <- txvp$Taxa
    s[[i]] <- txvp
  }
  sx <- suppressMessages(dplyr::bind_cols(s))
  tx.pv.lci <- apply(sx, 1, function(x){
    quantile(x,lower_conf,na.rm=TRUE)
  })
  tx.pv.uci <- apply(sx, 1, function(x){
    quantile(x,upper_conf,na.rm=TRUE)
  })

  cis <- data.frame(LowerCI = tx.pv.lci,
                    UpperCI = tx.pv.uci)
  return(cis)
}

#' @importFrom microbiome abundances
## Overall PV calculations
.taxa_pv <- function(x) {

  sp_tab_tax <- NULL
  # proportional variability index (PV)
  sp_tab_tax <- t(abundances(x))
  pv.d <- apply(sp_tab_tax, 2, .cal.pv)
  pv.d <- data.frame(pv = pv.d)
  return(pv.d)
}

# Modified PV Formula to account for sparsity
.cal.pv <- function(x){
  k = (1/100)*mean(x) # get low pseudocount
  pv_val_tax <- .pv.index(x + k)
  return(pv_val_tax)
}

# Main PV Formula
#' @importFrom utils combn
.pv.index <- function (Z){
  n = length(Z)
  pairs = combn(Z,2)
  min_z = apply(pairs,2, min)
  max_z = apply(pairs,2, max)
  z = 1- (min_z/max_z)
  PV=2*sum(z)/(n*(n-1))
  return(PV)
}
