% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxaPV.R
\name{taxaPV}
\alias{taxaPV}
\title{Calculate Proportional Variability For Taxa}
\usage{
taxaPV(
  x,
  subsample = TRUE,
  lower_conf = 0.025,
  upper_conf = 0.975,
  sampling_n = 9,
  round_pv_vals = 3
)
}
\arguments{
\item{x}{A phyloseq object}

\item{subsample}{Logical. Default=TRUE.}

\item{lower_conf}{Lower confidence. Default=0.025}

\item{upper_conf}{Lower confidence. Default=0.975}

\item{sampling_n}{Number of subsampling. Default=9}

\item{round_pv_vals}{Return PV values by rounding to certain value. Default=3}
}
\value{
A tibble with taxa ids, taxonomic information,
        two-group prevalence and fold change values.
}
\description{
Calculate Proportional Variability For Taxa
}
\details{
Proportional Variability of taxa using relative abundance.
         The Proportional Variability is a non-parametric measure of variability.
         \url{https://doi.org/10.1371/journal.pone.0084074}. Here, I have
         incorporated subsampling approach to quantify Proportional Variability
         in microbiome studies.
}
\examples{
library(biomeUtils)
library(dplyr)
library(microbiome)
data("FuentesIliGutData")
# Check for control samples how variable are core taxa
ps1 <- filterSampleData(FuentesIliGutData, ILI == "C") \%>\%
     microbiome::transform("compositional") \%>\%
     core(0.01, 75/100)
taxa_fd <- taxaPV(ps1)
taxa_fd

}
\references{
Heath JP, Borowski P. (2013). Quantifying Proportional Variability.
\emph{PLoS ONE} 8(12): e84074.
\url{https://doi.org/10.1371/journal.pone.0084074}

Heath J. (2006) Quantifying temporal variability in population abundances.
\emph{Oikos} 115, no. 3 (2006): 573-581.
\url{https://doi.org/10.1111/j.2006.0030-1299.15067.x}

Shetty SA (2021). Utilities for microbiome analytics.
\url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
}
\author{
Original Author: Leo Lahti. Adapted by Sudarshan A. Shetty
}
