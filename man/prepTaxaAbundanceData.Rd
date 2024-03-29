% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepTaxaAbundanceData.R
\name{prepTaxaAbundanceData}
\alias{prepTaxaAbundanceData}
\title{Prepares Taxa Abundance Data for Association analysis}
\usage{
prepTaxaAbundanceData(
  x,
  detection_threshold = 0.001,
  prevalence_threshold = 5/100
)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object}

\item{detection_threshold}{Remove taxa with less than this abundance in total dataset}

\item{prevalence_threshold}{Remove taxa with less than this prevalence in total dataset}
}
\value{
A tibble
}
\description{
Prepares Taxa Abundance Data for Association analysis
}
\details{
Prepares the \code{\link{phyloseq-class}} object to
         taxa abundance table for check associations
}
\examples{
library(biomeStats)
data("FuentesIliGutData")
ps <- FuentesIliGutData
gen_tab <- prepTaxaAbundanceData(ps,
                        detection_threshold = 0.01,
                        prevalence_threshold = 10/100)
head(gen_tab)

}
\references{
Ferreira JA, Fuentes S. (2020). Some comments on certain statistical
aspects of the study of the microbiome.
\emph{Briefings in bioinformatics} 21(4), pp.1487-1494.
\url{https://doi.org/10.1093/bib/bbz077}
}
\author{
Sudarshan A. Shetty
}
