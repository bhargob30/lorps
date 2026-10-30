\name{lorps-package}
\alias{lorps-package}
\alias{lorps}
\docType{package}
\title{Low-Rank Probit Stick-Breaking Multiple Imputation}
\description{
The main function \code{\link{lorps_mi}} fits a flexible normal-mixture model
to an incomplete numeric matrix, with information shared through latent row,
column, and mixture-component factors. The package compiles the sampler and
provides initialization, optional rank selection, and compact storage of
retained missing-entry draws. Use \code{\link{complete_lorps}} to construct
complete datasets and \code{\link{lorps_suggest_K}} to evaluate a sufficient
truncation level from the manuscript's prior bound.
}
\author{Bhargob Kakoty}
\keyword{package}
