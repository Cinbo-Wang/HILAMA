% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_hilama.R
\name{get_nie_pvalue_screen_JST}
\alias{get_nie_pvalue_screen_JST}
\title{Get the p-value threshold for the mediation effect test using Joint Significance Test (JST).}
\usage{
get_nie_pvalue_screen_JST(
  pvalue_Theta_mat,
  pvalue_beta,
  screen_N = NULL,
  alpha = 0.1
)
}
\arguments{
\item{pvalue_Theta_mat}{p-value matrix of exposures to mediators.}

\item{pvalue_beta}{p-value matrix of mediators to outcome.}

\item{screen_N}{Number of pairs used for multiple testing.}

\item{alpha}{The nominal False Discovery Rate level.}
}
\value{
\item{pvalue_screen_JST_mat}{p-value matrix of NIE using JST. The value of 1 indicates that the item has been screened out.}
\item{pvalue_cutoff}{p-value threshold that controls the FDR at `\alpha` level.}
}
\description{
Get the p-value threshold for the mediation effect test using Joint Significance Test (JST).
}
