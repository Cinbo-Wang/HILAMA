% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_hilama.R
\name{infer_ddlasso}
\alias{infer_ddlasso}
\title{Infer the coefficient of high dimensional linear model with latent confounder using Decorrelating and Debiasing approach.}
\usage{
infer_ddlasso(X, Y, ratio_singular_min = 1.5, Khat = NULL, score_ls = NULL)
}
\arguments{
\item{X}{A data matrix.}

\item{Y}{Continuous response vector.}

\item{ratio_singular_min}{Used to check the existence of hidden confounder. If \eqn{\frac{\lambda_k}{\lambda_{k+1}}} is always smaller than this value,
then no hidden confounder exists and only single debiased procedure is applied.}

\item{Khat}{The number of possible hidden confounders.}

\item{score_ls}{Object returned from esti_score_all() function.}
}
\description{
This function estimates and infers the coefficients of a high-dimensional linear model with a latent confounder using the ecorrelating and Debiasing approach proposed by Sun, Ma, and Xia (2023).
}
\references{
Sun, Y., Ma, L., & Xia, Y. (2023). A decorrelating and debiasing approach to simultaneous inference for high-dimensional confounded models. Journal of the American Statistical Association, 1-12.
}