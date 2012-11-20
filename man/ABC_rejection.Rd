\name{ABC_rejection}
\alias{ABC_rejection}
\title{
  Rejection sampling scheme for ABC}
\description{
  This function launches a series of \code{nb_simul} model simulations with model parameters drawn in the prior distribution specified in \code{prior_matrix}.
}
\usage{
  ABC_rejection(model, prior_matrix, nb_simul, use_seed = TRUE, seed_count = 0, n_cluster=1, progress_bar=FALSE)
}
\arguments{
  \item{model}{
    a \code{R} function implementing the model to be simulated. It must take as arguments an array containing a seed value (if \code{use_seed=TRUE}) and the model parameters and returning an array of summary statistics obtained at the end of the simulation.
    A tutorial is provided in the package's vignette to dynamically link a binary code to a \code{R} function.  Users may alternatively wish to wrap their binary executables using the provided functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}. The use of these functions is associated with slightly different constraints on the design of the binary code (see \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).}
  \item{prior_matrix}{
    a matrix or data frame containing the ranges of the prior distribution for each of the \eqn{p} model parameters. It has \eqn{p} lines and 2 columns. The first (second) column contains the lower (upper) bound of the prior distribution of the \eqn{p} parameters.}
  \item{nb_simul}{
    a positive integer equal to the desired number of simulations of the model.}
  \item{use_seed}{
    logical. If \code{TRUE} (default), \code{ABC_rejection} provides as input to the function \code{model} an array containing an integer seed value and the model parameters used for the simulation.
    The seed value should be used by \code{model} to initialize its pseudo-random number generators (if \code{model} is stochastic).
    If \code{FALSE} (default), \code{ABC_rejection} provides as input to the function \code{model} an array containing solely the model parameters used for the simulation.}
  \item{seed_count}{
    a positive integer, the initial seed value provided to the function \code{model} (if \code{use_seed=TRUE}). This value is incremented by 1 at each call of the function \code{model}.}
  \item{n_cluster}{
    a positive integer. If larger than 1 (the default value), \code{ABC_rejection} will launch \code{model} simulations in parallel on \code{n_cluster} cores of the computer.}
  \item{progress_bar}{
    logical, \code{FALSE} by default. If \code{TRUE}, \code{ABC_rejection} will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}
}
%\details{}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{weights}{
  	The weights of the different \code{model} simulations. In the standard rejection scheme, all \code{model} simulations have the same weights.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{computime}{
  	The computing time to perform the simulations.}
}
\references{
  Pritchard, J.K., and M.T. Seielstad and A. Perez-Lezaun and
  M.W. Feldman (1999) Population growth of human Y chromosomes: a study
  of Y chromosome microsatellites. \emph{Molecular Biology and
  Evolution}, \bold{16}, 1791--1798.
}
\author{Franck Jabot, Thierry Faure and Nicolas Dumoulin}
\seealso{
  \code{\link{binary_model}}, \code{\link{binary_model_cluster}}, \code{\link{ABC_sequential}}, \code{\link{ABC_mcmc}}
}
\examples{
  \dontrun{
    ## artificial example to show how to use the ABC_rejection function
    set.seed(1)
    n=10
    priormatrix=cbind(c(500,3,-2.3,1,-25,-0.7),c(500,5,1.6,1,125,3.2))
    priormatrix
    ABC_rej<-ABC_rejection(model=trait_model, prior_matrix=priormatrix, nb_simul=n, n_cluster=1, progress_bar=FALSE)
    ABC_rej

    ## To use the ABC_rejection function with multiple cores, simply change the value of n_cluster:
    ABC_rejb<-ABC_rejection(model=trait_model, prior_matrix=priormatrix, nb_simul=n, n_cluster=2, progress_bar=FALSE)
  }
}
\keyword{abc}
\keyword{model}
\keyword{inference}
