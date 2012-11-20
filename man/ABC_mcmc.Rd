\name{ABC_mcmc}
\alias{ABC_mcmc}
\title{
	Coupled to MCMC schemes for ABC
}
\description{
	This function implements three different algorithms to perform coupled to MCMC ABC.
}
\usage{
	ABC_mcmc(method, model, prior_matrix, n_obs, n_between_sampling, summary_stat_target, n_cluster = 1, ...)
}
\arguments{
  \item{method}{
	a character string indicating the ABC-MCMC algorithm to be used. Possible values are \code{"Marjoram_original"}, \code{"Marjoram"} and \code{"Wegmann"}.
	Note that the method \code{"Marjoram_original"} cannot be used with multiple cores.}
  \item{model}{
    	a \code{R} function implementing the model to be simulated. It must take as arguments an array containing a seed value (if \code{use_seed=TRUE}) and the model parameters and returning an array of summary statistics obtained at the end of the simulation.
	A tutorial is provided in the package's vignette to dynamically link a binary code to a \code{R} function.  Users may alternatively wish to wrap their binary executables using the provided functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}. The use of these functions is associated with slightly different constraints on the design of the binary code (see \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).}
  \item{prior_matrix}{
    	a matrix or data frame containing the ranges of the prior distribution for each of the \eqn{p} model parameters. It has \eqn{p} lines and 2 columns. The first (second) column contains the lower (upper) bound of the prior distribution of the \eqn{p} parameters.}
  \item{n_obs}{
    	a positive integer equal to the desired number of sampled points along the MCMC.}
  \item{n_between_sampling}{
    	a positive integer equal to the desired spacing between sampled points along the MCMC.}
  \item{summary_stat_target}{
	a vector containing the targeted (observed) summary statistics.}
  \item{n_cluster}{
    a positive integer. If larger than 1 (the default value), \code{ABC_mcmc} will launch \code{model} simulations in parallel on \code{n_cluster} cores of the computer.}
  \item{\dots}{Additional arguments can be passed depending on the choosen method (see below)}
}

\section{Additional paramaters}{
  Depending on the choosen method, you can specify the following arguments:
  \describe{
    \item{dist_max}{
	  a positive number, used when \code{method} is \code{"Marjoram_original"}.
	  This is the tolerance threshold used during the MCMC.}
    \item{tab_normalization}{
	  a vector of the same length as \code{summary_stat_target}, used when \code{method} is \code{"Marjoram_original"}.
	  Each element contains a positive number by which each summary statistics must be divided before the computation of the Euclidean distance between simulations and data.}
    \item{proposal_range}{
	  a vector of the same length as the number of model parameters, used when \code{method} is \code{"Marjoram_original"}.
	  Each element contains a positive number defining the range of MCMC jumps for each model parameter.}
    \item{n_calibration}{
	  a positive integer, used when \code{method} is \code{"Marjoram"} or \code{"Wegmann"}.
	  This is the number of simulations performed during the calibration step.
	  Default value is 10000.}
    \item{tolerance_quantile}{
	  a positive number between 0 and 1 (strictly), used when \code{method} is \code{"Marjoram"} or \code{"Wegmann"}.
	  This is the percentage of simulations retained during the calibration step to determine the tolerance threshold to be used during the MCMC.
	  Default value is 0.01.}
    \item{proposal_phi}{
	  a positive number, used when \code{method} is \code{"Marjoram"} or \code{"Wegmann"}.
	  This is a scaling factor defining the range of MCMC jumps.
	  Default value is 1.}
    \item{numcomp}{
	  a positive integer, used when \code{method} is \code{"Wegmann"}.
	  This is the number of components to be used for PLS transformations.
	  Default value is 0 which encodes that this number is equal to the number of summary statistics.}
    \item{use_seed}{
	  logical. If \code{TRUE} (default), \code{ABC_mcmc} provides as input to the function \code{model} an array containing an integer seed value and the model parameters used for the simulation.
	  The seed value should be used by \code{model} to initialize its pseudo-random number generators (if \code{model} is stochastic).
	  If \code{FALSE} (default), \code{ABC_mcmc} provides as input to the function \code{model} an array containing solely the model parameters used for the simulation.}
    \item{seed_count}{
	  a positive integer, the initial seed value provided to the function \code{model} (if \code{use_seed=TRUE}). This value is incremented by 1 at each call of the function \code{model}.}
    \item{verbose}{
	  logical. \code{FALSE} by default. If \code{TRUE}, \code{ABC_mcmc} writes in the current directory intermediary results at the end of the calibration step of the algorithm.}
    \item{progress_bar}{
      logical, \code{FALSE} by default. If \code{TRUE}, \code{ABC_mcmc} will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}
  }
}
\details{
  See the package's vignette for details on ABC-MCMC.
}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{dist}{
  	The distance of the simulations to the data.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations of the initial step.
	These values are used to normalize the summary statistics before the computation of the Euclidean distance between simulations and data.
	If \code{method} is \code{"Marjoram_original"}, this is equal to \code{tab_normalization}.
	If \code{method} is \code{"Wegmann"}, this is not provided.}
  \item{epsilon}{
  	The final maximal distance between simulations and data in the retained sample of particles.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{n_between_sampling}{
  	The spacing between two sampled points in the MCMC.}
  \item{computime}{
  	The computing time to perform the simulations.}
  \item{min_stats}{
  	The minimal values of each summary statistics during the calibration step, given when \code{method} is \code{"Wegmann"}.}
  \item{max_stats}{
  	The maximal values of each summary statistics during the calibration step, given when \code{method} is \code{"Wegmann"}.}
  \item{lambda}{
  	The lambda values of the Box-Cox transformation, given when \code{method} is \code{"Wegmann"}.}
  \item{geometric_mean}{
  	The geometric means, given when \code{method} is \code{"Wegmann"}.}
  \item{boxcox_mean}{
  	The means of Box-Cox transforms, given when \code{method} is \code{"Wegmann"}.}
  \item{boxcox_sd}{
  	The standard deviations of Box-Cox transforms, given when \code{method} is \code{"Wegmann"}.}
  \item{pls_transform}{
  	The matrix of PLS transformation, given when \code{method} is \code{"Wegmann"}.}
  \item{n_component}{
  	The number of used components for the PLS transformation, given when \code{method} is \code{"Wegmann"}.}
}
\references{
	Marjoram, P., Molitor, J., Plagnol, V. and Tavar\'e, S. (2003) Markov chain Monte Carlo without likelihoods. \emph{PNAS}, \bold{100}, 15324--15328.

	Wegmann, D., Leuenberger, C. and Excoffier, L. (2009) Efficient approximate Bayesian computation coupled with Markov chain Monte Carlo without likelihood. \emph{Genetics}, \bold{182}, 1207-1218.
}
\author{Franck Jabot, Thierry Faure and Nicolas Dumoulin}
\seealso{
  \code{\link{binary_model}}, \code{\link{binary_model_cluster}}, \code{\link{ABC_rejection}}, \code{\link{ABC_sequential}}
}
\examples{
  \dontrun{
    priormatrix=cbind(c(500,3,-2.3,1,-25,-0.7),c(500,5,1.6,1,125,3.2))
    sum_stat_obs=c(100,2.5,20,30000)
    n=10
    nbetweensampling=1
    distmax=8
    tabnormalization=c(50,1,20,10000)
    proposalrange=c(0,1,0.5,0,50,1)

    ## artificial example to perform the Marjoram et al. (2003)'s method.
    ##
    ABC_Marjoram_original<-ABC_mcmc(method="Marjoram_original", model=trait_model, prior_matrix=priormatrix,
    n_obs=n, n_between_sampling=nbetweensampling, summary_stat_target=sum_stat_obs, n_cluster=1, dist_max=distmax,
    tab_normalization=tabnormalization, proposal_range=proposalrange)
    ABC_Marjoram_original

    ## artificial example to perform the Marjoram et al. (2003)'s method, with modifications drawn from Wegmann et al. (2009) without Box-Cox and PLS transformations.
    ##
    ncalib=10
    tolquantile = 0.5
    ABC_Marjoram<-ABC_mcmc(method="Marjoram", model=trait_model, prior_matrix=priormatrix,
    n_obs=n, n_between_sampling=nbetweensampling, summary_stat_target=sum_stat_obs, n_cluster=1,
    n_calibration=ncalib, tolerance_quantile=tolquantile)
    ABC_Marjoram


    ## artificial example to perform the Wegmann et al. (2009)'s method.
    ##
    ABC_Wegmann<-ABC_mcmc(method="Wegmann", model=trait_model, prior_matrix=priormatrix,
    n_obs=n, n_between_sampling=nbetweensampling, summary_stat_target=sum_stat_obs, n_cluster=1,
    n_calibration=ncalib, tolerance_quantile=tolquantile)
    ABC_Wegmann
  }
}
\keyword{abc}
\keyword{model}
\keyword{inference}
\keyword{mcmc}
