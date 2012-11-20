\name{ABC_sequential}
\alias{ABC_sequential}
\title{
	Sequential sampling schemes for ABC}
\description{
	This function implements four different algorithms to perform sequential sampling schemes for ABC.
	Sequential sampling schemes consist in sampling initially model parameters in the prior distribution, just like in a standard rejection-based ABC, in order to obtain a rough posterior distribution of parameter values, and in subsequently sampling close to this rough posterior distribution to refine it.
	Sequential sampling schemes have been shown to be more efficient than standard rejection-based procedures.
}
\usage{
	ABC_sequential(method, model, prior_matrix, nb_simul, summary_stat_target, n_cluster = 1,...)
}

\arguments{
  \item{method}{
	a character string indicating the sequential algorithm to be used. Possible values are \code{"Beaumont"}, \code{"Drovandi"}, \code{"Delmoral"} and \code{"Lenormand"}.}
  \item{model}{
    	a \code{R} function implementing the model to be simulated. It must take as arguments an array containing a seed value (if \code{use_seed=TRUE}) and the model parameters and returning an array of summary statistics obtained at the end of the simulation.
	A tutorial is provided in the package's vignette to dynamically link a binary code to a \code{R} function.  Users may alternatively wish to wrap their binary executables using the provided functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}. The use of these functions is associated with slightly different constraints on the design of the binary code (see \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).}
  \item{prior_matrix}{
    	a matrix or data frame containing the ranges of the prior distribution for each of the \eqn{p} model parameters. It has \eqn{p} lines and 2 columns. The first (second) column contains the lower (upper) bound of the prior distribution of the \eqn{p} parameters.}
  \item{nb_simul}{
    	a positive integer equal to the desired number of simulations of the model below the tolerance threshold when \code{method} is \code{"Beaumont"}, \code{"Drovandi"} and \code{"Delmoral"}.
	When \code{method} is \code{"Lenormand"}, the number of simulations below the tolerance threshold is equal to \code{nb_simul * alpha}. See the package's vignette and Lenormand et al. (2012) for details.}
  \item{summary_stat_target}{
	a vector containing the targeted (observed) summary statistics.}
  \item{n_cluster}{
    a positive integer. If larger than 1 (the default value), \code{ABC_sequential} will launch \code{model} simulations in parallel on \code{n_cluster} cores of the computer.}
  \item{\dots}{Additional arguments can be passed depending on the choosen method (see below)}
}

\section{Additional paramaters}{
  Depending on the choosen method, you can specify the following arguments:
  \describe{
    \item{use_seed}{
	  logical. If \code{TRUE} (default), \code{ABC_sequential} provides as input to the function \code{model} an array containing an integer seed value and the model parameters used for the simulation.
	  The seed value should be used by \code{model} to initialize its pseudo-random number generators (if \code{model} is stochastic).
	  If \code{FALSE} (default), \code{ABC_sequential} provides as input to the function \code{model} an array containing solely the model parameters used for the simulation.}
    \item{seed_count}{
	  a positive integer, the initial seed value provided to the function \code{model} (if \code{use_seed=TRUE}). This value is incremented by 1 at each call of the function \code{model}.}
    \item{inside_prior}{
	  logical used when \code{method} is \code{"Beaumont"} or \code{"Lenormand"}. \code{TRUE} by default.
	  If \code{FALSE}, parameter sampling is not restricted to the initial ranges of the prior distribution during the subsequent algorithm steps.}
    \item{tolerance_tab}{
	  a vector containing the sequence of tolerance thresholds when \code{method} is \code{"Beaumont"}, or the targeted final tolerance threshold when \code{method} is \code{"Drovandi"}.}
    \item{alpha}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Drovandi"}, \code{"Delmoral"} or \code{"Lenormand"}.
	  \code{alpha} is the proportion of particles rejected at each step in the algorithm \code{"Drovandi"}.
	  This is the proportion of particles kept at each step in the algorithms \code{"Delmoral"} and \code{"Lenormand"}.
	  Default values are 0.5 when \code{method} is \code{"Drovandi"} or \code{"Lenormand"} and 0.9 for \code{"Delmoral"}.
	  See the package's vignette for details.}
    \item{c}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Drovandi"}.
	  This is the expected proportion of particles which are going to be duplicated at each step.
	  Default value is 0.01.
	  See the package's vignette and Drovandi and Pettitt (2011) for details.}
    \item{first_tolerance_level_auto}{
	  logical used when \code{method} is \code{"Drovandi"}.
	  Default value is \code{TRUE}. In this case, the first tolerance threshold is determined by the algorithm, by taking the 1-\code{alpha} quantile of the distances between the simulated and targeted summary statistics.
	  If \code{FALSE}, the initial tolerance threshold for the first step has to be provided as the first element of the vector \code{tolerance_tab}. In this case, the targeted final tolerance threshold is the second element of \code{tolerance_tab}.}
    \item{M}{
	  a positive integer used when \code{method} is \code{"Delmoral"}.
	  This is the number of \code{model} simulations performed for each parameter set.
	  Default value is 1. See the package's vignette and Del Moral et al. (2012) for details.}
    \item{nb_threshold}{
	  a positive integer used when \code{method} is \code{"Delmoral"}.
	  This is the minimal effective sample size below which a resampling step is launched. See the package's vignette and Del Moral et al. (2012) for details.}
    \item{tolerance_target}{
	  a positive number used when \code{method} is \code{"Delmoral"}.
	  This is the targeted final tolerance threshold.}
    \item{p_acc_min}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Lenormand"}.
	  This is the stopping criterion of the algorithm: a small number ensures a better convergence of the algorithm, but at a cost in computing time.
	  Default value is 0.05. See the package's vignette and Lenormand et al. (2012) for details.}
    \item{verbose}{
	  logical. \code{FALSE} by default. If \code{TRUE}, \code{ABC_sequential} writes in the current directory intermediary results at the end of each step of the algorithm.}
    \item{progress_bar}{
      logical, \code{FALSE} by default. If \code{TRUE}, \code{ABC_sequential} will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}
  }
}
\details{
  See the package's vignette for details on the four algorithms.
}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{weights}{
  	The weights of the different \code{model} simulations.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations of the initial step.
	These values are used to normalize the summary statistics before the computation of the Euclidean distance between simulations and data.}
  \item{epsilon}{
  	The final maximal distance between simulations and data in the retained sample of particles.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{computime}{
  	The computing time to perform the simulations.}
}
\references{
	Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009) Adaptive approximate Bayesian computation. \emph{Biometrika},\bold{96}, 983--990.

	Del Moral, P., Doucet, A., and Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}, 1009--1020.

	Drovandi, C. C. and Pettitt, A. N. (2011) Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. \emph{Biometrics}, \bold{67}, 225--233.

	Lenormand, M., Jabot, F., Deffuant G. (2012) Adaptive approximate Bayesian computation for complex models. http://arxiv.org/pdf/1111.1308.pdf
}
\author{Franck Jabot, Thierry Faure and Nicolas Dumoulin}
\seealso{
  \code{\link{binary_model}}, \code{\link{binary_model_cluster}}, \code{\link{ABC_rejection}}, \code{\link{ABC_mcmc}}
}
\examples{
  \dontrun{
    priormatrix=cbind(c(500,3,-2.3,1,-25,-0.7),c(500,5,1.6,1,125,3.2))
    sum_stat_obs=c(100,2.5,20,30000)

    ## artificial example to perform the Beaumont et al. (2009)'s method.
    ##
    ABC_Beaumont<-ABC_sequential(method="Beaumont", model=trait_model, prior_matrix=priormatrix,
    nb_simul=10, summary_stat_target=sum_stat_obs, n_cluster=1, tolerance_tab=c(8,5))
    ABC_Beaumont

    ## artificial example to perform the Drovandi and Pettitt (2011)'s method.
    ##
    ABC_Drovandi<-ABC_sequential(method="Drovandi", model=trait_model, prior_matrix=priormatrix,
    nb_simul=10, summary_stat_target=sum_stat_obs, n_cluster=1, tolerance_tab=3, c=0.7)
    ABC_Drovandi

    ## artificial example to perform the Del Moral et al. (2012)'s method.
    ##
    ABC_Delmoral<-ABC_sequential(method="Delmoral", model=trait_model, prior_matrix=priormatrix,
    nb_simul=10, summary_stat_target=sum_stat_obs, n_cluster=1, alpha=0.5, tolerance_target=3)
    ABC_Delmoral

    ## artificial example to perform the Lenormand et al. (2012)'s method.
    ##
    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=trait_model, prior_matrix=priormatrix,
    nb_simul=10, summary_stat_target=sum_stat_obs, n_cluster=1, p_acc_min=0.4)
    ABC_Lenormand
  }
}
\keyword{abc}
\keyword{model}
\keyword{inference}
\keyword{sequential monte carlo}
\keyword{population monte carlo}
