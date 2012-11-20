## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al. 2009)
##############################################################################
ABC_mcmc <-function(method,model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_cluster=1,...){
    ## checking errors in the inputs
    if(missing(method)) stop("'method' is missing")
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(n_obs)) stop("'n_obs' is missing")
    if(missing(n_between_sampling)) stop("'n_between_sampling' is missing")
    if(missing(summary_stat_target)) stop("'summary_stat_target' is missing")
    if(!any(method == c("Marjoram_original", "Marjoram", "Wegmann"))){
        stop("Method must be Marjoram_original, Marjoram or wegmann")
    }
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) stop("'prior_matrix' has to be a matrix or data.frame.")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns.")
    if(!is.vector(summary_stat_target)) stop("'summary_stat_target' has to be a vector.")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)


    	mcmc=NULL
	if (n_cluster==1){
		mcmc = .ABC_mcmc_internal(method,model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,...)
	}
	else{
		mcmc = .ABC_mcmc_cluster(method,model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_cluster,...)
	}
mcmc
}

