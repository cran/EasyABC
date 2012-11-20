## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################
ABC_rejection<-function(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0,n_cluster=1,progress_bar=FALSE){
    ## checking errors in the inputs
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) 			  stop("'prior_matrix' has to be a matrix or data.frame")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1")
    if(!is.logical(use_seed)) stop("'use_seed' has to be boolean")
    if(!is.vector(seed_count)) stop("'seed_count' has to be a number")
    if(length(seed_count)>1) stop("'seed_count' has to be a number")
    if (seed_count<0) stop ("'seed_count' has to be a positive number")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)

    nb_simul=floor(nb_simul)
    seed_count=floor(seed_count)
    rejection=NULL
    if (n_cluster==1){
    	rejection=.ABC_rejection(model,prior_matrix,nb_simul,use_seed,seed_count,progress_bar)
    }
    else{
	rejection=.ABC_rejection_cluster(model,prior_matrix,nb_simul,seed_count,n_cluster)
    }
list(param=rejection$param, stats=rejection$stats, weights=rejection$weights, stats_normalization=rejection$stats_normalization, nsim=rejection$nsim, computime=rejection$computime)
}


