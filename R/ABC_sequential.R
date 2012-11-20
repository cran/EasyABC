## FUNCTION ABC_sequential: Sequential ABC methods (Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
###################################################################################################################################

ABC_sequential <-function(method,model,prior_matrix,nb_simul,summary_stat_target,n_cluster=1,...){
    ## checking errors in the inputs
    if(missing(method)) stop("'method' is missing")
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(missing(summary_stat_target)) stop("'summary_stat_target' is missing")
    if(!any(method == c("Beaumont", "Drovandi", "Delmoral", "Lenormand"))) {
        stop("Method must be Beaumont, Drovandi, Delmoral or Lenormand")
    }
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) stop("'prior_matrix' has to be a matrix or data.frame.")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns.")
    if(!is.vector(nb_simul)) stop("'nb_simul' has to be a number.")
    if(length(nb_simul)>1) stop("'nb_simul' has to be a number.")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1.")
    nb_simul=floor(nb_simul)
    if(!is.vector(summary_stat_target)) stop("'summary_stat_target' has to be a vector.")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)

    	sequential=NULL
	if (n_cluster==1){
		sequential = .ABC_sequential(method,model,prior_matrix,nb_simul,summary_stat_target,...)
	}
	else{
		sequential = .ABC_sequential_cluster(method,model,prior_matrix,nb_simul,summary_stat_target,n_cluster,...)
	}
sequential
}



