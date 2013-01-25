## FUNCTION ABC_sequential: Sequential ABC methods (Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
###################################################################################################################################

ABC_sequential <-function(method,model,prior,nb_simul,summary_stat_target,n_cluster=1,use_seed=FALSE,verbose=FALSE,...){
    ## checking errors in the inputs
    if(missing(method)) stop("'method' is missing")
    if(missing(model)) stop("'model' is missing")
    if(missing(prior)) stop("'prior' is missing")
    if(!is.list(prior)) stop("'prior' has to be a list")
    l=length(prior)
    for (i in 1:l){
    	if(!any(prior[[i]][1] == c("unif", "normal", "lognormal", "exponential"))) {
        	stop("Prior distribution type must be unif, normal, lognormal or exponential")
    	}
	if (prior[[i]][1]=="exponential"){
		if (length(prior[[i]])<2){
			stop(paste("Incomplete prior information for parameter ",i,sep=""))
		}
	}
	else{
		if (length(prior[[i]])<3){
			stop(paste("Incomplete prior information for parameter ",i,sep=""))
		}
	}
    }
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(missing(summary_stat_target)) stop("'summary_stat_target' is missing")
    if(!any(method == c("Beaumont", "Drovandi", "Delmoral", "Lenormand"))) {
        stop("Method must be Beaumont, Drovandi, Delmoral or Lenormand")
    }
    if(!is.vector(nb_simul)) stop("'nb_simul' has to be a number.")
    if(length(nb_simul)>1) stop("'nb_simul' has to be a number.")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1.")
    nb_simul=floor(nb_simul)
    if(!is.vector(summary_stat_target)) stop("'summary_stat_target' has to be a vector.")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)
    if(!is.logical(use_seed)) stop("'use_seed' has to be boolean")
    if(!is.logical(verbose)) stop("'verbose' has to be boolean")

    	sequential=NULL
	if (n_cluster==1){
		sequential = .ABC_sequential(method,model,prior,nb_simul,summary_stat_target,use_seed,verbose,...)
	}
	else{
		sequential = .ABC_sequential_cluster(method,model,prior,nb_simul,summary_stat_target,n_cluster,use_seed,verbose,...)
	}
sequential
}



