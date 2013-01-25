#######################
## INTERNAL FUNCTIONS
#######################

## necessary libraries: -> imported thanks to the file NAMESPACE
# library(mnormt)
# library(lhs)
# library(pls)
# library(MASS)
# library(parallel)
# library(abc)

## function to compute a distance between a matrix of simulated statistics (row: different simulations, columns: different summary statistics) and the array of data summary statistics
#######################################################################################################################################################################################
.compute_dist<-function(summary_stat_target,simul,sd_simul){
  l=length(summary_stat_target)
  vartab=array(1,l)
  if (l>1){
   nsimul=dim(simul)[1]
   dist=array(0,nsimul)
   for (i in 1:l){
    vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
    dist=dist+vartab[i]*(simul[,i]-summary_stat_target[i])*(simul[,i]-summary_stat_target[i]) ## an euclidean distance is used
   }
  }
  else{
    vartab[1]=min(1,1/(sd_simul[1]*sd_simul[1])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
    dist=vartab[1]*(simul-summary_stat_target[1])*(simul-summary_stat_target[1]) ## an euclidean distance is used
  }
  dist
}

## same as .compute_dist when there is only one simulation
##########################################################
.compute_dist_single<-function(summary_stat_target,simul,sd_simul){
  l=length(summary_stat_target)
  dist=0
  vartab=array(1,l)
  if (l>1){
   for (i in 1:l){
    vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i]))
    dist=dist+vartab[i]*(simul[i]-summary_stat_target[i])*(simul[i]-summary_stat_target[i])
   }
  }
  else{
    vartab[1]=min(1,1/(sd_simul[1]*sd_simul[1])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
    dist=vartab[1]*(simul-summary_stat_target[1])*(simul-summary_stat_target[1]) 
  }
  dist
}

## function to select the simulations that are at a distance smaller than tol from the data
###########################################################################################
.selec_simul<-function(summary_stat_target,param,simul,sd_simul,tol){
  dist=.compute_dist(summary_stat_target,simul,sd_simul)
  dd=dim(param)[1]
  ll=length(dist[dist<tol])
  if (!is.null(dd)){
   if (ll>1){
    param2=param[dist<tol,]
   }
   else{
     if (ll==1){
    	param2=as.matrix(param[dist<tol,])
	dim(param2)=c(dim(param2)[2],dim(param2)[1])
     }
     else{
  	param2=NULL
     }
   }
  }
  else{
   if (ll>=1){
    param2=as.matrix(param[dist<tol])
   }
   else{
    param2=NULL
   }
  }

  dd=dim(simul)[1]
  if (!is.null(dd)){
   if (ll>1){
    simul2=simul[dist<tol,]
   }
   else{
     if (ll==1){
    	simul2=as.matrix(simul[dist<tol,])
	dim(simul2)=c(dim(simul2)[2],dim(simul2)[1])
     }
     else{
  	simul2=NULL
     }
   }
  }
  else{
   if (ll>=1){
    simul2=as.matrix(simul[dist<tol])
   }
   else{
    simul2=NULL
   }
  }
  cbind(param2,simul2)
}

## function to randomly pick a particle from a weighted array (of sum=1)
########################################################################
.particle_pick<-function(param,tab_weight){
  tab_weight2=tab_weight/sum(tab_weight)
  u=runif(1)
  weight_cum=cumsum(tab_weight2)
  pos=1:length(tab_weight)
  p=min(pos[weight_cum>u])
  res=NULL
  if (!is.null(dim(param)[1])){
	  res=param[p,]
  }
  else{
	res=param[p]
  }
res
}


## function to check whether moved parameters are still in the prior distribution
#################################################################################
.is_included<-function(res,prior){
  test=TRUE
  count=1
  for (i in 1:length(prior)){
	if (prior[[i]][1]=="unif"){
		if (as.numeric(prior[[i]][2])!=as.numeric(prior[[i]][3])){
   			if ((res[count]<as.numeric(prior[[i]][2]))||(res[count]>as.numeric(prior[[i]][3]))){
      				test=FALSE
   			}
		}
		else{
			count=count-1
		}
	}
	else{
		if (prior[[i]][1]=="exponential"){
			if (res[count]<0){
      				test=FALSE
   			}
		}
	}
	count=count+1
  }
test	
}

## function to check whether moved parameters are still in the prior distribution
#################################################################################
.is_included2<-function(res,prior){
  test=TRUE
  for (i in 1:length(prior)){
	if (prior[[i]][1]=="unif"){
   		if ((res[i]<as.numeric(prior[[i]][2]))||(res[i]>as.numeric(prior[[i]][3]))){
      				test=FALSE
  		}
	}
	else{
		if (prior[[i]][1]=="exponential"){
			if (res[i]<0){
      				test=FALSE
   			}
		}
	}
  }
test	
}


## function to move a particle
##############################
.move_particle<-function(param_picked,varcov_matrix){
  # with package mnormt
  #library(mnormt) 
  rmnorm(n = 1, mean = param_picked, varcov_matrix)
}

## function to move a particle
##############################
.move_particleb<-function(param_picked,varcov_matrix,prior){
  # with package mnormt
  #library(mnormt)
  test=FALSE
  counter=0
  while ((!test)&&(counter<100)){
    counter=counter+1
    res=rmnorm(n = 1, mean = param_picked, varcov_matrix)
    test=.is_included(res,prior)
  }
  if (counter==100){
    stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
  }
  res
}


.sample_prior_unif<-function(prior){
	runif(1,min=as.numeric(prior[2]),max=as.numeric(prior[3]))
}

.sample_prior_normal<-function(prior){
	rnorm(1,mean=as.numeric(prior[2]),sd=as.numeric(prior[3]))
}

.sample_prior_lognormal<-function(prior){
	rlnorm(1,meanlog=as.numeric(prior[2]),sdlog=as.numeric(prior[3]))
}

.sample_prior_exponential<-function(prior){
	rexp(1,rate=as.numeric(prior[2]))
}

.sample_prior<-function(prior){
	l=length(prior)
	param=NULL
	for (i in 1:l){
		param[i]=(switch(EXPR = prior[[i]][1],
		    "unif" = .sample_prior_unif(prior[[i]]),
		    "normal" = .sample_prior_normal(prior[[i]]),
		    "lognormal" = .sample_prior_lognormal(prior[[i]]),
		    "exponential" = .sample_prior_exponential(prior[[i]])))
	}
param
}


.ABC_rejection_internal<-function(model, prior, nb_simul, use_seed, seed_count, verbose=FALSE, progressbarwidth=0){
  options(scipen=50)
  tab_simul_summarystat=NULL
  tab_param=NULL
  pb = NULL
  if (progressbarwidth>0) {
    pb = .progressBar(width=progressbarwidth)
  }
  if (verbose){
     write.table(NULL,file="output",row.names=F,col.names=F,quote=F)
  }
  l=length(prior)
  start = Sys.time()
  for (i in 1:nb_simul){
    param=.sample_prior(prior)
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summarystat=model(param)
    tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(simul_summarystat))
    if (use_seed) {
      tab_param=rbind(tab_param,param[2:(l+1)])
      if (verbose){
	write.table(c(param[2:(l+1)],simul_summarystat),file="output",row.names=F,col.names=F,quote=F,append=T)
      }
    } else {
      tab_param=rbind(tab_param,param)
      if (verbose){
	write.table(c(param,simul_summarystat),file="output",row.names=F,col.names=F,quote=F,append=T)
      }
    }
    if (!is.null(pb)) {
      duration = difftime(Sys.time(), start, units="secs")
      text = "";
      if (i == nb_simul) {
        text = paste("Completed  in", format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
      } 
      else {
        text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/i*(nb_simul-i), tz="GMT"), "%H:%M:%S"));
      }
      .updateProgressBar(pb, i/nb_simul, text)
    }
  }
  if (!is.null(pb)) {
    close(pb)
  }
  options(scipen=0)
  
  list(param=as.matrix(tab_param), summarystat=as.matrix(tab_simul_summarystat), start=start)
}

## function to compute particle weights
#######################################  
.compute_weight_prior<-function(particle,prior){
    res=1
    count=1
    for (i in 1:length(prior)){
	if (prior[[i]][1]=="unif"){
		if (as.numeric(prior[[i]][2])!=as.numeric(prior[[i]][3])){
			res=res/(as.numeric(prior[[i]][3])-as.numeric(prior[[i]][2]))			
		}
		else{
			count=count-1
		}
	}
	else{
		if (prior[[i]][1]=="normal"){
			res=res*dnorm(particle[count],as.numeric(prior[[i]][2]),as.numeric(prior[[i]][3]))
		}
		else{
			if (prior[[i]][1]=="lognormal"){
				res=res*dlnorm(particle[count],as.numeric(prior[[i]][2]),as.numeric(prior[[i]][3]))
			}
			else{
				if (prior[[i]][1]=="exponential"){
					res=res*dexp(particle[count],as.numeric(prior[[i]][2]))
				}
			}
		}
	}
	count=count+1
    }
res
}


.compute_weight_prior_tab<-function(particle,prior){
	l=dim(particle)[1]
	res=array(0,l)
	for (i in 1:l){
		res[i]=.compute_weight_prior(particle[i,],prior)
	}
res
}

.compute_weight<-function(param_simulated,param_previous_step,tab_weight,prior){
  vmat=2*var(param_previous_step)
  n_particle=dim(param_previous_step)[1]
  n_new_particle=dim(param_simulated)[1]
  l=dim(param_previous_step)[2]
  if (is.null(n_particle)){
	n_particle=length(param_previous_step)
	n_new_particle=length(param_simulated)
	l=1
  }
  tab_weight_new=array(0,n_new_particle)
  invmat=0.5*solve(vmat)
  for (i in 1:n_particle){
    for (k in 1:n_new_particle){
      if (l>1){
       temp=as.numeric(param_simulated[k,])-param_previous_step[i,]
       tab_weight_new[k]=tab_weight_new[k]+tab_weight[i]*as.numeric(exp(- t(temp) %*% invmat %*% temp ))
      }
      else{
       temp=as.numeric(param_simulated[k])-param_previous_step[i]
       tab_weight_new[k]=tab_weight_new[k]+tab_weight[i]*as.numeric(exp(- t(temp) %*% invmat %*% temp ))
      }
    }
  }
  tab_weight_prior=.compute_weight_prior_tab(param_simulated,prior)
  tab_weight_new=tab_weight_prior/tab_weight_new
  tab_weight_new/sum(tab_weight_new)
}


## function to move a particle with a unidimensional normal jump
################################################################
.move_particle_uni<-function(param_picked,sd_array){
  res=param_picked
  for (i in 1:length(param_picked)){
    res[i]=rnorm(n = 1, mean = param_picked[i], sd_array[i])
  }
  res
}

## function to move a particle with a unidimensional normal jump
################################################################
.move_particleb_uni<-function(param_picked,sd_array,prior){
  test=FALSE
  res=param_picked
  counter=0
  while ((!test)&&(counter<100)){
    counter=counter+1
    for (i in 1:length(param_picked)){
      res[i]=rnorm(n = 1, mean = param_picked[i], sd_array[i])
    }
    test=.is_included(res,prior)
  }
  if (counter==100){
    stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
  }
  res
}


## function to compute particle weights with unidimensional jumps
#################################################################
.compute_weight_uni<-function(param_simulated,param_previous_step,tab_weight,prior){
  l=dim(param_previous_step)[2]
  if (!is.null(l)){
   n_particle=dim(param_previous_step)[1]
   n_new_particle=dim(param_simulated)[1]
   var_array=array(1,l)
   multi=(1/sqrt(2*pi))^l
   for (j in 1:l){
    var_array[j]=4*var(param_previous_step[,j])
    multi=multi*(1/sqrt(var_array[j]/2))
   }
  }
  else{
   l=1
   n_particle=length(param_previous_step)
   n_new_particle=length(param_simulated)
   multi=(1/sqrt(2*pi))
   var_array=4*var(param_previous_step)
   multi=multi*(1/sqrt(var_array/2))
  }
  var_array=as.numeric(var_array)
  tab_weight_new=array(0,n_new_particle)
  for (i in 1:n_particle){
    tab_temp=array(tab_weight[i]*multi,n_new_particle)
    if (l>1){
     for (k in 1:l){
      tab_temp=tab_temp*exp(-(as.numeric(param_simulated[,k])-as.numeric(param_previous_step[i,k]))*(as.numeric(param_simulated[,k])-as.numeric(param_previous_step[i,k]))/var_array[k])
     }
    }
    else{
     tab_temp=tab_temp*exp(-(as.numeric(param_simulated)-as.numeric(param_previous_step[i]))*(as.numeric(param_simulated)-as.numeric(param_previous_step[i]))/var_array)
    }
    tab_weight_new=tab_weight_new+tab_temp
  }
  tab_weight_prior=.compute_weight_prior_tab(param_simulated,prior)
  tab_weight_new=tab_weight_prior/tab_weight_new
  tab_weight_new/sum(tab_weight_new)
}


## function to perform ABC simulations from a non-uniform prior and with unidimensional jumps
#############################################################################################
.ABC_launcher_not_uniform_uni<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count,inside_prior){
  tab_simul_summarystat=NULL
  tab_param=NULL
  l=dim(param_previous_step)[2]
  lp=1
  if (is.null(l)){
	lp=0
	l=length(param_previous_step)
	covmat=2*cov.wt(as.matrix(param_previous_step),as.vector(tab_weight))$cov # computation of a WEIGHTED variance
  }
  else{
  	covmat=2*cov.wt(as.matrix(param_previous_step[,tab_unfixed_param]),as.vector(tab_weight))$cov # computation of a WEIGHTED variance
  }
  l_array=dim(param_previous_step[,tab_unfixed_param])[2]
  if (is.null(l_array)){
	l_array=1
  }
  sd_array=array(1,l_array)
  for (j in 1:l_array){
    sd_array[j]=sqrt(covmat[j,j])
    
  }
  for (i in 1:nb_simul){
    if (!inside_prior){
      # pick a particle
      if (lp==0){
	      param_picked=.particle_pick(as.matrix(param_previous_step),tab_weight)
      }
      else{
	      param_picked=.particle_pick(as.matrix(param_previous_step[,tab_unfixed_param]),tab_weight)
      }
      # move it
      if (lp==0){
     	param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
      }
      else{
     	param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
      }
    }
    else{
      test=FALSE
      counter=0
      while ((!test)&&(counter<100)){
        counter=counter+1
        # pick a particle
	if (lp==0){
	        param_picked=.particle_pick(param_previous_step,tab_weight)
        }
	else{
	        param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
        }
        # move it
	if (lp==0){
	        param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
         	test=.is_included(param_moved,prior)
 	}
	else{
	        param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
        	test=.is_included(param_moved,prior)
	}
      }
      if (counter==100){
        stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
      }
    }
    if (lp==0){
	    param=param_previous_step[1]
    }
    else{
	    param=param_previous_step[1,]
    }
    param[tab_unfixed_param]=param_moved
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summarystat=model(param)
    tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
    if (use_seed) {
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    else{
      tab_param=rbind(tab_param,param)
    }
  }
  cbind(tab_param,tab_simul_summarystat)
}


## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################
.ABC_rejection<-function(model,prior,nb_simul,use_seed,seed_count,verbose,progress_bar){
    nb_simul=floor(nb_simul)
    seed_count=floor(seed_count)
    pgwidth=0
    if (progress_bar){
	pgwidth=50
    }

    rejection = .ABC_rejection_internal(model,prior,nb_simul,use_seed,seed_count,verbose,progressbarwidth=pgwidth)

    sd_simul=sapply(as.data.frame(rejection$summarystat), sd)
    
list(param=rejection$param, stats=as.matrix(rejection$summarystat), weights=array(1/nb_simul,nb_simul), stats_normalization=as.numeric(sd_simul), nsim=nb_simul, computime=as.numeric(difftime(Sys.time(), rejection$start, units="secs")))
}


## PMC ABC algorithm: Beaumont et al. Biometrika 2009
#####################################################
.ABC_PMC<-function(model,prior,nb_simul,summary_stat_target,use_seed,verbose,seed_count=0,inside_prior=TRUE,tolerance_tab=-1,progress_bar=FALSE){
  ## checking errors in the inputs
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(inside_prior)) stop("'inside_prior' has to be boolean.")
  if(!is.vector(tolerance_tab)) stop("'tolerance_tab' has to be a vector.")
  if(tolerance_tab[1]==-1) stop("'tolerance_tab' is missing")
  if (min(tolerance_tab)<=0) stop ("tolerance values have to be strictly positive.")
  lll=length(tolerance_tab)
  if (lll<=1) stop("at least two tolerance values need to be provided.")
  if (min(tolerance_tab[1:(lll-1)]-tolerance_tab[2:lll])<=0) stop ("tolerance values have to decrease.")
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")
  
  if (progress_bar){
	  print("    ------ Beaumont et al. (2009)'s algorithm ------")
  }
 

  start = Sys.time()
  seed_count_ini=seed_count
  T=length(tolerance_tab)
  nparam=length(prior)
  if (is.null(nparam)){
	nparam=1
  }
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }

  
  ## step 1
  nb_simul_step=nb_simul
  simul_below_tol=NULL
  while (nb_simul_step>0){
    if (nb_simul_step>1){
      # classic ABC step
      tab_ini=.ABC_rejection_internal(model,prior,nb_simul_step,use_seed,seed_count)
      if (nb_simul_step==nb_simul){     
        sd_simul=sapply(as.data.frame(tab_ini$summarystat),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
      }
      seed_count=seed_count+nb_simul_step
      # selection of simulations below the first tolerance level
      simul_below_tol=rbind(simul_below_tol,.selec_simul(summary_stat_target,tab_ini$param,tab_ini$summarystat,sd_simul,tolerance_tab[1]))
      if (length(simul_below_tol)>0){
        nb_simul_step=nb_simul-dim(simul_below_tol)[1]
      }
    }
    else{
      tab_ini=.ABC_rejection_internal(model,prior,nb_simul_step,use_seed,seed_count)
      seed_count=seed_count+nb_simul_step
      if (.compute_dist_single(summary_stat_target,tab_ini$summarystat,sd_simul)<tolerance_tab[1]){
        simul_below_tol=rbind(simul_below_tol,cbind(tab_ini$param, tab_ini$summarystat))
        nb_simul_step=0
      }
    }
  } # until we get nb_simul simulations below the first tolerance threshold
  # initially, weights are equal
  tab_weight=array(1/nb_simul,nb_simul)
  if (verbose==TRUE){
    write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  ## steps 2 to T
  for (it in 2:T){
    nb_simul_step=nb_simul
    simul_below_tol2=NULL
    while (nb_simul_step>0){
      if (nb_simul_step>1){
        # Sampling of parameters around the previous particles
	tab_ini=.ABC_launcher_not_uniform_uni(model,prior,as.matrix(simul_below_tol[,1:nparam]),tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
        seed_count=seed_count+nb_simul_step
        simul_below_tol2=rbind(simul_below_tol2,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[it]))
        if (length(simul_below_tol2)>0){
          nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
        }
      }
      else{
        tab_ini=.ABC_launcher_not_uniform_uni(model,prior,as.matrix(simul_below_tol[,1:nparam]),tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
        seed_count=seed_count+nb_simul_step
        if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[it]){
          simul_below_tol2=rbind(simul_below_tol2,tab_ini)
          nb_simul_step=0
        }
      }
    } # until we get nb_simul simulations below the it-th tolerance threshold
    # update of particle weights
    tab_weight2=.compute_weight_uni(as.matrix(as.matrix(simul_below_tol2[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(simul_below_tol[,1:nparam])[,tab_unfixed_param]),tab_weight,prior)
    # update of the set of particles and of the associated weights for the next ABC sequence
    tab_weight=tab_weight2
    simul_below_tol=matrix(0,nb_simul,(nparam+nstat))
    for (i1 in 1:nb_simul){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    if (verbose==TRUE){
      write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
	    print(paste("step ",it," completed",sep=""))
    }
  }
  list(param=as.matrix(simul_below_tol[,1:nparam]),stats=as.matrix(simul_below_tol[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(simul_below_tol[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}

## function to select the alpha quantile closest simulations
############################################################
.selec_simul_alpha<-function(summary_stat_target,param,simul,sd_simul,alpha){
  dist=.compute_dist(summary_stat_target,simul,sd_simul)
  n_alpha=ceiling(alpha*length(dist))
  tol=sort(dist)[n_alpha]
  ll=length(dist[dist<tol])
  if (ll==0){
	param2=NULL
	simul2=NULL
  }
  else{
    if (!is.null(dim(param)[1])){
     if (ll>1){
      param2=param[dist<=tol,]
     }
     else{
      param2=as.matrix(param[dist<=tol,])
      dim(param2)=c(dim(param2)[2],dim(param2)[1])
     }
    }
    else{
     param2=as.matrix(param[dist<=tol])
    }
    if (!is.null(dim(simul)[1])){
     if (ll>1){
      simul2=simul[dist<=tol,]
     }
     else{
      simul2=as.matrix(simul[dist<=tol,])
      dim(simul2)=c(dim(simul2)[2],dim(simul2)[1])
     }
    }
    else{
     simul2=as.matrix(simul[dist<=tol])
    }
  }
  cbind(param2,simul2)
}

## function to select the simulations that are at a distance smaller are equal to tol from the data
###################################################################################################
.selec_simulb<-function(summary_stat_target,param,simul,sd_simul,tol){
  dist=.compute_dist(summary_stat_target,simul,sd_simul)
  ll=length(dist[dist<=tol])
  if (ll==0){
	param2=NULL
	simul2=NULL
  }
  else{
   if (!is.null(dim(param)[1])){
    if (ll>1){
     param2=param[dist<=tol,]
    }
    else{
     param2=as.matrix(param[dist<=tol,])
     dim(param2)=c(dim(param2)[2],dim(param2)[1])
    }
   }
   else{
    param2=as.matrix(param[dist<=tol])
   }
   if (!is.null(dim(simul)[1])){
    if (ll>1){
     simul2=simul[dist<=tol,]
    }
    else{
     simul2=as.matrix(simul[dist<=tol,])
     dim(simul2)=c(dim(simul2)[2],dim(simul2)[1])
    }
   }
   else{
    simul2=as.matrix(simul[dist<=tol])
   }
  }
  cbind(param2,simul2)
}

## sequential algorithm of Drovandi & Pettitt 2011 - the proposal used is a multivariate normal (cf paragraph 2.2 - p. 227 in Drovandi & Pettitt 2011)
######################################################################################################################################################
.ABC_Drovandi<-function(model,prior,nb_simul,summary_stat_target,use_seed,verbose,tolerance_tab=-1,alpha=0.5,c=0.01,first_tolerance_level_auto=TRUE,seed_count=0,progress_bar=FALSE){
  ## checking errors in the inputs
  if(!is.vector(tolerance_tab)) stop("'tolerance_tab' has to be a vector.")
  if(tolerance_tab[1]==-1) stop("'tolerance_tab' is missing")
  if (min(tolerance_tab)<=0) stop ("'tolerance values have to be strictly positive.")
  lll=length(tolerance_tab)
  if (lll>1){
    if (min(tolerance_tab[1:(lll-1)]-tolerance_tab[2:lll])<=0) stop ("'tolerance values have to decrease.")
  }
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.")
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.")
  if(!is.vector(c)) stop("'c' has to be a vector.")
  if(length(c)>1) stop("'c' has to be a number.")
  if (c<=0) stop ("'c' has to be between 0 and 1.")
  if (c>=1) stop ("'c' has to be between 0 and 1.")
  if(!is.logical(first_tolerance_level_auto)) stop("'first_tolerance_level_auto' has to be boolean.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")
  
  start = Sys.time()
  
  progressbarwidth=0
  if (progress_bar){
  	print("    ------ Drovandi & Pettitt (2011)'s algorithm ------")
	progressbarwidth=50
  } 
  
  seed_count_ini=seed_count
  n_alpha=ceiling(nb_simul*alpha)
  nparam=length(prior)
  if (is.null(nparam)){
	nparam=1
  }
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  if (first_tolerance_level_auto){
    tol_end=tolerance_tab[1]
  }
  else{
    tol_end=tolerance_tab[2]
  }
  
  if (progress_bar){
	print(paste("targetted tolerance = ",tol_end,sep=""))
  } 
  
  ## step 1
  nb_simul_step=ceiling(nb_simul/(1-alpha))
  simul_below_tol=NULL
  if (first_tolerance_level_auto){
    # classic ABC step
    tab_ini=.ABC_rejection_internal(model,prior,nb_simul_step,use_seed,seed_count,progressbarwidth)
    sd_simul=sapply(as.data.frame(tab_ini$summarystat),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
    seed_count=seed_count+nb_simul_step
    # selection of simulations below the first tolerance level
    simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,tab_ini$param,tab_ini$summarystat,sd_simul,(1-alpha)))
    simul_below_tol=as.matrix(simul_below_tol[1:nb_simul,]) # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile
  }
  else{
    nb_simul_step=nb_simul
    while (nb_simul_step>0){
      if (nb_simul_step>1){
        # classic ABC step
        tab_ini=.ABC_rejection_internal(model,prior,nb_simul_step,use_seed,seed_count)
        if (nb_simul_step==nb_simul){
          sd_simul=sapply(as.data.frame(tab_ini$summarystat),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
        }
        seed_count=seed_count+nb_simul_step
        # selection of simulations below the first tolerance level
        simul_below_tol=rbind(simul_below_tol,.selec_simulb(summary_stat_target,tab_ini$param,tab_ini$summarystat,sd_simul,tolerance_tab[1])) # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        if (length(simul_below_tol)>0){
          nb_simul_step=nb_simul-dim(simul_below_tol)[1]
        }
      }
      else{
        tab_ini=.ABC_rejection_internal(model,prior,nb_simul_step,use_seed,seed_count)
        seed_count=seed_count+nb_simul_step
        if (.compute_dist_single(summary_stat_target,tab_ini$summarystat,sd_simul)<=tolerance_tab[1]){
          simul_below_tol=rbind(simul_below_tol,cbind(tab_ini$param, tab_ini$summarystat))
          nb_simul_step=0
        }
      }
    } # until we get nb_simul simulations below the first tolerance threshold
  }
  # initially, weights are equal
  tab_weight=array(1/nb_simul,nb_simul)
  if (verbose==TRUE){
    write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  
  ## following steps until tol_end is reached
  tol_next=tolerance_tab[1]
  if (first_tolerance_level_auto){
    tol_next=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul))
  } 
  R=1
  l=dim(simul_below_tol)[2]
  it=1
  while (tol_next>tol_end){
    it=it+1
    i_acc=0
    nb_simul_step=n_alpha
    # compute epsilon_next
    tol_next=sort(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul))[(nb_simul-n_alpha)]
    # drop the n_alpha poorest particles
    simul_below_tol2=.selec_simulb(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul,tol_next) # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
    simul_below_tol=matrix(0,(nb_simul-n_alpha),(nparam+nstat))
    for (i1 in 1:(nb_simul-n_alpha)){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    
    simul_below_tol2=NULL
    startb = Sys.time()
    # progress bar
    pb=NULL
    if (progress_bar){
	    pb <- .progressBar(width=50)
    }
    duration = 0;
    for (i in 1:nb_simul_step){
      # pick a particle
      simul_picked=.particle_pick(simul_below_tol,tab_weight[1:(nb_simul-n_alpha)])
      for (j in 1:R){
        # move it
        param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(simul_below_tol[,1:nparam])[,tab_unfixed_param])),prior)
        param=simul_picked[1:nparam]
        param[tab_unfixed_param]=param_moved
        if (use_seed) {
          param=c((seed_count+j),param)
        }
        # perform a simulation
        new_simul=c(param,model(param))
        if (use_seed) {
          new_simul=new_simul[2:(l+1)]
        }
        # check whether it is below tol_next and undo the move if it is not
        if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
          simul_picked=as.numeric(new_simul)
          i_acc=i_acc+1
        }
      }
      seed_count=seed_count+R
      simul_below_tol2=rbind(simul_below_tol2,simul_picked)
      # for progressbar message and time evaluation
	if (progress_bar){
      	duration = difftime(Sys.time(), startb, units="secs")
      	text="";
      	if (i==nb_simul_step) {
        		text = paste("Step ",it," completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
      	} 
	      else {
      	  	text = paste("Time elapsed during step ",it,":",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining for step ",it,":",format(.POSIXct(duration/i*(nb_simul_step-i), tz="GMT"), "%H:%M:%S"));
      	}
      	.updateProgressBar(pb, i/nb_simul_step, text)
	}
    }
    if (progress_bar){
	    close(pb)
    }
    simul_below_tol3=matrix(0,nb_simul,(nparam+nstat))
    for (i1 in 1:(nb_simul-n_alpha)){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol3[i1,i2]=as.numeric(simul_below_tol[i1,i2])
      }
    }
    for (i1 in (nb_simul-n_alpha+1):nb_simul){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol3[i1,i2]=as.numeric(simul_below_tol2[(i1-nb_simul+n_alpha),i2])
      }
    }
    simul_below_tol=simul_below_tol3
    p_acc=max(1,i_acc)/(nb_simul_step*R) # to have a strictly positive p_acc
    Rp=R
    R=ceiling(log(c)/log(1-p_acc))
    if (verbose==TRUE){
      write.table(as.numeric(Rp),file=paste("R_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(tol_next),file=paste("tolerance_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    tol_next=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul))
    if (progress_bar){
	    print(paste("step ",it," completed - R used = ",Rp," - tol = ",tol_next," - next R used will be ",R,sep=""))
    }
  }
  
  ## final step to diversify the n_alpha particles
  simul_below_tol2=NULL
  for (i in 1:nb_simul){
    simul_picked=simul_below_tol[i,]
    for (j in 1:R){
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(simul_below_tol[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      if (use_seed) {
        param=c((seed_count+j),param)
      }
      # perform a simulation
      new_simul=c(param,model(param))
      if (use_seed) {
        new_simul=new_simul[2:(l+1)]
      }
      # check whether it is below tol_next and undo the move if it is not
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        simul_picked=as.numeric(new_simul)
      }
    }
    seed_count=seed_count+R
    simul_below_tol2=rbind(simul_below_tol2,as.numeric(simul_picked))
  }
  list(param=as.matrix(simul_below_tol2[,1:nparam]),stats=as.matrix(simul_below_tol2[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol2)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}


## rejection algorithm with M simulations per parameter set
############################################################
.ABC_rejection_M<-function(model,prior,nb_simul,M,use_seed,seed_count){
  tab_simul_summarystat=NULL
  tab_param=NULL
  l=length(prior)
  if (is.null(l)){
	l=1
  }
  for (i in 1:nb_simul){
    param=.sample_prior(prior)
    for (k in 1:M){
      if (use_seed) {
        param=c((seed_count+1),param)	
      }
      seed_count=seed_count+1
      simul_summarystat=model(param)
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(simul_summarystat))
      if (use_seed) {
        param=param[2:(l+1)]
      }
      tab_param=rbind(tab_param,as.numeric(param))
    }
  }
  cbind(tab_param,tab_simul_summarystat)
}

## function to compute a distance between a matrix of simulated statistics and the array of data summary statistics - for M replicates simulations
##################################################################################################################################################
.compute_dist_M<-function(M,summary_stat_target,simul,sd_simul){
  l=length(summary_stat_target)
  nsimul=dim(simul)[1]
  if (is.null(nsimul)){
	nsimul=length(simul)
  }
  vartab=array(1,l)
  dist=array(0,nsimul)
  if (l>1){
   for (i in 1:l){
    vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
    dist=dist+vartab[i]*(simul[,i]-summary_stat_target[i])*(simul[,i]-summary_stat_target[i]) ## an euclidean distance is used
   }
  }
  else{
    vartab=min(1,1/(sd_simul*sd_simul)) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
    dist=dist+vartab*(simul-summary_stat_target[i])*(simul-summary_stat_target[i]) ## an euclidean distance is used
  }
  distb=matrix(0,nsimul/M,M)
  for (i in 1:(nsimul/M)){
    for (j in 1:M){
      distb[i,j]=dist[((i-1)*M+j)]
    }
  }
  distb
}

## function to compute the updated weights, given a new tolerance value
#######################################################################
.compute_weight_delmoral<-function(particle_dist_mat,tolerance){
  n_particle=dim(particle_dist_mat)[1]
  new_weight=array(0,n_particle)
  for (i in 1:n_particle){
    new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
  }
  new_weight/sum(new_weight)
}

## function to compute the updated ESS, given a new tolerance value
###################################################################
.compute_ESS<-function(particle_dist_mat,tolerance){
  n_particle=dim(particle_dist_mat)[1]
  new_weight=array(0,n_particle)
  for (i in 1:n_particle){
    new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
  }
  new_weight=new_weight/sum(new_weight)
  1/(sum(new_weight*new_weight))
}

## function to compute the number of simul below a new tolerance value
######################################################################
.compute_below<-function(particle_dist_mat,tolerance){
  n_particle=dim(particle_dist_mat)[1]
  new_weight=array(0,n_particle)
  for (i in 1:n_particle){
    new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
  }
  new_weight
}

## function to randomly pick a particle from a weighted array (of sum=1) for the Del Moral algorithm
####################################################################################################
.particle_pick_delmoral<-function(simul_below_tol,tab_weight,M){
  u=runif(1)
  tab_weight2=tab_weight/sum(tab_weight)
  weight_cum=cumsum(tab_weight2)
  pos=1:length(tab_weight)
  p=min(pos[weight_cum>u])
  simul_below_tol[((1:M)+(p-1)*M),]
}

## function to replicate each cell of tab_weight M times
########################################################
.replicate_tab<-function(tab_weight,M){
  l=length(tab_weight)
  tab_weight2=array(0,M*l)
  for (i in 1:l){
    tab_weight2[((i-1)*M+(1:M))]=tab_weight[i]
  }
  tab_weight2
}


## sequential algorithm of Del Moral et al. 2012 - the proposal used is a normal in each dimension (cf paragraph 3.2 in Del Moral et al. 2012)
##############################################################################################################################################
.ABC_Delmoral<-function(model,prior,nb_simul,summary_stat_target,use_seed,verbose,alpha=0.9,M=1,nb_threshold=floor(nb_simul/2),tolerance_target=-1,seed_count=0,progress_bar=FALSE){
  
  ## checking errors in the inputs
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.")
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.")
  if(!is.vector(M)) stop("'M' has to be a number.")
  if(length(M)>1) stop("'M' has to be a number.")
  if (M<1) stop ("'M' has to be a positive integer.")
  M=floor(M)
  if(!is.vector(nb_threshold)) stop("'nb_threshold' has to be a number.")
  if(length(nb_threshold)>1) stop("'nb_threshold' has to be a number.")
  if (nb_threshold<1) stop ("'nb_threshold' has to be a positive integer.")
  nb_threshold=floor(nb_threshold)
  if(!is.vector(tolerance_target)) stop("'tolerance_target' has to be a number.")
  if(length(tolerance_target)>1) stop("'tolerance_target' has to be a number.")
  if (tolerance_target<=0) stop ("'tolerance_target' has to be positive.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")
  
  start = Sys.time()

  if (progress_bar){
	  print("    ------ Delmoral et al. (2012)'s algorithm ------") 
  }
  
  seed_count_ini=seed_count
  nparam=length(prior)
  if (is.null(nparam)){
	nparam=1
  }
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  # step 1
  # classic ABC step
  simul_below_tol=.ABC_rejection_M(model,prior,nb_simul,M,use_seed,seed_count)
  seed_count=seed_count+M*nb_simul
  tab_weight=rep(1/nb_simul,nb_simul)
  ESS=nb_simul
  uu=(1:nb_simul)*M  # to compute sd_simul with only one simulation per parameter set
  sd_simul=sapply(as.data.frame(simul_below_tol[uu,(nparam+1):(nparam+nstat)]),sd)  # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
  l=dim(simul_below_tol)[2]
  if (M>1){
    particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  }
  else{
    particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  }
  dim(particle_dist_mat)<-c(nb_simul,M)
  new_tolerance=max(particle_dist_mat)
  
  tab_weight2=.replicate_tab(tab_weight,M)
  if (verbose==TRUE){
    write.table(as.matrix(cbind(tab_weight2,simul_below_tol)),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  
  # following steps
  kstep=1
  while(new_tolerance>tolerance_target){	
    kstep=kstep+1
    # determination of the new tolerance
    ESS_target=alpha*ESS
    
    tolerance_list=sort(as.numeric(names(table(particle_dist_mat))),decreasing=TRUE)
    i=1
    test=FALSE
    while((!test)&&(i<length(tolerance_list))){
      i=i+1
      # computation of new ESS with the new tolerance value
      new_ESS=.compute_ESS(particle_dist_mat,tolerance_list[i])
      # check whether this value is below ESS_targ
      if (new_ESS<ESS_target){
        new_tolerance=tolerance_list[(i-1)]
        test=TRUE
      }
    }
    
    # if effective sample size is too small, resampling of particles
    ESS=.compute_ESS(particle_dist_mat,new_tolerance)
    tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
    tab_below=.compute_below(particle_dist_mat,new_tolerance)
    particles=matrix(0,(nb_simul*M),(nparam+nstat))
    if (ESS<nb_threshold){
      # sample nb_simul particles 
      for (i in 1:nb_simul){
        particles[((1:M)+(i-1)*M),]=as.matrix(.particle_pick_delmoral(simul_below_tol,tab_weight,M))
      }
      simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))
      for (i1 in 1:(nb_simul*M)){
        for (i2 in 1:(nparam+nstat)){
          simul_below_tol[i1,i2]=as.numeric(particles[i1,i2])
        }
      }
      particles=as.matrix(particles[uu,1:nparam])
      if (M>1){
        particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
      }
      else{
        particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
      }
      dim(particle_dist_mat)<-c(nb_simul,M)
      tab_below=.compute_below(particle_dist_mat,new_tolerance)
      # reset their weight to 1/nb_simul
      tab_weight=rep(1/nb_simul,nb_simul)
      ESS=nb_simul
    }
    else{
      particles=as.matrix(simul_below_tol[,1:nparam])
    }
    
    # MCMC move
    covmat=2*cov.wt(as.matrix(as.matrix(particles[uu,tab_unfixed_param])[tab_weight>0,]),as.vector(tab_weight[tab_weight>0]))$cov
    l_array=dim(particles[,tab_unfixed_param])[2]
    if (is.null(l_array)){
	l_array=1
    }
    sd_array=array(1,l_array)
    for (j in 1:l_array){
      sd_array[j]=sqrt(covmat[j,j])
      
    }
    simul_below_tol2=simul_below_tol
    simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))
    
    startb = Sys.time()
    # progress bar
    if (progress_bar){
	    pb <- .progressBar(width=50)
    }
    duration = 0;
    for (i in 1:nb_simul){
      if (tab_weight[i]>0){
        tab_new_simul=NULL
        # move it
        param_moved=.move_particleb_uni(as.numeric(particles[(i*M),tab_unfixed_param]),sd_array,prior)
        param=particles[(i*M),]
        param[tab_unfixed_param]=param_moved
        if (use_seed) {
          param=c((seed_count+1),param)
        }
        # perform M simulations
        for (j in 1:M){
          new_simul=c(param,model(param))
	  if (use_seed) {
	    param[1]=param[1]+1
	  }
          seed_count=seed_count+1
          if (use_seed) {
            new_simul=new_simul[2:(l+1)]
          }
          tab_new_simul=rbind(tab_new_simul,new_simul)
        }
        if (M>1){
          tab_new_simul2=matrix(0,M,(nparam+nstat))
          for (i1 in 1:M){
            for (i2 in 1:(nparam+nstat)){
              tab_new_simul2[i1,i2]=as.numeric(tab_new_simul[i1,i2])
            }
          }
        }
        else{
          tab_new_simul2=as.numeric(tab_new_simul)
        }
        dim(tab_new_simul2)<-c(M,(nparam+nstat))
        # check whether the move is accepted
        n_acc=1
        if (M>1){
          new_dist=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(tab_new_simul2)[,(nparam+1):(nparam+nstat)]),sd_simul)
          n_acc=length(new_dist[new_dist<new_tolerance])
        }
        else{
          new_dist=.compute_dist(summary_stat_target,rbind(tab_new_simul2[(nparam+1):(nparam+nstat)],tab_new_simul2[(nparam+1):(nparam+nstat)]),sd_simul)
          if (new_dist[1]>new_tolerance){
            n_acc=0
          }
        }
        MH=min(1,(n_acc/tab_below[i]))
        uuu=runif(1)
        if (uuu<=MH){
          for (i1 in 1:M){
            for (i2 in 1:(nparam+nstat)){
              simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(tab_new_simul2[i1,i2])
            }
          }
          
        }
        else{
          for (i1 in 1:M){
            for (i2 in 1:(nparam+nstat)){
              simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
            }
          }
        }
      }
      else{
        for (i1 in 1:M){
          for (i2 in 1:(nparam+nstat)){
            simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
          }
        }
      }
    # for progressbar message and time evaluation
    if (progress_bar){
      duration = difftime(Sys.time(), startb, units="secs")
      text="";
      if (i==nb_simul) {
	text = paste("Step ",kstep," completed in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
      } else {
	text = paste("Time elapsed during step ",kstep,":",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining for step ",kstep,":",format(.POSIXct(duration/i*(nb_simul-i), tz="GMT"), "%H:%M:%S"));
      }
	.updateProgressBar(pb, i/nb_simul, text)
      }
    }
    if (progress_bar){
      close(pb)
    }
    if (M>1){
      particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
    }
    else{
      particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
    }
    dim(particle_dist_mat)<-c(nb_simul,M)
    tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
    tab_weight2=.replicate_tab(tab_weight,M)
    if (verbose==TRUE){
      write.table(as.numeric(new_tolerance),file=paste("tolerance_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.matrix(cbind(tab_weight2,simul_below_tol)),file=paste("output_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
	    print(paste("step ",kstep," completed - tol =",new_tolerance,sep=""))
    }
  }
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),weights=tab_weight2/sum(tab_weight2),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}

.all_unif<-function(prior){
	res=TRUE
	for (i in 1:length(prior)){
		res=res&&(prior[[i]][1]=="unif")
	}
res
}

## function to sample in the prior distributions using a Latin Hypercube sample
###############################################################################
.ABC_rejection_lhs<-function(model,prior,nb_simul,tab_unfixed_param,use_seed,seed_count){
  #library(lhs)
  tab_simul_summarystat=NULL
  tab_param=NULL
  l=length(prior)
  nparam=length(tab_unfixed_param[tab_unfixed_param])
  random_tab=NULL
  all_unif_prior=.all_unif(prior)
  if (all_unif_prior){
	  random_tab=randomLHS(nb_simul,nparam)
  }
  
  for (i in 1:nb_simul){
    param=array(0,l)
    if (!all_unif_prior){
	param=.sample_prior(prior)
    }
    else{
      count=1
      for (j in 1:l){
	if (tab_unfixed_param[j]){
        	param[j]=as.numeric(prior[[j]][2])+(as.numeric(prior[[j]][3])-as.numeric(prior[[j]][2]))*random_tab[i,count]
		count=count+1
	}
	else{
		param[j]=as.numeric(prior[[j]][2])
	}
      }
    }
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summarystat=model(param)
    tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
    if (use_seed) {
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    else{
      tab_param=rbind(tab_param,param)
    }
  }
  cbind(tab_param,tab_simul_summarystat)
}

## function to compute particle weights without normalizing to 1
################################################################
.compute_weightb<-function(param_simulated,param_previous_step,tab_weight2,prior){
  tab_weight=tab_weight2/sum(tab_weight2)
  vmat=2*var(param_previous_step)
  n_particle=dim(param_previous_step)[1]
  n_new_particle=dim(param_simulated)[1]
  tab_weight_new=array(0,n_new_particle)
  
  l=dim(param_previous_step)[2]
  multi=exp(-0.5*l*log(2*pi))/sqrt(abs(det(vmat)))
  invmat=0.5*solve(vmat)
  for (i in 1:n_particle){
    for (k in 1:n_new_particle){
      temp=param_simulated[k,]-param_previous_step[i,]
      tab_weight_new[k]=tab_weight_new[k]+tab_weight[i]*as.numeric(exp(- t(temp) %*% invmat %*% temp ))
    }
  }
  prior_density=.compute_weight_prior_tab(param_simulated,prior)
  tab_weight_new=prior_density/(multi*tab_weight_new)
  tab_weight_new
}

## function to compute particle weights with unidimensional jumps without normalizing to 1
##########################################################################################
.compute_weightb_uni<-function(param_simulated,param_previous_step,tab_weight2,prior){
  tab_weight=tab_weight2/sum(tab_weight2)
  l=dim(param_previous_step)[2]
  var_array=array(1,l)
  multi=(1/sqrt(2*pi))^l
  for (j in 1:l){
    var_array[j]=4*var(param_previous_step[,j])
    multi=multi*(1/sqrt(var_array[j]/2))
  }
  var_array=as.numeric(var_array)
  n_particle=dim(param_previous_step)[1]
  n_new_particle=dim(param_simulated)[1]
  tab_weight_new=array(0,n_new_particle)
  for (i in 1:n_particle){
    tab_temp=array(tab_weight[i]*multi,n_new_particle)
    for (k in 1:l){
      tab_temp=tab_temp*exp(-(as.numeric(param_simulated[,k])-as.numeric(param_previous_step[i,k]))*(as.numeric(param_simulated[,k])-as.numeric(param_previous_step[i,k]))/var_array[k])
    }
    tab_weight_new=tab_weight_new+tab_temp
  }
  prior_density=.compute_weight_prior_tab(param_simulated,prior)
  tab_weight_new=prior_density/tab_weight_new
  tab_weight_new
}



## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniformc<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count,inside_prior,progress_bar){
  tab_simul_summarystat=NULL
  tab_param=NULL
  k_acc=0
  
  startb = Sys.time()
  # progress bar
  if (progress_bar){
    pb <- .progressBar(width=50)
  }
  duration = 0;
  for (i in 1:nb_simul){
    l=dim(param_previous_step)[2]
    if (!inside_prior){
      k_acc=k_acc+1
      # pick a particle
      param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
      # move it
      param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
    }
    else{
      test=FALSE
      counter=0
      while ((!test)&&(counter<100)){
        counter=counter+1
        k_acc=k_acc+1
        # pick a particle
        param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
        # move it
        param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
        test=.is_included(param_moved,prior)
      }
      if (counter==100){
        stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
      }
      
    }
    param=param_previous_step[1,]
    param[tab_unfixed_param]=param_moved
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summarystat=model(param)
    tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
    if (use_seed) {
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    else{
      tab_param=rbind(tab_param,param)
    }
    # for progressbar message and time evaluation
    if (progress_bar){
    	duration = difftime(Sys.time(), startb, units="secs")
    	text="";
    	if (i==nb_simul) {
      	text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
    	} 
    	else {
     	 text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/i*(nb_simul-i), tz="GMT"), "%H:%M:%S"));
    	}
    	.updateProgressBar(pb, i/nb_simul, text)
    }
  }
  if (progress_bar){
	  close(pb)
  }
  list(cbind(tab_param,tab_simul_summarystat),nb_simul/k_acc)
}


## sequential algorithm of Lenormand et al. 2012 
################################################
.ABC_Lenormand<-function(model,prior,nb_simul,summary_stat_target,use_seed,verbose,alpha=0.5,p_acc_min=0.05,seed_count=0,inside_prior=TRUE,progress_bar=FALSE){  
  
  ## checking errors in the inputs
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.")
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.")
  if(!is.vector(p_acc_min)) stop("'p_acc_min' has to be a number.")
  if(length(p_acc_min)>1) stop("'p_acc_min' has to be a number.")
  if (p_acc_min<=0) stop ("'p_acc_min' has to be between 0 and 1.")
  if (p_acc_min>=1) stop ("'p_acc_min' has to be between 0 and 1.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(inside_prior)) stop("'inside_prior' has to be boolean.")
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")
  
  start = Sys.time()
  if (progress_bar){
	  print("    ------ Lenormand et al. (2012)'s algorithm ------") 
  }
  
  
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
     if (prior[[i]][1]!="unif"){
	stop("Prior distributions must be uniform to use the Lenormand et al. (2012)'s algorithm.")
     }
  }
  n_alpha=ceiling(nb_simul*alpha)
  
  ## step 1
  # ABC rejection step with LHS
  tab_ini=.ABC_rejection_lhs(model,prior,nb_simul,tab_unfixed_param,use_seed,seed_count)
  seed_count=seed_count+nb_simul
  sd_simul=sapply(as.data.frame(tab_ini[,(nparam+1):(nparam+nstat)]),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
  
  # selection of the alpha quantile closest simulations
  simul_below_tol=NULL
  simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,as.matrix(as.matrix(tab_ini)[,1:nparam]),as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul,alpha))
  simul_below_tol=simul_below_tol[1:n_alpha,] # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile
  
  # initially, weights are equal
  tab_weight=array(1,n_alpha)
  
  tab_dist=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  tol_next=max(tab_dist)
  if (verbose==TRUE){
    write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(tol_next),file="tolerance_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	print("step 1 completed")
  }
  
  ## following steps
  p_acc=p_acc_min+1
  nb_simul_step=nb_simul-n_alpha
  it=1
  while (p_acc>p_acc_min){
    it=it+1
    simul_below_tol2=NULL
    tab_inic=.ABC_launcher_not_uniformc(model,prior,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),tab_unfixed_param,tab_weight/sum(tab_weight),nb_simul_step,use_seed,seed_count,inside_prior,progress_bar)
    tab_ini=as.matrix(tab_inic[[1]])
    tab_ini=as.numeric(tab_ini)
    dim(tab_ini)=c(nb_simul_step,(nparam+nstat))
    seed_count=seed_count+nb_simul_step
    if (!inside_prior){
      tab_weight2=.compute_weightb(as.matrix(as.matrix(tab_ini[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(simul_below_tol[,1:nparam])[,tab_unfixed_param]),tab_weight/sum(tab_weight),prior)
    }
    else{
      tab_weight2=tab_inic[[2]]*(.compute_weightb(as.matrix(as.matrix(tab_ini[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(simul_below_tol[,1:nparam])[,tab_unfixed_param]),tab_weight/sum(tab_weight),prior))
    }
    simul_below_tol2=rbind(as.matrix(simul_below_tol),as.matrix(tab_ini))
    tab_weight=c(tab_weight,tab_weight2)
    tab_dist2=.compute_dist(summary_stat_target,as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul)
    p_acc=length(tab_dist2[tab_dist2<=tol_next])/nb_simul_step
    tab_dist=c(tab_dist,tab_dist2)
    tol_next=sort(tab_dist)[n_alpha]
    simul_below_tol2=simul_below_tol2[tab_dist<=tol_next,]
    tab_weight=tab_weight[tab_dist<=tol_next]
    tab_weight=tab_weight[1:n_alpha]
    tab_dist=tab_dist[tab_dist<=tol_next]
    tab_dist=tab_dist[1:n_alpha]
    simul_below_tol=matrix(0,n_alpha,(nparam+nstat))
    for (i1 in 1:n_alpha){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    if (verbose==TRUE){
      write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(p_acc),file=paste("p_acc_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(tol_next),file=paste("tolerance_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
    	print(paste("step ",it," completed - p_acc = ",p_acc,sep=""))
    }
  }
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}



## FUNCTION ABC_sequential: Sequential ABC methods (Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
###################################################################################################################################

.ABC_sequential <-function(method,model,prior,nb_simul,summary_stat_target,use_seed,verbose,...){
    options(scipen=50)

    ## general function regrouping the different sequential algorithms     
    ## [Beaumont et al., 2009] Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika,96(4):983-990.
    ## [Drovandi & Pettitt 2011] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. Biometrics, 67(1):225-233.
    ## [Del Moral et al. 2012] Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation, Statistics and Computing., 22(5):1009-1020.
    ## [Lenormand et al. 2012] Lenormand, M., Jabot, F., Deffuant G. (2012). Adaptive approximate Bayesian computation for complex models, submitted to Comput. Stat. )
    return(switch(EXPR = method,
	    Beaumont = .ABC_PMC(model,prior,nb_simul,summary_stat_target,use_seed,verbose,...),
	    Drovandi = .ABC_Drovandi(model,prior,nb_simul,summary_stat_target,use_seed,verbose,...),
	    Delmoral = .ABC_Delmoral(model,prior,nb_simul,summary_stat_target,use_seed,verbose,...),
	    Lenormand = .ABC_Lenormand(model,prior,nb_simul,summary_stat_target,use_seed,verbose,...)))

    options(scipen=0)
}



## function to move a particle with a unidimensional uniform jump
#################################################################
.move_particle_uni_uniform<-function(param_picked,sd_array,prior){
  test=FALSE
  res=param_picked
  counter=0
  while ((!test)&&(counter<100)){
    counter=counter+1
    for (i in 1:length(param_picked)){
      res[i]=runif(n = 1, min = param_picked[i]-sd_array[i],max=param_picked[i]+sd_array[i])
    }
    test=.is_included2(res,prior)
  }
  if (counter==100){
    stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
  }
  res
}

.make_proposal_range<-function(prior){
	res=NULL
	for (i in 1:length(prior)){
		if (prior[[i]][1]=="unif"){
			res[i]=(as.numeric(prior[[i]][3])-as.numeric(prior[[i]][2]))/50
		}
		if (prior[[i]][1]=="normal"){
			res[i]=as.numeric(prior[[i]][3])/20
		}
		if (prior[[i]][1]=="lognormal"){
			res[i]=exp(as.numeric(prior[[i]][3])*as.numeric(prior[[i]][3]))/20
		}
		if (prior[[i]][1]=="exponential"){
			res[i]=1/(20*as.numeric(prior[[i]][2]))
		}
	}
res
}

## ABC-MCMC algorithm of Marjoram et al. 2003
#############################################
.ABC_MCMC<-function(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,dist_max=0,tab_normalization=summary_stat_target,proposal_range=vector(mode = "numeric", length = length(prior)),seed_count=0,progress_bar=FALSE){
  
  ## checking errors in the inputs
  if(!is.vector(dist_max)) stop("'dist_max' has to be a number.")
  if(length(dist_max)>1) stop("'dist_max' has to be a number.")
  if (dist_max<0) stop ("'dist_max' has to be positive.")
  if(!is.vector(tab_normalization)) stop("'tab_normalization' has to be a vector.")
  if(length(tab_normalization)!=length(summary_stat_target)) stop("'tab_normalization' must have the same length as 'summary_stat_target'.")
  if(!is.vector(proposal_range)) stop("'proposal_range' has to be a vector.")
  if(length(proposal_range)!=length(prior)) stop("'proposal_range' must have the same length as the number of model parameters.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")

  start = Sys.time()
  if (progress_bar){
	  print("    ------ Marjoram et al. (2003)'s algorithm ------") 
  }
 
  if (sum(abs(tab_normalization-summary_stat_target))==0){
	print("Warning: summary statistics are normalized by default through a division by the target summary statistics - it may not be appropriate to your case.")
	print("Consider providing normalization constants for each summary statistics in the option 'tab_normalization' or using the method 'Marjoram' which automatically determines these constants.")
  }
  for (ii in 1:length(tab_normalization)){
	if (tab_normalization[ii]==0){
		tab_normalization[ii]=1
	}
  }

  if (sum(proposal_range)==0){
	proposal_range=.make_proposal_range(prior)
	print("Warning: default values for proposal distributions are used - they may not be appropriate to your case.")
	print("Consider providing proposal range constants for each parameter in the option 'proposal_range' or using the method 'Marjoram' which automatically determines these constants.")
  }



  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_simul_summary_stat=NULL
  tab_param=NULL
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  
  # initial draw of a particle below the tolerance dist_max
  test=FALSE
  dist_simul=NULL

  

  if (dist_max==0){
    param=.sample_prior(prior)
    if (use_seed) {
      param=c((seed_count+1),param)
    }
    simul_summary_stat=model(param)
    dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),tab_normalization)
    dist_max=dist_simul/2
    seed_count=seed_count+1
    print("Warning: a default value for the tolerance has been computed - it may not be appropriate to your case.")
    print("Consider providing a tolerance value in the option 'dist_max' or using the method 'Marjoram' which automatically determines this value.")
  }

  while (!test){
    param=.sample_prior(prior)
    if (use_seed) {
      param=c((seed_count+1),param)
    }
    simul_summary_stat=model(param)
    dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),tab_normalization)
    if (dist_simul<dist_max){
      test=TRUE
    }
    seed_count=seed_count+1
  }
  tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
  tab_param=rbind(tab_param,param)
  if (use_seed) {
    tab_param=tab_param[,2:(nparam+1)]
  }
  tab_simul_ini=as.numeric(simul_summary_stat)
  param_ini=tab_param
  dist_ini=dist_simul
  if (verbose==TRUE){
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(NULL,file="output_mcmc",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("initial draw performed ")
  }
  
  # chain run
  # progress bar
  if (progress_bar){
	  pb <- .progressBar(width=50)
	  duration = 0;
  }
  
  tab_param=param_ini
  tab_simul_summary_stat=tab_simul_ini
  tab_dist=as.numeric(dist_ini)
  for (is in 2:n_obs){
    for (i in 1:n_between_sampling){
      param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior)
      if (use_seed) {
        param=c(seed_count,param)
      }
      simul_summary_stat=model(param)
      if (use_seed) {
        param=param[2:(nparam+1)]
      }
      dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),tab_normalization)
      if (dist_simul<dist_max){
        param_ini=param
        tab_simul_ini=as.numeric(simul_summary_stat)
        dist_ini=dist_simul
      }
      seed_count=seed_count+1
    }
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
    tab_param=rbind(tab_param,as.numeric(param_ini))
    tab_dist=rbind(tab_dist,as.numeric(dist_ini))
    if (verbose==TRUE){
    	write.table(c(as.numeric(param_ini),tab_simul_ini,as.numeric(dist_ini)),file="output_mcmc",row.names=F,col.names=F,quote=F,append=T)
    }
    if (progress_bar){
	    # for progressbar message and time evaluation
	    duration = difftime(Sys.time(), start, units="secs")
	    text="";
	    if (is==n_obs) {
	      text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
	    } 
	    else {
	      text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/is*(n_obs-is), tz="GMT"), "%H:%M:%S"));
	    }
	    .updateProgressBar(pb, is/n_obs, text)
    }
  }
  if (progress_bar){
	  close(pb)
  }
  tab_param2=matrix(0,dim(tab_param)[1],dim(tab_param)[2])
  for (i in 1:dim(tab_param)[1]){
    for (j in 1:dim(tab_param)[2]){
      tab_param2[i,j]=tab_param[i,j]
    }
  }
  tab_simul_summary_stat2=matrix(0,dim(tab_simul_summary_stat)[1],dim(tab_simul_summary_stat)[2])
  for (i in 1:dim(tab_simul_summary_stat)[1]){
    for (j in 1:dim(tab_simul_summary_stat)[2]){
      tab_simul_summary_stat2[i,j]=tab_simul_summary_stat[i,j]
    }
  }
  tab_dist2=array(0,length(tab_dist))
  for (i in 1:length(tab_dist)){
    tab_dist2[i]=tab_dist[i]
  }
  list(param=as.matrix(tab_param2),stats=as.matrix(tab_simul_summary_stat2),dist=tab_dist2,stats_normalization=as.numeric(tab_normalization),epsilon=max(tab_dist),nsim=(seed_count-seed_count_ini),n_between_sampling=n_between_sampling,computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}

## ABC-MCMC2 algorithm of Marjoram et al. 2003 with automatic determination of the tolerance and proposal range following Wegmann et al. 2009
############################################################################################################################################
.ABC_MCMC2<-function(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,seed_count=0,progress_bar=FALSE){
  
  
  ## checking errors in the inputs
  if(!is.vector(n_calibration)) stop("'n_calibration' has to be a number.")
  if(length(n_calibration)>1) stop("'n_calibration' has to be a number.")
  if (n_calibration<1) stop ("'n_calibration' has to be positive.")
  n_calibration=floor(n_calibration)
  if(!is.vector(tolerance_quantile)) stop("'tolerance_quantile' has to be a number.")
  if(length(tolerance_quantile)>1) stop("'tolerance_quantile' has to be a number.")
  if (tolerance_quantile<=0) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if (tolerance_quantile>=1) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if(!is.vector(proposal_phi)) stop("'proposal_phi' has to be a number.")
  if(length(proposal_phi)>1) stop("'proposal_phi' has to be a number.")
  if (proposal_phi<=0) stop ("'proposal_phi' has to be positive.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")

  
  start = Sys.time()
  if (progress_bar){
	  print("    ------ Marjoram et al. (2003)'s algorithm with modifications drawn from Wegmann et al. (2009) related to automatization ------") 
  }
  
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_simul_summary_stat=NULL
  tab_param=NULL
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  
  # initial draw of a particle
  for (i in 1:(n_calibration)){
    param=.sample_prior(prior)
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summary_stat=model(param)
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
    tab_param=rbind(tab_param,param)
  }
  seed_count=seed_count+n_calibration
  if (use_seed) {
    tab_param=tab_param[,2:(nparam+1)]
  }
  sd_simul=array(0,nstat)
  for (i in 1:nstat){
    sd_simul[i]=sd(tab_simul_summary_stat[,i])
  }
  simuldist=.compute_dist(summary_stat_target,tab_simul_summary_stat,sd_simul)
  ord_sim=order(simuldist,decreasing=F)
  nmax=ceiling(tolerance_quantile*n_calibration)
  dist_max=simuldist[(ord_sim[nmax])]
  tab_param=tab_param[(ord_sim[1:nmax]),]
  proposal_range=vector(mode = "numeric", length = nparam)
  for (i in 1:nparam){
    proposal_range[i]=sd(as.matrix(tab_param)[,i])*proposal_phi/2
  }
  n_ini=sample(nmax,1)
  tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
  dist_ini=simuldist[(ord_sim[n_ini])]
  param_ini=as.matrix(tab_param)[n_ini,]
  if (verbose==TRUE){
    write.table((seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(NULL,file="output_mcmc",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("initial calibration performed ")
  }
  
  # chain run
  # progress bar
  startb = Sys.time()
    if (progress_bar){
	pb <- .progressBar(width=50)
	duration = 0;
    }
  
  tab_param=param_ini
  tab_simul_summary_stat=tab_simul_ini
  tab_dist=as.numeric(dist_ini)
  seed_count=seed_count+1
  for (is in 2:n_obs){
    for (i in 1:n_between_sampling){
      param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior)
      if (use_seed) {
        param=c(seed_count,param)
      }	
      simul_summary_stat=model(param)
      if (use_seed) {
        param=param[2:(nparam+1)]
      }
      dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),sd_simul)
      if (dist_simul<dist_max){
        param_ini=param
        tab_simul_ini=as.numeric(simul_summary_stat)
        dist_ini=dist_simul
      }
      seed_count=seed_count+1
    }
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
    tab_param=rbind(tab_param,as.numeric(param_ini))
    tab_dist=rbind(tab_dist,as.numeric(dist_ini))
    if (verbose==TRUE){
    	write.table(c(as.numeric(param_ini),tab_simul_ini,as.numeric(dist_ini)),file="output_mcmc",row.names=F,col.names=F,quote=F,append=T)
    }
    if (progress_bar){
	    # for progressbar message and time evaluation
		    duration = difftime(Sys.time(), startb, units="secs")
	    text="";
	    if (is==n_obs) {
	      text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
	    } 
	    else {
	      text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/is*(n_obs-is), tz="GMT"), "%H:%M:%S"));
	    }
	    .updateProgressBar(pb, is/n_obs, text)
    }
  }
  if (progress_bar){
	  close(pb)
  }
  tab_param2=matrix(0,dim(tab_param)[1],dim(tab_param)[2])
  for (i in 1:dim(tab_param)[1]){
    for (j in 1:dim(tab_param)[2]){
      tab_param2[i,j]=tab_param[i,j]
    }
  }
  tab_simul_summary_stat2=matrix(0,dim(tab_simul_summary_stat)[1],dim(tab_simul_summary_stat)[2])
  for (i in 1:dim(tab_simul_summary_stat)[1]){
    for (j in 1:dim(tab_simul_summary_stat)[2]){
      tab_simul_summary_stat2[i,j]=tab_simul_summary_stat[i,j]
    }
  }
  tab_dist2=array(0,length(tab_dist))
  for (i in 1:length(tab_dist)){
    tab_dist2[i]=tab_dist[i]
  }
  list(param=as.matrix(tab_param2),stats=as.matrix(tab_simul_summary_stat2),dist=tab_dist2,stats_normalization=as.numeric(sd_simul),epsilon=max(tab_dist),nsim=(seed_count-seed_count_ini),n_between_sampling=n_between_sampling,computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}


## ABC-MCMC3 algorithm of Wegmann et al. 2009 - the PLS step is drawn from the manual of ABCtoolbox (figure 9) - NB: for consistency with ABCtoolbox, AM11-12 are not implemented in the algorithm
#################################################################################################################################################################################################
.ABC_MCMC3<-function(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,numcomp=0,seed_count=0,progress_bar=FALSE){  
  
  ## checking errors in the inputs
  if(!is.vector(n_calibration)) stop("'n_calibration' has to be a number.")
  if(length(n_calibration)>1) stop("'n_calibration' has to be a number.")
  if (n_calibration<1) stop ("'n_calibration' has to be positive.")
  n_calibration=floor(n_calibration)
  if(!is.vector(tolerance_quantile)) stop("'tolerance_quantile' has to be a number.")
  if(length(tolerance_quantile)>1) stop("'tolerance_quantile' has to be a number.")
  if (tolerance_quantile<=0) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if (tolerance_quantile>=1) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if(!is.vector(proposal_phi)) stop("'proposal_phi' has to be a number.")
  if(length(proposal_phi)>1) stop("'proposal_phi' has to be a number.")
  if (proposal_phi<=0) stop ("'proposal_phi' has to be positive.")
  
  if(!is.vector(numcomp)) stop("'numcomp' has to be a number.")
  if(length(numcomp)>1) stop("'numcomp' has to be a number.")
  if (numcomp<0) stop ("'numcomp' has to be positive.")
  if (numcomp>length(summary_stat_target)) stop ("'numcomp' has to smaller or equal to the number of summary statistics.")
  numcomp=floor(numcomp)
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean.")

  #library(pls)
  #library(MASS)
  
  start = Sys.time()
  if (progress_bar){
	print("    ------ Wegmann et al. (2009)'s algorithm ------") 
  }
  
  
  ##AM1
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  if (nstat<=1){
	stop("A single summary statistic is used, use the method 'Marjoram' instead")
  }
  if (numcomp==0){
    numcomp=nstat
  }
  tab_simul_summary_stat=NULL
  tab_param=NULL
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  if (length(tab_unfixed_param[tab_unfixed_param])<=1){
	stop("A single parameter is varying, use the method 'Marjoram' instead")
  }
  
  # initial draw of a particle
  for (i in 1:(n_calibration)){
    param=.sample_prior(prior)
    if (use_seed) {
      param=c((seed_count+i),param)
    }
    simul_summary_stat=model(param)
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
    tab_param=rbind(tab_param,param)
  }
  if (use_seed) {
    tab_param=tab_param[,2:(nparam+1)]
  }
  seed_count=seed_count+n_calibration
  if (verbose){
  	write.table((seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    	write.table(NULL,file="output_mcmc",row.names=F,col.names=F,quote=F)
  }
  
  ## AM2: PLS step
  #print("AM2 ")
  #standardize the params
  sparam=as.matrix(tab_param[,tab_unfixed_param])
  ls=dim(sparam)[2]
  if (is.null(ls)){
	ls=1
  }
  for(i in 1:ls){
    sparam[,i]=(sparam[,i]-mean(sparam[,i]))/sd(sparam[,i])
  }
  #force stat in [1,2]
  myMax<-c()
  myMin<-c()
  lambda<-c()
  myGM<-c()
  stat=tab_simul_summary_stat
  #print("stat 1 ")
  #print(stat)
  summary_stat_targ=summary_stat_target
  for (i in 1:nstat){
    myMax<-c(myMax,max(stat[,i]))
    myMin<-c(myMin,min(stat[,i]))
    stat[,i]=1+(stat[,i]-myMin[i])/(myMax[i]-myMin[i])
    summary_stat_targ[i]=1+(summary_stat_targ[i]-myMin[i])/(myMax[i]-myMin[i])
  }
  #print("stat 2 ")
  #print(stat)
  #transform statistics via boxcox
  dmat=matrix(0,n_calibration,(ls+1))
  for(i in 1:nstat){
    d=cbind(as.vector(as.numeric(stat[,i])),as.matrix(sparam))
    for (i1 in 1:n_calibration){
      for (i2 in 1:(ls+1)){
        dmat[i1,i2]=as.numeric(d[i1,i2])
      }
    }
    
    save(dmat,file="dmat.RData")
    load("dmat.RData",.GlobalEnv)
    file.remove("dmat.RData")
    
    mylm<-lm(as.formula(as.data.frame(dmat)),data=as.data.frame(dmat))
    #mylm<-lm(as.formula(as.data.frame(dmat)))
    #mylm<-lm(stat[,i]~as.matrix(sparam))
    myboxcox<-boxcox(mylm,lambda=seq(-20,100,1/10),interp=T,eps=1/50,plotit=FALSE)
    lambda<-c(lambda,myboxcox$x[myboxcox$y==max(myboxcox$y)])
    myGM<-c(myGM, exp(mean(log(stat[,i]))))
  }
  #standardize the BC-stat
  myBCMeans<-c()
  myBCSDs<-c()
  for(i in 1:nstat){
    stat[,i]<-((stat[,i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)))
    summary_stat_targ[i]<-((summary_stat_targ[i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)))
    myBCSDs<-c(myBCSDs, sd(stat[,i]))
    myBCMeans<-c(myBCMeans, mean(stat[,i]))
    stat[,i]<-(stat[,i]-myBCMeans[i])/myBCSDs[i]
    summary_stat_targ[i]<-(summary_stat_targ[i]-myBCMeans[i])/myBCSDs[i]
  }
  #perform pls
  myPlsr<-plsr(as.matrix(sparam)~as.matrix(stat), scale=F,ncomp=numcomp,validation='LOO')
  pls_transformation=matrix(0,numcomp,nstat)
  for (i in 1:numcomp){
    pls_transformation[i,]=as.numeric(myPlsr$loadings[,i])
  }
  
  ## AM3
  #print("AM3 ")
  summary_stat_targ=t(pls_transformation %*% as.vector(summary_stat_targ))
  stat_pls=t(pls_transformation %*% t(stat))
  simuldist=.compute_dist(summary_stat_targ,stat_pls,rep(1,numcomp))
  
  ## AM4
  #print("AM4 ")
  ord_sim=order(simuldist,decreasing=F)
  nmax=ceiling(tolerance_quantile*n_calibration)
  dist_max=simuldist[(ord_sim[nmax])]
  
  tab_param=tab_param[(ord_sim[1:nmax]),]
  proposal_range=vector(mode = "numeric", length = nparam)
  for (i in 1:nparam){
    proposal_range[i]=sd(tab_param[,i])*proposal_phi/2
  }
  if (progress_bar){
	  print("initial calibration performed ")
  }
  
  ## AM5: chain run
  # progress bar
  startb = Sys.time()
  if (progress_bar){
	  pb <- .progressBar(width=50)
	  duration = 0;
  }
  
  #print("AM5 ")
  n_ini=sample(nmax,1)
  tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
  param_ini=tab_param[n_ini,]
  tab_param=param_ini
  tab_simul_summary_stat=tab_simul_ini
  dist_ini=simuldist[(ord_sim[n_ini])]
  tab_dist=as.numeric(dist_ini)
  seed_count=seed_count+1
  for (is in 2:n_obs){
    for (i in 1:n_between_sampling){
      ## AM6
      #print("AM6 ")
      param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior)
      if (use_seed) {
        param=c(seed_count,param)
      }
      ## AM7	
      #print("AM7 ")
      simul_summary_stat=model(param)
      simul_summary_stat_output=simul_summary_stat
      if (use_seed) {
        param=param[2:(nparam+1)]
      }
      for (ii in 1:nstat){
        simul_summary_stat[ii]=1+(simul_summary_stat[ii]-myMin[ii])/(myMax[ii]-myMin[ii])
      }
      for(ii in 1:nstat){
        simul_summary_stat[ii]<-((simul_summary_stat[ii]^lambda[ii]) - 1)/(lambda[ii]*(myGM[ii]^(lambda[ii]-1)))
        simul_summary_stat[ii]<-(simul_summary_stat[ii]-myBCMeans[ii])/myBCSDs[ii]
      }
      simul_summary_stat=as.matrix(simul_summary_stat)
      dim(simul_summary_stat)<-c(nstat,1)
      simul_summary_stat=t(pls_transformation %*% simul_summary_stat)
      dist_simul=.compute_dist_single(summary_stat_targ,as.numeric(simul_summary_stat),rep(1,numcomp))
      ## AM8-9
      #print("AM8-9 ")
      if (dist_simul<dist_max){
        param_ini=param
        tab_simul_ini=as.numeric(simul_summary_stat_output)
        dist_ini=dist_simul
      }
      seed_count=seed_count+1
    }
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
    tab_param=rbind(tab_param,as.numeric(param_ini))
    tab_dist=rbind(tab_dist,as.numeric(dist_ini))
    if (verbose==TRUE){
    	write.table(c(as.numeric(param_ini),tab_simul_ini,as.numeric(dist_ini)),file="output_mcmc",row.names=F,col.names=F,quote=F,append=T)
    }
    if (progress_bar){
	    # for progressbar message and time evaluation
	    duration = difftime(Sys.time(), startb, units="secs")
	    text="";
	    if (is==n_obs) {
	      text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
	    } 
	    else {
	      text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/is*(n_obs-is), tz="GMT"), "%H:%M:%S"));
	    }
	    .updateProgressBar(pb, is/n_obs, text)
    }
  }
  if (progress_bar){
	  close(pb)
  }
  tab_param2=matrix(0,dim(tab_param)[1],dim(tab_param)[2])
  for (i in 1:dim(tab_param)[1]){
    for (j in 1:dim(tab_param)[2]){
      tab_param2[i,j]=tab_param[i,j]
    }
  }
  tab_simul_summary_stat2=matrix(0,dim(tab_simul_summary_stat)[1],dim(tab_simul_summary_stat)[2])
  for (i in 1:dim(tab_simul_summary_stat)[1]){
    for (j in 1:dim(tab_simul_summary_stat)[2]){
      tab_simul_summary_stat2[i,j]=tab_simul_summary_stat[i,j]
    }
  }
  tab_dist2=array(0,length(tab_dist))
  for (i in 1:length(tab_dist)){
    tab_dist2[i]=tab_dist[i]
  }
  list(param=as.matrix(tab_param2),stats=as.matrix(tab_simul_summary_stat2),dist=tab_dist2,epsilon=max(tab_dist),nsim=(seed_count-seed_count_ini),n_between_sampling=n_between_sampling,min_stats=myMin,max_stats=myMax,lambda=lambda,geometric_mean=myGM,boxcox_mean=myBCMeans,boxcox_sd=myBCSDs,pls_transform=pls_transformation,n_component=numcomp,computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}

## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al. 2009)
##############################################################################
.ABC_mcmc_internal <-function(method,model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,...){
	options(scipen=50)

	    return(switch(EXPR = method,
	       Marjoram_original = .ABC_MCMC(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,...),
	       Marjoram = .ABC_MCMC2(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,,...),
	       Wegmann = .ABC_MCMC3(model,prior,n_obs,n_between_sampling,summary_stat_target,use_seed,verbose,...)))

	options(scipen=0)
}





###################### parallel functions ###############


.ABC_rejection_internal_cluster<-function(model,prior,nb_simul,seed_count=0,n_cluster=1){
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  
  tab_simul_summarystat=NULL
  tab_param=NULL
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  if (npar>0){
   for (irun in 1:npar){
    for (i in 1:(100*n_cluster)){
      l=length(prior)
      param=.sample_prior(prior)
      #if (use_seed) { # NB: we force the value use_seed=TRUE
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+100*n_cluster
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:(100*n_cluster)){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
   }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      l=length(prior)
      param=.sample_prior(prior)
      #if (use_seed) { # NB: we force the value use_seed=TRUE
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  cbind(tab_param,tab_simul_summarystat)
}




## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################

.ABC_rejection_cluster <-function(model,prior,nb_simul,seed_count=0,n_cluster=1,verbose){
    if (verbose){
	write.table(NULL,file="output",row.names=F,col.names=F,quote=F)
    }
     
	#library(parallel)

	start = Sys.time()
	
	options(scipen=50)
	cl <- makeCluster(getOption("cl.cores", n_cluster))

	tab_simul_summarystat=NULL
	tab_param=NULL
	list_param=list(NULL)
	npar=floor(nb_simul/(100*n_cluster))
	n_end=nb_simul-(npar*100*n_cluster)
	if (npar>0){
	 for (irun in 1:npar){
	  paramtemp=NULL
	  simultemp=NULL
	  for (i in 1:(100*n_cluster)){
		l=length(prior)
		param=.sample_prior(prior)
		#if (use_seed) { # NB: we force the value use_seed=TRUE
		param=c((seed_count+i),param)
		list_param[[i]]=param
		tab_param=rbind(tab_param,param[2:(l+1)])
		paramtemp=rbind(paramtemp,param[2:(l+1)])
	  }
	  seed_count=seed_count+n_cluster
	  list_simul_summarystat=parLapplyLB(cl,list_param,model)
	  for (i in 1:(100*n_cluster)){
		tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
		simultemp=rbind(simultemp,as.numeric(list_simul_summarystat[[i]]))
	  }
	  if (verbose){
		write.table(cbind(as.matrix(paramtemp),as.matrix(simultemp)),file="output",row.names=F,col.names=F,quote=F,append=T)
	  }
         }
	}
	if (n_end>0){
	  #stopCluster(cl)
	  #cl <- makeCluster(getOption("cl.cores", 1))
	  list_param=list(NULL)
	  paramtemp=NULL
	  simultemp=NULL
	  for (i in 1:n_end){
		l=length(prior)[1]
		param=.sample_prior(prior)
		#if (use_seed) { # NB: we force the value use_seed=TRUE
		param=c((seed_count+i),param)
		list_param[[i]]=param
		tab_param=rbind(tab_param,param[2:(l+1)])
		paramtemp=rbind(paramtemp,param[2:(l+1)])
	  }
	  seed_count=seed_count+n_end
	  list_simul_summarystat=parLapplyLB(cl,list_param,model)
	  for (i in 1:n_end){
		tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
		simultemp=rbind(simultemp,as.numeric(list_simul_summarystat[[i]]))
	  }
 	  if (verbose){
		write.table(cbind(as.matrix(paramtemp),as.matrix(simultemp)),file="output",row.names=F,col.names=F,quote=F,append=T)
	  }
    	  stopCluster(cl)
	}
	else{
	  stopCluster(cl)
	}
	options(scipen=0)
	sd_simul=sapply(as.data.frame(tab_simul_summarystat),sd)
list(param=as.matrix(tab_param),stats=as.matrix(tab_simul_summarystat),weights=array(1/nb_simul,nb_simul),stats_normalization=as.numeric(sd_simul),nsim=nb_simul,computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}




## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniform_cluster<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,seed_count,inside_prior,n_cluster){
  tab_simul_summarystat=NULL
  tab_param=NULL
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  if (npar>0){
    for (irun in 1:npar){
      for (i in 1:(100*n_cluster)){
        l=dim(param_previous_step)[2]
        if (!inside_prior){
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
        }
        else{
          test=FALSE
          counter=0
          while ((!test)&&(counter<100)){
            counter=counter+1
            # pick a particle
            param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
            # move it
            param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
            test=.is_included(param_moved,prior)
          }
          if (counter==100){
            stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
          }
        }
        param=param_previous_step[1,]
        param[tab_unfixed_param]=param_moved
        param=c((seed_count+i),param)
        list_param[[i]]=param
        tab_param=rbind(tab_param,param[2:(l+1)])
      }
      seed_count=seed_count+n_cluster
      list_simul_summarystat=parLapplyLB(cl,list_param,model)
      for (i in 1:(100*n_cluster)){
        tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
      }
    }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      l=dim(param_previous_step)[2]
      if (!inside_prior){
        # pick a particle
        param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
        # move it
        param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
      }
      else{
        test=FALSE
        counter=0
        while ((!test)&&(counter<100)){
          counter=counter+1
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
          test=.is_included(param_moved,prior)
        }
        if (counter==100){
          stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
        }
      }
      param=param_previous_step[1,]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  cbind(tab_param,tab_simul_summarystat)
}


## function to perform ABC simulations from a non-uniform prior and with unidimensional jumps
#############################################################################################
.ABC_launcher_not_uniform_uni_cluster<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,seed_count,inside_prior,n_cluster){
  tab_simul_summarystat=NULL
  tab_param=NULL
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  l=dim(param_previous_step)[2]
  l_array=dim(param_previous_step[,tab_unfixed_param])[2]
  if (is.null(l_array)){
	l_array=1
  }
  sd_array=array(1,l_array)
  covmat=2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov # computation of a WEIGHTED variance
  for (j in 1:l_array){
    sd_array[j]=sqrt(covmat[j,j])
    
  }
  if (npar>0){
    for (irun in 1:npar){
      for (i in 1:(100*n_cluster)){
        if (!inside_prior){
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
        }
        else{
          test=FALSE
          counter=0
          while ((!test)&&(counter<100)){
            counter=counter+1
            # pick a particle
            param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
            # move it
            param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
            test=.is_included(param_moved,prior)
          }
          if (counter==100){
            stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
          }
        }
        param=param_previous_step[1,]
        param[tab_unfixed_param]=param_moved
        param=c((seed_count+i),param)
        list_param[[i]]=param
        tab_param=rbind(tab_param,param[2:(l+1)])
      }
      seed_count=seed_count+n_cluster
      list_simul_summarystat=parLapplyLB(cl,list_param,model)
      for (i in 1:(100*n_cluster)){
        tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
      }
    }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      if (!inside_prior){
        # pick a particle
        param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
        # move it
        param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
      }
      else{
        test=FALSE
        counter=0
        while ((!test)&&(counter<100)){
          counter=counter+1
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
          test=.is_included(param_moved,prior)
        }
        if (counter==100){
          stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
        }
      }
      param=param_previous_step[1,]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  cbind(tab_param,tab_simul_summarystat)
}

## PMC ABC algorithm: Beaumont et al. Biometrika 2009
#####################################################
.ABC_PMC_cluster<-function(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,seed_count=0,inside_prior=TRUE,tolerance_tab=-1,progress_bar=FALSE){
  
  ## checking errors in the inputs
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(inside_prior)) stop("'inside_prior' has to be boolean.")
  if(!is.vector(tolerance_tab)) stop("'tolerance_tab' has to be a vector.")
  if(tolerance_tab[1]==-1) stop("'tolerance_tab' is missing")
  if (min(tolerance_tab)<=0) stop ("tolerance values have to be strictly positive.")
  lll=length(tolerance_tab)
  if (lll<=1) stop("at least two tolerance values need to be provided.")
  if (min(tolerance_tab[1:(lll-1)]-tolerance_tab[2:lll])<=0) stop ("tolerance values have to decrease.")
  
  if (progress_bar){
	print("    ------ Beaumont et al. (2009)'s algorithm ------")
  }
  start = Sys.time()
  
  seed_count_ini=seed_count
  T=length(tolerance_tab)
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  
  ## step 1
  nb_simul_step=nb_simul
  simul_below_tol=NULL
  while (nb_simul_step>0){
    if (nb_simul_step>1){
      # classic ABC step
      tab_ini=.ABC_rejection_internal_cluster(model,prior,nb_simul_step,seed_count,n_cluster)
      if (nb_simul_step==nb_simul){
        sd_simul=sapply(as.data.frame(tab_ini[,(nparam+1):(nparam+nstat)]),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
      }
      seed_count=seed_count+nb_simul_step
      # selection of simulations below the first tolerance level
      simul_below_tol=rbind(simul_below_tol,.selec_simul(summary_stat_target,as.matrix(as.matrix(tab_ini)[,1:nparam]),as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul,tolerance_tab[1]))
      if (length(simul_below_tol)>0){
        nb_simul_step=nb_simul-dim(simul_below_tol)[1]
      }
    }
    else{
      tab_ini=.ABC_rejection_internal_cluster(model,prior,nb_simul_step,seed_count,n_cluster)
      seed_count=seed_count+nb_simul_step
      if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[1]){
        simul_below_tol=rbind(simul_below_tol,tab_ini)
        nb_simul_step=0
      }
    }
  } # until we get nb_simul simulations below the first tolerance threshold
  # initially, weights are equal
  tab_weight=array(1/nb_simul,nb_simul)
  if (verbose==TRUE){
    write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table((seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  ## steps 2 to T
  for (it in 2:T){
    nb_simul_step=nb_simul
    simul_below_tol2=NULL
    while (nb_simul_step>0){
      if (nb_simul_step>1){
        # Sampling of parameters around the previous particles
        tab_ini=.ABC_launcher_not_uniform_uni_cluster(model,prior,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),tab_unfixed_param,tab_weight,nb_simul_step,seed_count,inside_prior,n_cluster)
        seed_count=seed_count+nb_simul_step
        simul_below_tol2=rbind(simul_below_tol2,.selec_simul(summary_stat_target,as.matrix(as.matrix(tab_ini)[,1:nparam]),as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul,tolerance_tab[it]))
        if (length(simul_below_tol2)>0){
          nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
        }
      }
      else{
        tab_ini=.ABC_launcher_not_uniform_uni_cluster(model,prior,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),tab_unfixed_param,tab_weight,nb_simul_step,seed_count,inside_prior,n_cluster)
        seed_count=seed_count+nb_simul_step
        if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[it]){
          simul_below_tol2=rbind(simul_below_tol2,tab_ini)
          nb_simul_step=0
        }
      }
    } # until we get nb_simul simulations below the it-th tolerance threshold
    # update of particle weights
    tab_weight2=.compute_weight_uni(as.matrix(as.matrix(as.matrix(simul_below_tol2)[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(as.matrix(simul_below_tol)[,1:nparam])[,tab_unfixed_param]),tab_weight,prior)
    # update of the set of particles and of the associated weights for the next ABC sequence
    tab_weight=tab_weight2
    simul_below_tol=matrix(0,nb_simul,(nparam+nstat))
    for (i1 in 1:nb_simul){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    if (verbose==TRUE){
      write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table((seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
   	 print(paste("step ",it," completed",sep=""))
    }
  }
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(as.matrix(simul_below_tol))[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs")))
}


.move_drovandi_ini_cluster<-function(nb_simul_step,simul_below_tol,tab_weight,nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul){
  i_acc=0
  res=NULL
  npar=floor(nb_simul_step/(100*n_cluster))
  n_end=nb_simul_step-(npar*100*n_cluster)
  l=length(prior)
  list_param=list(NULL)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  if (npar>0){
   for (irun in 1:npar){
    tab_param=NULL
    tab_picked=NULL
    for (i in 1:(100*n_cluster)){
      # pick a particle
      simul_picked=.particle_pick(simul_below_tol,tab_weight)
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(simul_below_tol)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_cluster
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:(100*n_cluster)){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
   }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    tab_picked=NULL
    tab_param=NULL
    for (i in 1:n_end){
      # pick a particle
      simul_picked=.particle_pick(simul_below_tol,tab_weight)
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(simul_below_tol)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_end
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  list(res,i_acc)
}

.move_drovandi_end_cluster<-function(nb_simul_step,new_particles,nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul){
  i_acc=0
  res=NULL
  npar=floor(nb_simul_step/(100*n_cluster))
  n_end=nb_simul_step-(npar*100*n_cluster)
  l=length(prior)
  list_param=list(NULL)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  if (npar>0){
   for (irun in 1:npar){
    tab_param=NULL
    tab_picked=NULL
    for (i in 1:(100*n_cluster)){
      # pick a particle
      simul_picked=new_particles[((irun-1)*n_cluster+i),]
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(new_particles)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_cluster
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:(100*n_cluster)){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
   }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    tab_picked=NULL
    tab_param=NULL
    for (i in 1:n_end){
      # pick a particle
      simul_picked=new_particles[(npar*n_cluster+i),]
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(new_particles)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_end
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  list(res,i_acc)
}

.move_drovandi_diversify_cluster<-function(nb_simul_step,new_particles,nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul){
  i_acc=0
  res=NULL
  npar=floor(nb_simul_step/(100*n_cluster))
  n_end=nb_simul_step-(npar*100*n_cluster)
  l=length(prior)
  list_param=list(NULL)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  if (npar>0){
   for (irun in 1:npar){
    tab_param=NULL
    tab_picked=NULL
    for (i in 1:(100*n_cluster)){
      # pick a particle
      simul_picked=new_particles[((irun-1)*n_cluster+i),]
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(new_particles)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_cluster
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:(100*n_cluster)){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
   }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    tab_picked=NULL
    tab_param=NULL
    for (i in 1:n_end){
      # pick a particle
      simul_picked=new_particles[(npar*n_cluster+i),]
      # move it
      param_moved=.move_particleb(simul_picked[1:nparam][tab_unfixed_param],2*var(as.matrix(as.matrix(as.matrix(new_particles)[,1:nparam])[,tab_unfixed_param])),prior)
      param=simul_picked[1:nparam]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      tab_param=rbind(tab_param,param[2:(l+1)])
      list_param[[i]]=param
      tab_picked=rbind(tab_picked,as.numeric(simul_picked))
    }
    seed_count=seed_count+n_end
    # perform simulations
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      # check whether it is below tol_next and undo the move if it is not
      new_simul=c(as.numeric(tab_param[i,]),as.numeric(list_simul_summarystat[[i]]))
      if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        tab_picked[i,]=as.numeric(new_simul)
        i_acc=i_acc+1
      }
    }
    res=rbind(res,tab_picked)
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  res
}


## sequential algorithm of Drovandi & Pettitt 2011 - the proposal used is a multivariate normal (cf paragraph 2.2 - p. 227 in Drovandi & Pettitt 2011)
######################################################################################################################################################
.ABC_Drovandi_cluster<-function(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,seed_count=0,tolerance_tab=-1,alpha=0.5,c=0.01,first_tolerance_level_auto=TRUE,progress_bar=FALSE){
  ## checking errors in the inputs
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.vector(tolerance_tab)) stop("'tolerance_tab' has to be a vector.")
  if(tolerance_tab[1]==-1) stop("'tolerance_tab' is missing")
  if (min(tolerance_tab)<=0) stop ("tolerance values have to be strictly positive.")
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.") 
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.") 
  if(!is.vector(c)) stop("'c' has to be a vector.") 
  if(length(c)>1) stop("'c' has to be a number.") 
  if (c<=0) stop ("'c' has to be between 0 and 1.") 
  if (c>=1) stop ("'c' has to be between 0 and 1.") 
  if(!is.logical(first_tolerance_level_auto)) stop("'first_tolerance_level_auto' has to be boolean.")
  
  if (progress_bar){
	  print("    ------ Drovandi & Pettitt (2011)'s algorithm ------")
  }
  start = Sys.time()
  
  seed_count_ini=seed_count
  n_alpha=ceiling(nb_simul*alpha)
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  if (first_tolerance_level_auto){
    tol_end=tolerance_tab[1]
  }
  else{
    tol_end=tolerance_tab[2]
  }
  
  ## step 1
  nb_simul_step=ceiling(nb_simul/(1-alpha))
  simul_below_tol=NULL
  if (first_tolerance_level_auto){
    # classic ABC step
    tab_ini=.ABC_rejection_internal_cluster(model,prior,nb_simul_step,seed_count,n_cluster)
    sd_simul=sapply(as.data.frame(tab_ini[,(nparam+1):(nparam+nstat)]),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
    seed_count=seed_count+nb_simul_step
    # selection of simulations below the first tolerance level
    simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,as.matrix(as.matrix(tab_ini)[,1:nparam]),as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul,(1-alpha)))
    simul_below_tol=simul_below_tol[1:nb_simul,] # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile
  }
  else{
    nb_simul_step=nb_simul
    while (nb_simul_step>0){
      if (nb_simul_step>1){
        # classic ABC step
        tab_ini=.ABC_rejection_internal_cluster(model,prior,nb_simul_step,seed_count,n_cluster)
        if (nb_simul_step==nb_simul){
          sd_simul=sapply(as.data.frame(tab_ini[,(nparam+1):(nparam+nstat)]),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
        }
        seed_count=seed_count+nb_simul_step
        # selection of simulations below the first tolerance level
        simul_below_tol=rbind(simul_below_tol,.selec_simulb(summary_stat_target,as.matrix(as.matrix(tab_ini)[,1:nparam]),as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul,tolerance_tab[1])) # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
        if (length(simul_below_tol)>0){
          nb_simul_step=nb_simul-dim(simul_below_tol)[1]
        }
      }
      else{
        tab_ini=.ABC_rejection_internal_cluster(model,prior,nb_simul_step,seed_count,n_cluster)
        seed_count=seed_count+nb_simul_step
        if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<=tolerance_tab[1]){
          simul_below_tol=rbind(simul_below_tol,tab_ini)
          nb_simul_step=0
        }
      }
    } # until we get nb_simul simulations below the first tolerance threshold
  }
  # initially, weights are equal
  tab_weight=array(1/nb_simul,nb_simul)
  if (verbose==TRUE){
    write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  
  ## following steps until tol_end is reached
  tol_next=tolerance_tab[1]
  if (first_tolerance_level_auto){
    tol_next=max(.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul))
  } 
  R=1
  l=dim(simul_below_tol)[2]
  it=1
  while (tol_next>tol_end){
    it=it+1
    i_acc=0
    nb_simul_step=n_alpha
    # compute epsilon_next
    tol_next=sort(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul))[(nb_simul-n_alpha)]
    # drop the n_alpha poorest particles
    simul_below_tol2=.selec_simulb(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul,tol_next) # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
    simul_below_tol=matrix(0,(nb_simul-n_alpha),(nparam+nstat))
    for (i1 in 1:(nb_simul-n_alpha)){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    
    md=.move_drovandi_ini_cluster(nb_simul_step,simul_below_tol,tab_weight[1:(nb_simul-n_alpha)],nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul)
    new_particles=md[[1]]
    i_acc=i_acc+md[[2]]
    seed_count=seed_count+nb_simul_step
    if (R>1){
      for (j in 2:R){
        md=.move_drovandi_end_cluster(nb_simul_step,new_particles,nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul)
        new_particles=md[[1]]
        i_acc=i_acc+md[[2]]
        seed_count=seed_count+nb_simul_step
      }
    }
    simul_below_tol3=matrix(0,nb_simul,(nparam+nstat))
    for (i1 in 1:(nb_simul-n_alpha)){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol3[i1,i2]=as.numeric(simul_below_tol[i1,i2])
      }
    }
    for (i1 in (nb_simul-n_alpha+1):nb_simul){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol3[i1,i2]=as.numeric(new_particles[(i1-nb_simul+n_alpha),i2])
      }
    }
    simul_below_tol=simul_below_tol3
    p_acc=max(1,i_acc)/(nb_simul_step*R) # to have a strictly positive p_acc
    Rp=R
    R=ceiling(log(c)/log(1-p_acc))
    if (verbose==TRUE){
      write.table(as.numeric(tol_next),file=paste("tolerance_step",it,sep=""),row.names=F,col.names=F,quote=F) 
      write.table(as.numeric(Rp),file=paste("R_step",it,sep=""),row.names=F,col.names=F,quote=F) 
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(cbind(tab_weight,simul_below_tol),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    tol_next=max(.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul))
    if (progress_bar){
    	print(paste("step ",it," completed - R used = ",Rp," - tol = ",tol_next," - next R used will be ",R,sep=""))
    }
  }
  
  ## final step to diversify the n_alpha particles
  simul_below_tol2=NULL
  for (j in 1:R){
    simul_below_tol=.move_drovandi_diversify_cluster(nb_simul,simul_below_tol,nparam,nstat,tab_unfixed_param,prior,summary_stat_target,tol_next,seed_count,n_cluster,model,sd_simul)
    seed_count=seed_count+nb_simul
  }
  
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs"))) 
}


## rejection algorithm with M simulations per parameter set
############################################################
.ABC_rejection_M_cluster<-function(model,prior,nb_simul,M,seed_count,n_cluster){
  tab_simul_summarystat=NULL
  tab_param=NULL
  list_param=list(NULL)
  npar=floor(nb_simul*M/(100*n_cluster))
  n_end=nb_simul*M-(npar*100*n_cluster)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  l=length(prior)
  param2=array(0,(l+1))
  if (npar>0){
    for (irun in 1:npar){
      for (irun2 in 1:(100*n_cluster)){
	ii=(100*n_cluster*(irun-1)+irun2)%%M
	if ((ii==1)||(M==1)){
		param=.sample_prior(prior)
		param2=c((seed_count+1),param)
		seed_count=seed_count+1
	}
	else{
		param2[1]=param2[1]+1
		seed_count=seed_count+1
  	}
	list_param[[irun2]]=param2
        tab_param=rbind(tab_param,param2[2:(l+1)])
      }
      list_simul_summarystat=parLapplyLB(cl,list_param,model)
      for (i in 1:(100*n_cluster)){
          tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
      }
    }
  }
  if (n_end>0){
    list_param=list(NULL)
    for (irun2 in 1:n_end){
	ii=(100*n_cluster*npar+irun2)%%M
        if ((ii==1)||(M==1)){
		param=.sample_prior(prior)
		param2=c((seed_count+1),param)
		seed_count=seed_count+1
	}
	else{
		param2[1]=param2[1]+1
		seed_count=seed_count+1
  	}
	list_param[[irun2]]=param2
        tab_param=rbind(tab_param,param2[2:(l+1)])
    }
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
        tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
     }
  }
  stopCluster(cl)
  cbind(tab_param,tab_simul_summarystat)
}

## sequential algorithm of Del Moral et al. 2012 - the proposal used is a normal in each dimension (cf paragraph 3.2 in Del Moral et al. 2012)
##############################################################################################################################################
.ABC_Delmoral_cluster<-function(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,alpha=0.9,M=1,nb_threshold=floor(nb_simul/2),tolerance_target=-1,seed_count=0,progress_bar=FALSE){
  ## checking errors in the inputs
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.")
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.")
  if(!is.vector(M)) stop("'M' has to be a number.")
  if(length(M)>1) stop("'M' has to be a number.")
  if (M<1) stop ("'M' has to be a positive integer.")
  M=floor(M)
  if(!is.vector(nb_threshold)) stop("'nb_threshold' has to be a number.")
  if(length(nb_threshold)>1) stop("'nb_threshold' has to be a number.")
  if (nb_threshold<1) stop ("'nb_threshold' has to be a positive integer.")
  nb_threshold=floor(nb_threshold)
  if(!is.vector(tolerance_target)) stop("'tolerance_target' has to be a number.")
  if(length(tolerance_target)>1) stop("'tolerance_target' has to be a number.")
  if (tolerance_target<=0) stop ("'tolerance_target' has to be positive.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  
  start = Sys.time()
  if (progress_bar){
	  print("    ------ Delmoral et al. (2012)'s algorithm ------") 
  }
  

  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  
  # step 1
  # classic ABC step
  simul_below_tol=.ABC_rejection_M_cluster(model,prior,nb_simul,M,seed_count,n_cluster)
  seed_count=seed_count+M*nb_simul
  tab_weight=rep(1/nb_simul,nb_simul)
  ESS=nb_simul
  uu=(1:nb_simul)*M  # to compute sd_simul with only one simulation per parameter set
  sd_simul=sapply(as.data.frame(simul_below_tol[uu,(nparam+1):(nparam+nstat)]),sd)  # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
  l=dim(simul_below_tol)[2]
  if (M>1){
    particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  }
  else{
    particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  }
  dim(particle_dist_mat)<-c(nb_simul,M)
  new_tolerance=max(particle_dist_mat)
  
  tab_weight2=.replicate_tab(tab_weight,M)
  if (verbose==TRUE){ 
    write.table(cbind(tab_weight2,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  
  # following steps
  kstep=1
  l=length(prior)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  while(new_tolerance>tolerance_target){	
    kstep=kstep+1
    # determination of the new tolerance
    ESS_target=alpha*ESS
    
    tolerance_list=sort(as.numeric(names(table(particle_dist_mat))),decreasing=TRUE)
    i=1
    test=FALSE
    while((!test)&&(i<length(tolerance_list))){
      i=i+1
      # computation of new ESS with the new tolerance value
      new_ESS=.compute_ESS(particle_dist_mat,tolerance_list[i])
      # check whether this value is below ESS_targ
      if (new_ESS<ESS_target){
        new_tolerance=tolerance_list[(i-1)]
        test=TRUE
      }
    }
    
    # if effective sample size is too small, resampling of particles
    ESS=.compute_ESS(particle_dist_mat,new_tolerance)
    tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
    tab_below=.compute_below(particle_dist_mat,new_tolerance)
    particles=matrix(0,(nb_simul*M),(nparam+nstat))
    if (ESS<nb_threshold){
      # sample nb_simul particles 
      for (i in 1:nb_simul){
        particles[((1:M)+(i-1)*M),]=as.matrix(.particle_pick_delmoral(simul_below_tol,tab_weight,M))
      }
      simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))
      for (i1 in 1:(nb_simul*M)){
        for (i2 in 1:(nparam+nstat)){
          simul_below_tol[i1,i2]=as.numeric(particles[i1,i2])
        }
      }
      particles=as.matrix(particles[uu,1:nparam])
      if (M>1){
        particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
      }
      else{
        particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
      }
      dim(particle_dist_mat)<-c(nb_simul,M)
      tab_below=.compute_below(particle_dist_mat,new_tolerance)
      # reset their weight to 1/nb_simul
      tab_weight=rep(1/nb_simul,nb_simul)
      ESS=nb_simul
    }
    else{
      particles=as.matrix(simul_below_tol[,1:nparam])
    }
    
    # MCMC move
    npar=floor(nb_simul*M/(100*n_cluster))
    n_end=nb_simul*M-(npar*100*n_cluster)
    covmat=2*cov.wt(as.matrix(as.matrix(particles[uu,tab_unfixed_param])[tab_weight>0,]),as.vector(tab_weight[tab_weight>0]))$cov
    l_array=dim(particles[,tab_unfixed_param])[2]
    if (is.null(l_array)){
     l_array=1
    }
    sd_array=array(1,l_array)
    for (j in 1:l_array){
      sd_array[j]=sqrt(covmat[j,j])
      
    }
    simul_below_tol2=simul_below_tol
    simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))

    param2=array(0,(l+1))
    tab_param=NULL
    tab_simul_summarystat=NULL
    if (npar>0){
      for (irun in 1:npar){
        list_param=list(NULL)
	ic=1
	for (irun2 in 1:(100*n_cluster)){
		ii=(100*n_cluster*(irun-1)+irun2)%%M
		ii2=ceiling((100*n_cluster*(irun-1)+irun2)/M)
		if (tab_weight[ii2]>0){
			if ((ii==1)||(M==1)){
				param_moved=.move_particleb_uni(as.numeric(particles[(ii2*M),tab_unfixed_param]),sd_array,prior)
        			param=particles[(ii2*M),]
        			param[tab_unfixed_param]=param_moved
        			param2=c((seed_count+1),param)
				seed_count=seed_count+1
			}
			else{
				param2[1]=param2[1]+1
				seed_count=seed_count+1
  			}
			list_param[[ic]]=param2
			ic=ic+1
        		tab_param=rbind(tab_param,param2[2:(l+1)])
		}
	}
	list_simul_summarystat=parLapplyLB(cl,list_param,model)
	ic=ic-1
	for (ii in 1:ic){
          tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[ii]]))
      	}
      }
    }
    if (n_end>0){
	list_param=list(NULL)
 	ic=1
	for (irun2 in 1:n_end){
		ii=(100*n_cluster*npar+irun2)%%M
		ii2=ceiling((100*n_cluster*npar+irun2)/M)
		if (tab_weight[ii2]>0){
			if ((ii==1)||(M==1)){
				param_moved=.move_particleb_uni(as.numeric(particles[(ii2*M),tab_unfixed_param]),sd_array,prior)
        			param=particles[(ii2*M),]
        			param[tab_unfixed_param]=param_moved
        			param2=c((seed_count+1),param)
				seed_count=seed_count+1
			}
			else{
				param2[1]=param2[1]+1
				seed_count=seed_count+1
  			}
			list_param[[ic]]=param2
			ic=ic+1
        		tab_param=rbind(tab_param,param2[2:(l+1)])
		}
	}
	list_simul_summarystat=parLapplyLB(cl,list_param,model)
	ic=ic-1
	for (ii in 1:ic){
          tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[ii]]))
      	}
    }
    tab_new_simul=cbind(tab_param,tab_simul_summarystat)
    tab_new_simul2=matrix(0,dim(tab_new_simul)[1],dim(tab_new_simul)[2])
    for (i1 in 1:(dim(tab_new_simul)[1])){
   	for (i2 in 1:(nparam+nstat)){
        	tab_new_simul2[i1,i2]=as.numeric(tab_new_simul[i1,i2])
        }
    }
    iii=0
    for (i in 1:nb_simul){
	if (tab_weight[i]>0){
		iii=iii+1
		n_acc=1
		if (M>1){
			new_dist=.compute_dist_M(M,summary_stat_target,tab_new_simul2[(((iii-1)*M+1):(iii*M)),(nparam+1):(nparam+nstat)],sd_simul)
			n_acc=length(new_dist[new_dist<new_tolerance])
		}
		else{
			new_dist=.compute_dist(summary_stat_target,rbind(tab_new_simul2[iii,(nparam+1):(nparam+nstat)],tab_new_simul2[iii,(nparam+1):(nparam+nstat)]),sd_simul)
			if (new_dist[1]>new_tolerance){
           			n_acc=0
		        }
		}
		MH=min(1,(n_acc/tab_below[i]))
		uuu=runif(1)
		if (uuu<=MH){
          		for (i1 in 1:M){
           			for (i2 in 1:(nparam+nstat)){
              				simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(tab_new_simul2[((iii-1)*M+i1),i2])
            			}
          		}
        	}
		else{
          		for (i1 in 1:M){
            			for (i2 in 1:(nparam+nstat)){
              				simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
            			}
          		}
        	}
	}
	else{
        	for (i1 in 1:M){
          		for (i2 in 1:(nparam+nstat)){
            			simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
          		}
        	}
      	}
    }
    if (M>1){
      particle_dist_mat=.compute_dist_M(M,summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
    }
    else{
      particle_dist_mat=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
    }
    dim(particle_dist_mat)<-c(nb_simul,M)
    tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
    tab_weight2=.replicate_tab(tab_weight,M)
    if (verbose==TRUE){ 
      write.table(as.numeric(new_tolerance),file=paste("tolerance_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
      write.table(cbind(tab_weight2,simul_below_tol),file=paste("output_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
      print(paste("step ",kstep," completed - tol =",new_tolerance,sep=""))
    }
  }
  stopCluster(cl)	
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),weights=tab_weight2/sum(tab_weight2),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs"))) 
}

## function to sample in the prior distributions using a Latin Hypercube sample
###############################################################################
.ABC_rejection_lhs_cluster<-function(model,prior,nb_simul,tab_unfixed_param,seed_count,n_cluster){
  #library(lhs)
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  tab_simul_summarystat=NULL
  tab_param=NULL
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  nparam=length(tab_unfixed_param[tab_unfixed_param])
  l=length(prior)
  random_tab=NULL
  all_unif_prior=.all_unif(prior)
  if (all_unif_prior){
	  random_tab=randomLHS(nb_simul,nparam)
  }
  
  if (npar>0){
   for (irun in 1:npar){
    for (i in 1:(100*n_cluster)){
      param=array(0,l)
      if (!all_unif_prior){
	param=.sample_prior(prior)
      }
      else{
        count=1
        for (j in 1:l){
	  if (tab_unfixed_param[j]){
        	param[j]=as.numeric(prior[[j]][2])+(as.numeric(prior[[j]][3])-as.numeric(prior[[j]][2]))*random_tab[((irun-1)*100*n_cluster+i),count]
		count=count+1
	  }
	  else{
		param[j]=as.numeric(prior[[j]][2])
	  }
        }
      }
      #if (use_seed) { # NB: we force the value use_seed=TRUE
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_cluster
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:(100*n_cluster)){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
   }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      param=array(0,l)
      if (!all_unif_prior){
	param=.sample_prior(prior)
      }
      else{
        count=1
        for (j in 1:l){
	  if (tab_unfixed_param[j]){
        	param[j]=as.numeric(prior[[j]][2])+(as.numeric(prior[[j]][3])-as.numeric(prior[[j]][2]))*random_tab[(npar*100*n_cluster+i),count]
		count=count+1
	  }
	  else{
		param[j]=as.numeric(prior[[j]][2])
	  }
        }
      }
      #if (use_seed) { # NB: we force the value use_seed=TRUE
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  options(scipen=0)
  cbind(tab_param,tab_simul_summarystat)
}



## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniformc_cluster<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,seed_count,inside_prior,n_cluster){
  tab_simul_summarystat=NULL
  tab_param=NULL
  k_acc=0
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  if (npar>0){
    for (irun in 1:npar){
      for (i in 1:(100*n_cluster)){
        l=dim(param_previous_step)[2]
        if (!inside_prior){
          k_acc=k_acc+1
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
        }
        else{
          test=FALSE
          counter=0
          while ((!test)&&(counter<100)){
            counter=counter+1
            k_acc=k_acc+1
            # pick a particle
            param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
            # move it
            param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
            test=.is_included(param_moved,prior)
          }
          if (counter==100){
            stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
          }
        }
        param=param_previous_step[1,]
        param[tab_unfixed_param]=param_moved
        param=c((seed_count+i),param)
        list_param[[i]]=param
        tab_param=rbind(tab_param,param[2:(l+1)])
      }
      seed_count=seed_count+n_cluster
      list_simul_summarystat=parLapplyLB(cl,list_param,model)
      for (i in 1:(100*n_cluster)){
        tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
      }
    }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      l=dim(param_previous_step)[2]
      if (!inside_prior){
        k_acc=k_acc+1
        # pick a particle
        param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
        # move it
        param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
      }
      else{
        test=FALSE
        counter=0
        while ((!test)&&(counter<100)){
          k_acc=k_acc+1
          counter=counter+1
          # pick a particle
          param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
          # move it
          param_moved=.move_particle(as.numeric(param_picked),2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov) # only variable parameters are moved, computation of a WEIGHTED variance
          test=.is_included(param_moved,prior)
        }
        if (counter==100){
          stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
        }
      }
      param=param_previous_step[1,]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  list(cbind(tab_param,tab_simul_summarystat),nb_simul/k_acc)
}

## function to perform ABC simulations from a non-uniform prior and with unidimensional jumps
#############################################################################################
.ABC_launcher_not_uniformc_uni_cluster<-function(model,prior,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,seed_count,inside_prior,n_cluster){
  tab_simul_summarystat=NULL
  tab_param=NULL
  k_acc=0
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  list_param=list(NULL)
  npar=floor(nb_simul/(100*n_cluster))
  n_end=nb_simul-(npar*100*n_cluster)
  l=dim(param_previous_step)[2]
  l_array=dim(param_previous_step[,tab_unfixed_param])[2]
  if (is.null(l_array)){
	l_array=1
  }
  sd_array=array(1,l_array)
  covmat=2*cov.wt(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),as.vector(tab_weight))$cov # computation of a WEIGHTED variance
  for (j in 1:l_array){
    sd_array[j]=sqrt(covmat[j,j])
    
  }
  if (npar>0){
    for (irun in 1:npar){
      for (i in 1:(100*n_cluster)){
        if (!inside_prior){
          k_acc=k_acc+1
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
        }
        else{
          test=FALSE
          counter=0
          while ((!test)&&(counter<100)){
            k_acc=k_acc+1
            counter=counter+1
            # pick a particle
            param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
            # move it
            param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
            test=.is_included(param_moved,prior)
          }
          if (counter==100){
            stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
          }
        }
        param=param_previous_step[1,]
        param[tab_unfixed_param]=param_moved
        param=c((seed_count+i),param)
        list_param[[i]]=param
        tab_param=rbind(tab_param,param[2:(l+1)])
      }
      seed_count=seed_count+n_cluster
      list_simul_summarystat=parLapplyLB(cl,list_param,model)
      for (i in 1:(100*n_cluster)){
        tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
      }
    }
  }
  if (n_end>0){
    #stopCluster(cl)
    #cl <- makeCluster(getOption("cl.cores", 1))
    list_param=list(NULL)
    for (i in 1:n_end){
      if (!inside_prior){
        k_acc=k_acc+1
        # pick a particle
        param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
        # move it
        param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
      }
      else{
        test=FALSE
        counter=0
        while ((!test)&&(counter<100)){
          counter=counter+1
          k_acc=k_acc+1
          # pick a particle
          param_picked=.particle_pick(as.matrix(as.matrix(param_previous_step)[,tab_unfixed_param]),tab_weight)
          # move it
          param_moved=.move_particle_uni(as.numeric(param_picked),sd_array) # only variable parameters are moved
          test=.is_included(param_moved,prior)
        }
        if (counter==100){
          stop("The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution")
        }
      }
      param=param_previous_step[1,]
      param[tab_unfixed_param]=param_moved
      param=c((seed_count+i),param)
      list_param[[i]]=param
      tab_param=rbind(tab_param,param[2:(l+1)])
    }
    seed_count=seed_count+n_end
    list_simul_summarystat=parLapplyLB(cl,list_param,model)
    for (i in 1:n_end){
      tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
    }
    stopCluster(cl)
  }
  else{
    stopCluster(cl)
  }
  list(cbind(tab_param,tab_simul_summarystat),nb_simul/k_acc)
}



## sequential algorithm of Lenormand et al. 2012 
################################################
.ABC_Lenormand_cluster<-function(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,alpha=0.5,p_acc_min=0.05,seed_count=0,inside_prior=TRUE,progress_bar=FALSE){
  ## checking errors in the inputs
  if(!is.vector(alpha)) stop("'alpha' has to be a number.")
  if(length(alpha)>1) stop("'alpha' has to be a number.")
  if (alpha<=0) stop ("'alpha' has to be between 0 and 1.")
  if (alpha>=1) stop ("'alpha' has to be between 0 and 1.")
  if(!is.vector(p_acc_min)) stop("'p_acc_min' has to be a number.")
  if(length(p_acc_min)>1) stop("'p_acc_min' has to be a number.")
  if (p_acc_min<=0) stop ("'p_acc_min' has to be between 0 and 1.")
  if (p_acc_min>=1) stop ("'p_acc_min' has to be between 0 and 1.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  if(!is.logical(inside_prior)) stop("'inside_prior' has to be boolean.")
  
  start = Sys.time()
  if (progress_bar){
	  print("    ------ Lenormand et al. (2012)'s algorithm ------") 
  }
  
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
     if (prior[[i]][1]!="unif"){
	stop("Prior distributions must be uniform to use the Lenormand et al. (2012)'s algorithm.")
     }
  }
  n_alpha=ceiling(nb_simul*alpha)
  
  ## step 1
  # ABC rejection step with LHS
  tab_ini=.ABC_rejection_lhs_cluster(model,prior,nb_simul,tab_unfixed_param,seed_count,n_cluster)
  seed_count=seed_count+nb_simul
  sd_simul=sapply(as.data.frame(tab_ini[,(nparam+1):(nparam+nstat)]),sd) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
  
  # selection of the alpha quantile closest simulations
  simul_below_tol=NULL
  simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,alpha))
  simul_below_tol=simul_below_tol[1:n_alpha,] # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile
  
  # initially, weights are equal
  tab_weight=array(1,n_alpha)
  
  tab_dist=.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)
  tol_next=max(tab_dist)
  if (verbose==TRUE){ 
    write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(as.numeric(tol_next),file="tolerance_step1",row.names=F,col.names=F,quote=F) 
  }
  if (progress_bar){
	  print("step 1 completed")
  }
  
  ## following steps
  p_acc=p_acc_min+1
  nb_simul_step=nb_simul-n_alpha
  it=1
  while (p_acc>p_acc_min){
    it=it+1
    simul_below_tol2=NULL
    tab_inic=.ABC_launcher_not_uniformc_cluster(model,prior,as.matrix(as.matrix(simul_below_tol)[,1:nparam]),tab_unfixed_param,tab_weight/sum(tab_weight),nb_simul_step,seed_count,inside_prior,n_cluster)
    tab_ini=as.matrix(tab_inic[[1]])
    tab_ini=as.numeric(tab_ini)
    dim(tab_ini)=c(nb_simul_step,(nparam+nstat))
    seed_count=seed_count+nb_simul_step
    if (!inside_prior){
      tab_weight2=.compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(as.matrix(simul_below_tol)[,1:nparam])[,tab_unfixed_param]),tab_weight/sum(tab_weight),prior)
    }
    else{
      tab_weight2=tab_inic[[2]]*(.compute_weightb(as.matrix(as.matrix(as.matrix(tab_ini)[,1:nparam])[,tab_unfixed_param]),as.matrix(as.matrix(as.matrix(simul_below_tol)[,1:nparam])[,tab_unfixed_param]),tab_weight/sum(tab_weight),prior))
    }
    simul_below_tol2=rbind(as.matrix(simul_below_tol),as.matrix(tab_ini))
    tab_weight=c(tab_weight,tab_weight2)
    tab_dist2=.compute_dist(summary_stat_target,as.matrix(as.matrix(tab_ini)[,(nparam+1):(nparam+nstat)]),sd_simul)
    p_acc=length(tab_dist2[tab_dist2<=tol_next])/nb_simul_step
    tab_dist=c(tab_dist,tab_dist2)
    tol_next=sort(tab_dist)[n_alpha]
    simul_below_tol2=simul_below_tol2[tab_dist<=tol_next,]
    tab_weight=tab_weight[tab_dist<=tol_next]
    tab_weight=tab_weight[1:n_alpha]
    tab_dist=tab_dist[tab_dist<=tol_next]
    tab_dist=tab_dist[1:n_alpha]
    simul_below_tol=matrix(0,n_alpha,(nparam+nstat))
    for (i1 in 1:n_alpha){
      for (i2 in 1:(nparam+nstat)){
        simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
      }
    }
    if (verbose==TRUE){ 
      write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(seed_count-seed_count_ini),file=paste("n_simul_tot_step",it,sep=""),row.names=F,col.names=F,quote=F)
      write.table(as.numeric(p_acc),file=paste("p_acc_step",it,sep=""),row.names=F,col.names=F,quote=F) 
      write.table(as.numeric(tol_next),file=paste("tolerance_step",it,sep=""),row.names=F,col.names=F,quote=F)
    }
    if (progress_bar){
    		print(paste("step ",it," completed - p_acc = ",p_acc,sep=""))
    }
  }
  list(param=as.matrix(as.matrix(simul_below_tol)[,1:nparam]),stats=as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),weights=tab_weight/sum(tab_weight),stats_normalization=as.numeric(sd_simul),epsilon=max(.compute_dist(summary_stat_target,as.matrix(as.matrix(simul_below_tol)[,(nparam+1):(nparam+nstat)]),sd_simul)),nsim=(seed_count-seed_count_ini),computime=as.numeric(difftime(Sys.time(), start, units="secs"))) 
}


.ABC_sequential_cluster <-
function(method,model,prior,nb_simul,summary_stat_target,n_cluster,use_seed,verbose,...){
    if (use_seed==FALSE){
		stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
    }


	options(scipen=50)
	#library(parallel)

	       ## general function regrouping the different sequential algorithms     
	  ## [Beaumont et al., 2009] Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika,96(4):983â??990.
	  ## [Drovandi & Pettitt 2011] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. Biometrics, 67(1):225â??233.
	  ## [Del Moral et al. 2012] Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation, Statistics and Computing., 22(5):1009-1020.
	  ## [Lenormand et al. 2012] Lenormand, M., Jabot, F., Deffuant G. (2012). Adaptive approximate Bayesian computation for complex models, submitted to Comput. Stat. )
	  return(switch(EXPR = method,
         	 Beaumont = .ABC_PMC_cluster(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,,...),
	         Drovandi = .ABC_Drovandi_cluster(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,...),
	         Delmoral = .ABC_Delmoral_cluster(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,...),
	         Lenormand = .ABC_Lenormand_cluster(model,prior,nb_simul,summary_stat_target,n_cluster,verbose,...)))

	options(scipen=0)

}



## ABC-MCMC algorithm of Marjoram et al. 2003 with automatic determination of the tolerance and proposal range following Wegmann et al. 2009
############################################################################################################################################
.ABC_MCMC2_cluster<-function(model,prior,n_obs,n_between_sampling,summary_stat_target,n_cluster,verbose,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,seed_count=0){
  ## checking errors in the inputs
  if(!is.vector(n_calibration)) stop("'n_calibration' has to be a number.")
  if(length(n_calibration)>1) stop("'n_calibration' has to be a number.")
  if (n_calibration<1) stop ("'n_calibration' has to be positive.")
  n_calibration=floor(n_calibration)
  if(!is.vector(tolerance_quantile)) stop("'tolerance_quantile' has to be a number.")
  if(length(tolerance_quantile)>1) stop("'tolerance_quantile' has to be a number.")
  if (tolerance_quantile<=0) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if (tolerance_quantile>=1) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if(!is.vector(proposal_phi)) stop("'proposal_phi' has to be a number.")
  if(length(proposal_phi)>1) stop("'proposal_phi' has to be a number.")
  if (proposal_phi<=0) stop ("'proposal_phi' has to be positive.")
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  
  start = Sys.time()
  if (verbose){
  	print("    ------ Marjoram et al. (2003)'s algorithm with modifications drawn from Wegmann et al. (2009) related to automatization ------")
  }
  
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  tab_simul_summary_stat=NULL
  tab_param=NULL
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  
  # initial draw of a particle
  initial=.ABC_rejection_internal_cluster(model,prior,n_calibration,seed_count,n_cluster)
  seed_count=seed_count+n_calibration
  tab_param=as.matrix(as.matrix(initial)[,1:nparam])
  tab_simul_summary_stat=as.matrix(initial[,(nparam+1):(nparam+nstat)])
  
  sd_simul=array(0,nstat)
  for (i in 1:nstat){
    sd_simul[i]=sd(tab_simul_summary_stat[,i])
  }
  simuldist=.compute_dist(summary_stat_target,tab_simul_summary_stat,sd_simul)
  ord_sim=order(simuldist,decreasing=F)
  nmax=ceiling(tolerance_quantile*n_calibration)
  dist_max=simuldist[(ord_sim[nmax])]
  tab_param=as.matrix(tab_param)[(ord_sim[1:nmax]),]
  proposal_range=vector(mode = "numeric", length = nparam)
  for (i in 1:nparam){
    proposal_range[i]=sd(as.matrix(tab_param)[,i])*proposal_phi/2
  }
  n_ini=sample(nmax,1)
  tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
  dist_ini=simuldist[(ord_sim[n_ini])]
  param_ini=as.matrix(tab_param)[n_ini,]
  if (verbose==TRUE){ 
    write.table((seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(NULL,file="output_mcmc",row.names=F,col.names=F,quote=F)
    print("initial calibration performed ")
  }

  # chain run
  tab_param=param_ini
  tab_simul_summary_stat=tab_simul_ini
  tab_dist=as.numeric(dist_ini)
  for (is in 2:n_obs){
    for (i in 1:n_between_sampling){
      param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior)
      param=c((seed_count+i),param)
      simul_summary_stat=model(param)
      param=param[2:(nparam+1)]
      dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),sd_simul)
      if (dist_simul<dist_max){
        param_ini=param
        tab_simul_ini=as.numeric(simul_summary_stat)
        dist_ini=dist_simul
      } 
    }
    seed_count=seed_count+n_between_sampling
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
    tab_param=rbind(tab_param,as.numeric(param_ini))
    tab_dist=rbind(tab_dist,as.numeric(dist_ini))
    if (verbose==TRUE){
    	write.table(c(as.numeric(param_ini),tab_simul_ini,as.numeric(dist_ini)),file="output_mcmc",row.names=F,col.names=F,quote=F,append=T)
     	if (is%%100==0){
      		print(paste(is," ",sep=""))
     	}
    }
  }
  tab_param2=matrix(0,dim(tab_param)[1],dim(tab_param)[2])
  for (i in 1:dim(tab_param)[1]){
    for (j in 1:dim(tab_param)[2]){
      tab_param2[i,j]=tab_param[i,j]
    }
  }
  tab_simul_summary_stat2=matrix(0,dim(tab_simul_summary_stat)[1],dim(tab_simul_summary_stat)[2])
  for (i in 1:dim(tab_simul_summary_stat)[1]){
    for (j in 1:dim(tab_simul_summary_stat)[2]){
      tab_simul_summary_stat2[i,j]=tab_simul_summary_stat[i,j]
    }
  }
  tab_dist2=array(0,length(tab_dist))
  for (i in 1:length(tab_dist)){
    tab_dist2[i]=tab_dist[i]
  } 		
  list(param=as.matrix(tab_param2),stats=as.matrix(tab_simul_summary_stat2),dist=tab_dist2,stats_normalization=as.numeric(sd_simul),epsilon=max(tab_dist),nsim=(seed_count-seed_count_ini),n_between_sampling=n_between_sampling,computime=as.numeric(difftime(Sys.time(), start, units="secs"))) 
}


## ABC-MCMC algorithm of Wegmann et al. 2009 - the PLS step is drawn from the manual of ABCtoolbox (figure 9) - NB: for consistency with ABCtoolbox, AM11-12 are not implemented in the algorithm
#################################################################################################################################################################################################
.ABC_MCMC3_cluster<-function(model,prior,n_obs,n_between_sampling,summary_stat_target,n_cluster,verbose,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,numcomp=0,seed_count=0){
  ## checking errors in the inputs
  if(!is.vector(n_calibration)) stop("'n_calibration' has to be a number.")
  if(length(n_calibration)>1) stop("'n_calibration' has to be a number.")
  if (n_calibration<1) stop ("'n_calibration' has to be positive.")
  n_calibration=floor(n_calibration)
  if(!is.vector(tolerance_quantile)) stop("'tolerance_quantile' has to be a number.")
  if(length(tolerance_quantile)>1) stop("'tolerance_quantile' has to be a number.")
  if (tolerance_quantile<=0) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if (tolerance_quantile>=1) stop ("'tolerance_quantile' has to be between 0 and 1.")
  if(!is.vector(proposal_phi)) stop("'proposal_phi' has to be a number.")
  if(length(proposal_phi)>1) stop("'proposal_phi' has to be a number.")
  if (proposal_phi<=0) stop ("'proposal_phi' has to be positive.")
  
  if(!is.vector(numcomp)) stop("'numcomp' has to be a number.")
  if(length(numcomp)>1) stop("'numcomp' has to be a number.")
  if (numcomp<0) stop ("'numcomp' has to be positive.")
  if (numcomp>length(summary_stat_target)) stop ("'numcomp' has to smaller or equal to the number of summary statistics.")
  numcomp=floor(numcomp)
  if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
  if(length(seed_count)>1) stop("'seed_count' has to be a number.")
  if (seed_count<0) stop ("'seed_count' has to be a positive number.")
  seed_count=floor(seed_count)
  
  start = Sys.time()
  if (verbose){
	  print("    ------ Wegmann et al. (2009)'s algorithm ------")
  }
  #library(pls)
  #library(MASS)
  
  ##AM1
  seed_count_ini=seed_count
  nparam=length(prior)
  nstat=length(summary_stat_target)
  if (nstat<=1){
	stop("A single summary statistic is used, use the method 'Marjoram' instead")
  }
  if (numcomp==0){
    numcomp=nstat
  }
  tab_simul_summary_stat=NULL
  tab_param=NULL
  tab_unfixed_param=array(TRUE,nparam)
  for (i in 1:nparam){
     tab_unfixed_param[i]=!((prior[[i]][1]=="unif")&&(as.numeric(prior[[i]][2])==as.numeric(prior[[i]][3])))
  }
  if (length(tab_unfixed_param[tab_unfixed_param])<=1){
	stop("A single parameter is varying, use the method 'Marjoram' instead")
  }
  
  # initial draw of a particle
  initial=.ABC_rejection_internal_cluster(model,prior,nb_simul=n_calibration,seed_count,n_cluster)
  seed_count=seed_count+n_calibration
  tab_param=as.matrix(as.matrix(initial)[,1:nparam])
  tab_simul_summary_stat=as.matrix(initial[,(nparam+1):(nparam+nstat)])
  if (verbose==TRUE){ 
    write.table((seed_count-seed_count_ini),file="n_simul_tot_step1",row.names=F,col.names=F,quote=F)
    write.table(NULL,file="output_mcmc",row.names=F,col.names=F,quote=F)
  }
  
  ## AM2: PLS step
  #print("AM2 ")
  #standardize the params
  sparam=as.matrix(tab_param)[,tab_unfixed_param]
  ls=dim(sparam)[2]
  for(i in 1:ls){
    sparam[,i]=(sparam[,i]-mean(sparam[,i]))/sd(sparam[,i])
  }
  #force stat in [1,2]
  myMax<-c()
  myMin<-c()
  lambda<-c()
  myGM<-c()
  stat=tab_simul_summary_stat
  #print("stat 1 ")
  #print(stat)
  summary_stat_targ=summary_stat_target
  for (i in 1:nstat){
    myMax<-c(myMax,max(stat[,i]))
    myMin<-c(myMin,min(stat[,i]))
    stat[,i]=1+(stat[,i]-myMin[i])/(myMax[i]-myMin[i])
    summary_stat_targ[i]=1+(summary_stat_targ[i]-myMin[i])/(myMax[i]-myMin[i])
  }
  #print("stat 2 ")
  #print(stat)
  #transform statistics via boxcox
  dmat=matrix(0,n_calibration,(ls+1))
  for(i in 1:nstat){
    d=cbind(as.vector(as.numeric(stat[,i])),as.matrix(sparam))
    for (i1 in 1:n_calibration){
      for (i2 in 1:(ls+1)){
        dmat[i1,i2]=as.numeric(d[i1,i2])
      }
    }
    
    save(dmat,file="dmat.RData")
    load("dmat.RData",.GlobalEnv)
    file.remove("dmat.RData")
    
    mylm<-lm(as.formula(as.data.frame(dmat)),data=as.data.frame(dmat))
    #mylm<-lm(as.formula(as.data.frame(dmat)))
    #mylm<-lm(stat[,i]~as.matrix(sparam))
    myboxcox<-boxcox(mylm,lambda=seq(-20,100,1/10),interp=T,eps=1/50,plotit=FALSE)
    lambda<-c(lambda,myboxcox$x[myboxcox$y==max(myboxcox$y)])
    myGM<-c(myGM, exp(mean(log(stat[,i]))))
  }
  #standardize the BC-stat
  myBCMeans<-c()
  myBCSDs<-c()
  for(i in 1:nstat){
    stat[,i]<-((stat[,i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)))
    summary_stat_targ[i]<-((summary_stat_targ[i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)))
    myBCSDs<-c(myBCSDs, sd(stat[,i]))
    myBCMeans<-c(myBCMeans, mean(stat[,i]))
    stat[,i]<-(stat[,i]-myBCMeans[i])/myBCSDs[i]
    summary_stat_targ[i]<-(summary_stat_targ[i]-myBCMeans[i])/myBCSDs[i]
  }
  #perform pls
  myPlsr<-plsr(as.matrix(sparam)~as.matrix(stat), scale=F,ncomp=numcomp,validation='LOO')
  pls_transformation=matrix(0,numcomp,nstat)
  for (i in 1:numcomp){
    pls_transformation[i,]=as.numeric(myPlsr$loadings[,i])
  }
  
  ## AM3
  #print("AM3 ")
  summary_stat_targ=t(pls_transformation %*% as.vector(summary_stat_targ))
  stat_pls=t(pls_transformation %*% t(stat))
  simuldist=.compute_dist(summary_stat_targ,stat_pls,rep(1,numcomp))
  
  ## AM4
  #print("AM4 ")
  ord_sim=order(simuldist,decreasing=F)
  nmax=ceiling(tolerance_quantile*n_calibration)
  dist_max=simuldist[(ord_sim[nmax])]
  
  tab_param=as.matrix(tab_param)[(ord_sim[1:nmax]),]
  proposal_range=vector(mode = "numeric", length = nparam)
  for (i in 1:nparam){
    proposal_range[i]=sd(tab_param[,i])*proposal_phi/2
  }
  if (verbose){
	  print("initial calibration performed ")
  }
  
  ## AM5: chain run
  #print("AM5 ")
  n_ini=sample(nmax,1)
  tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
  param_ini=tab_param[n_ini,]
  tab_param=param_ini
  tab_simul_summary_stat=tab_simul_ini
  dist_ini=simuldist[(ord_sim[n_ini])]
  tab_dist=as.numeric(dist_ini)
  for (is in 2:n_obs){
    for (i in 1:n_between_sampling){
      ## AM6
      #print("AM6 ")
      param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior)
      param=c((seed_count+i),param)
      ## AM7	
      #print("AM7 ")
      simul_summary_stat=model(param)
      param=param[2:(nparam+1)]
      simul_summary_stat_output=simul_summary_stat
      for (ii in 1:nstat){
        simul_summary_stat[ii]=1+(simul_summary_stat[ii]-myMin[ii])/(myMax[ii]-myMin[ii])
      }
      for(ii in 1:nstat){
        simul_summary_stat[ii]<-((simul_summary_stat[ii]^lambda[ii]) - 1)/(lambda[ii]*(myGM[ii]^(lambda[ii]-1)))
        simul_summary_stat[ii]<-(simul_summary_stat[ii]-myBCMeans[ii])/myBCSDs[ii]
      }
      simul_summary_stat=as.matrix(simul_summary_stat)
      dim(simul_summary_stat)<-c(nstat,1)
      simul_summary_stat=t(pls_transformation %*% simul_summary_stat)
      dist_simul=.compute_dist_single(summary_stat_targ,as.numeric(simul_summary_stat),rep(1,numcomp))
      ## AM8-9
      #print("AM8-9 ")
      if (dist_simul<dist_max){
        param_ini=param
        tab_simul_ini=as.numeric(simul_summary_stat_output)
        dist_ini=dist_simul
      }
    }
    seed_count=seed_count+n_between_sampling
    tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
    tab_param=rbind(tab_param,as.numeric(param_ini))
    tab_dist=rbind(tab_dist,as.numeric(dist_ini))
    if (verbose==TRUE){
    	write.table(c(as.numeric(param_ini),tab_simul_ini,as.numeric(dist_ini)),file="output_mcmc",row.names=F,col.names=F,quote=F,append=T)
     	if (is%%100==0){
      		print(paste(is," ",sep=""))
     	}
    }
  }	
  tab_param2=matrix(0,dim(tab_param)[1],dim(tab_param)[2])
  for (i in 1:dim(tab_param)[1]){
    for (j in 1:dim(tab_param)[2]){
      tab_param2[i,j]=tab_param[i,j]
    }
  }
  tab_simul_summary_stat2=matrix(0,dim(tab_simul_summary_stat)[1],dim(tab_simul_summary_stat)[2])
  for (i in 1:dim(tab_simul_summary_stat)[1]){
    for (j in 1:dim(tab_simul_summary_stat)[2]){
      tab_simul_summary_stat2[i,j]=tab_simul_summary_stat[i,j]
    }
  }
  tab_dist2=array(0,length(tab_dist))
  for (i in 1:length(tab_dist)){
    tab_dist2[i]=tab_dist[i]
  }
  list(param=as.matrix(tab_param2),stats=as.matrix(tab_simul_summary_stat2),dist=tab_dist2,epsilon=max(tab_dist),nsim=(seed_count-seed_count_ini),n_between_sampling=n_between_sampling,min_stats=myMin,max_stats=myMax,lambda=lambda,geometric_mean=myGM,boxcox_mean=myBCMeans,boxcox_sd=myBCSDs,pls_transform=pls_transformation,n_component=numcomp,computime=as.numeric(difftime(Sys.time(), start, units="secs"))) 
}

## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al. 2009)
##############################################################################
.ABC_mcmc_cluster <-function(method,model,prior,n_obs,n_between_sampling,summary_stat_target,n_cluster=1,use_seed,verbose,...){
    if (use_seed==FALSE){
		stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
    }


	options(scipen=50)
	#library(parallel)
## Note that we do not consider the original Marjoram's algortithm, which is not prone to parallel computing. (no calibration step)
	    return(switch(EXPR = method,
	       Marjoram = .ABC_MCMC2_cluster(model,prior,n_obs,n_between_sampling,summary_stat_target,n_cluster,verbose,...),
	       Wegmann = .ABC_MCMC3_cluster(model,prior,n_obs,n_between_sampling,summary_stat_target,n_cluster,verbose,...)))

	options(scipen=0)
}



######################### Utilities functions
## Progress Bar
###############

.progressBar <- function (min = 0, max = 1, initial = 0, text = "", char = "=", width = NA, 
                          title, label, style = 1, file = "") 
{
  if (!identical(file, "") && !(inherits(file, "connection") && 
    isOpen(file))) 
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have max > min")
  up <- function(value, text) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste(c("\r  |", rep.int(" ", nw * width + 6)), collapse = ""), file = file)
    cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
                                                    nw * (width - nb)), sprintf("| %3d%% %s", pc,text)), collapse = ""), file = file)
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up(initial, text)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

.updateProgressBar <- function (pb, value, text = "") 
{
  if (!inherits(pb, "txtProgressBar")) 
    stop("'pb' is not from class \"txtProgressBar\"")
  oldval <- pb$getVal()
  pb$up(value,text)
  invisible(oldval)
}







