R# library nécessaire : mnormt


rm(list = ls())
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}



# pour eviter les effets de bords
sourceDir('R')



# brute-force ABC
# .ABC_PMC2(.binary_model("./parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix


initAndClear <- function(path,...){
model=trait_model
sourceDir(path)
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix

}

###### END INIT

tmps.ini=Sys.time()         #top chrono
initAndClear('R')
ABC_rejection(model,prior_matrix,10,TRUE)
cat("temps mis : ",difftime(Sys.time(),tmps.ini,units="sec")," s\n") #top fin chrono


# sequential algorithm
# VALID : OK
# .ABC_PMC(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
#initial version
.ABC_PMC(model,prior_matrix,20,c(50,2.5),use_seed=TRUE,inside_prior=TRUE,tolerance_tab=c(0.8,0.6,0.4))
# easyABC version
initAndClear('R')
ABC_sequential('Beaumont',model,prior_matrix,20,c(50,2.5),tolerance_tab=c(0.8,0.6,0.4))

#.ABC_Drovandi<-function (model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha=0.5,c=0.01,first_tolerance_level_auto=TRUE,use_seed=TRUE,seed_count=0,...)
#.ABC_Drovandi(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha,c,first_tolerance_level_auto,use_seed,seed_count,tolerance_tab=1.0)
# easyABC version
initAndClear('R')
ABC_sequential('Drovandi',model,prior_matrix,20,c(50,2.5),alpha=0.5,c=0.01,tolerance_tab=1.0)

# VALID : NO
#Delmoral = .ABC_Delmoral(model,prior_matrix,nb_simul,alpha,M,nb_threshold,tolerance_target,summary_stat_target,use_seed,seed_count),
#.ABC_Delmoral(.binary_model("./parthy"),prior_matrix,20,0.5,1,5,c(50,2.5),tolerance_target=0.8)

initAndClear('R')
.ABC_Delmoral(model,prior_matrix,10,c(50,2.5),0.5,1,5,5)

initAndClear('R')
ABC_sequential('Delmoral',model,prior_matrix,10,c(50,2.5),0.5,1,5,5)
# ABC_sequential('Delmoral',model,prior_matrix,20,alpha=0.5,M=1,nb_threshold=c(1,0.8),tolerance_target=c(50,2.5))
#[1] "    ------ Delmoral algoritm ------"
#Erreur dans length(summary_stat_target) : 'summary_stat_target' est manquant

#Lenormand = .ABC_Lenormand
initAndClear('R')
.ABC_Lenormand(model,prior_matrix,10,c(50,2.5),p_acc_min=0.4)

initAndClear('R')
ABC_sequential('Lenormand',model,prior_matrix,10,c(50,2.5),p_acc_min=0.4)


#### Validation of MCMC algorithm
initAndClear('R')
.ABC_MCMC(model,prior_matrix,10,1,c(50,2.5),8,c(50,1,20,10000),c(0,1,0.5,0,50,1),use_seed=TRUE,seed_count=0)

