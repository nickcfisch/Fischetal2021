#############################################################
#Logistics of reading each .dat file, writing reports, etc. 
#############################################################

#Fitting ADMB models by changing .dat file 
require(R2admb)
setup_admb("/apps/gcc/5.2.0/admb/12.0")
#system("admb")

Assess_fit<-function(Operating_Mod, Likelihood, N, verbose, Francis_wt, Output_name, rnd_comp, jfactor){
  for (i in N){
#For supercomputer

############################################################################################################################################################################################
#Writing data file that rerouts to correct model
  if (rnd_comp==FALSE){
   write(x=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/dat_file.dat"), file=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.dat"))
  } else if (rnd_comp==TRUE){
   write(x=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/dat_file_rnd.dat"), file=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.dat"))
  }
#############################################################################################################################################################################################
  if (Francis_wt==TRUE){
   Bio_sample_yr<-readRDS(paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/Bio_sample_yr.rds"))
   FIM_biosample<-readRDS(paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/FIM_biosample.rds"))
   Start_SS<-rowSums(Bio_sample_yr[51:150,])
   Start_survss<-c(rep(0,40),apply(FIM_biosample,1,FUN = sum))
 
   ESS_func<-function(location){
     
     #ESS Calcs 
     results<-dget(paste0(location,"/SimSCAA.rdat"))
     
     #Sams parameterization
     #Function to help get rid of zeros during calculations.  m is for matrix, v is for (Ntilde) vector
     any0 <- function(m, v=NULL){ #get rid of zeros
       rs <- apply(m, 1, sum)
       if(is.null(v)) v <- rep(1, nrow(m))
       out <- rs!=0 & v!=0
       note <- ifelse(any(rs==0), 1, 0)
       return(list(zvec=out, note=note))
     }
     
     # zerof <- function(m){
     # rs <- apply(m, 1, sum)
     # return(rs==0)
     # }
     
     #standardize b/c some props coming out of the assessment aren't exactly 0<p<1
     stdf <- function(m){
       m2 <- sweep(m, 1, apply(m, 1, sum), '/')
       m2[is.nan(m2)] <- NA
       return(m2)
     }
     
     get_v <- function(x, E){
       Ebar <- apply(E, 1, function(z) sum(z*x))
       v=c()
       for(i in 1:nrow(E))
         v[i]=sum(x^2 * E[i,]) - Ebar[i]^2
       return(v)
     }
     
     TA1.8 <- function(Ntilde, E, O, x, type=NULL){ #do something about removing the NA
       # TA1.8 <- function(E, O, x, type=NULL){ #do something about removing the NA
       zeros <- any0(O, Ntilde)
       # zeros <- any0(O)#any0(O, Ntilde)
       Ntildef <- Ntilde[zeros$zvec]
       Ef <- stdf(E[zeros$zvec,])
       Of <- stdf(O[zeros$zvec,])
       v <- get_v(x, Ef)
       Obar <- apply(Of, 1, function(z) sum(z*x))
       Ebar <- apply(Ef, 1, function(z) sum(z*x))
       w <- 1 / var( (Obar - Ebar) / (v / Ntildef)^0.5 )
       # ESS <- 1 / var( (Obar - Ebar) / (v)^0.5 )
       
       wR <- rep(NA, length(zeros$zvec))
       wR[zeros$zvec] <- w
       wR[is.na(wR)] <- 0
       
       # ESSR <- rep(NA, length(zeros$zvec))
       # ESSR[zeros$zvec] <- ESS
       # ESSR[is.na(ESSR)] <- 0
       
       # return(list(ESS=ESSR, note=zeros$note, zvec=zeros$zvec))
       return(list(w=wR, note=zeros$note, zvec=zeros$zvec))
     }
     
     ESS_f<-results$SS_fishery*TA1.8(results$SS_fishery,results$Pred_comp,results$Obs_comp,0:20)$w
     ESS_surv<-results$SS_FIM*TA1.8(results$SS_FIM,results$PredFIM_comp,results$Obs_FIM_comp,0:20)$w
     
     return(list(ESS_F=ESS_f, ESS_surv=ESS_surv))
   }
   
   pin_file<-list(Start_SS, Start_survss, -5, -14.265, rep(0,120), 0.99, 17.32035,-2.995732, -1.203973, -1.38629, -1.84093, 2, 2, 0.4, 2, 2, 2, 0.75, rep(-2,100))
   
   unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
   lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
   
   #Running tpl, must be compiled first
   setwd(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))
   run_admb(fn="SCAA_SnapperModel_SR", verbose = verbose)   #running code  
   warnings()
   #Extracting results
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   
   devs<-c(results$SS_fishery-round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_F),
           results$SS_FIM[41:100]-round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_surv[41:100]))
   counter<-0
   while (mean(devs)>5 & counter < 10){
     #Writing PIN File
     pin_file<-list(round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_F),
                    round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_surv),
                    -5, -14.265, rep(0,120), 0.99, 17.32035, -2.995732, -1.203973, -1.38629, -1.84093, 2, 2, 0.4, 2, 2, 2, 0.75, rep(-2,100))
     
     unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
     lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
     
     #running admb, have to set wd first
     setwd(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))
     run_admb(fn="SCAA_SnapperModel_SR", verbose = verbose)   #running code
     warnings()
     results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
     devs<-c(results$SS_fishery-round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_F),
             results$SS_FIM[41:100]-round(ESS_func(location=paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))$ESS_surv[41:100]))
     counter<-sum(counter,1)
   }
   results$counter<-counter
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))
  } else if (Francis_wt==FALSE) {
################################################################################################################################################################ 
 #Setting working directory to where tpl is, R2admb needs this
  setwd(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood))
#################################################################################################################################################################
   if (Likelihood=="DM_Linear" | Likelihood=="DM_Asymptotic"){
    counter<-0
    warning_index<-1
    while (warning_index==1 & counter < 10){
  
     func<-function(){
     tryCatch(run_admb(fn="SCAA_SnapperModel_SR", verbose = FALSE), warning = function(w) { 
     pin_file<-list(jitter(-5, factor=jfactor) ,jitter(-14.265, factor=jfactor), rep(0,120),0.99, 17.32035,-2.995732,-1.203973, -1.38629, -1.84093, jitter(2, factor=jfactor), jitter(2, factor=jfactor), jitter(-1, factor=jfactor),jitter(-1, factor=jfactor), rep(1,21), 0.4, 2, 2, 2, 0.75, rep(-2,100))
      unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
      lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
      print(w)  })}
      
      func_run<-func()
      if(length(func_run)>1){   #if there is a warning, the length will be 2... i think 
       warning_index<-1       #indicator that keeps the while loop going if there is a warning
      } else {
       warning_index<-0        #warning index 0 stops the while loop 
      }
     counter<-sum(counter,1)
   }
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   results$counter<-counter   #saving the counter to see if it got to the end
#Saving results
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))
#######################################################################################################################################################################################
  } else if (Likelihood== "LN_AR1"){
    counter<-0
    warning_index<-1
    while (warning_index==1 & counter < 10){
  
     func<-function(){
     tryCatch(run_admb(fn="SCAA_SnapperModel_SR", verbose = FALSE), warning = function(w) {
      pin_file<-list(jitter(-5, factor=jfactor),jitter(-14.265, factor=jfactor), rep(0,120), 0.99, 17.32035,
               -2.995732, -1.203973, -1.38629, -1.84093, jitter(2, factor=jfactor), jitter(2, factor=jfactor), jitter(0.5, factor=jfactor), jitter(-0.6931472, factor=jfactor), jitter(0.2, factor=jfactor), jitter(-1, factor=jfactor), 0.4 , 2, 2, 2, 0.75, rep(-2,100))

      unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
      lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
      print(w)  })}
      
      func_run<-func()
      if(length(func_run)>1){   #if there is a warning, the length will be 2... i think 
       warning_index<-1       #indicator that keeps the while loop going if there is a warning
      } else {
       warning_index<-0        #warning index 0 stops the while loop 
      }
     counter<-sum(counter,1)
   }
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   results$counter<-counter   #saving the counter to see if it got to the end
#Saving results
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))
############################################################################################################################################################################################
  } else if (Likelihood=="LN_AR2" | Likelihood=="LN_ARMA"){
    counter<-0
    warning_index<-1
    while (warning_index==1 & counter < 10){
  
     func<-function(){
     tryCatch(run_admb(fn="SCAA_SnapperModel_SR", verbose = FALSE), warning = function(w) {
      pin_file<-list(jitter(-5, factor=jfactor),jitter(-14.265, factor=jfactor), rep(0,120), 0.99, 17.32035,
               -2.995732, -1.203973, -1.38629, -1.84093, jitter(2, factor=jfactor), jitter(2, factor=jfactor), jitter(0.5, factor=jfactor), jitter(0.2, factor=jfactor), jitter(-0.6931472, factor=jfactor), jitter(0.2, factor=jfactor), jitter(0., factor=jfactor), jitter(-1, factor=jfactor), 0.4 , 2, 2, 2, 0.75, rep(-2,100))

      unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
      lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
      print(w)  })}
      
      func_run<-func()
      if(length(func_run)>1){   #if there is a warning, the length will be 2... i think 
       warning_index<-1       #indicator that keeps the while loop going if there is a warning
      } else {
       warning_index<-0        #warning index 0 stops the while loop 
      }
     counter<-sum(counter,1)
   }
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   results$counter<-counter   #saving the counter to see if it got to the end
#Saving results
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))
############################################################################################################################################################################################
  }  else if (Likelihood== "Dirichlet"){
    counter<-0
    warning_index<-1
    while (warning_index==1 & counter < 10){
  
     func<-function(){
     tryCatch(run_admb(fn="SCAA_SnapperModel_SR", verbose = FALSE), warning = function(w) {
      pin_file<-list(jitter(-5, factor=jfactor),jitter(-14.265, factor=jfactor), rep(0,120), 0.99, 17.32035,
               -2.995732, -1.203973, -1.38629, -1.84093, jitter(2, factor=jfactor), jitter(2, factor=jfactor), jitter(0, factor=jfactor), jitter(0, factor=jfactor), 0.4, 2, 2, 2, 0.75, rep(-2,100))

      unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
      lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
      print(w)  })}
      
      func_run<-func()
      if(length(func_run)>1){   #if there is a warning, the length will be 2... i think 
       warning_index<-1       #indicator that keeps the while loop going if there is a warning
      } else {
       warning_index<-0        #warning index 0 stops the while loop 
      }
     counter<-sum(counter,1)
   }
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   results$counter<-counter   #saving the counter to see if it got to the end
#Saving results
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))
############################################################################################################################################################################################
  } else {
    counter<-0
    warning_index<-1
    while (warning_index==1 & counter < 10){
  
     func<-function(){
     tryCatch(run_admb(fn="SCAA_SnapperModel_SR", verbose = FALSE), warning = function(w) {
      pin_file<-list(jitter(-5, factor=jfactor), jitter(-14.265, factor=jfactor), rep(0,120), 0.99, 17.32035,
               -2.995732, -1.203973, -1.38629, -1.84093, jitter(2, factor=jfactor), jitter(2, factor=jfactor), 0.4, 2, 2, 2, 0.75, rep(-2,100))

      unlink(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), recursive = TRUE)
      lapply(pin_file, write, paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SCAA_SnapperModel_SR.pin"), append=TRUE, ncolumns=50)
      print(w)  })}
      
      func_run<-func()
      if(length(func_run)>1){   #if there is a warning, the length will be 2... i think 
       warning_index<-1       #indicator that keeps the while loop going if there is a warning
      } else {
       warning_index<-0        #warning index 0 stops the while loop 
      }
     counter<-sum(counter,1)
   }
   results<-dget(paste0("/blue/edvcamp/nfisch/Spatial_Model/EM_CPUE/",Likelihood,"/SimSCAA.rdat")) 
   results$counter<-counter   #saving the counter to see if it got to the end
#Saving results
   saveRDS(results, file=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",Operating_Mod,i,"/",Likelihood,"_",Output_name,".rds"))

   }
  }
 }
}

#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="Multinomial", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="Multinomial_Robust", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="DM_Linear", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="DM_Asymptotic", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="Dirichlet", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="LN_AR1", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="LN_AR2", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="LN_ARMA", N=1:1000, Francis_wt=FALSE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)

#Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="Multinomial_FrWt", N=1:1000, Francis_wt=TRUE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)
Assess_fit(Operating_Mod="SQ_PE_GE_00062580_FQ_", Likelihood="Multinomial_Robust_FrWt", N=1:1000, Francis_wt=TRUE, verbose=FALSE, Output_name="results", rnd_comp=FALSE, jfactor=10)


warnings()