
###############################################
#This code gets OM info and creates new data
###############################################
rm(list=ls(all=TRUE))

Spatial_Model<-function(save_wd,               #Working directory for saving output and data
                        extract_wd, 
                        seed,                  #random number seed
                        PE,                    #Process Error (turned on or off via T/F)
                        nyear,                 #Number of years in simulation, first 50 years will be unfished 
                        nyear_init,            #Number of unfished years to start simulation    
                        num_ages,              #Number of ages
                        num_cells,
                        num_ports,
                        lam_data,              #Whether or not fish movement cost is based on Data
                        lam_moveCost,          #The distance cost function for fish movement
                        DD_rate,               #DD rate for negative exponential
                        DD_Thresh_quant,       #Quantile of unfished cells used for DD threshold
                        sig_logR,              #Recruitment variation
                        q,                     #Fishery catchability
                        sel_grow,              #Fishery Selectivity (logistic growth rate)
                        sel_midpt,             #Fishery Selectivity (logistic midpoint)
                        lam_Costdist,          #Negative exponential parameter for fishing distance cost function 
                        Abunpref_grate,        #Logistic growth rate of fisher abundance preference
                        eff_scalar,            #Scalar for effort relationship
                        cv_totaleff,           #CV of total effort (puts spread around logistic) 
                        eff_midpt,             #midpoint of effort logistic
                        eff_grate,             #Logistic growth rate for effort logistic
                        cv_harv,               #CV of sampled harvest
                        cv_cpue,               #CV of CPUE
                        cv_effort,             #CV of sampled effort
                        Perc_eff_smpld,        #Percentage of effort sampled
                        Prop_sample,           #Are we proportionally sampling fish from the fishery catches or not
                        Perc_sample_pb,        #Percentage of samples taken per boat 
                        Sample_size_pb,        #Number of fish sampled from each trip
                        Perc_yrs_FIM,          #Percent of years with FIM data
                        FIM_q,                 #FIM q
                        num_FIM_samples,       #right now saying they take 20 trips a year
                        w_dist,                #power weight for fish movement distance 
                        w_depth,               #power weight for depth preference
                        w_substrate,           #power weight for substrate preference
                        w_DD,                  #power weight for DD preference
                        w_cost,                #power weight for fisher cost function
                        w_profit,              #power weight for fisher profit function
                        Super_Comp)            #Is this being run on the super computer?
{
  set.seed(seed)
  #Von-Bert
  Linf<-85.64   #L-Infinity (cm)
  k<-0.19       #Brody Growth Coefficient
  tnot<--0.39   #T-not
  Lt<-Linf*(1-exp(-k*(0:20-tnot)))
  
  #W-L Relationship
  a<-1.7E-5
  b<-3
  Wt<-a*Lt^b

  setwd(paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd))
  N_wSpace_postM<-readRDS("N_wSpace_postM.rds")
  Effort_midpoints<-readRDS("Effort_midpoints.rds")
  Sel_FIM<-readRDS("Sel_FIM.rds")
  Catch_bio_space<-readRDS("Catch_bio_space.rds")
  Catch_obs<-readRDS("Catch_obs.rds")
  CPUE_obs<-readRDS("CPUE_obs.rds")
  FIM_index_data<-readRDS("FIM_index_data.rds")
  FIM_biodata<-readRDS("FIM_biodata.rds")
  FIM_biosample<-readRDS("FIM_biosample.rds")
  Effort_space<-readRDS("Effort_space.rds")
  Catch_numage_space<-readRDS("Catch_numage_space.rds")
  Sel<-readRDS("Selectivity.rds")
  County_Distance<-as.matrix(read.table("/blue/edvcamp/nfisch/Spatial_Model/County_Distance.txt"))
  M_vec<-c(2,1.2,0.19,0.15,0.129,0.115,0.106,0.099,0.095,0.091,0.088,0.086,0.085,0.083,0.082,0.081,0.081,0.08,0.08,0.079,0.078)
  Fec<-data.frame("Age"=0:20, "Fecundity"=c(0,0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6))

  P_Fish<-array(0,dim=c(num_ports,num_cells,nyear))      #Matrix describing the probability of fishing a cell 
  ExploitBiomass_space<-matrix(0, nrow=nyear+1, ncol=num_cells)  #Exploitable Biomass (for fishing preference)
  ExploitBiomass_space[1,]<-colSums(Sel*N_wSpace_postM[1,,]*Wt)     #Exploitable Biomass, vectorized
  for (i in 2:(nyear+1)){ #Pop Loop (i year)     
      for (p in 1:num_ports){ #Each port
        #There is a probability of fishing matrix each year because abundance changes
        P_Fish[p,,i-1]<-exp(-lam_Costdist*County_Distance[p,])^w_cost * (1/(1+exp(-Abunpref_grate*(ExploitBiomass_space[i-1,]-median(ExploitBiomass_space[i-1,])))))^w_profit
        P_Fish[p,,i-1]<-P_Fish[p,,i-1]/sum(P_Fish[p,,i-1]) #Standardizing so each row sums to 1
      }
      ExploitBiomass_space[i,]<-colSums(Sel*N_wSpace_postM[i,,]*Wt)     #Exploitable Biomass
    }
  #############################
  #Getting Data for dat file
  #############################
  
  #Composition sampling
  #Effort_sample<-Effort_midpoints*Perc_eff_smpld
  Effort_sample<-matrix(0, nrow=nyear,ncol=num_ports)
  for (i in (nyear_init+1):nyear){
    if(sum(Effort_space[i,])>0){
     Effort_sample[i,]<-rmultinom(n=1, size=sum(Effort_space[i,])*Perc_eff_smpld, prob=Effort_midpoints[i,])  #adding in error to which ports are sampled
    }
  }
  
  #Choosing which cells to sample (matrix is filled in with cell indicator to be sampled with 1 unit of effort)
  Cells_sampled<-array(0, dim=c(nyear, num_ports, ceiling(max(Effort_sample))))
  #Bio_Samples contains year, ports, sample, and ages
  Bio_samples<-array(NA,dim=c(nyear,num_ports,ceiling(max(Effort_sample)),num_ages))
  #Sampling the Fishery
  for (i in 1:nyear){
    for (p in 1:num_ports){
      if (Effort_sample[i,p]!=0){  #If there are actual units of effort sampled from that port, then which cells do they go to
        Cells_sampled[i,p,1:ceiling(Effort_sample[i,p])]<-sample(1:num_cells,size=ceiling(Effort_sample[i,p]),prob=P_Fish[p,,i], replace=TRUE)
        for (k in 1:ceiling(Effort_sample[i,p])){
          if (Prop_sample==TRUE){
            if (Effort_space[i,Cells_sampled[i,p,k]]>0){ #If the effort going into that cell was higher than zero, can sample it
              #The number of samples taken number of fish caught per unit of effort that went to that cell, times the proportion sampled per unit of effort  
              if (sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])>0){ #if there is no catch then can't sample
                Bio_samples[i,p,k,]<-rmultinom(n=1, size=round(sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])/Effort_space[i,Cells_sampled[i,p,k]]*Perc_sample_pb), prob=Catch_numage_space[i,,Cells_sampled[i,p,k]]/sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]))
              }
            }
          }
        }
      }
    }
  }
  
  #Total sample for each year
  Bio_sample_yr<-apply(X=Bio_samples, MARGIN=c(1,4), FUN=sum, na.rm=T)
  
  #Now getting the data that would feed into assessment, the above but in composition form
  Age_Comp_Data<-Bio_sample_yr/rowSums(Bio_sample_yr)    #Data
  Age_Comp_Data<-ifelse(Age_Comp_Data==0,Age_Comp_Data+1E-5,Age_Comp_Data)  #Suppressing zeroes by adding a small constant
  Age_Comp_Data<-Age_Comp_Data/rowSums(Age_Comp_Data) #renormalizing
    
  #RANDOMLY sampling catch at age in each year
  catch_age<-rowSums(Catch_numage_space, dims=2)
  catch_age_smpld<-matrix(NA, nrow=100, ncol=21)
  for (j in 51:150){
    catch_age_smpld[j-50,]<-rmultinom(n=1, size=rowSums(Bio_sample_yr)[j], prob=catch_age[j,]) 
  }
  rnd_age_comp<-catch_age_smpld/ifelse(rowSums(catch_age_smpld)==0,1,rowSums(catch_age_smpld))
  rnd_age_comp<-ifelse(rnd_age_comp==0,rnd_age_comp+1E-5,rnd_age_comp)  #Suppressing zeroes by adding a small constant
  rnd_age_comp<-rnd_age_comp/rowSums(rnd_age_comp) #renormalizing
  
  #############################################
  #Adding FIM index CPUE and Catch Composition
  #############################################
  
  #FIM CPUE (only for latter 3/5ths of time series)
  FIM_years<-(nyear-nyear_init)*Perc_yrs_FIM
      
  #RANDOMLY sampling survey in each year
  surv_cage<-t(t(rowSums(N_wSpace_postM, dims=2))*Sel[,1])
  surv_cage_smpl<-matrix(NA, nrow=60, ncol=21)
  for (k in 91:150){
    surv_cage_smpl[k-90,]<-rmultinom(n=1, size=apply(FIM_biosample,1,FUN = sum)[k-90], prob=surv_cage[k,])
  }
  rnd_surv_comp<-surv_cage_smpl/rowSums(surv_cage_smpl)
  rnd_surv_comp<-ifelse(rnd_surv_comp==0,rnd_surv_comp+1E-5,rnd_surv_comp)  #Suppressing zeroes by adding a small constant
  rnd_surv_comp<-rnd_surv_comp/rowSums(rnd_surv_comp) #renormalizing
  
  dir.create(save_wd)
  
  write(c(seed,PE,nyear,nyear_init,num_ages,lam_data,lam_moveCost,DD_rate,DD_Thresh_quant,sig_logR,q, sel_grow, sel_midpt, lam_Costdist,
          Abunpref_grate, eff_scalar, eff_midpt, eff_grate, cv_harv, cv_cpue, cv_effort, Perc_eff_smpld, Prop_sample, Perc_sample_pb,
          Sample_size_pb, Perc_yrs_FIM,FIM_q, num_FIM_samples, w_dist, w_depth, w_DD, w_cost, w_profit),
        file = paste0(save_wd, "/Parameters.txt"), ncolumns=32)

  #Writing Dat File
    dat_file<-list(1, nyear-nyear_init, 0, num_ages-1, M_vec, 
                 Catch_obs[(nyear_init+1):nyear],                 #Observed Catch
                 CPUE_obs,                                        #Observed CPUE
                 t(Age_Comp_Data[(nyear_init+1):nyear,]),             #Fishery Age Composition data
                 rowSums(Bio_sample_yr[(nyear_init+1):nyear,]),    #SS for Fishery composition
                 c(rep(0,nyear-nyear_init-(nyear-nyear_init)*Perc_yrs_FIM),FIM_index_data),        #Survey Index
                 t(rbind(matrix(0,nrow=(nyear-nyear_init)*(1-Perc_yrs_FIM), ncol=num_ages),FIM_biodata)),    #Survey Composition
                 c(rep(0,nyear-nyear_init-(nyear-nyear_init)*Perc_yrs_FIM),apply(FIM_biosample,1,FUN = sum)),  #Sample size of FIM Composition
                 Fec$Fecundity,                                    #Fecundity
                 Wt,                                               #Wt-at-age
                 c(11,22,33))                                      #test Vector

  unlink(paste0(save_wd,"/dat_file.dat"), recursive = TRUE)
  lapply(dat_file, write, paste0(save_wd,"/dat_file.dat"), append=TRUE, ncolumns=21)
  
  dat_file_rnd<-list(1, nyear-nyear_init, 0, num_ages-1, M_vec, 
                Catch_obs[(nyear_init+1):nyear],                 #Observed Catch
                CPUE_obs,                                        #Observed CPUE
                t(rnd_age_comp[1:100,]),             #Fishery Age Composition data
                rowSums(Bio_sample_yr[51:150,]),    #SS for Fishery composition
                c(rep(0,40),FIM_index_data),        #Survey Index
                t(rbind(matrix(0,nrow=40, ncol=21),rnd_surv_comp)),    #Survey Composition
                c(rep(0,40),apply(FIM_biosample,1,FUN = sum)),  #Sample size of FIM Composition
                Fec$Fecundity,                                    #Fecundity
                Wt,                                               #Wt-at-age
                c(11,22,33))                                      #test Vector
                 
  unlink(paste0(save_wd,"/dat_file_rnd.dat"), recursive = TRUE)
  lapply(dat_file_rnd, write, paste0(save_wd,"/dat_file_rnd.dat"), append=TRUE, ncolumns=21)
  
  #Saving Ouput from model
  #Operating Model
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/N_wSpace_postM.rds" ),to=paste0(save_wd, "/N_wSpace_postM.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/SSB_space.rds" ),to=paste0(save_wd, "/SSB_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Catch_bio_space.rds" ),to=paste0(save_wd, "/Catch_bio_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Catch_num_space.rds" ),to=paste0(save_wd, "/Catch_num_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Catch_numage_space.rds" ),to=paste0(save_wd, "/Catch_numage_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/F_space.rds" ),to=paste0(save_wd, "/F_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Effort_space.rds" ),to=paste0(save_wd, "/Effort_space.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Selectivity.rds" ),to=paste0(save_wd, "/Selectivity.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/Sel_FIM.rds" ),to=paste0(save_wd, "/Sel_FIM.rds" ))
  file.copy(from=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/",extract_wd,"/FIM_sample.rds" ),to=paste0(save_wd, "/FIM_sample.rds" ))
  
  #Data section below
  saveRDS(Catch_obs, file=paste0(save_wd, "/Catch_obs.rds" ))
  saveRDS(CPUE_obs, file=paste0(save_wd, "/CPUE_obs.rds" ))
  saveRDS(Cells_sampled, file=paste0(save_wd, "/Cells_sampled.rds" ))
  saveRDS(Bio_samples, file=paste0(save_wd, "/Bio_samples.rds" ))
  saveRDS(Bio_sample_yr, file=paste0(save_wd, "/Bio_sample_yr.rds" ))
  saveRDS(Age_Comp_Data, file=paste0(save_wd, "/Age_Comp_Data.rds" ))
  saveRDS(FIM_biosample, file=paste0(save_wd, "/FIM_biosample.rds" ))
  saveRDS(FIM_index_data, file=paste0(save_wd, "/FIM_index_data.rds" ))
  saveRDS(FIM_biodata, file=paste0(save_wd, "/FIM_biodata.rds" ))
}

#Status Quo Model
for (i in 1:1000){
Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/SQ_PE_GE_00062580_FQ_",i), extract_wd=paste0("SQ_PE_GE_0210_FQ_",i),
              seed=i,PE=TRUE,nyear=150,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0.03, Abunpref_grate=0.00025, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.000625, Prop_sample=TRUE, Perc_sample_pb=0.80, Sample_size_pb=20, Perc_yrs_FIM=0.6,
              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE)          
}

#RF Model
#for (i in 1:1000){
#Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/RF_PE_GE_00062580_FQ_",i), extract_wd=paste0("RF_PE_GE_0210_FQ_",i),
#              seed=i,PE=TRUE,nyear=150,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.000625, Prop_sample=TRUE, Perc_sample_pb=0.80, Sample_size_pb=20, Perc_yrs_FIM=0.6,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE)          
#}

#Hybrid Model
#for (i in 1:1000){
#Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/RF_Hybrid_GE_1020_FQ_",i), extract_wd=paste0("RF_Hybrid_GE_0210_FQ_",i),
#              seed=i,PE="Hybrid",nyear=150,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.10, Prop_sample=TRUE, Perc_sample_pb=0.20, Sample_size_pb=20, Perc_yrs_FIM=0.6,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE)          
#}
warnings()
