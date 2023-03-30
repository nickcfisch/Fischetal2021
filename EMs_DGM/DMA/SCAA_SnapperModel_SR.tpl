//Red Snapper Assessment Model Fit to Data Generated from Operating Model
//Dirichlet-Multinomial Saturating Parameterization

DATA_SECTION
  init_adstring datafile
!!ad_comm::change_datafile_name(datafile);

  init_int fyear                                //first year
  init_int lyear                                //last year
  init_int fage                                 //first age
  init_int lage                                 //last age
  init_vector M_vec(fage,lage)                 //Natural mortality vector at age from assessment
  init_vector obs_harv(fyear,lyear)	//observed total fishery harvest
  init_vector obs_fishery_cpue(fyear,lyear)	//observed fishery catch-per-unit-effort
  init_matrix obs_fishery_comp(fyear,lyear,fage,lage) //Age composition from the fishery catch
  init_vector SS_fishery(fyear,lyear)	//sample size for biodata

  init_vector obs_FIM_CPUE(fyear,lyear)                   //FIM CPUE
  init_matrix obs_FIM_comp(fyear,lyear,fage,lage)     //FIM composition
  init_vector SS_FIM(fyear,lyear)                 //FIM sample size    

  init_vector Fecund_aa(fage,lage)   //Fecundity at-age
  init_vector Wtage(fage,lage)      //Mean Weight at-age
  
  init_vector test(1,3)				//test vector for ensuring data are read in correctly
  int i                                         //index for year loop
  int j                                         //index for age loop

  matrix M(fyear,lyear,fage,lage)               //Natural mortality matrix

//  !!cout << test << endl;
//  !!exit(88);

INITIALIZATION_SECTION
  //Have it in a pin file

PARAMETER_SECTION
  //Declaration of parameters to be estimated
  init_bounded_number log_q(-20.,1.,1)            //Fishery Catchability 
  init_bounded_number log_q_FIM(-20.,1.,1)            //FIM Catchability 

  init_bounded_dev_vector log_recruit_devs(fage,lyear+lage-1,-10.,10.,1)  //Recruitment Deviations
  init_bounded_number steepness(0.,1.,-2)             //Steepness for recruitment
  init_bounded_number log_R0_FLA(10.,25.,1)          //Unfished recruitment for Florida cells

  init_bounded_number log_cv_fishery(-5.,2.,-1)           //log cv for fishery catch
  init_bounded_number log_sigma_rec(-5.,2.,-1)            //log sd for recruitment
  init_bounded_number log_cv_fishery_CPUE(-5.,2.,-1)          //log cv for Fishery CPUE
  init_bounded_number log_cv_FIM_CPUE(-5.,2.,1)          //log cv for FIM CPUE

  init_bounded_number fishery_sellogis_k(-2.,5.,2)        //logistic selectivity parameter 1 for Fishery 
  init_bounded_number fishery_sellogis_midpt(0.,20.,2)    //logistic selectivity parameter 2 for Fishery

  init_bounded_number FIM_sellogis_k(-2.,5.,2)            //logistic selectivity parameter 1 for Survey 
  init_bounded_number FIM_sellogis_midpt(0.,20.,2)        //logistic selectivity parameter 2 for Survey 

  init_bounded_number ltheta(-10.,20.,1)            //theta param for D-M
  init_bounded_number ltheta_FIM(-10.,20.,1)        //theta param for D-M

  init_vector ones(fage,lage,-1)

  init_bounded_vector log_fint(fyear,lyear,-20.,0.,1)   //Fishing intensity for each year

 //Derived variables
  vector fishery_sel(fage,lage)                  //Fishery selectivity
  vector FIM_sel(fage,lage)                      //FIM selectivity

  matrix F(fyear,lyear,fage,lage)                //Instantaneous fishing mortality matrix
  matrix Z(fyear,lyear,fage,lage)                //Total mortality matrix
  matrix S(fyear,lyear,fage,lage)                //Annual survival matrix

  matrix N(fyear,lyear+1,fage,lage)	         //Predicted abundance at age, +1 is to incorporate final catch 
  matrix Nbar(fyear,lyear,fage,lage)	         //Nbar

  matrix pred_caa(fyear,lyear,fage,lage)         //Predicted fishery catch at age
  matrix pred_fishery_comp(fyear,lyear,fage,lage)//Predicted age composition from the fishery catch
  vector pred_harv(fyear,lyear)                 //Predicted Harvest weight (kg) 
  vector pred_totcatch_fishery(fyear,lyear)      //Predicted Fishery total catch, used as denominator to get age composition

  matrix pred_FIM_caa(fyear,lyear,fage,lage)     //Predicted FIM catch at age
  matrix pred_FIM_comp(fyear,lyear,fage,lage)    //Predicted age composition from the survey
  vector pred_totcatch_FIM(fyear,lyear)          //Predicted FIM total catch, used as denominator to get age composition
  vector pred_fishery_cpue(fyear,lyear)          //Predicted fishery cpue

  vector log_rec_devs(fage,lyear+lage)          //rec devs vector that adds in the current year

  vector spbiomass(fyear,lyear+1)      //Spawning biomass
  vector N0_FLA_age(fage,lage)               //Unfished numbers at age
  vector lxo(fage,lage)               //Unfished numbers at age
  number SSB0_FLA                             //Unfished spawning biomass
  
///Sigmas
  number cv_fishery                   //sigma for fishery catch
  number sigma_rec                       //sigma for recruit deviations
  number cv_FIM_CPUE                  //sigma for survey CPUE
  number cv_fishery_CPUE                  //cv for fishery CPUE

  vector Beta(fyear,lyear)
  number theta
  number theta_FIM
  vector Beta_FIM(fyear,lyear)

//Likelihood components
  number L1
  number L2
  number L3
  number L4
  number L5
  number L6
  objective_function_value NLL;

  //!!cout << FIM_sellogis_midpt << endl;
  //!!exit(88);

PRELIMINARY_CALCS_SECTION

PROCEDURE_SECTION
  get_selectivity();
  get_mortality();
  get_population();
  get_catch();
  get_FIM();
  get_objective();

FUNCTION get_selectivity
//FISHERY SELECTIVITY
//  /*
  for (j=fage;j<=lage;j++){
  fishery_sel(j)=1/(1+mfexp(-fishery_sellogis_k*(int(j)-fishery_sellogis_midpt)));  //fishery selectivity as logistic function
  }
  fishery_sel=fishery_sel/max(fishery_sel);
//  */
  
  //cout<< fishery_sel <<endl;
 // exit(1);

FUNCTION get_mortality
  //MORTALITY
  for(i=fyear;i<=lyear;i++)
   {
    for(j=fage;j<=lage;j++)
     {
     F(i,j) = fishery_sel(j)*mfexp(log_fint(i));  //setting year and age specific fishing mortality as the product of selectivity and year specific fishing intensity
     }
   }
   
//Filling in natural mortality matrix
  for(i=fyear;i<=lyear;i++){
   M(i)=M_vec;
  }

  Z = M + F;             //Total instantaneous mortality
  S = mfexp(-1.0*Z);     //Annual survival

//  cout<< fishery_sel <<endl;
//  cout<< F(1) <<endl;
//  cout<< M(1) <<endl;
//  exit(1);
FUNCTION get_population
  //Unfished Spawning biomass calcs
  lxo(fage)=1;
  for(j=fage+1;j<=lage;j++){
   lxo(j)=lxo(j-1)*mfexp(-M_vec(j-1));   //cumulative survival
  }
  lxo(lage)=lxo(lage)/(1-mfexp(-M_vec(lage)));
  N0_FLA_age=mfexp(log_R0_FLA)*lxo;
  SSB0_FLA=sum(elem_prod(N0_FLA_age,Fecund_aa));   //Numbers at age * Fecundity

  //Filling in log_rec_devs
  log_rec_devs(fage,lyear+lage-1)=log_recruit_devs(fage,lyear+lage-1);
  log_rec_devs(lyear+lage)=0;

  N.initialize(); 
  //Abundance at age in the first year
  N(fyear)(fage,lage)=mfexp(log_R0_FLA)*elem_prod(mfexp(log_recruit_devs(fage,lage)),lxo(fage,lage)); //set abundance in the fist year for age 2 and older fish as unfished * natural mortality 

  //Spawning Biomass in the first year
  spbiomass(fyear)=sum(elem_prod(N(fyear),Fecund_aa));  // Getting spawning biomass for the first year (acounting for rec dev, unlike SSB0_FLA)

//Population loop
  for(i=fyear+1;i<=lyear+1;i++)
  {
    for(j=fage+1;j<=lage;j++)
     {
      N(i,j)=N(i-1,j-1)*S(i-1,j-1);
     }
   N(i,lage)+=N(i-1,lage)*S(i-1,lage);   //Plus group
   spbiomass(i)=sum(elem_prod(N(i),Fecund_aa)); //Spawning Biomass
   N(i,fage)=((4.*steepness*mfexp(log_R0_FLA)*spbiomass(i))/(SSB0_FLA*(1.-steepness)+spbiomass(i)*(5.*steepness-1.)))*mfexp(log_rec_devs(i+lage-1));  //Recruitment
  }

//  cout<< spbiomass <<endl;
//  cout<< N(fyear) <<endl;
//  cout<<N0_FLA_age<<endl;
//  exit(1);
FUNCTION get_catch
//Catch
  pred_caa=elem_prod(elem_div(F,Z),elem_prod(1.0-S,N));	//Baranov catch equation for predicting catch at age

  for(i=fyear;i<=lyear;i++){
   pred_totcatch_fishery(i)=sum(pred_caa(i));          // total predicted catch by year
   pred_harv(i)=sum(elem_prod(pred_caa(i),Wtage));  //Predicted total harvest weight each year
   pred_fishery_comp(i)=pred_caa(i)/(pred_totcatch_fishery(i)+0.0001);  //calculating predicted catch age composition
   Nbar(i)=elem_prod(elem_prod(N(i),(1.0-mfexp(-1.0*Z(i)))),1.0/Z(i));
   pred_fishery_cpue(i)=exp(log_q)*sum(elem_prod(elem_prod(Nbar(i),Wtage),fishery_sel));
  }

  //cout<<gammln(pred_fishery_comp(5)+3)<<endl;
  //exit(1);
FUNCTION get_FIM
 //FIM Survey
  for (j=fage;j<=lage;j++){
  FIM_sel(j)=1/(1+mfexp(-FIM_sellogis_k*(int(j)-FIM_sellogis_midpt)));  //Logistic Selectivity
  }
  FIM_sel=FIM_sel/max(FIM_sel);

  for (i=fyear;i<=lyear;i++){
  pred_FIM_caa(i)=mfexp(log_q_FIM)*elem_prod(FIM_sel,N(i));	//Predicted survey catch-at-age each year 
  pred_totcatch_FIM(i)=sum(pred_FIM_caa(i));          // total predicted survey cpue
  pred_FIM_comp(i)=pred_FIM_caa(i)/(pred_totcatch_FIM(i)+0.0001);  //calculating predicted survey catch age composition
  }

  //cout<<pred_FIM_caa(50)<<endl;
  //exit(1);

FUNCTION get_objective  //working on sigmas, look into variance ratios

  cv_fishery=mfexp(log_cv_fishery);    //sd for fishery harvest, sometimes treated as CV
  cv_FIM_CPUE=mfexp(log_cv_FIM_CPUE);  //Survey CPUE sd
  cv_fishery_CPUE=mfexp(log_cv_fishery_CPUE);  //Fishery CPUE sd
  sigma_rec=mfexp(log_sigma_rec);   //sd for recruitments

//LIKELIHOODS
  L1.initialize();
  // Normal by CV
  for(i=fyear;i<=lyear;i++){
  L1 += log(cv_fishery*pred_harv(i)+0.0001)+0.5*square((obs_harv(i)-pred_harv(i))/(cv_fishery*pred_harv(i)+0.0001));   //Normal likelihood using CV
  } 

//Dirichlet-Multinomial
  theta=mfexp(ltheta);
  Beta=theta*SS_fishery;
  L2.initialize();
  for(i=fyear;i<=lyear;i++){
  //Saturating Version from Thorson
  if(SS_fishery(i)>0){
  L2 += -1*(gammln(SS_fishery(i)+1)-sum(gammln(SS_fishery(i)*obs_fishery_comp(i)+ones))+gammln(theta)-gammln(SS_fishery(i)+theta)+sum(gammln(SS_fishery(i)*obs_fishery_comp(i)+theta*pred_fishery_comp(i))-gammln(theta*pred_fishery_comp(i)))); //Simplified
  }}

//Fishery CPUE
  L3.initialize();
  for(i=fyear;i<=lyear;i++){
  L3 += log(cv_fishery_CPUE*pred_fishery_cpue(i)+0.0001)+0.5*square((obs_fishery_cpue(i)-pred_fishery_cpue(i))/(cv_fishery_CPUE*pred_fishery_cpue(i)+0.0001)); //Likelihood for Fishery CPUE
  }

  theta_FIM=mfexp(ltheta_FIM);
  Beta_FIM=theta_FIM*SS_FIM;
  L4.initialize();
  L5.initialize();
  for(i=lyear*2/5+fyear;i<=lyear;i++){
  //Thorson Saturating
  if(SS_FIM(i)>0){
  L4 += -1*(gammln(SS_FIM(i)+1)-sum(gammln(SS_FIM(i)*obs_FIM_comp(i)+ones))+gammln(theta_FIM)-gammln(SS_FIM(i)+theta_FIM)+sum(gammln(SS_FIM(i)*obs_FIM_comp(i)+theta_FIM*pred_FIM_comp(i))-gammln(theta_FIM*pred_FIM_comp(i)))); //Simplified
  }
  L5 += log(cv_FIM_CPUE*pred_totcatch_FIM(i)+0.0001)+0.5*square((obs_FIM_CPUE(i)-pred_totcatch_FIM(i))/(cv_FIM_CPUE*pred_totcatch_FIM(i)+0.0001)); //Likelihood for FIM CPUE
  }

  L6 = double(size_count(log_recruit_devs))*log_sigma_rec+(1.0/(2.0*square(sigma_rec)))*norm2(log_recruit_devs); //Recruitment deviations


  NLL=L1+L2+L3+L4+L5+L6;

  /*
  cout<<NLL<<endl;
  cout<<L1<<endl;
  cout<<pred_harv<<endl;
  cout<<obs_harv<<endl;
  cout<<L2<<endl;
  cout<<L3<<endl;
  cout<<L4<<endl;
  cout<<L5<<endl;
  cout<<obs_FIM_CPUE<<endl;
  cout<<pred_totcatch_FIM<<endl;
  cout<<L6<<endl;
  exit(1);
  */

GLOBALS_SECTION
  #include "admodel.h"
  #include "admb2r.cpp"

REPORT_SECTION
  report<<"Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
  report << " " << endl;
  report << "Total Abundance" << endl;
  for (i=fyear;i<=lyear;i++)
  {  report << sum(N(i)) << endl;
  }
  report << "  " << endl;
  report << "Spawning Biomass" << endl;
  report << spbiomass << endl;
  report << " " << endl;
  report << "Recruitment" << endl;
  for (i=fyear;i<=lyear;i++)
  {  report << N(i,fage) << endl;
  }
  report << "   "  << endl;
  report << "Recruitment deviations" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << mfexp(log_recruit_devs(i)) << endl;
  }
  report << " " << endl;
  report << "Abundance at Age" << endl;
  report << N << endl;
  report << " " << endl;
  report << "Fishing Mortality" << endl;
  report << F << endl;
  report << " " << endl;
  report << "Natural Mortality" << endl;
  report << M << endl;
  report << " " << endl;
  report << "Fishery Selectivity" << endl;
  report << fishery_sel << endl;
  report << " " << endl;
  report << "Survey Selectivity" << endl;
  report << FIM_sel << endl;
  report << " " << endl;
  report << "Survey catchability" << endl;
  report << mfexp(log_q_FIM) << endl;
  report << " " << endl;
  report << "Observed Harvest" << " " << "Predicted Harvest" << endl;
  for (i=fyear;i<=lyear;i++)
      {
          report << obs_harv(i) << "                 " << pred_harv(i) << endl;
      }
  report << " " << endl;
  report << "Fishery Catch CV" << endl;  
  report << mfexp(log_cv_fishery) << endl;
  report << " " << endl;
  report << "FIM CPUE CV" << endl;  
  report << mfexp(log_cv_FIM_CPUE) << endl;
  report << " " << endl;
  report << "Observed Survey CPUE" << "     " << "Predicted Survey CPUE" << endl;
  for (i=fyear;i<=lyear;i++)
      {
          report << obs_FIM_CPUE(i) << "                 " << pred_totcatch_FIM(i) << endl;
      }
  report << " " << endl;

//Writing results into R
  open_r_file("SimSCAA.rdat");
   open_r_vector("Gradient_NLL");
    wrt_r_item("Gradient",objective_function_value::pobjfun->gmax);
    wrt_r_item("NLL",NLL);
    wrt_r_item("L1",L1);
    wrt_r_item("L2",L2);
    wrt_r_item("L3",L3);
    wrt_r_item("L4",L4);
    wrt_r_item("L5",L5);
    wrt_r_item("L6",L6);
   close_r_vector();
    wrt_r_complete_vector("Predharv", pred_harv);
    wrt_r_complete_vector("Obsharv", obs_harv);
    wrt_r_complete_vector("Pred_CPUE", pred_fishery_cpue);
    wrt_r_complete_vector("obs_CPUE", obs_fishery_cpue);
    wrt_r_complete_vector("PredFIM_CPUE", pred_totcatch_FIM);
    wrt_r_complete_vector("obsFIM_CPUE", obs_FIM_CPUE);
    open_r_matrix("Obs_comp");
      wrt_r_matrix(obs_fishery_comp);
    close_r_matrix();
    open_r_matrix("Pred_comp");
      wrt_r_matrix(pred_fishery_comp);
    close_r_matrix();
    open_r_matrix("Obs_FIM_comp");
      wrt_r_matrix(obs_FIM_comp);
    close_r_matrix();
    open_r_matrix("PredFIM_comp");
      wrt_r_matrix(pred_FIM_comp);
    close_r_matrix();
    open_r_matrix("N");
      wrt_r_matrix(N);
    close_r_matrix();
     wrt_r_complete_vector("Sel" ,fishery_sel);
     wrt_r_complete_vector("Wt", Wtage);
     wrt_r_complete_vector("FIM_Sel" ,FIM_sel);
     wrt_r_complete_vector("Catch_num",pred_totcatch_fishery);
    open_r_matrix("F");
      wrt_r_matrix(F);
    close_r_matrix();
     wrt_r_complete_vector("M",M_vec);
     wrt_r_complete_vector("log_recruit_devs",log_recruit_devs);  //Recruitment Deviations
     wrt_r_complete_vector("Fint",log_fint);     //log catchability deviations
    open_r_vector("params");
      wrt_r_item("log_q",log_q);          //Fishery Catchability 
      wrt_r_item("log_q_FIM",log_q_FIM);          //FIM Catchability 
      wrt_r_item("steepness",steepness);      //Steepness for recruitment
      wrt_r_item("log_R0_FLA",log_R0_FLA);          //Unfished recruitment for Florida cells
      wrt_r_item("log_cv_fishery",log_cv_fishery);           //log cv for fishery catch
      wrt_r_item("log_sigma_rec",log_sigma_rec);            //log sd for recruitment
      wrt_r_item("log_cv_fishery_CPUE",log_cv_fishery_CPUE);          //log cv for fishery CPUE
      wrt_r_item("log_cv_FIM_CPUE",log_cv_FIM_CPUE);          //log cv for FIM CPUE
      wrt_r_item("fishery_sellogis_k",fishery_sellogis_k);        //logistic selectivity parameter 1 for Fishery 
      wrt_r_item("fishery_sellogis_midpt",fishery_sellogis_midpt);    //logistic selectivity parameter 2 for Fishery
      wrt_r_item("FIM_sellogis_k",FIM_sellogis_k);            //logistic selectivity parameter 1 for Survey 
      wrt_r_item("FIM_sellogis_midpt",FIM_sellogis_midpt);        //logistic selectivity parameter 2 for Survey 
      wrt_r_item("ltheta",ltheta);               //Logistic-normal rho parameter for fishery comp
      wrt_r_item("ltheta_FIM",ltheta_FIM);          //Logistic-normal sigma parameter for fishery comp
     close_r_vector();
    close_r_file();
