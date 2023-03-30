# Fischetal2021
Code for Fisch et al., (2021) - Fisheries Research

SpatialModel_Function.R includes the spatially explicit simulation operating model

Get_Data.R includes the Sampling model

Other files in the main directory include data necessary for the spatial model


EMs includes estimation models with double logistic selectivity 
  Dirichlet - EMs fit using Dirichlet likelihood for composition data 
  DMA - EMs fit using Dirichlet-multinomial (asymptotic formulation) likelihood for composition data 
  DML - EMs fit using Dirichlet-multinomial (linear formulation) likelihood for composition data 
  LNAR1 - EMs fit using Logistic-normal with an AR(1) parameterization of the variance-covariance matrix likelihood for composition data 
  LNAR2 - EMs fit using Logistic-normal with an AR(2) parameterization of the variance-covariance matrix likelihood for composition data 
  LNARMA - EMs fit using Logistic-normal with an autoregressive moving average parameterization of the variance-covariance matrix likelihood for composition data 
  MN - EMs fit using Multinomial likelihood for composition data 
  MNFr - EMs fit using Multinomial likelihood for composition data, iteratively weighted using Francis (2011) TA1.8
  MNR - EMs fit using Robust Multinomial likelihood for composition data 
  MNRFr - EMs fit using Robust Multinomial likelihood for composition data, iteratively weighted using Francis (2011) TA1.8

EMs_DGM includes estimation models with logistic selectivity
  with the same acronyms as above 







