### Functions for importing the parameters of the PBPK model -------------------------------------------------------------------

# Author: Nicola Melillo
# Date: 04/03/2021
# Version: v0.1

# example on how to use these functions is at the end of the file

### output of getPBPKParam
# The output is a list param.PBPK containing 14 elements:
# - organ_bf     [L/h]    organs blood flow
# - organ_w      [Kg]     organs weight
# - organ_d      [Kg/L]   organs density
# - organ_v      [L]      organs volumes
# - ACAT_v_lum   [L]      volumes of the ACAT compartments
# - ACAT_l       [cm]     length of the various ACAT sections
# - ACAT_d       [cm]     diameter of the various ACAT sections
# - ACAT_pH      []       pH in a given ACAT section
# - ACAT_v_ent   [L]      enterocytes volumes of a given ACAT section
# - ACAT_f_CO    []       fraction of the cardiac output (CO) directed to the given enterocytes section
# - general_p             - "weight"    [Kg]  weight of the subjects
#                         - "Htc"       []    haematocrit
#                         - "water_po"  [L]   volume of water drinked together with the oral dose (standard 250 mL)
# - param_drug            drug related parameters specified by the user before calling the function
# - comp_PBPK_names       names of PBPK compartments
# - comp_ACAT_names       names of the ACAT compartments
# - part_coeff   []       partition coefficients
# - physical_param

### PBPK organs for which we have the various parameters (in order)
# 1  lungs
# 2  brain
# 3  heart
# 4  kidneys
# 5  bone
# 6  muscle
# 7  stomach
# 8  spleen
# 9  liver
# 10 gut
# 11 pancreas
# 12 skin
# 13 fat
# 14 arterial_blood
# 15 venous_blood

### ACAT compartments
# 1 stomach
# 2 duodenum
# 3 jejunum1
# 4 jejunum2
# 5 ileum1
# 6 ileum2
# 7 ileum3
# 8 ileum4


### general notes
# - as now only a mean human adult female and male subject of 70kg is supported
# - as now only Poulin & Theil (PT) and Berezhkovsky (bere) partition coefficients are supported

### todos
# - RR partition coefficients

### function get partition coefficients ----------------------------------------------------------------------------------------
# PT: Poulin & Theil (https://doi.org/10.1002/jps.10005)
# bere: Berezhkovsky (https://doi.org/10.1002/jps.20073)
getPartitionCoeff <- function(type, param_drug, param_part_coeff_RR, param_part_coeff_PT) {
  
  # derive drug parameters
  Pow    <- param_drug["Pow"]
  Dvow_s <- param_drug["Dvow_s"]
  fup    <- param_drug["fup"]
  fut    <- param_drug["fut"]
  BP     <- param_drug["BP"]
  
  if (type == "PT" || type == "bere"){ # Poulin & Theil or Berezhkovsky
    
    lo <- length(param_part_coeff_PT$organ_name)
    part_coeff_na <- structure(numeric(lo-1), names=param_part_coeff_PT$organ_name[1:lo-1]) # nonadipose partition coefficient
    part_coeff_a  <- structure(numeric(lo-1), names=param_part_coeff_PT$organ_name[1:lo-1]) # adipose partition coefficient
    part_coeff    <- structure(numeric(lo-1), names=param_part_coeff_PT$organ_name[1:lo-1]) # nonadipose partition coefficient
    
    # convert to named vectors all the various fractions
    Vw_PT  <- structure(param_part_coeff_PT$Vw, names=param_part_coeff_PT$organ_name)
    Vnl_PT <- structure(param_part_coeff_PT$Vnl, names=param_part_coeff_PT$organ_name)
    Vph_PT <- structure(param_part_coeff_PT$Vph, names=param_part_coeff_PT$organ_name)

    #f_n_p <- structure(param_part_coeff_RR$f_n_p, names=param_part_coeff_RR$organ_name)
    #f_n_l <- structure(param_part_coeff_RR$f_n_l, names=param_part_coeff_RR$organ_name)
    #f_ew  <- structure(param_part_coeff_RR$f_ew,  names=param_part_coeff_RR$organ_name)
    #f_iw  <- structure(param_part_coeff_RR$f_iw,  names=param_part_coeff_RR$organ_name)
    #AR    <- structure(param_part_coeff_RR$AR,  names=param_part_coeff_RR$organ_name)
    #LR    <- structure(param_part_coeff_RR$LR,  names=param_part_coeff_RR$organ_name)
    
    # rearranging the parameters...
    Vw_t  <- Vw_PT[-lo]
    Vw_p  <- Vw_PT[lo]
    Vnl_t <- Vnl_PT[-lo]
    Vnl_p <- Vnl_PT[lo]
    Vph_t <- Vph_PT[-lo]
    Vph_p <- Vph_PT[lo]
    #f_n_p_t <- f_n_p[-lo]
    #f_n_l_t <- f_n_l[-lo]
    #f_n_p_p <- f_n_p[lo]
    #f_n_l_p <- f_n_l[lo]
    
    # define some useful parameters for RR
    HCT <- 0.45 #hematocrit
    Kpu_bc <- (HCT - 1 + BP)/(HCT)
    
    if(type == 0){ # netruals
      X <- 0
      Y <- 0
      Z <- 1
    }else if(type == 1){ # monoprotic acids
      X <- 10^(pH_IW-pKa)
      Y <- 10^(pH_P-pKa)
      Z <- 1
    }else if(type ==2){ # monoprotic bases
      X <- 10^(pKa-pH_IW)
      Y <- 10^(pKa-pH_P)
      Z <- 10^(pKa-pH_RBC)
    }
    
    # derive the partition coefficients
    if (type=="PT"){                    # Poulin & Theil
      
      for(i in 1:lo-1){
        
        Vw_ti  <- Vw_PT[param_part_coeff_PT$organ_name[i]]
        Vnl_ti <- Vnl_PT[param_part_coeff_PT$organ_name[i]]
        Vph_ti <- Vph_PT[param_part_coeff_PT$organ_name[i]]
        
        # nonadipose partition coefficients
        part_coeff_na[param_part_coeff_PT$organ_name[i]] <- (  Pow*(Vnl_ti + 0.3*Vph_ti) + (Vw_ti + 0.7*Vph_ti) )/( Pow*(Vnl_p + 0.3*Vph_p) + (Vw_p + 0.7*Vph_p) )*fup/fut
        
        # adipose partition coefficients
        part_coeff_a[param_part_coeff_PT$organ_name[i]]  <- (  Dvow_s*(Vnl_ti + 0.3*Vph_ti) + (Vw_ti + 0.7*Vph_ti) )/( Dvow_s*(Vnl_p + 0.3*Vph_p) + (Vw_p + 0.7*Vph_p) )*fup
      }
      
    }else if(type=="bere"){             # Berezhkovsky
      
      for(i in 1:lo-1){
        
        Vw_ti  <- Vw_PT[param_part_coeff_PT$organ_name[i]]
        Vnl_ti <- Vnl_PT[param_part_coeff_PT$organ_name[i]]
        Vph_ti <- Vph_PT[param_part_coeff_PT$organ_name[i]]
        
        # nonadipose partition coefficients
        part_coeff_na[param_part_coeff_PT$organ_name[i]] <- (  Pow*(Vnl_ti + 0.3*Vph_ti) + (Vw_ti/fut + 0.7*Vph_ti) )/( Pow*(Vnl_p + 0.3*Vph_p) + (Vw_p/fup + 0.7*Vph_p) );
        
        # adipose partition coefficients
        part_coeff_a[param_part_coeff_PT$organ_name[i]]  <- (  Dvow_s*(Vnl_ti + 0.3*Vph_ti) + (Vw_ti/fut + 0.7*Vph_ti) )/( Dvow_s*(Vnl_p + 0.3*Vph_p) + (Vw_p/fup + 0.7*Vph_p) );
      }
      
    }#else if(type=="RR_neutrals"){
      
      #for(i in 1:lo-1){
        
      #  Vw_p  <- Vw_PT[lo]
      #  Vnl_p <- Vnl_PT[lo]
      #  Vph_p <- Vph_PT[lo]
        
      #  Vw_ti  <- Vw_PT[param_part_coeff_PT$organ_name[i]]
      #  Vnl_ti <- Vnl_PT[param_part_coeff_PT$organ_name[i]]
      #  Vph_ti <- Vph_PT[param_part_coeff_PT$organ_name[i]]
        
      #  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
      #  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
        
      #  Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
      #  Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*fup  #lipid
     # }
      
    #}
    
    # set the nonadpose partition coefficients for all the organs except the fat, where the adipose partition coefficients are set
    part_coeff <- part_coeff_na
    part_coeff["fat"] <- part_coeff_a["fat"]
    
  }
  
  return(part_coeff)
  
}


### function get PBPK model parameters -----------------------------------------------------------------------------------------
getPBPKParam <- function(filename, param_drug, type_part_coeff, specie){
  
  ### notes:
  # - filename must be formatted as the file I origially used
  # - to read the excel files I used read_excel from readxl package. This return tibbles. I want as output named vectors, thus, I need to convert them.
  # - all the vectors relative to organs magnitude MUST have the same organs order
  
  library(readxl)
  
  ### read the files (sheet names hard coded - modify for rats)
  sheet_param_organs          <- "parameters_organs"
  sheet_param_general         <- "parameters_general"
  sheet_param_acat            <- "parameters_ACAT"
  sheet_param_part_coeff_RR   <- "parameters_partition_coeff_RR"
  sheet_param_part_coeff_PT   <- "parameters_partition_coeff_PT"
  sheet_param_dissolution     <- "physical_param_dissolution"
  
  param_organs        <- read_excel(filename, sheet=sheet_param_organs)
  param_general       <- read_excel(filename, sheet=sheet_param_general)
  param_ACAT          <- read_excel(filename, sheet=sheet_param_acat)
  param_part_coeff_RR <- read_excel(filename, sheet=sheet_param_part_coeff_RR)
  param_part_coeff_PT <- read_excel(filename, sheet=sheet_param_part_coeff_PT)
  param_diss          <- read_excel(filename, sheet=sheet_param_dissolution)
  
  ### get partition coefficients
  part_coeff <- getPartitionCoeff(type_part_coeff, param_drug, param_part_coeff_RR, param_part_coeff_PT)
  
  ### readapt ORGANS parameters
  organ_bf <- structure(param_organs$blood_flow, names=param_organs$organ_name)  # [L/h]
  organ_d  <- structure(param_organs$density, names=param_organs$organ_name)     # [Kg/L]
  organ_v  <- structure(param_organs$volume, names=param_organs$organ_name)        # [Kg]
  organ_w  <- organ_v*organ_d                                                    # [L]
  
  
  ### readapt ACAT parameters
  ACAT_v_lum <- structure(param_ACAT$volume_lum, names=param_ACAT$compartment)   # [L]
  ACAT_l     <- structure(param_ACAT$length, names=param_ACAT$compartment)       # [cm]
  ACAT_d     <- structure(param_ACAT$diameter, names=param_ACAT$compartment)     # [cm]
  ACAT_pH    <- structure(param_ACAT$pH, names=param_ACAT$compartment)
  ACAT_v_ent <- structure(param_ACAT$volume_ent, names=param_ACAT$compartment)   # [L]
  ACAT_f_CO  <- structure(param_ACAT$fraction_CO, names=param_ACAT$compartment)
  
  ### readapt and define some dissolution parameters
  physical_param <- structure(param_diss$value, names=param_diss$parameter)
  
  ### readapt GENERAL parameters
  # weight [Kg]
  # Htc    
  # water_po
  general_p <- structure(param_general$value, names=param_general$parameter)
  
  ### get PBPK compartments names (add sink compartments)
  comp_PBPK_names <- names(organ_v)
  comp_PBPK_names <- c(comp_PBPK_names, "sink_CLh", "sink_CLr")
  
  ### define some parameter for the dissolution
  # for details of all these equations see supplementary materials of https://doi.org/10.1007/s10928-018-9615-8 
  kb    <- physical_param["kb"]            # [J/K]              Boltzmann constant
  an    <- physical_param["avogadro_num"]  # [mol^-1]           Avogadro number
  eta_w <- physical_param["eta_w"]         # [Pa*s]=[s*J/m^3]   dynamic viscosity of water at 37 °C
  temp  <- physical_param["temp"]          # [K]                absolute temperature 
  
  # derive diffusivity parameter
  # diffusion coefficient calculated with the Stockes Einstein equation
  # hypothesis is that the drug molecule is spherical in shape
  Rh   <- ( (3*param_drug["mw"])/(4*pi*an*param_drug["rho"]) )^(1/3)*10^-1  # [m]                hydrodynamic radius of diffusing drug
  Diff <- (kb*temp/(6*pi*eta_w*Rh))*10^4*60*60;                             # [m^2/s]->[cm^2/h]  diffusion coefficient
  
  # calculate the effective thickness of the hydrodynamic diffusion layer
  # Hintz and Johnson model
  if(param_drug["r"]>30){
    ht <- 30               # [um]
  }else{
    ht <- param_drug["r"]  # [um]
  }
  
  # add Diff and h to the drug related parameters
  param_drug["ht"]   <- ht
  param_drug["Diff"] <- Diff
  
  ### define ACAT compartment names
  comp_ACAT_names_or <- names(ACAT_v_lum)
  comp_ACAT_names_s  <- comp_ACAT_names_or                                  # solid names
  comp_ACAT_names_d  <- comp_ACAT_names_or                                  # dissolved names
  comp_ACAT_names_e  <- comp_ACAT_names_or[2:(length(comp_ACAT_names_or))]  # enterocytes names
  comp_ACAT_names_ea <- comp_ACAT_names_e                                   # drug absorbed
  comp_ACAT_names_ec <- comp_ACAT_names_e                                   # drug cleared by the hepatocytes
  
  for(i in 1:length(comp_ACAT_names_or)){
    comp_ACAT_names_s[i] <- paste(comp_ACAT_names_s[i], "s", sep="_")
    comp_ACAT_names_d[i] <- paste(comp_ACAT_names_d[i], "d", sep="_")
  }
  for(i in 2:length(comp_ACAT_names_or)){
    comp_ACAT_names_e[i-1]  <- paste(comp_ACAT_names_e[i-1], "e", sep="_")
    comp_ACAT_names_ea[i-1] <- paste(comp_ACAT_names_ea[i-1], "ea", sep="_")
    comp_ACAT_names_ec[i-1] <- paste(comp_ACAT_names_ec[i-1], "ec", sep="_")
  }
  comp_ACAT_names <- c(comp_ACAT_names_s, comp_ACAT_names_d, comp_ACAT_names_e, comp_ACAT_names_ea, comp_ACAT_names_ec, "sink_s", "sink_d")
  
  ### define output list
  output_list <- list(organ_bf = organ_bf,
                      organ_w = organ_w,
                      organ_d = organ_d,
                      organ_v = organ_v,
                      ACAT_v_lum = ACAT_v_lum,
                      ACAT_l = ACAT_l,
                      ACAT_d = ACAT_d,
                      ACAT_pH = ACAT_pH,
                      ACAT_v_ent = ACAT_v_ent,
                      ACAT_f_CO = ACAT_f_CO,
                      general_p = general_p,
                      param_drug = param_drug,
                      comp_PBPK_names = comp_PBPK_names,
                      comp_ACAT_names = comp_ACAT_names,
                      part_coeff = part_coeff,
                      physical_param = physical_param
                      )
  
  return(output_list)
  
}






### reorganize parameters for RxODE --------------------------------------------------------------------------------------------

reorganizeParam.rxode <- function(param.PBPK){

  # get names compartments and remove from param.PBPK list
  comp_PBPK_names <- param.PBPK$comp_PBPK_names
  comp_ACAT_names <- param.PBPK$comp_ACAT_names
  param.PBPK$comp_PBPK_names <- NULL
  param.PBPK$comp_ACAT_names <- NULL
  
  names.list <- names(param.PBPK)
  param.PBPK.v <- c()
  
  for(i in 1:length(names.list)){
    vect.i <- unname(param.PBPK[[names.list[i]]])
    names(vect.i) <- paste(names.list[i],names(param.PBPK[[names.list[i]]]),sep="_")
    param.PBPK.v <- c(param.PBPK.v, vect.i)
  }
  
  return(param.PBPK.v)
  
  }




