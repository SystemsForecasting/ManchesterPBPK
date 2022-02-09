### MAIN - simulated PBPK-ACAT coupled model

### units
#   mass    [mg]
#   volumes [L]
#   time    [h]

### notes
# - the model structure is hard coded, therefore, it is not straightforward changing compartments and compartment order...


#setwd("")

### load libraries & functions -------------------------------------------------------------------------------------

# libraries
library(readxl)
library(RxODE)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# my functions
source("./functions/import_param.R")
source("./functions/PBPK_model_rxode.R")
source("./functions/functions_plot3.R")


### define drug parameters -----------------------------------------------------------------------------------------
# notes:
#     - diffusion coefficient (Diff) and effective thickness of the hydrodynamic diffusion layer are caculated in getPBPKParam function
#     - types: 0 neutral | 1 acid | 2 base

### parameters related to the molecule
# mw     [g/mol]   molecular weight
# type   flag      0 neutral | 1 acid | 2 base
# Pow    [adim]    octanole to water partition coefficient
# Dvow   [adim]    olive oil to water partition coefficient      - generally this parameter is derived from the logP
# Dvow_s [adim]    corrected oil to water partition coefficient  - generally it is derived from logDov
# fup    [adim]    fraction unbound in plasma
# fut    [adim]    fraction unbound in tissues                   - used to calculate the partition coefficients, derived from fup
# BP     [adim]    blood to plasma ratio
# pKa    [adim]
mw     <- 300
type   <- 0
logPow <- 1
fup    <- 0.1
BP     <- 0.5
pKa    <- 7

# derive the other molecular parameters (*1)
# Handerson Hasselback equation
fut     <- 1/( 1 + 0.5*(1-fup)/fup )
logDvow <- 1.115 * logPow - 1.35
pH.tiss <- 7.4                         # HP, see (*1)
if(type==1){
  logDvow_s <- logDvow - log10(1 + 10^(pH.tiss - pKa))
}else if(type==2){
  logDvow_s <- logDvow - log10(1 + 10^(-pH.tiss + pKa))
}else{
  logDvow_s <- logDvow
}

### formulation related parameters
# r      [um]      radius of the particle size of the formulation
# rho    [g/L]     density of the formulation
# Csint  [mg/L]    intrinsic solubility
# Peff   [cm/h]    effective permeability across gut wall
r     <- 2.5
rho   <- 1000
Csint <- 100

# if you have the water solubility and the pH of the solvent, you can derive the intrinsic solubility
Csw <- 100
pHw <- 6
if(type==1){
  Csint_w <- Csw / (1 + 10^(pHw - pKa))
}else if(type==2){
  Csint_w <- Csw / (1 + 10^(-pHw + pKa))
}else{
  Csint_w <- Csw
}
# Csint <- Csint_w  # uncomment if you have water solubility!



# FOR humans: from Papp to Peff, regression!! must be in 10^-4 cm/s (*2)
# FOR mice:   consider to use directly the caco2 permeability... must be in 10^-4 cm/s 
Peff_caco2  <- 1000  # [10^-4 cm/s]
logPeff     <- 0.4926 * log10(Peff_caco2) - 0.1454
Peff        <- 10^(logPeff) * 10^-4 * 3600

### clearances
# CLh    [L/h]     intrinsic hepatic clearance
# CLr    [L/h]     intrinsic renal clearance
# CLent  [L/h]     enterocyte clerance
CLh   <- 10
CLr   <- 10
CLent <- 0

### choose partition coefficients
# "PT"    - Poulin & Theil
# "bere"  - Berezhkovsky
type_part_coeff <- "PT"

### build parameters vector
param_drug = c(Pow    = 10^(logPow),
               Dvow   = 10^(logDvow),
               Dvow_s = 10^(logDvow_s),
               fup    = fup,
               fut    = fut,
               BP     = BP,
               CLh    = CLh,
               CLr    = CLr,
               CLent  = CLent,
               Peff   = Peff,
               r      = r,
               mw     = mw,
               rho    = rho,
               Csint  = Csint,
               pKa    = pKa,
               type   = type
)

# load PBPK parameters
specie          <- "human"
sex             <- "female"
if(specie=="human"){
  if(sex=="female"){
    filename <- "./data/PBPK_parameters/2021_02_27_pbpk_parameters_human_female.xlsx"
  }else{
    filename <- "./data/PBPK_parameters/2021_02_27_pbpk_parameters_human_male.xlsx"
  }
}else if(specie=="mouse"){
  filename <- "./data/PBPK_parameters/2021_03_23_pbpk_parameters_mices.xlsx"
}else if(specie=="beagle"){
  filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_beagles.xlsx"
}else if(specie=="dog"){
  filename <- "./data/PBPK_parameters/2021_04_16_pbpk_parameters_dogs.xlsx"
}

param.PBPK       <- getPBPKParam(filename, param_drug, type_part_coeff, specie)
param.PBPK.rxode <- reorganizeParam.rxode(param.PBPK)
comp.names       <- c(param.PBPK$comp_PBPK_names, param.PBPK$comp_ACAT_names)
names.PBPK       <- param.PBPK$comp_PBPK_names
names.ACAT       <- param.PBPK$comp_ACAT_names
lo               <- length(comp.names)

param.rxode <- c(param.PBPK.rxode, param_drug)

# add glomerular fraction rate (1 yes | 0 no)
param.rxode["GFR_flag"] <- 1

### define changes for multiple simulations ------------------------------------------------------------------------

inits <- c()

ev <- eventTable(amount.units="mg", time.units="hr") %>%
  add.dosing(dose=as.double(10), dosing.to="stomach_s", nbr.doses=1, dosing.interval=12) %>%
  add.sampling(seq(0,24,by=0.01))

ev

# if you want to change some param_drug you need to re-run function for param.PBPK! Especially for formulation properties and parameters used to derive the partition coefficients
# param.PBPK  <- getPBPKParam(filename, param_drug, type_part_coeff)
sim1 <- list(param.rxode = param.rxode, ev = ev, inits = inits)


sim2 <- sim1
sim2$param.rxode["Csint"] <- 0

sim3 <- sim1
sim3$param.rxode["Peff"] <- sim3$param.rxode["Peff"]*0


# here you can add any simulation param that you want...
list.sim <- list(sim1, sim2, sim3)
l.param.set <- length(list.sim)

names.sim <- c("baseline", "Csint*0", "Peff*0")

### simulate the model & plot the system ---------------------------------------------------------------------------
system.out.list <- list()
for(i in 1:l.param.set){
  simi <- list.sim[[i]]
  system.out.list[[i]] <- PBPK.ACAT %>% rxSolve(simi$param.rxode, simi$ev, simi$inits)
}


# plot
flag_log <- 0
data.PK <- list()
axis.limits <- list(xaxis = NA,
                    yaxis = NA)
p.tot <- plotPBPK(system.out.list, list.sim, names.sim = names.sim, names.PBPK, names.ACAT, flag_log, data.PK, axis.limits)

# plot all the organs and tissues (except absorbed and sink enterocytes compartments)
#do.call("grid.arrange", c(p.tot$p.pbpk, ncol=5, nrow=4))      # plot all organs mass PK
#do.call("grid.arrange", c(p.tot$p.acat.1, ncol=6, nrow=4))    # plot ACAT compartments mass PK
p.tot$p.f.excr    # plot fraction excreted
p.tot$p.f.abs     # plot fraction absorbed
p.tot$p.plasma    # plot plasma concentration PK

### some references ------------------------------------------------------------------------------------------------
# (*1) Poulin & Theil 2001 https://doi.org/10.1002/jps.10005
# (*2) Sun 2002 https://doi.org/10.1023/a:1020483911355




