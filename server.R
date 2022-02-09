# Server PBPK model in R
# Nicola Melillo, Hitesh Mistry, 19/06/2021

#setwd("C:/Users/nicol/Documents/POSTDOC/projects/systemsforcasting/PBPK/codes/2021_05_19_PBPK_app_v08")

# libraries
library(shiny)
library(shinyBS)
library(shinyjs)
library(readxl)
library(writexl)
library(RxODE)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(PKNCA)
library(shinybusy)

# my functions
source("./functions/import_param.R")
source("./functions/PBPK_model_rxode.R")
source("./functions/functions_plot3.R")



shinyServer(function(input, output, session) {
  
  
  values <- reactiveValues(system.out.list = list(),
                           list.sim = list(),
                           names.PBPK = c(),
                           names.ACAT = c(),
                           names.sim = c(),
                           ev_list = c(),
                           NCA = list(),
                           ptot = list(),
                           type_sim = list(),
                           count = 1)  
  flag.clear <- reactiveValues(flag=0)
  flag.simulated <- reactiveValues(flag=0)
  param.drug.library <- reactiveValues(param.drug = list(), 
                                       PK.data = list(), 
                                       PK.data.sel = list())
  axis.limits <- reactiveValues(xaxis = NA,
                                yaxis = NA)
  
  timeShiftIVBolus <- reactiveValues(val = 2/60)
  
  ### other functions --------------------------------------------
  
  # remove initial 2 minutes from all model states following an IV bolus
  removeInitialDataBolusIV <- function(system.out.list, ev_list){
    
    # derive information about the schedule
    idx_sel <- !is.na(ev_list[[1]]$amt)
    if(is.na(ev_list[[1]]$addl[idx_sel])){
      n_doses <- 1
      ii      <- 0
    }else{
      n_doses <- ev_list[[1]]$addl[idx_sel] + 1
      ii      <- as.numeric(ev_list[[1]]$ii[idx_sel])
    }
    
    # remove first 2 minutes for each bolus
    system.out    <- system.out.list[[1]]
    time_dose_ith <- 0
    time_2min     <- timeShiftIVBolus$val # [h]
    variables_PK  <- colnames(system.out)
    
    for(i in 1:n_doses){
      
      idx_sel_i <- (as.numeric(system.out$time) >= time_dose_ith) & (as.numeric(system.out$time) < (time_dose_ith + time_2min))
      system.out <- system.out[!idx_sel_i,]
      time_dose_ith <- time_dose_ith + ii
      
    }
    
    system.out.list.ret <- list(system.out)
    return(system.out.list.ret)
    
  }
  
  ### run the model ------------------------------------
  #runmodel <- eventReactive(input$runButton, {
  observeEvent(input$runButton, {
    withProgress(message = 'Simulating the model', value = 0, {
      mw     <- input$MW
      type   <- as.numeric(input$type)
      logPow <- input$logPow
      fup    <- input$fup
      BP     <- input$BP
      pKa    <- input$pKa
      
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
      r     <- input$r
      rho   <- input$rho
      Csint <- input$Csint
      
      # if you have the water solubility and the pH of the solvent, you can derive the intrinsic solubility
      #Csw <- 100
      #pHw <- 6
      #if(type==1){
      #  Csint_w <- Csw / (1 + 10^(pHw - pKa))
      #}else if(type==2){
      #  Csint_w <- Csw / (1 + 10^(-pHw + pKa))
      #}else{
      #  Csint_w <- logDvow
      #}
      # Csint <- Csint_w  # uncomment if you have water solubility!
      
      
      # FOR humans: from Papp to Peff, regression!! must be in 10^-4 cm/s (*2)
      # FOR mice:   consider to use directly the caco2 permeability... must be in 10^-4 cm/s 
      Peff_caco2  <- input$Peff_caco2 * 3600 / 10000 # [10^-4 cm/s] -> [10^-4 cm/h]
      #logPeff     <- 0.4926 * log10(Peff_caco2) - 0.1454
      Peff        <- Peff_caco2 # 10^(logPeff) * 10^-4 * 3600
      
      ### clearances
      # CLh    [L/h]     intrinsic hepatic clearance
      # CLr    [L/h]     intrinsic renal clearance
      # CLent  [L/h]     enterocyte clerance
      CLh   <- input$Clh
      CLr   <- input$Clr
      CLent <- input$Clent
      
      
      ### choose partition coefficients
      # "PT"    - Poulin & Theil
      # "bere"  - Berezhkovsky
      type_part_coeff <- input$PCM
      
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
      specie          <- input$Species
      if(specie=="human"){
        if(input$Sex=="female"){
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
      }else if(specie=="rat"){
        filename <- "./data/PBPK_parameters/2021_07_16_pbpk_parameters_rats.xlsx"
      }
      
      param.PBPK       <- getPBPKParam(filename, param_drug, type_part_coeff, specie)
      param.PBPK.rxode <- reorganizeParam.rxode(param.PBPK)
      comp.names       <- c(param.PBPK$comp_PBPK_names, param.PBPK$comp_ACAT_names)
      names.PBPK       <- param.PBPK$comp_PBPK_names
      names.ACAT       <- param.PBPK$comp_ACAT_names
      lo               <- length(comp.names)
      
      param.rxode <- c(param.PBPK.rxode, param_drug)
      
      # add flag to include GFR
      if(input$GFR_flag){
        param.rxode["GFR_flag"] <- 1
      }else{
        param.rxode["GFR_flag"] <- 0
      }
      
      # define schedule
      n_rep <- input$daily_admin * input$days
      time_interval <- 24 / input$daily_admin
      time_interval_end <- 24
      
      
      # select compartment for initial conditions
      if (input$Route==1 || input$Route == 4){
        comp_dose = "venous_blood"
      }else if(input$Route==2){
        comp_dose <- "stomach_s"
      }else if(input$Route==3){
        comp_dose <- "stomach_d"
      }
      
      # dose in mg or in mg/kg
      if(input$dose_unit==1){ # mg
        dose <- input$dose
      }else{                  # mg/kg
        dose <- input$dose * param.PBPK$general_p["weight"]
      }
      
      if(input$Route == 4){
        ev <- eventTable(amount.units="mg", time.units="hr") %>%
          add.dosing(dose=as.double(dose), dosing.to=comp_dose, nbr.doses=n_rep, dosing.interval=time_interval, dur=input$inf_dur) %>%
          add.sampling(seq(0,n_rep * time_interval + time_interval_end,by=0.01))
      }else{
        ev <- eventTable(amount.units="mg", time.units="hr") %>%
          add.dosing(dose=as.double(dose), dosing.to=comp_dose, nbr.doses=n_rep, dosing.interval=time_interval) %>%
          add.sampling(seq(0,n_rep * time_interval + time_interval_end,by=0.01))
      }
      
      inits <- c()
      sim1 <- list(param.rxode = param.rxode, ev = ev, inits = inits)
      list.sim <- list(sim1)
      l.param.set <- length(list.sim)
      
      #incProgress(1/2)
      
      ### simulate the model & plot the system 
      system.out.list <- list()
      system.out.list[[1]] <- PBPK.ACAT %>% rxSolve(sim1$param.rxode, sim1$ev, sim1$inits)
      ev_list <- list()
      ev_list[[1]] <- sim1$ev
      
      # remove first 2min following a IV bolus
      #if(input$Route == 1){
      #  system.out.list <- removeInitialDataBolusIV(system.out.list,ev_list)
      #}
      
      
      if(input$keepPlots && flag.clear$flag==0){
        names.sim <- paste("sim",values$count,sep="")
        values$system.out.list <- c(values$system.out.list, system.out.list)
        values$list.sim        <- c(values$list.sim, list.sim)
        values$names.PBPK      <- names.PBPK
        values$names.ACAT      <- names.ACAT
        values$names.sim       <- c(values$names.sim, names.sim)
        values$ev_list         <- c(values$ev_list, ev_list)
        values$type_sim        <- c(values$type_sim, list(input$Route))
        values$count           <- values$count + 1
        values$NCA             <- list()
      }else{
        names.sim <- paste("sim","1",sep="")
        values$system.out.list <- c(system.out.list)
        values$list.sim        <- c(list.sim)
        values$names.PBPK      <- names.PBPK
        values$names.ACAT      <- names.ACAT
        values$names.sim       <- c(names.sim)
        values$ev_list         <- c(ev_list)
        values$type_sim        <- c(list(input$Route))
        values$count           <- 2
        values$NCA             <- list()
        flag.clear$flag <- 0
      }
      
      axis.limits$xaxis <- NA
      axis.limits$yaxis <- NA
      
      
      list.out <- list(system.out.list, list.sim, names.sim, names.PBPK, names.ACAT)
      
      flag.simulated$flag = 1
      incProgress(1/2)
      
      
    })
    
  })
  
  
  ### plot --------------------------------------------------------------------------
  output$PK<-renderPlot({
    withProgress(message = 'Plotting plasma PK', value = 0, {
      if(input$logscale){
        logscale <- 1
      }else{
        logscale <- 0
      }
      
      
      if(flag.clear$flag==0 && flag.simulated$flag == 1){
        
        # remove first 2 min for plotting when there is a IV bolus
        system.out.plot <- values$system.out.list
        for(i in 1:length(values$type_sim)){
          if(as.numeric(values$type_sim[i]) == 1){
            a <- 1
            system.out.plot[i] <- removeInitialDataBolusIV(system.out.plot[i],values$ev_list[i])
          }
        }
        
        # plot
        p.tot <- plotPBPK(system.out.plot, values$list.sim, names.sim = values$names.sim, values$names.PBPK, values$names.ACAT,logscale,param.drug.library$PK.data.sel,axis.limits)
        incProgress(1/4)
        ptlist <- list(p.tot$p.plasma, p.tot$p.f.excr, p.tot$p.f.abs)
        do.call("grid.arrange", c(ptlist, ncol=1, nrow=3))
        flag.simulated$flag = 1
        values$ptot <- p.tot
        incProgress(3/4)
      }
    })
  })
  
  output$PK_comp_PBPK<-renderPlot({
    
    if(flag.clear$flag==0 && flag.simulated$flag==1 && input$plotOrgansPK){
      withProgress(message = 'Plotting organs PK', value = 0, {
        if(input$logscale){
          logscale <- 1
        }else{
          logscale <- 0
        }
        
        # remove first 2 min for plotting when there is a IV bolus
        system.out.plot <- values$system.out.list
        for(i in 1:length(values$type_sim)){
          if(as.numeric(values$type_sim[i]) == 1){
            a <- 1
            system.out.plot[i] <- removeInitialDataBolusIV(system.out.plot[i],values$ev_list[i])
          }
        }
        
        p.tot <- plotPBPK(system.out.plot, values$list.sim, names.sim = values$names.sim, values$names.PBPK, values$names.ACAT,logscale,param.drug.library$PK.data.sel,axis.limits)
        incProgress(1/4)
        do.call("grid.arrange", c(c(p.tot$p.pbpk,p.tot$p.acat.1), ncol=3, nrow=13))
        values$ptot <- p.tot
        incProgress(3/4)
      })
    }
    
  })
  
  
  ### axis appearance --------------------------------------------------------------------------------------------------
  
  # xaxis slider
  output$UI_xaxis_slider <- renderUI({
    
    # find max value
    if(flag.simulated$flag == 1){
      n.sim <- length(values$system.out.list)
      value_axis_max <- 0
      for(i in 1:n.sim){
        l.sim <- length(values$system.out.list[[i]]$time)
        if(value_axis_max<as.numeric(values$system.out.list[[i]]$time[l.sim])){
          value_axis_max <- as.numeric(values$system.out.list[[i]]$time[l.sim])
        }
      }
    }else{
      value_axis_max <- 48
    }
    
    # min value set to default 0
    value_axis_min <- 0
    
    # define slider
    digits.signif <- 4
    sliderInput("xaxis_slider", label = "x axis limits", min = signif(value_axis_min, digits=digits.signif), 
                max = signif(value_axis_max, digits=digits.signif), value = c(signif(value_axis_min, digits=digits.signif), signif(value_axis_max, digits=digits.signif)))
    

  })
  
  # yaxis slider
  output$UI_yaxis_slider <- renderUI({
    
    # find max value
    if(flag.simulated$flag == 1){
      n.sim <- length(values$system.out.list)
      value_axis_max <- 0
      for(i in 1:n.sim){
        if(value_axis_max<max(as.numeric(values$system.out.list[[i]]$plasma_conc))){
          value_axis_max <- max(as.numeric(values$system.out.list[[i]]$plasma_conc))
        }
      }
    }else{
      value_axis_max <- 10
    }
    
    # min set to default 0
    value_axis_min <- 0
    
    # define slider
    digits.signif <- 4
    if(input$logscale){
      value_axis_min <- log10(value_axis_max/1e+6)
      value_axis_max <- log10(value_axis_max)
    }
    sliderInput("yaxis_slider", label = "y axis limits", min = signif(value_axis_min, digits=digits.signif), 
                max = signif(value_axis_max, digits=digits.signif), value = c(signif(value_axis_min, digits=digits.signif), signif(value_axis_max, digits=digits.signif)))
    
  })
  
  # set axis limits values according when relative button is pressed
  observeEvent(input$scalex, {
    axis.limits$xaxis <- input$xaxis_slider
  })
  
  observeEvent(input$scaley, {
    axis.limits$yaxis <- input$yaxis_slider
  })
  
  ### functions for activating/deactivating options according to given events -------------------------------------------
  observeEvent(input$Route, {
    if(input$Route == 4){
      shinyjs::enable("inf_dur")
    }else{
      shinyjs::disable("inf_dur")
    }
  })
  
  observeEvent(input$Species, {
    if(input$Species == "human"){
      shinyjs::enable("Sex")
    }else{
      shinyjs::disable("Sex")
    }
  })
  
  observeEvent(input$clearPlot, {
    #hide("PK")
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    param.drug.library$PK.data.sel <- list()
    axis.limits$xaxis <- NA
    axis.limits$yaxis <- NA
    values$type_sim <- list()
  })
  
  observeEvent(input$download_menu, {
    if(flag.simulated$flag == 1){
      shinyjs::enable("downloadPK")
      shinyjs::enable("radio_download_PK")
    }else{
      shinyjs::disable("downloadPK")
      shinyjs::disable("radio_download_PK")
    }
  })
  
  observeEvent(input$download_menu, {
    if(length(values$NCA) == 0){
      shinyjs::disable("downloadNCA")
    }else if(flag.simulated$flag == 0){
      shinyjs::disable("downloadNCA")
      values$NCA <- list()
    }else{
      shinyjs::enable("downloadNCA")
    }
  })
  
  observeEvent(input$download_menu, {
    if(flag.simulated$flag == 1){
      shinyjs::enable("downloadPlots")
      shinyjs::enable("checkPlots")
    }else{
      shinyjs::disable("downloadPlots")
      shinyjs::disable("checkPlots")
    }
  })
  
  observeEvent(input$checkPlots, {
    if(is.null(input$checkPlots)){
      shinyjs::disable("downloadPlots")
    }else{
      shinyjs::enable("downloadPlots")
    }
  }, ignoreNULL = FALSE)
  
  
  ### upload parameters values & PK from library ------------------------------------------------------
  output$libraryDrugs <- renderUI({
    param_drug <- read_excel("./data/library_drugs/drugsParam.xlsx", sheet="drug_param")
    param.drug.library$param.drug <- param_drug
    mydata <- param_drug$drug
    selectInput('selectedDrug', 'Select drug', c(Choose='', mydata), selectize=FALSE)
  })
  
  # upload the parameters from library
  observeEvent(input$uploadDrugParam, {
    
    # remove plots & old stored data for other drugs
    param.drug.library$PK.data.sel <- list()
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    
    if(input$selectedDrug!=""){
      
      # select drug parameters
      library.param <- param.drug.library$param.drug[param.drug.library$param.drug$drug==input$selectedDrug,]
      
      # molecular related parameters
      updateNumericInput(session,"MW", value = library.param$mw)
      if(library.param$type=="neutral"){
        updateRadioButtons(session,"type",selected=0)
        Csint <- library.param$Cs_w
      }else if(library.param$type=="acid"){
        updateRadioButtons(session,"type",selected=1)
        updateNumericInput(session,"pKa", value = library.param$pKa)
        Csint <- library.param$Cs_w / (1 + 10^(-library.param$pKa + library.param$pH_ref))
      }else{
        updateRadioButtons(session,"type",selected=2)
        updateNumericInput(session,"pKa", value = library.param$pKa)
        Csint <- library.param$Cs_w / (1 + 10^( library.param$pKa - library.param$pH_ref))
      }
      updateNumericInput(session,"logPow", value = library.param$logPow)
      updateNumericInput(session,"fup", value = library.param$fup)
      updateNumericInput(session,"BP", value = library.param$BP)
      
      # dissolution and absorption parameters
      updateNumericInput(session,"r", value = library.param$r)
      updateNumericInput(session,"rho", value = library.param$rho)
      updateNumericInput(session,"Csint", value = Csint)
      updateNumericInput(session,"Peff_caco2", value = library.param$Peff)
      
      # clearance parameters
      updateNumericInput(session,"Clh", value = library.param$Clh)
      updateNumericInput(session,"Clr", value = library.param$Clr)
      updateNumericInput(session,"Clent", value = library.param$Clent)
    }
    
  })
  
  # define selectInput with drug specific data
  output$PKData <- renderUI({
    req(input$selectedDrug)
    PK.data <- read_excel("./data/library_drugs/dataPK.xlsx", sheet=input$selectedDrug)
    param.drug.library$PK.data <- PK.data
    mydata <- unique(PK.data$schedule_ID)
    selectInput('selectPKData', 'Select PK data', c(Choose='', mydata), selectize=FALSE)
  })
  
  # show action button only if a drug is chosen
  output$UploadPKData <- renderUI({
    req(input$selectedDrug)
    actionButton("uploadDrugPK", "Upload PK", width='100px')
  })
  
  # upload the chosen PK data
  observeEvent(input$uploadDrugPK,{
    param.drug.library$PK.data.sel <- param.drug.library$PK.data[param.drug.library$PK.data$schedule_ID==input$selectPKData,]
    
    # update schedule
    updateNumericInput(session,"dose", value = param.drug.library$PK.data.sel$dose[1])
    updateRadioButtons(session,"Species",selected="human")
    updateRadioButtons(session,"Sex",selected="female")
    if(param.drug.library$PK.data.sel$route[1]=="PO"){
      updateRadioButtons(session,"Route",selected=2)
    }else if(param.drug.library$PK.data.sel$route[1]=="IV"){
      
      if(is.na(param.drug.library$PK.data.sel$inf_dur[1])){
        updateRadioButtons(session,"Route",selected=1)
      }else{
        updateRadioButtons(session,"Route",selected=4)
        updateNumericInput(session,"inf_dur", value = as.numeric(param.drug.library$PK.data.sel$inf_dur[1])/60) # infusion duration in PK dataset is always in minutes 
      }
        
    }
  })
  
  ### reset button --------------------------------------------------------------------------------------
  observeEvent(input$resetButton,{
    
    # remove plots & old stored data for other drugs
    param.drug.library$PK.data <- list()
    param.drug.library$PK.data.sel <- list()
    values$type_sim <- list()
    flag.clear$flag <- 1
    flag.simulated$flag <- 0
    updateSelectInput(session, "selectedDrug",selected="")
    axis.limits$xaxis <- NA
    axis.limits$yaxis <- NA
    
    ### default parameters
    
    # schedule
    updateRadioButtons(session,"Route",selected=2) 
    updateNumericInput(session,"inf_dur", value = 0.5)
    updateNumericInput(session,"dose", value = 10)
    updateRadioButtons(session,"dose_unit", selected = 1)
    updateNumericInput(session,"days", value = 1)
    updateSliderInput(session,"daily_admin", value=1)
    
    # molecular related parameters
    updateNumericInput(session,"MW", value = 500)
    updateRadioButtons(session,"type",selected=0)
    updateNumericInput(session,"logPow", value = 2)
    updateNumericInput(session,"fup", value = 0.8)
    updateNumericInput(session,"BP", value = 0.8)
    
    # dissolution and absorption parameters
    updateNumericInput(session,"r", value = 25)
    updateNumericInput(session,"rho", value = 1000)
    updateNumericInput(session,"Csint", value = 100)
    updateNumericInput(session,"Peff_caco2", value = 2)
    
    # clearance parameters
    updateNumericInput(session,"Clh", value = 10)
    updateNumericInput(session,"Clr", value = 10)
    updateNumericInput(session,"Clent", value = 0)
    updateRadioButtons(session,"PCM",selected="PT")
    
    # specie and sex
    updateRadioButtons(session,"Species",selected="human")
    updateRadioButtons(session,"Sex",selected="female")
    
    
  })
  
  ### noncompartmental analysis ----------------------------------------------------------------------------------
  # NCA is performed with the package PKNCA
  # for multiple doses I select only the last dose!!
  output$table_nca <- renderDataTable({
    
    if(flag.simulated$flag==1){
      show_modal_spinner(spin="circle", text="performing NCA, it can take some time...")
      
      df.conc <- data.frame(matrix(ncol = 3, nrow = 0))
      colnames(df.conc) <- c("Time", "conc", "Subject")
      df.dose <- data.frame(matrix(ncol=5, nrow=0))
      colnames.dose <- c("Time", "Dose", "Subject", "route", "duration")
      colnames(df.dose) <- colnames.dose
      
      for(i in 1:length(values$system.out.list)){
        
        # create concentration dataset
        conc    <- values$system.out.list[[i]]$plasma_conc
        Time    <- as.numeric(values$system.out.list[[i]]$time)
        Subject <- numeric(length(time)) + i
        df.conc <- rbind(df.conc,data.frame(Time,conc,Subject))
        
        # create dose dataset (choose last administration for multiple doses)
        idx_sel   <- !is.na(values$ev_list[[i]]$amt)
        if(!is.na(values$ev_list[[i]]$addl[idx_sel])){
          Time.dose <- as.numeric(values$ev_list[[i]]$addl[idx_sel] * values$ev_list[[i]]$ii[idx_sel])
        }else{
          Time.dose <- as.numeric(values$ev_list[[i]]$time[idx_sel])
        }
        
        if(values$ev_list[[i]]$cmt[idx_sel]=="venous_blood"){
          route <- "intravascular"
        }else{
          route <- "extravascular"
        }
        df.dose.i <- data.frame(Time.dose,
                                as.numeric(values$ev_list[[i]]$amt[idx_sel]),
                                i,
                                route,
                                as.numeric(values$ev_list[[i]]$dur[idx_sel]))
        colnames(df.dose.i) <- colnames.dose
        df.dose <- rbind(df.dose, df.dose.i)
        
      }
      
      # create NCA objects & perform NCA
      conc_obj <- PKNCAconc(df.conc, conc~Time|Subject)
      dose_obj <- PKNCAdose(df.dose, Dose~Time|Subject, route = df.dose$route, duration = df.dose$duration)
      
      data_obj <- PKNCAdata(conc_obj, dose_obj)
      
      # add cl.all and vss.obs
      n.subj <- length(values$system.out.list)
      data_obj$intervals[seq(from=2, to=2*n.subj, by=2),"cl.all"]  <- TRUE
      data_obj$intervals[seq(from=1, to=(2*n.subj-1), by=2),"vss.obs"] <- TRUE
      
      # calculate NCA
      results_obj <- pk.nca(data_obj)
      
      # create dataset containing NCA information to print
      df.nca.out <- data.frame(matrix(ncol=8, nrow=0))
      colnames.nca <- c("simulation","AUClast [h*mg/L]", "Cmax [mg/L]", "tmax [h]", "half_life [h]","AUCinf [h*mg/L]","plasma CL/F [L/h]","Vss obs [L]")
      colnames(df.nca.out) <- colnames.nca
      for(i in 1:length(values$system.out.list)){
        idx_sel_i <- results_obj$result$Subject==i
        result_i  <- results_obj$result[idx_sel_i,]
        
        # find parameters
        if(values$type_sim[[i]]!=1){
          df.nca.i <- data.frame(i,
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "auclast"],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "cmax"],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "tmax"][1],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "half.life"][1],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "aucinf.obs"][2],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "cl.all"],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "vss.obs"],4))
        }else{
          df.nca.i <- data.frame(i,
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "auclast"],4),
                                 NA,
                                 NA,
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "half.life"][1],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "aucinf.obs"][2],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "cl.all"],4),
                                 signif(result_i$PPORRES[result_i$PPTESTCD == "vss.obs"],4))
        }
        
        colnames(df.nca.i) <- colnames.nca
        df.nca.out <- rbind(df.nca.out, df.nca.i)
      }
      
      remove_modal_spinner()
      
    }else{
      df.nca.out <- data.frame()
    }
    values$NCA <- df.nca.out
    return(df.nca.out)
  })
  
  ### download -------------------------------------------------------------------------------------------------
  
  output$downloadNCA <- downloadHandler(
    filename = function() { "NCA.xlsx"},
    content = function(file) {
      
      # find information regarding dose etc for each simulation
      length.sim <- length(values$system.out.list)
      df.nca.out <- data.frame(matrix(ncol=5, nrow=0))
      colnames.df <- c("amt [mg]", "cmt", "interval [h]", "additional doses", "infusion duration [h]")
      colnames(df.nca.out) <- colnames.df
      
      for(i in 1:length.sim){
        
        system.out <- values$system.out.list[[i]]
        length.so <- length(system.out$time)
        
        # derive information about the schedule
        idx_sel <- !is.na(as.numeric(values$ev_list[[i]]$amt))
        amt_v   <- as.numeric(values$ev_list[[i]]$amt[idx_sel])
        cmt_v   <- values$ev_list[[i]]$cmt[idx_sel]
        ii_v    <- as.numeric(values$ev_list[[i]]$ii[idx_sel])
        addl_v  <- as.numeric(values$ev_list[[i]]$addl[idx_sel])
        dur_v   <- as.numeric(values$ev_list[[i]]$dur[idx_sel])
        
        df.i <- data.frame(amt_v, cmt_v, ii_v, addl_v, dur_v)
        colnames(df.i) <- colnames.df
        df.nca.out <- rbind(df.nca.out, df.i)
        
      }
      
      df.out <- cbind(values$NCA, df.nca.out)
      
      write_xlsx(df.out, path = file)
      }
  ) 
  
  organizeTablePK <- function(){
    length.sim <- length(values$system.out.list)
    
    # initialize df.out
    if(input$radio_download_PK==1){
      names.col <- c("simulation","dose","cmt","interval [h]","additional doses", "infusion duration [h]", "time [h]", "plasma conc [mg/L]")
      df.out <- data.frame(matrix(ncol = length(names.col), nrow = 0))
      colnames(df.out) <- names.col
    }else{
      names.col <- c("simulation","dose","cmt","interval [h]","additional doses", "infusion duration [h]", "time [h]", "plasma conc [mg/L]", paste(values$names.PBPK, "[mg]"), paste(values$names.ACAT, "[mg]"))
      df.out <- data.frame(matrix(ncol = length(names.col), nrow = 0))
      colnames(df.out) <- names.col
    }
    
    
    # fill df.out
    for(i in 1:length.sim){
      
      system.out <- values$system.out.list[[i]]
      length.so <- length(system.out$time)
      
      # derive information about the schedule
      idx_sel <- !is.na(as.numeric(values$ev_list[[i]]$amt))
      amt_v   <- rep(as.numeric(values$ev_list[[i]]$amt[idx_sel]), length.so)
      cmt_v   <- rep((values$ev_list[[i]]$cmt[idx_sel]), length.so)
      ii_v    <- rep(as.numeric(values$ev_list[[i]]$ii[idx_sel]), length.so)
      addl_v  <- rep(as.numeric(values$ev_list[[i]]$addl[idx_sel]), length.so)
      dur_v   <- rep(as.numeric(values$ev_list[[i]]$dur[idx_sel]), length.so)
      subj_v  <- rep(i,length.so)
      
      if(input$radio_download_PK==1){
        df.out.i <- data.frame(subj_v, amt_v, cmt_v, ii_v, addl_v, dur_v, system.out$time, system.out$plasma_conc)
      }else{
        df.out.i <- data.frame(subj_v, amt_v, cmt_v, ii_v, addl_v, dur_v, system.out$time, system.out$plasma_conc, system.out[values$names.PBPK], system.out[values$names.ACAT])
      }
      
      colnames(df.out.i) <- names.col
      df.out <- rbind(df.out, df.out.i)
    }
    
    return(df.out)
  }
  
  
  output$downloadPK <- downloadHandler(
    filename = function() {"PK_data.xlsx"},
    content = function(file) {
      write_xlsx(organizeTablePK(), path = file)
    }
  )
  
  output$downloadPlots <- downloadHandler(
    filename = function() {"plots.zip"},
    
    content = function(file) {
      
      show_modal_spinner(spin="circle", text="building plots...")
      
      dir.original <- getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
      names.field <- c("p.plasma","p.f.excr","p.f.abs","p.pbpk","p.acat.1")
      filename    <- c("plasma_conc", "fraction_eliminated", "fraction_absorbed", "organs_PK", "ACAT_PK")
      fs <- c()
      
      for (i in 1:length(input$checkPlots)) {

        dpi.i <- 600
        scale.i <- 1
        filename_i <- paste(filename[ as.numeric(input$checkPlots[i]) ],".png",sep="")
        
        if(as.numeric(input$checkPlots[i])<4){
          plot_i <- values$ptot[[ names.field[ as.numeric(input$checkPlots[i]) ] ]]
        }else if(as.numeric(input$checkPlots[i])==4){
          plot_i <- do.call("grid.arrange", c(c(values$ptot$p.pbpk), ncol=5, nrow=4))
        }else{
          plot_i <- do.call("grid.arrange", c(c(values$ptot$p.acat.1), ncol=5, nrow=5))
        }
        
        
        ggsave(
          filename_i,
          device = "png",
          plot = plot_i,
          path = "./",
          scale = 1,
          units = c("cm"),
          dpi = dpi.i,
          width = 25,
          height = 25
        )
        
        fs <- c(fs, paste(filename_i,sep=""))
        
      }
      
      remove_modal_spinner()
      
      zip(zipfile=file, files=fs)
      setwd(dir.original)
    }
  )
  
})

