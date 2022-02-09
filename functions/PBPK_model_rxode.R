### functions describing the PBPK model

### Units
# x        - [mg]
# volumes  - [L]
# time     - [h]

### PBPK organs for which we have the various parameters (in order)

# PBPK distribution
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

# PBPK sink
# 16 sink liver
# 17 sink kidney

# ACAT solid
# 18 stomach_s
# 19 duodenum_s
# 20 jejunum1_s
# 21 jejunum2_s
# 22 ileum1_s
# 23 ileum2_s
# 24 ileum3_s

# ACAT dissolved
# 25 stomach_d
# 26 duodenum_d
# 27 jejunum1_d
# 28 jejunum2_d
# 29 ileum1_d
# 30 ileum2_d
# 31 ileum3_d

# ACAT enterocytes
# 32 duodenum_e
# 33 jejunum1_e
# 34 jejunum2_e
# 35 ileum1_e
# 36 ileum2_e
# 37 ileum3_e

# ACAT_absorbed
# 38 duodenum_ea
# 39 jejunum1_ea
# 40 jejunum2_ea
# 41 ileum1_ea
# 42 ileum2_ea
# 43 ileum3_ea

# ACAT cleared
# 44 duodenum_ec
# 45 jejunum1_ec
# 46 jejunum2_ec
# 47 ileum1_ec
# 48 ileum2_ec
# 49 ileum3_ec

# ACAT sink
# 50 sink_s
# 51 sink_d


PBPK.ACAT <- RxODE({
  
  ### ACAT model ---------------------------------------------------------------------------------
  
  input_GI <- 0
  
  # term in front of the dissolution equation
  Kd1 <- 10^5*(3*param_drug_Diff/(rho*param_drug_ht*r)) # [L/(h*mg)] 
  
  # define Henderson-Hasselback term for dissolution
  if(type==0){           # neutral
    Ys  <- 1
    Yd  <- 1
    Yj1 <- 1
    Yj2 <- 1
    Yi1 <- 1
    Yi2 <- 1
    Yi3 <- 1
  }else if(type==1){     # acid
    Ys  <- 1 + 10^( ACAT_pH_stomach  - pKa )
    Yd  <- 1 + 10^( ACAT_pH_duodenum - pKa )
    Yj1 <- 1 + 10^( ACAT_pH_jejunum1 - pKa )
    Yj2 <- 1 + 10^( ACAT_pH_jejunum2 - pKa )
    Yi1 <- 1 + 10^( ACAT_pH_ileum1   - pKa )
    Yi2 <- 1 + 10^( ACAT_pH_ileum2   - pKa )
    Yi3 <- 1 + 10^( ACAT_pH_ileum3   - pKa )
  } else if (type==2){   # base
    Ys  <- 1 + 10^( -ACAT_pH_stomach  + pKa )
    Yd  <- 1 + 10^( -ACAT_pH_duodenum + pKa )
    Yj1 <- 1 + 10^( -ACAT_pH_jejunum1 + pKa )
    Yj2 <- 1 + 10^( -ACAT_pH_jejunum2 + pKa )
    Yi1 <- 1 + 10^( -ACAT_pH_ileum1   + pKa )
    Yi2 <- 1 + 10^( -ACAT_pH_ileum2   + pKa )
    Yi3 <- 1 + 10^( -ACAT_pH_ileum3   + pKa ) 
  }
  
  ### stomach equations
  
  # define parameters and useful indices
  Cs_st <- Csint * Ys
  Kd_st <- Kd1 * ( Cs_st - stomach_d/(ACAT_v_lum_stomach + general_p_water_po))
  kst   <- 1/general_p_GET # [1/h]
  
  # define differential equations
  output_s_s   <- kst * stomach_s
  output_d_s   <- kst * stomach_d
  d/dt(stomach_s) <- -output_s_s - Kd_st * stomach_s
  d/dt(stomach_d) <- -output_d_s + Kd_st * stomach_s
  
  
  ###### equations for all the small intestine transit compartments
  
  length_small_int <- ACAT_l_duodenum + ACAT_l_jejunum1 + ACAT_l_jejunum2 + ACAT_l_ileum1 + ACAT_l_ileum2 + ACAT_l_ileum3
  CO <- organ_bf_venous_blood  # cardiac output [L/h]
  
  ### duodenum 
  Cs_it <- Csint * Yd
  Kd_it <- Kd1 * ( Cs_it - duodenum_d/ACAT_v_lum_duodenum)
  kit   <- 1/(general_p_SITT * ACAT_l_duodenum/length_small_int)       # [1/h]
  ka    <- 2 * Peff / (ACAT_d_duodenum/2)                      # [1/h]
  Qe_i  <- ACAT_f_CO_duodenum * CO                         # [L/h]
  
  # define differential equations
  d/dt(duodenum_s)  <- +output_s_s - Kd_it * duodenum_s - kit * duodenum_s
  d/dt(duodenum_d)  <- +output_d_s + Kd_it * duodenum_s - kit * duodenum_d - ka * duodenum_d
  d/dt(duodenum_e)  <- ka * duodenum_d - Qe_i * duodenum_e/ACAT_v_ent_duodenum - CLent * duodenum_e/ACAT_v_ent_duodenum
  d/dt(duodenum_ea) <- Qe_i  * duodenum_e/ACAT_v_ent_duodenum
  d/dt(duodenum_ec) <- CLent * duodenum_e/ACAT_v_ent_duodenum
  
  output_s_d   <- kit * duodenum_s
  output_d_d   <- kit * duodenum_d
  input_GI     <- input_GI + Qe_i  * duodenum_e/ACAT_v_ent_duodenum
  
  
  ### jejunum 1 
  Cs_it_j1 <- Csint * Yj1
  Kd_it_j1 <- Kd1 * ( Cs_it_j1 - jejunum1_d/ACAT_v_lum_jejunum1)
  kit      <- 1/(general_p_SITT * ACAT_l_jejunum1/length_small_int)       # [1/h]
  ka       <- 2 * Peff / (ACAT_d_jejunum1/2)                      # [1/h]
  Qe_i     <- ACAT_f_CO_jejunum1 * CO                         # [L/h]
  
  # define differential equations
  d/dt(jejunum1_s)  <- +output_s_d - Kd_it_j1 * jejunum1_s - kit * jejunum1_s
  d/dt(jejunum1_d)  <- +output_d_d + Kd_it_j1 * jejunum1_s - kit * jejunum1_d - ka * jejunum1_d
  d/dt(jejunum1_e)  <- ka * jejunum1_d - Qe_i * jejunum1_e/ACAT_v_ent_jejunum1 - CLent * jejunum1_e/ACAT_v_ent_jejunum1
  d/dt(jejunum1_ea) <- Qe_i  * jejunum1_e/ACAT_v_ent_jejunum1
  d/dt(jejunum1_ec) <- CLent * jejunum1_e/ACAT_v_ent_jejunum1
  
  output_s_j1   <- kit * jejunum1_s
  output_d_j1   <- kit * jejunum1_d
  input_GI      <- input_GI + Qe_i  * jejunum1_e/ACAT_v_ent_jejunum1
  
  ### jejunum 2 
  Cs_it_j2 <- Csint * Yj2
  Kd_it_j2 <- Kd1 * ( Cs_it_j2 - jejunum2_d/ACAT_v_lum_jejunum2)
  kit      <- 1/(general_p_SITT * ACAT_l_jejunum2/length_small_int)      # [1/h]
  ka       <- 2 * Peff / (ACAT_d_jejunum2/2)                      # [1/h]
  Qe_i     <- ACAT_f_CO_jejunum2 * CO                         # [L/h]
  
  # define differential equations
  d/dt(jejunum2_s)  <- +output_s_j1 - Kd_it_j2 * jejunum2_s - kit * jejunum2_s
  d/dt(jejunum2_d)  <- +output_d_j1 + Kd_it_j2 * jejunum2_s - kit * jejunum2_d - ka * jejunum2_d
  d/dt(jejunum2_e)  <- ka * jejunum2_d - Qe_i * jejunum2_e/ACAT_v_ent_jejunum2 - CLent * jejunum2_e/ACAT_v_ent_jejunum2
  d/dt(jejunum2_ea) <- Qe_i  * jejunum2_e/ACAT_v_ent_jejunum2
  d/dt(jejunum2_ec) <- CLent * jejunum2_e/ACAT_v_ent_jejunum2
  
  output_s_j2   <- kit * jejunum2_s
  output_d_j2   <- kit * jejunum2_d
  input_GI      <- input_GI +  Qe_i * jejunum2_e/ACAT_v_ent_jejunum2
  
  ### ileum 1
  Cs_it_i1 <- Csint * Yi1
  Kd_it_i1 <- Kd1 * ( Cs_it_i1 - ileum1_d/ACAT_v_lum_ileum1)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum1/length_small_int)       # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum1/2)                      # [1/h]
  Qe_i     <- ACAT_f_CO_ileum1 * CO                         # [L/h]
  
  # define differential equations
  d/dt(ileum1_s)  <- +output_s_j2 - Kd_it_i1 * ileum1_s - kit * ileum1_s
  d/dt(ileum1_d)  <- +output_d_j2 + Kd_it_i1 * ileum1_s - kit * ileum1_d - ka * ileum1_d
  d/dt(ileum1_e)  <- ka * ileum1_d - Qe_i * ileum1_e/ACAT_v_ent_ileum1 - CLent * ileum1_e/ACAT_v_ent_ileum1
  d/dt(ileum1_ea) <- Qe_i  * ileum1_e/ACAT_v_ent_ileum1
  d/dt(ileum1_ec) <- CLent * ileum1_e/ACAT_v_ent_ileum1
  
  output_s_i1   <- kit * ileum1_s
  output_d_i1   <- kit * ileum1_d
  input_GI      <- input_GI + Qe_i  * ileum1_e/ACAT_v_ent_ileum1
  
  
  ### ileum 2
  Cs_it_i2 <- Csint * Yi2
  Kd_it_i2 <- Kd1 * ( Cs_it_i2 - ileum2_d/ACAT_v_lum_ileum2)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum2/length_small_int)      # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum2/2)                      # [1/h]
  Qe_i     <- ACAT_f_CO_ileum2 * CO                         # [L/h]
  
  # define differential equations
  d/dt(ileum2_s)  <- +output_s_i1 - Kd_it_i2 * ileum2_s - kit * ileum2_s
  d/dt(ileum2_d)  <- +output_d_i1 + Kd_it_i2 * ileum2_s - kit * ileum2_d - ka * ileum2_d
  d/dt(ileum2_e)  <- ka * ileum2_d - Qe_i * ileum2_e/ACAT_v_ent_ileum2 - CLent * ileum2_e/ACAT_v_ent_ileum2
  d/dt(ileum2_ea) <- Qe_i  * ileum2_e/ACAT_v_ent_ileum2
  d/dt(ileum2_ec) <- CLent * ileum2_e/ACAT_v_ent_ileum2
  
  output_s_i2   <- kit * ileum2_s
  output_d_i2   <- kit * ileum2_d
  input_GI      <- input_GI + Qe_i  * ileum2_e/ACAT_v_ent_ileum2
  
  ### ileum 3
  Cs_it_i3 <- Csint * Yi3
  Kd_it_i3 <- Kd1 * ( Cs_it_i3 - ileum3_d/ACAT_v_lum_ileum3)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum3/length_small_int )      # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum3/2)                  # [1/h]
  Qe_i     <- ACAT_f_CO_ileum3 * CO                         # [L/h]
  
  # define differential equations
  d/dt(ileum3_s)  <- +output_s_i2 - Kd_it_i3 * ileum3_s - kit * ileum3_s
  d/dt(ileum3_d)  <- +output_d_i2 + Kd_it_i3 * ileum3_s - kit * ileum3_d - ka * ileum3_d
  d/dt(ileum3_e)  <- ka * ileum3_d - Qe_i * ileum3_e/ACAT_v_ent_ileum3 - CLent * ileum3_e/ACAT_v_ent_ileum3
  d/dt(ileum3_ea) <- Qe_i  * ileum3_e/ACAT_v_ent_ileum3
  d/dt(ileum3_ec) <- CLent * ileum3_e/ACAT_v_ent_ileum3
  
  output_s_i3   <- kit * ileum3_s
  output_d_i3   <- kit * ileum3_d
  input_GI      <- input_GI + Qe_i  * ileum3_e/ACAT_v_ent_ileum3
  
  ### equations for the sink compartments & model output
  
  # equations for the sink compartments
  d/dt(sink_s) <- output_s_i3
  d/dt(sink_d) <- output_d_i3
  
  
  
  ### PBPK model ------------------------------------------------------------------------------------------------
  
  # preallocation of some support variable...
  bf_tot      <- 0
  outflow_tot <- 0
  outflow_spl <- 0
  
  ### parallel tissues equations...  c("brain","heart","kidneys","bone","muscle","skin","fat")
  
  # brain
  brain_outflow <- organ_bf_brain * (brain/organ_v_brain)/(part_coeff_brain/BP)
  d/dt(brain)   <- organ_bf_brain * arterial_blood/organ_v_arterial_blood - brain_outflow
  bf_tot        <- bf_tot + organ_bf_brain
  outflow_tot   <- outflow_tot + brain_outflow
  
  # heart
  heart_outflow <- organ_bf_heart * (heart/organ_v_heart)/(part_coeff_heart/BP)
  d/dt(heart)   <- organ_bf_heart * arterial_blood/organ_v_arterial_blood - heart_outflow
  bf_tot        <- bf_tot + organ_bf_heart
  outflow_tot   <- outflow_tot + heart_outflow
  
  # kidneys
  clear_kid       <- CLr * fup * (kidneys/organ_v_kidneys)/(part_coeff_kidneys) + general_p_GFR * fup * arterial_blood/organ_v_arterial_blood/BP * GFR_flag
  kidneys_outflow <- organ_bf_kidneys * (kidneys/organ_v_kidneys)/(part_coeff_kidneys/BP)
  d/dt(kidneys)   <- organ_bf_kidneys * arterial_blood/organ_v_arterial_blood - kidneys_outflow - clear_kid
  bf_tot          <- bf_tot + organ_bf_kidneys
  outflow_tot     <- outflow_tot + kidneys_outflow
  
  # bone
  bone_outflow  <- organ_bf_bone * (bone/organ_v_bone)/(part_coeff_bone/BP)
  d/dt(bone)    <- organ_bf_bone * arterial_blood/organ_v_arterial_blood - bone_outflow
  bf_tot        <- bf_tot + organ_bf_bone
  outflow_tot   <- outflow_tot + bone_outflow
  
  # muscle
  muscle_outflow <- organ_bf_muscle * (muscle/organ_v_muscle)/(part_coeff_muscle/BP)
  d/dt(muscle)   <- organ_bf_muscle * arterial_blood/organ_v_arterial_blood - muscle_outflow
  bf_tot         <- bf_tot + organ_bf_muscle
  outflow_tot    <- outflow_tot + muscle_outflow
  
  # skin
  skin_outflow <- organ_bf_skin * (skin/organ_v_skin)/(part_coeff_skin/BP)
  d/dt(skin)   <- organ_bf_skin * arterial_blood/organ_v_arterial_blood - skin_outflow
  bf_tot       <- bf_tot + organ_bf_skin
  outflow_tot  <- outflow_tot + skin_outflow
  
  # fat
  fat_outflow <- organ_bf_fat * (fat/organ_v_fat)/(part_coeff_fat/BP)
  d/dt(fat)   <- organ_bf_fat * arterial_blood/organ_v_arterial_blood - fat_outflow
  bf_tot      <- bf_tot + organ_bf_fat
  outflow_tot <- outflow_tot + fat_outflow
  
  
  ### splanchnic organs equations c("stomach", "spleen", "pancreas", "gut")
  
  # stomach
  stomach_outflow <- organ_bf_stomach * (stomach/organ_v_stomach)/(part_coeff_stomach/BP)
  d/dt(stomach)   <- organ_bf_stomach * arterial_blood/organ_v_arterial_blood - stomach_outflow
  bf_tot          <- bf_tot + organ_bf_stomach
  outflow_spl     <- outflow_spl + stomach_outflow
  
  # spleen
  spleen_outflow <- organ_bf_spleen * (spleen/organ_v_spleen)/(part_coeff_spleen/BP)
  d/dt(spleen)   <- organ_bf_spleen * arterial_blood/organ_v_arterial_blood - spleen_outflow
  bf_tot          <- bf_tot + organ_bf_spleen
  outflow_spl     <- outflow_spl + spleen_outflow
  
  # pancreas
  pancreas_outflow <- organ_bf_pancreas * (pancreas/organ_v_pancreas)/(part_coeff_pancreas/BP)
  d/dt(pancreas)   <- organ_bf_pancreas * arterial_blood/organ_v_arterial_blood - pancreas_outflow
  bf_tot           <- bf_tot + organ_bf_pancreas
  outflow_spl      <- outflow_spl + pancreas_outflow
  
  # gut
  gut_outflow <- organ_bf_gut * (gut/organ_v_gut)/(part_coeff_gut/BP)
  d/dt(gut)   <- organ_bf_gut * arterial_blood/organ_v_arterial_blood - gut_outflow
  bf_tot      <- bf_tot + organ_bf_gut
  outflow_spl <- outflow_spl + gut_outflow
  
  
  ### liver equations
  organ_bf_liv_tot <- organ_bf_liver + organ_bf_gut + organ_bf_pancreas + organ_bf_spleen + organ_bf_stomach
  outflow_liv      <- organ_bf_liv_tot * (liver/organ_v_liver)/(part_coeff_liver/BP)
  clear_liv        <- CLh * fup * (liver/organ_v_liver)/(part_coeff_liver)
  d/dt(liver)      <- outflow_spl + organ_bf_liver * arterial_blood/organ_v_arterial_blood - outflow_liv - clear_liv + input_GI
  bf_tot           <- bf_tot + organ_bf_liver
  outflow_tot      <- outflow_tot + outflow_liv
  
  ### venous blood equation
  outflow_ven <- bf_tot * venous_blood/organ_v_venous_blood
  d/dt(venous_blood) <- outflow_tot - outflow_ven
  
  ### lung equation
  outflow_lung <- bf_tot * ( lungs/organ_v_lungs ) / ( part_coeff_lungs/BP )
  d/dt(lungs) <- outflow_ven - outflow_lung
  
  # arterial blood equation
  d/dt(arterial_blood) <- outflow_lung - bf_tot * arterial_blood/organ_v_arterial_blood
  
  # sink compartments
  d/dt(sink_CLh) <- clear_liv # liver sink
  d/dt(sink_CLr) <- clear_kid # kidney sink
  
  plasma_conc <- venous_blood/BP/organ_v_venous_blood
  
  
})








