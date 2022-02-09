### plot functions


### PBPK states ---------------------------------------------------

## PBPK distribution
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

## PBPK sink
# 16 sink liver
# 17 sink kidney

## ACAT solid
# 18 stomach_s
# 19 duodenum_s
# 20 jejunum1_s
# 21 jejunum2_s
# 22 ileum1_s
# 23 ileum2_s
# 24 ileum3_s

## ACAT dissolved
# 25 stomach_d
# 26 duodenum_d
# 27 jejunum1_d
# 28 jejunum2_d
# 29 ileum1_d
# 30 ileum2_d
# 31 ileum3_d

## ACAT enterocytes
# 32 duodenum_e
# 33 jejunum1_e
# 34 jejunum2_e
# 35 ileum1_e
# 36 ileum2_e
# 37 ileum3_e

## ACAT_absorbed
# 38 duodenum_ea
# 39 jejunum1_ea
# 40 jejunum2_ea
# 41 ileum1_ea
# 42 ileum2_ea
# 43 ileum3_ea

## ACAT cleared
# 44 duodenum_ea
# 45 jejunum1_ea
# 46 jejunum2_ea
# 47 ileum1_ea
# 48 ileum2_ea
# 49 ileum3_ea

## ACAT sink
# 50 sink_s
# 51 sink_d

### other functions ------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### plot PBPK distribution model ---------------------------------------------------------------------------------
plotPBPK <- function(system.out.list, list.sim, names.sim, names.PBPK, names.ACAT, flag_log, data.PK, axis.limits){
  
  # initialize some useful variable outside the loop
  l.s <- length(list.sim)
  col <- gg_color_hue(l.s)
  
  lo.pbpk <- length(names.PBPK)
  df.pbpk.list <- list()
  p.pbpk       <- list()
  for(i in 1:lo.pbpk){
    df.pbpk.list[[i]] <- data.frame()
    p.pbpk[[i]]       <- ggplot()
  }
  
  lo.acat.1 <- length(names.ACAT)
  x.names.1 <- names.ACAT[c(1:20, 33, 34)]
  lo.acat.2 <- length(x.names.1)
  df.acat.list <- list()
  p.acat.1     <- list()
  for(i in 1:lo.acat.2){
    df.acat.list[[i]] <- data.frame()
    p.acat.1[[i]]     <- ggplot()
  }
  
  dfp <- data.frame()
  
  # derive dataframes for 
  for(jj in 1:l.s){
    
    system.out <- system.out.list[[jj]]
    
    x.pbpk   <- system.out[names.PBPK]
    n.points <- length(system.out$time) 
    
    # plot all the states dynamics
    for(i in 1:lo.pbpk){
      dfi.pbpk <- data.frame(time=as.vector(system.out$time), mass=as.vector( x.pbpk[[names.PBPK[i]]]), n.sim = rep(names.sim[jj],n.points))
      df.pbpk.list[[i]] <- rbind(df.pbpk.list[[i]], dfi.pbpk)
    }
    
    # plot venous blood
    c.plasma <- system.out$plasma_conc
    dfp_jj  <- data.frame(time=as.vector(system.out$time), mass=c.plasma, n.sim = rep(names.sim[jj],n.points))
    dfp     <- rbind(dfp, dfp_jj)
    
    ### plot ACAT distribution ---------------------------------------------------------------
    
    x.acat  <- system.out[names.ACAT]
    
    x.acat.1  <- x.acat[c(1:20, 33, 34)]
    
    for(i in 1:lo.acat.2){
      dfi.acat <- data.frame(time=as.vector(system.out$time), mass=as.vector( x.acat.1[[x.names.1[i]]]), n.sim = rep(names.sim[jj],n.points))
      df.acat.list[[i]] <- rbind(df.acat.list[[i]], dfi.acat)
    }
    
  }
  
  ### plot time-mass/conc
  
  # plasma concentration
  p.plasma <- ggplot() + geom_line(data=dfp, aes(x=time, y=mass, color=n.sim), lwd=1) + 
    labs(y='Conc [mg/L]') + 
    labs(x='Time (hours)')+ theme_bw(base_size=10) + 
    theme(text = element_text(size=16)) +
    ggtitle("plasma concentration")
   
  if(length(data.PK)!=0){
    p.plasma <- p.plasma + geom_point(data=data.PK, aes(x=time, y=PK, shape=study), color="gray38")
  }
  
  if(!is.na(axis.limits$xaxis[1])){
    p.plasma <- p.plasma + xlim(axis.limits$xaxis)
  }
  
  if(flag_log == 1){
    if(!is.na(axis.limits$yaxis[1])){
      axis.limits_log <- 10^(axis.limits$yaxis)
      p.plasma <- p.plasma + scale_y_log10(limits=axis.limits_log)
    }else{
      p.plasma <- p.plasma + scale_y_log10()
    }
  }else{
    if(!is.na(axis.limits$yaxis[1])){
      p.plasma <- p.plasma + ylim(axis.limits$yaxis)
    }
  }
  
  # PBPK states dynamics
  for(i in 1:lo.pbpk){
    dfi <- df.pbpk.list[[i]]
    p.pbpk[[i]] <- ggplot() + geom_line(data=dfi, aes(x=time, y=mass, color=n.sim)) + 
      labs(y='Mass (mg)') + 
      labs(x='Time (hours)')+ theme_bw(base_size=8) + 
      theme(legend.position="none", text = element_text(size=12)) + 
      ggtitle(names.PBPK[i])
    if(flag_log == 1){
      p.pbpk[[i]] <- p.pbpk[[i]] + scale_y_log10()
    }
    if(!is.na(axis.limits$xaxis[1])){
      p.pbpk[[i]] <- p.pbpk[[i]] + xlim(axis.limits$xaxis)
    }
  }
  #do.call("grid.arrange", c(p.pbpk, ncol=5, nrow=4))
  
  # ACAT states dynamics
  for(i in 1:lo.acat.2){
    dfi <- df.acat.list[[i]]
    p.acat.1[[i]] <- ggplot() + geom_line(data=dfi, aes(x=time, y=mass, color = n.sim)) + 
      labs(y='Mass (mg)') + 
      labs(x='Time (hours)')+ theme_bw(base_size=8) + 
      theme(legend.position="none", text = element_text(size=12)) + 
      ggtitle(x.names.1[i])
    if(!is.na(axis.limits$xaxis[1])){
      p.acat.1[[i]] <- p.acat.1[[i]] + xlim(axis.limits$xaxis)
    }
  }
  #do.call("grid.arrange", c(p.acat.1, ncol=6, nrow=4))
  
  ### plot fractions excreted & absorbed --------------------------------------------------------------
  
  names     <- c()
  fract     <- c()
  idx       <- c()
  ord       <- c()
  f_abs     <- c()
  names_abs <- c()
  idx_abs   <- c()
  ord_abs   <- c()
  
  # gut clearance
  gut_cl_name <- c("duodenum_ec",
                   "jejunum1_ec",
                   "jejunum2_ec",
                   "ileum1_ec",
                   "ileum2_ec",
                   "ileum3_ec")
  
  gut_abs_name <- c("duodenum_ea",
                   "jejunum1_ea",
                   "jejunum2_ea",
                   "ileum1_ea",
                   "ileum2_ea",
                   "ileum3_ea")
  
  
  # derive the vectors to make the dataframes
  for(jj in 1:l.s){
    
    ### excreted fraction
    system.out <- system.out.list[[jj]]
    
    dose_tot <- sum( as.vector(list.sim[[jj]]$ev$amt) * (as.vector(list.sim[[jj]]$ev$addl)+ 1), na.rm=TRUE)
    
    end.time <- length(system.out$time)
    f_s <- system.out$sink_s[end.time]/dose_tot
    f_d <- system.out$sink_d[end.time]/dose_tot
    f_h <- system.out$sink_CLh[end.time]/dose_tot
    f_r <- system.out$sink_CLr[end.time]/dose_tot
    f_g <- sum(system.out[end.time,gut_cl_name])/dose_tot
    f_tot <- f_s + f_d + f_h + f_r + f_g
    
    names <- c(names, "GI solid", "GI diss", "CL ent", "CL hep", "CL ren", "total")
    fract <- c(fract, f_s, f_d, f_g, f_h, f_r, f_tot)
    idx   <- c(idx, rep(names.sim[jj],6))
    ord   <- c(ord, 1:6)
    
    ### absorbed fraction
    f_abs_s   <- unlist(system.out[end.time, gut_abs_name]/dose_tot)
    f_abs_t   <- sum(f_abs_s)
    f_abs     <- c(f_abs, f_abs_s, f_abs_t)
    names_abs <- c("duodenum","jejunum1","jejunum2","ileum1","ileum2","ileum3","total")
    idx_abs   <- c(idx_abs, rep(names.sim[jj],7))
    ord_abs   <- c(ord_abs, 1:7)
  }
  
  ### plot fraction excreted
  df.bar <- data.frame(names = names, fract = fract, ord = ord, idx = idx)
  
  p.bar <- ggplot(data=df.bar, aes(x=reorder(names,ord), y=fract, fill=factor(idx))) +
    geom_bar(stat="identity", position="dodge")+
    theme(axis.text.x=element_text(angle=45, hjust=1), text = element_text(size=16))+
    ggtitle("fraction eliminated")+ylab("fraction")
  #p.bar
  
  ### plot fraction absorbed 
  df.abs <- data.frame(names = names_abs, val = f_abs, ord = ord_abs, idx=idx_abs)
  
  p.bar2 <- ggplot(data=df.abs, aes(x=reorder(names,ord), y=val, fill=factor(idx))) +
    geom_bar(stat="identity", position="dodge")+
    theme(axis.text.x=element_text(angle=45, hjust=1), text = element_text(size=16)) +
    ggtitle("fraction absorbed")+ylab("fraction")
  #p.bar2

  return(list(p.pbpk = p.pbpk, p.plasma=p.plasma, p.acat.1=p.acat.1, p.f.excr=p.bar, p.f.abs=p.bar2))
  
}
