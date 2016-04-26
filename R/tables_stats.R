
get_variables <- function() {
  list(
    vars = c("lmlf_mst", "lalf_mst", "llma","lalf_astba2", "lastba2_mst"),
    labels = c("$M_F/M_S$","$A_F/M_S$", "$M_F/A_F$", "$A_F/A_S$", "$A_S/M_S$")
    )
}

# Independent effects analysis (IEA)
do_hierpart <- function(dep, indep, data, ...){
  
  df <- data[,c(dep,indep)]
  df <- df[complete.cases(df),]
  
  y <- df[,dep]
  m <- df[,indep]
  
  h <- hier.part(y, m, ...) 
  
  h$N <- nrow(df)
  
  return(h)
}

make_table_hierpart <- function(dataset){
  
  # interaction
  dataset$MATMAP <- with(dataset, MAT*MAP)
  
  vars <- get_variables()[["vars"]]
  indep_vars <- c("lh.t","pft","MAT","MAP","MATMAP")
  hp <- lapply(vars, function(x){
    do_hierpart(x, indep_vars, dataset, gof="Rsqu", barplot=FALSE)
  })
  
  hp_maketablerow <- function(x){
    iperc <- c(x$I.perc[c("lh.t","pft"),], sum(x$I.perc[c("MAP","MAT","MATMAP"),]))
    r2 <- max(x$gfs)
    
    ipercrel <- r2* iperc / 100
    
    c(ipercrel,r2)
  }
  
  tab <- as.data.frame(t(sapply(hp, hp_maketablerow)))
  
  tab <- cbind(get_variables()[["labels"]], tab)
  colnames(tab) <- c("Variable","$H_T$","PFT","Climate","$R^2$ total")
  tab <- as.data.frame(tab)
  return(tab)
}



# Explained variance with mixed-effects models with various combinations of predictors 
mixedr2 <- function(data, returnfit=FALSE){
  
  runmixmodels <- function(yvar, dat){
    
    f <- list()
    dat <- dat[!is.na(dat$lh.t) & !is.na(dat$pft) & !is.na(dat$MAP) & !is.na(dat$MAT),]
    
    f[[1]] <- as.formula(paste(yvar,"~ lh.t*I(lh.t^2) + (1|Group)"))
    f[[2]] <- as.formula(paste(yvar,"~ lh.t*I(lh.t^2)*pft + (1|Group)"))
    f[[3]] <- as.formula(paste(yvar,"~ lh.t*I(lh.t^2)*pft*MAP*MAT + (1|Group)"))
                               
    g <- lapply(f, function(x)lmer(formula=x, data=dat))
    
    return(g)
  }
  
  vars <- get_variables()[["vars"]]
  varlabel <- get_variables()[["labels"]]
  
  mods <- lapply(vars, runmixmodels, dat=data)
  if(returnfit){
    mods <- unlist(mods)
    indepvar <- c("H","H-PFT","H-PFT-CLIM")
    e <- expand.grid(indepvar, vars)
    names(mods) <- paste(e[[2]], e[[1]], sep="-")
    
    return(mods)
  }
  
  # r.squared.merMod
  xr2 <- function(x)unlist(sapply(x,  r.squaredGLMM)["R2m",])
  r2g <- as.data.frame(t(sapply(mods, xr2)))
  
  tabg <- cbind(as.data.frame(varlabel), r2g)
  names(tabg) <- c("Variable","H","H, PFT","H, PFT, MAP, MAT")
  
return(tabg)
}



# Explained variance with GAMs with various combinations of predictors
gamr2 <- function(data, ranef=FALSE, climvar1="MI", climvar2="mgdd0", kgam=4){

  testmapmatgam2 <- function(yvar){

    f <- list()

    f[[1]] <- as.formula(paste(yvar,"~ te(lh.t)"))
    f[[2]] <- as.formula(paste(yvar,"~ pft"))
    f[[3]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft)"))
    f[[4]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(",climvar1,", k=",kgam,")"))
    f[[5]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(",climvar2,", k=",kgam,")"))
    f[[6]] <- as.formula(paste(yvar,"~ pft + te(",climvar1,", k=",kgam,")"))
    f[[7]] <- as.formula(paste(yvar,"~ pft + te(",climvar2,", k=",kgam,")"))

    if(!ranef)
      g <- lapply(f, function(x)gam(formula=x, data=data))
    else
      g <- lapply(f, function(x)gamm(formula=x, random=list(Group=~1), data=data))

    return(g)
  }

  vars <- get_variables()[["vars"]]
  gams <- lapply(vars, testmapmatgam2)

  r2g <- do.call(rbind,lapply(1:length(vars),
                              function(i)unlist(sapply(gams[[i]],function(x){
                                if(ranef)summary(x$gam)$r.sq else summary(x)$r.sq
                              }))))

  tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
  names(tabg) <- c("Variable","H","PFT","H, PFT",paste0("H, PFT, ",climvar1),paste0("H, PFT, ",climvar2),
                   paste0("PFT, ",climvar1), paste0("PFT, ",climvar2))

  list(r2table=tabg, fits=gams)
}

# GAM explained variance with various combinations of predictors
# - this version is for SuppInfo methods comparisons.
gamr2old <- function(data, ranef=FALSE, climvar1="MI", climvar2="mgdd0", kgam=4){
  
  testmapmatgam2 <- function(yvar){
    
    f <- list()
    
    f[[1]] <- as.formula(paste(yvar,"~ te(lh.t)"))
    f[[2]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft)"))
    f[[3]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(",climvar1,", k=",kgam,
                               ") + te(",climvar2,", k=",kgam,")"))

    if(!ranef)
      g <- lapply(f, function(x)gam(formula=x, data=data))
    else
      g <- lapply(f, function(x)gamm(formula=x, random=list(Group=~1), data=data))
    
    return(g)
  }
  
  vars <- get_variables()[["vars"]]
  gams <- lapply(vars, testmapmatgam2)
  
  r2g <- do.call(rbind,lapply(1:length(vars),
                              function(i)unlist(sapply(gams[[i]],function(x){
                                if(ranef)summary(x$gam)$r.sq else summary(x)$r.sq
                              }))))
  
  tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
  names(tabg) <- c("Variable","H","H, PFT",paste0("H, PFT, ",climvar1,climvar2))
  
  list(r2table=tabg, fits=gams)
}

make_table_gamr2MATMAP_old <- function(dataset) {
  g0 <- gamr2old(dataset, kgam=4, climvar1="MAT", climvar2="MAP")$r2table
  g0$Variable <- get_variables()[["labels"]]
  
  g0
}


make_table_gamr2MATARID <- function(dataset) {
  g0 <- gamr2(dataset, kgam=4, ranef=TRUE, climvar1="MAT", climvar2="aridity")
  
  r2table <- g0$r2table
  r2table$Variable <- get_variables()[["labels"]]

  r2table
}



# Count number of observations.
tabFun <- function(x, vars, pftvar="pft"){
  
  x <- x[complete.cases(x[,vars]),]
  x$pft <- x[,pftvar]
  
  xt <- addmargins(xtabs( ~ pft, data=x))
  names(dimnames(xt)) <- pftvar
  
  x <- x[!duplicated(x$species,x$pft),]
  xs <- addmargins(xtabs( ~ pft, data=x))
  
  m <- paste0(xt, " (", xs, ")")
  
  m[m == "0 (0)"] <- NA
  
  return(m)
}

make_samplesize_table <- function(dataset){
  tb <- t(sapply(get_variables()[["vars"]],function(x)tabFun(dataset,x)))
  rownames(tb) <- get_variables()[["labels"]]
  colnames(tb) <- c("Deciduous Angiosperm","Evergreen Angiosperm","Evergreen Gymnosperm","Total")
  tb
}




# For SuppInfo Figs. S5 and S6
reichstyle_climate_pvals <- function(dataset){

  dataset$lhtclass <- quantcut(dataset$lh.t, 5)
  dataset$pft2 <- ifelse(dataset$pft == "EG", "Gymnosperm", "Angiosperm")
  
  lmesp <- function(PFT){

    df <- subset(dataset, pft2 == PFT)
    dfsp <- split(df, df$lhtclass)
    
    
    pv <- c()
    for(i in 1:length(dfsp)){
      
      assign("dat", dfsp[[i]], envir=.GlobalEnv)
      fit <- lme(lmlf_mst ~ MAT, random=~1|Group, data=dat, na.action=na.omit)
      pv[i] <- Anova(fit, test.statistic="F")[[3]]
    }
    
    return(pv)
  }
  
  data.frame(lhtclass=levels(dataset$lhtclass),
             pval_angio = lmesp("Angiosperm"),
             pval_gymno = lmesp("Gymnosperm"))

}

get_gamr2 <- function(x)summary(x$gam)$r.sq

af_as_stat <- function(dataset){
  m0 <- gamm(lalf ~ te(a.stba2), random=list(Group=~1), data=dataset)
  m1 <- gamm(lalf ~ pft + te(a.stba2, by=pft), random=list(Group=~1), data=dataset)
  lik <- anova(m0$lme, m1$lme)
  list(LRT=lik, R2a=get_gamr2(m0), R2b=get_gamr2(m1))
}

ms_as_stat <- function(dataset){
  m0 <- gamm(lmst ~ te(a.stba2), random=list(Group=~1), data=dataset)
  m1 <- gamm(lmst ~ pft + te(a.stba2, by=pft), random=list(Group=~1), data=dataset)
  lik <- anova(m0$lme, m1$lme)
  list(LRT=lik, R2a=get_gamr2(m0), R2b=get_gamr2(m1))
}

# as discussed on http://stackoverflow.com/questions/14530770/calculating-r2-for-a-nonlinear-model
nlsr2 <- function(mod){
  cor(fitted(mod), resid(mod)+fitted(mod))^2
}





