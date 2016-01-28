
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
  
  vars <- c("lmlf_mst","lalf_mst", "lalf_astba2","llma","lastba2_mst")
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
  
  tab <- cbind(c("$M_F/M_S$","$A_F/M_S$", "$A_F/A_S$","$M_F/A_F$","$M_S/A_S$"), tab)
  colnames(tab) <- c("Variable","$H_T$","PFT","Climate","$R^2$ total")
  tab <- as.data.frame(tab)
  return(tab)
}




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
  
  vars <- c("lmlf_mst","lalf_mst","lalf_astba2","llma","lastba2_mst")
  varlabel <- c("$M_F/M_S$","$A_F/M_S$","$A_F/A_S$","$M_F/A_F$","$A_S/M_S$")
  
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





# GAM explained variance with mgdd0 and MI
gamr2 <- function(data, ranef=FALSE, climvar1="MI", climvar2="mgdd0", kgam=4){

  testmapmatgam2 <- function(yvar, mgdd0=TRUE){

    f <- list()

    f[[1]] <- as.formula(paste(yvar,"~ te(lh.t)"))
    f[[2]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft)"))
    f[[3]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(",climvar1,", k=",kgam,")",
                               if(mgdd0)" + te(",climvar2,", k=",kgam,")"))

    if(!ranef)
      g <- lapply(f, function(x)gam(formula=x, data=data))
    else
      g <- lapply(f, function(x)gamm(formula=x, random=list(Group=~1), data=data))

    return(g)
  }

  vars <- c("lmlf_mst","lalf_mst","lalf_astba2","llma","lastba2_mst")
  gams <- lapply(vars, testmapmatgam2, mgdd0=TRUE)

  r2g <- do.call(rbind,lapply(1:length(vars),
                              function(i)unlist(sapply(gams[[i]],function(x){
                                if(ranef)summary(x$gam)$r.sq else summary(x)$r.sq
                              }))))

  tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
  names(tabg) <- c("Variable","H","H, PFT",paste0("H, PFT, ",climvar1,", ",climvar2))

  list(r2table=tabg, fits=gams)
}

make_table_gamr2MATMAP <- function(dataset) {
  g0 <- gamr2(dataset, kgam=4, climvar1="MAT", climvar2="MAP")$r2table
  g0$Variable <- c("$M_F/M_S$","$A_F/M_S$", "$A_F/A_S$","$M_F/A_F$","$M_S/A_S$")

  
  g0
}







# Count number of observations.
tabFun <- function(x, vars, pftvar="pft", vegvar="bortemptrop"){
  
  x <- x[complete.cases(x[,vars]),]
  x$pft <- x[,pftvar]
  x$veg <- x[,vegvar]
  
  xt <- addmargins(xtabs( ~ pft + veg, data=x))
  names(dimnames(xt)) <- c(pftvar,vegvar)
  
  x <- x[!duplicated(x$species,x$pft,x$veg),]
  xs <- addmargins(xtabs( ~ pft + veg, data=x))
  
  m <- matrix(paste0(xt, " (", xs, ")"), ncol=ncol(xt))
  dimnames(m) <- dimnames(xt)
  
  m[m == "0 (0)"] <- NA
  
  return(m)
}



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







