
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
  
  vars <- c("lmlf_mst","lalf_astba2","llma","lastba2_mst")
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
  tab <- cbind(c("$M_F/M_S$","$A_F/A_S$","$M_F/A_F$","$M_S/A_S$"), tab)
  colnames(tab) <- c("Variable","$H_T$","PFT","Climate","$R^2$ total")
  tab <- as.data.frame(tab)
  return(tab)
}


make_r2_table <- function(data,  variable_name ){
  data.frame(Variable = c(variable_name, rep(NA, length(data)-1)),
    Predictors=names(data),
    r2=sapply(data,function(x) x[,"Marginal"]))
}


make_table_pipemodel_varpart <- function(dataset2) {

  # Datasets with NAs removed, of key variables.
  assign("dat_mlf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2))), envir = .GlobalEnv)
  assign("dat_alf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2))), envir = .GlobalEnv)

  mlf_formulas <- list()
  mlf_formulas[["H"]]         <- lmlf_astba2 ~ log10(h.t) + (1|Group)
  mlf_formulas[["H, PFT"]]    <- lmlf_astba2 ~ log10(h.t)*pft + (1|Group)
  mlf_formulas[["H, PFT, B"]] <- lmlf_astba2 ~ log10(h.t)*pft*bortemptrop + (1|Group)
  
  mlf_lmer_fits <- lapply(mlf_formulas, function(x) lmer(x, data=dat_mlf))
  mlf_lmer_r2 <- lapply( mlf_lmer_fits, r.squared.merMod)

  alf_formulas <- gsub("lmlf_astba2", "lalf_astba2", mlf_formulas)
  names(alf_formulas) <-   names(mlf_formulas)
  alf_lmer_fits <- lapply(alf_formulas, function(x) lmer(x, data=dat_alf))
  alf_lmer_r2 <- lapply( alf_lmer_fits, r.squared.merMod)

  rm(dat_mlf, dat_alf, envir = .GlobalEnv)
  rbind(make_r2_table(mlf_lmer_r2, "$M_F/A_S$"),
        make_r2_table(alf_lmer_r2, "$A_F/A_S$"))
}





make_table_lmaforpft <- function(dataset2){
  
  assign("dat_mlf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2))), envir = .GlobalEnv)
  assign("dat_alf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2))), envir = .GlobalEnv)

  
  model1 <- lmer(lmlf_astba2 ~ pft + log10(h.t) + pft:log10(h.t) + I(log10(h.t)^2):pft + (1|Group),
                 data=dat_mlf)
  model2 <- lmer(lmlf_astba2 ~ llma + log10(h.t) + llma:log10(h.t) + I(log10(h.t)^2):llma + (1|Group),
                 data=dat_mlf)
  
  model3 <- lmer(lmlf_mso ~ pft + log10(h.t) + pft:log10(h.t) + I(log10(h.t)^2):pft + (1|Group),
                 data=dat_mlf)
  model4 <- lmer(lmlf_mso ~ llma + log10(h.t) + llma:log10(h.t) + I(log10(h.t)^2):llma + (1|Group),
                 data=dat_mlf)
  
  
  r2 <- do.call(rbind,lapply(list(model1, model2, model3, model4),r.squared.merMod))
  rownames(r2) <- c("mfas_pft","mfas_lma","lmf_pft","lmf_lma")
return(r2)
}







make_table_LMFLAR_varpart <- function(dataset2) {

  assign("dat_mlf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2))), envir = .GlobalEnv)
  assign("dat_alf", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2))), envir = .GlobalEnv)

  assign("dat_mlfmso", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_mso))), envir = .GlobalEnv)
  assign("dat_alfmso", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso))), envir = .GlobalEnv)


  mlf_formulas <- list()
  mlf_formulas[["H"]]         <- log10(m.lf/m.so) ~ log10(h.t) + I(log10(h.t)^2)+ (1|Group)
  mlf_formulas[["H, PFT"]]    <- log10(m.lf/m.so) ~ log10(h.t)*pft + pft:I(log10(h.t)^2)+ (1|Group)
  mlf_formulas[["H, PFT, B"]] <- log10(m.lf/m.so) ~ pft*log10(h.t)*bortemptrop + pft:I(log10(h.t)^2) + (1|Group)

  mlf_lmer_fits <- lapply(mlf_formulas, function(x) lmer(x, data=dat_mlfmso))
  mlf_lmer_r2 <- lapply( mlf_lmer_fits, r.squared.merMod)

  alf_formulas <- gsub("m.lf/m.so", "a.lf/m.so", mlf_formulas)
  names(alf_formulas) <-   names(mlf_formulas)
  alf_lmer_fits <- lapply(alf_formulas, function(x) lmer(x, data=dat_alfmso))
  alf_lmer_r2 <- lapply( alf_lmer_fits, r.squared.merMod)

  rm(dat_mlfmso, dat_alfmso, envir = .GlobalEnv)
  rbind(make_r2_table( mlf_lmer_r2 ,"$M_F/M_T$"),
        make_r2_table(alf_lmer_r2  , "$A_F/M_T$"))
}


make_table_LMA_varpart <- function(dataset2) {
  
  assign("dat_lma", droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(llma))), envir = .GlobalEnv)
  
  mlf_formulas <- list()
  mlf_formulas[["H"]]         <- llma ~ log10(h.t) + I(log10(h.t)^2)+ (1|Group)
  mlf_formulas[["H, PFT"]]    <- llma ~ log10(h.t)*pft + pft:I(log10(h.t)^2)+ (1|Group)
  mlf_formulas[["H, PFT, B"]] <- llma ~ pft*log10(h.t)*bortemptrop + pft:I(log10(h.t)^2) + (1|Group)
  
  lma_lmer_fits <- lapply(mlf_formulas, function(x) lmer(x, data=dat_lma))
  lma_lmer_r2 <- lapply( lma_lmer_fits, r.squared.merMod)
  
  rm(dat_lma, envir = .GlobalEnv)
  df <- make_r2_table( lma_lmer_r2 ,"$M_F/A_F$")
        
return(df)
}


make_table_varpart2 <- function(dataset2){
  
  tab1 <- make_table_LMFLAR_varpart(dataset2)
  tab1$Variable <- as.character(tab1$Variable)
  tab2 <- make_table_pipemodel_varpart(dataset2)
  tab2$Variable <- as.character(tab2$Variable)
  tab3 <- make_table_LMA_varpart(dataset2)
  tab3$Variable <- as.character(tab3$Variable)
  
  
  # Reshape
  f <- function(x){
    
    if(any(duplicated(x$Predictors))){
      n <- nrow(x)/2
      x$Variable <- as.character(c(rep(x$Variable[1],n), rep(x$Variable[n+1],n)))
    } else {
      x$Variable <- as.character(rep(x$Variable[1],nrow(x)))
    }
    y <- reshape(x, direction="wide", timevar="Predictors", idvar="Variable")
    
  return(y)
  }
  nm <- rownames(tab1)[1:(nrow(tab1)/2)]
  
  df <- rbind(f(tab1), f(tab2), f(tab3))
  names(df)[2:ncol(df)] <- nm
  

  df <- cbind(data.frame(Description = c("Leaf mass fraction","Leaf area ratio",
                                         "Leaf mass / stem basal area","Leaf area / stem basal area",
                                         "Leaf mass per area")), df)
  
  
return(df)
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

make_table_mlfastba <- function(x) tabFun(x, c("lmlf_astba2","h.t"))
make_table_alfastba  <- function(x) tabFun(x, c("lalf_astba2","h.t"))
make_table_lmf  <- function(x) tabFun(x, c("lmlf_mso","h.t"))
make_table_lar  <- function(x) tabFun(x, c("lalf_mso","h.t"))


# #-----------------------------------------------------------------------------------------#

# # Test of effect of MAT and MAP on leaf - stem scaling

# # Data subset
# d <- droplevels(subset(dataset2, !is.na(m.st) & !is.na(m.lf) & !is.na(MAT)))

# mlfmst_lme0 <- lme(log10(m.lf) ~ log10(m.st)*pft, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)
# mlfmst_lme1 <- lme(log10(m.lf) ~ log10(m.st)*pft*MAT, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)
# mlfmst_lme2 <- lme(log10(m.lf) ~ log10(m.st)*pft*MAP, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)

# save(d, mlfmst_lme0,mlfmst_lme1,mlfmst_lme2,
#      file="manuscript/tables/Fits_lme_mlfmst_MAPMAT.RData")


# #------------------------------------------------------------------------------------#

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

  vars <- c("lmlf_mso","lalf_mso","lmlf_astba2","lalf_astba2","llma")
  gams <- lapply(vars, testmapmatgam2, mgdd0=TRUE)

  r2g <- do.call(rbind,lapply(1:length(vars),
                              function(i)unlist(sapply(gams[[i]],function(x){
                                if(ranef)summary(x$gam)$r.sq else summary(x)$r.sq
                              }))))

  tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
  names(tabg) <- c("Variable","H","H, PFT",paste0("H, PFT, ",climvar1,", ",climvar2))

  list(r2table=tabg, fits=gams)
}

make_table_gamr2MIgdd0 <- function(dataset) {
  g0 <- gamr2(dataset, kgam=4)$r2table
  g0$Variable <- c("$M_F/M_T$","$A_F/M_T$","$M_F/A_S$",
                   "$A_F/A_S$","$M_F/A_F$")
  
  
  g0 <- cbind(data.frame(Description = c("Leaf mass fraction","Leaf area ratio",
                                         "Leaf mass / stem basal area","Leaf area / stem basal area",
                                         "Leaf mass per area")), g0)
  g0
}

make_table_gamr2MATMAP <- function(dataset) {
  g0 <- gamr2(dataset, kgam=4, climvar1="MAT", climvar2="MAP")$r2table
  g0$Variable <- c("$M_F/M_T$","$A_F/M_T$","$M_F/A_S$",
                   "$A_F/A_S$","$M_F/A_F$")
  
  
  g0 <- cbind(data.frame(Description = c("Leaf mass fraction","Leaf area ratio",
                                         "Leaf mass / stem basal area","Leaf area / stem basal area",
                                         "Leaf mass per area")), g0)
  
  g0
}

