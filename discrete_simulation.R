# PACKAGES AND PRE-REQUISITES
require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
require(ggforce)
source("solarized.R")

Rcpp::sourceCpp("src/perm_test.cpp")
source("discrete_fn.R")

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- TRUE # save newly created data

#SIMULATION PARAMETERS
beta <- -0.3 # max absolute intervention effect
del_e <- max_dist <- 0.2
del_i <- 0.1 # ignore if using one-way model
n_locs <- 8 # number of intervention sites
radius <- 0.2 # radius of intervention areas
n_seed <- 10 # number of seed locations for observations
n_child <- 150 # number of children per seed
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
max_del <- max_de(n_locs) # upper bound on del_e
adjust <- TRUE # whether to adjust the analysis
oneway <- TRUE # whether a two-way or one-way spillover model, one-way is external only


# GENERATE BASE DATA INCL. SAMPLE POINTS AND LATENT SURFACE

# if prior data exists
if(use_data){
  if(!file.exists(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))){
    stop("Data for this combination of parameters does not exist")
  } else {
    dfp <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))
    dists_i <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_dists_sp.RDS"))
    dfi <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_dfi_sp.RDS"))
  }
} else {
  # create some fake locations
  # use a point process type method with children
  
  df <- data.frame(x = runif(n_seed,-1,1), y = runif(n_seed,-1,1))
  for(i in 1:n_seed){
    dfnew <- data.frame(
      x = rnorm(n_child,df$x[i],0.25),
      y = rnorm(n_child,df$y[i],0.25)
    )
    df <- rbind(df,dfnew)
  }
  df <- df[abs(df$x) < 1 & abs(df$y) < 1, ]
  
  df$t <- 1
  dfp <- rts2::create_points(df,pos_vars = c('x','y'), t_var = "t")
  dfp_coord <- as.data.frame(st_coordinates(dfp))
  dfp <- dfp[!duplicated(paste0(dfp_coord$X,dfp_coord$Y)),]
  dfp_coord <- as.data.frame(st_coordinates(dfp))
  dfp <- cbind(dfp,dfp_coord)
  dfp$t <- 1
  
  # function to generate new
  
  # generate a set of circular areas, calculate the potential distances for all
  # 1. generate potential centroids with min distance apart
  dfi <- data.frame(x = rep(NA,n_locs*2), y = rep(NA,n_locs*2))
  dfi[1,1] <- runif(1,-1,1)
  dfi[1,2] <- runif(1,-1,1)
  
  while(any(is.na(dfi$x))){
    new_i <- runif(2,-1,1)
    min_dist <- 10
    for(i in 1:nrow(dfi[!is.na(dfi$x),])){
      new_min_dist <- sqrt((new_i[1] - dfi$x[i])^2 + (new_i[2] - dfi$y[i])^2)
      if(new_min_dist < min_dist) min_dist <- new_min_dist
    }
    if(min_dist > max_dist){
      dfi[is.na(dfi$x),][1,1] <- new_i[1]
      dfi[is.na(dfi$y),][1,2] <- new_i[2]
    }
  }
  
  # 2. create potential distances
  df_coords <- st_coordinates(dfp)
  dists_i <- matrix(NA,nrow=nrow(dfp),ncol=nrow(dfi))
  
  for(i in 1:nrow(dfp)){
    for(j in 1:nrow(dfi)){
      dists_i[i,j] <- sqrt((df_coords[i,1] - dfi$x[j])^2 + (df_coords[i,2] - dfi$y[j])^2) - radius
    }
  }
  
  dfp$distance_potential <- apply(dists_i,1,min)
  
  #generate indicators
  dfp$d1 <- 0
  dfp$d1[dfp$distance_potential >= -0.2 & dfp$distance_potential < -0.1] <- 1
  dfp$d2 <- 0
  dfp$d2[dfp$distance_potential >= -0.1 & dfp$distance_potential < 0.0] <- 1
  dfp$d3 <- 0
  dfp$d3[dfp$distance_potential >= 0.0 & dfp$distance_potential < 0.1] <- 1
  dfp$d4 <- 0
  dfp$d4[dfp$distance_potential >= 0.1 & dfp$distance_potential < 0.2] <- 1
  
  
  if(save_data){
    saveRDS(dfp,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))
    saveRDS(dists_i,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_dists_sp.RDS"))
    saveRDS(dfi,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_dfi_sp.RDS"))
  } 
}

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  theme_bw()+
  ggtitle("Simulated locations"); p_loc


# GENERATE RELEVANT MATRICES 

if(adjust){
  model2b <- Model$new(
    ~ d2 + (1|fexp(X,Y)),
    data=as.data.frame(dfp),
    covariance = cov_pars,
    family = gaussian()
  )
  
  form <- ifelse(oneway,"~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*((distance)/(del_e))) + exp((-25)*(1+sign0(distance))))))^(4))^(8)) + d2","~ twoway2(distance,8,4,50) + d2")
  form2 <- "sim_y ~ fn + d2"
} else {
  model2b <- Model$new(
    ~ (1|fexp(X,Y)),
    data=as.data.frame(dfp),
    covariance = cov_pars,
    mean = 0.01,
    family = gaussian()
  )
  
  form <- ifelse(oneway,"~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*((distance)/(del_e))) + exp((-25)*(1+sign0(distance))))))^(4))^(8))","~ twoway2(distance,8,4,50)")
  form2 <- "sim_y ~ fn"
}

# this uses the full exponential GP model, so is a little slow, but only needs to be run once
S <- model2b$Sigma()
L <- t(chol(S))
Li <- solve(L)
rm(model2b)

#test
dfp <- generate_intervention(dfp, del_e, ifelse(oneway,0,del_i), beta, n_locs, dfi, TRUE)
dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]

#####################################################
# PERMUTATION TEST
####################################################

beta <- 0
pvals <- c()
pvals_ml <- c()

for(zz in 1:1000){
  cat("\nITER: ",zz,"\n")
  dfp <- generate_intervention(dfp, del_e, del_i, beta, n_locs, dfi, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),solve(L),dfanal$distance,dists_i,c(0.01),c(0.4),ifelse(oneway,0.0,-0.2),c(0.0),50,4,8,dfanal$t)
  
  # starting values
  mean_pars <- c(0.01,-0.01, 0.10)
  if(!oneway) mean_pars <- c(mean_pars, 0.10)
  if(adjust) mean_pars <- c(mean_pars,rep(0.05,1))
  
  # # maximum likelihood model
  model2 <- Model$new(
    as.formula(paste0(form," + (1|hsgp_fexp(X,Y))")),
    data=dfanal,
    covariance = cov_pars,
    mean = mean_pars,
    family = gaussian()
  )

  # model2$set_trace(1)
  model2$covariance$hsgp(m = c(10,10), L = c(1.05,1.05))
  model2$update_parameters(cov.pars = cov_pars)

  # bounds
  lbound <- c(-10,-10)
  ubound <- c(10,10)
  if(!oneway){
    lbound <- c(lbound, 0.01,0.01)
    ubound <- c(ubound, 10.0, 1.0)
  } else {
    lbound <- c(lbound, 0.01)
    ubound <- c(ubound, 1.0)
  }
  if(adjust){
    lbound <- c(lbound, rep(-10,1))
    ubound <- c(ubound, rep(10,1))
  }

  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
                               reml = FALSE, # not compatible with HSGP approximation
                               lower.bound = lbound,
                               upper.bound = ubound),
                   error = function(i)return(list()))
  
  model2b <-  Model$new(
    as.formula(paste0(form," + (1|fexp(X,Y))")),
    data=dfanal,
    covariance = model2$covariance$parameters, 
    mean = model2$mean$parameters,
    family = gaussian()
  )
  M <- model2b$information_matrix()
  se <- sqrt(diag(solve(M)))
  
  pvals_ml <- c(pvals_ml, 2*(1-pnorm(abs(fit2$coefficients$est[2]/se[2]))))

  pval_new <- permute_p_value_b(pt,n_locs,0.3,200,c(0.2,0.1))
  pvals <- c(pvals, pval_new)
  rm(pt)
  
  if(save_data & zz %% 10){
    saveRDS(pvals,paste0(getwd(),"/results/pvals_b0",ifelse(adjust,"adj",""),"_dis.RDS"))
    saveRDS(pvals_ml,paste0(getwd(),"/results/pvals_ml_b0",ifelse(adjust,"adj",""),"_dis.RDS"))
  }
}

mean(pvals < 0.05)

############################################
#  CONFIDENCE INTERVALS
##############################################

dfci <- data.frame(iter = 1:1000, lower = NA, upper = NA, lower2 = NA, upper2 = NA, lower_ml = NA, upper_ml = NA, b= NA, bp = NA, pval = NA)
dfcid <- data.frame(iter = 1:1000, lower_e = NA, upper_e = NA, lower_e2 = NA, upper_e2 = NA, lower_ml_e = NA, upper_ml_e = NA, lower_i = NA, upper_i = NA, lower_i2 = NA, upper_i2 = NA, lower_ml_i = NA, upper_ml_i = NA, d_e = NA, dp_e =NA, d_i = NA, dp_i =NA)

cl <- parallel::makeCluster(6)
parallel::clusterExport(cl,c('fn2','L','genrep'))

# reset beta

beta <- -0.3

for(i in 1:500){
  cat("\nITER: ",i,"\n")
  dfp <- generate_intervention(dfp, del_e, ifelse(oneway,0,del_i), beta, n_locs, dfi, i%%10 == 0)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),as.matrix(Li),dfanal$distance,dists_i,
                   c(0.01),c(0.4),ifelse(oneway,0.0,-2*del_i),c(0.0),50,4,8,dfanal$t)
  
  #starting values
  mean_pars <- c(0.01,-0.3, 0.10)
  if(!oneway) mean_pars <- c(mean_pars, c(0.10))
  if(adjust) mean_pars <- c(mean_pars,rep(0.05,1))
  
  # # maximum likelihood model  
  model2 <- Model$new(
    as.formula(paste0(form, "+ (1|hsgp_fexp(X,Y))")),
    data=dfanal,
    covariance = cov_pars,
    mean = mean_pars,
    family = gaussian()
  )
  
  # model2 <- Model$new(
  #   ~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*(((distance) + exp(del_i))/(exp(del_e) + exp(del_i)))) + exp((-25)*(1+sign0(distance))))))^(4))^(8)) + b1*d1 + (1|hsgp_fexp(X,Y)),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01,0.1,log(0.1),log(0.1),0.0),
  #   family = gaussian()
  # )
  
  # model2 <- Model$new(
  #   ~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*((distance)/(del_e))) + exp((-25)*(1+sign0(distance))))))^(4))^(8)) + (1|hsgp_fexp(X,Y)),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01,-0.1,0.1),
  #   family = gaussian()
  # )
  
  # model2$mcmc_options$samps <- 200
  
  model2$covariance$hsgp(m = c(10,10), L = c(1.05,1.05))
  model2$update_parameters(cov.pars = cov_pars)

  # bounds
  lbound <- c(-10,-10)
  ubound <- c(10,10)
  if(!oneway){
    lbound <- c(lbound, 1e-4,1e-4)
    ubound <- c(ubound, 10.0, 1.0)
  } else {
    lbound <- c(lbound, 1e-4)
    ubound <- c(ubound, 1.0)
  }
  if(adjust){
    lbound <- c(lbound, rep(-10,1))
    ubound <- c(ubound, rep(10,1))
  }

  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
                               reml = FALSE,
                               max.iter = 10,
                               lower.bound = lbound,
                               upper.bound = ubound),
                   error = function(i)return(list()))

  if(is(fit2,"mcml")){
    b0 <- fit2$coefficients$est[2]
    if(oneway){
      di0 <- 0
      de0 <- fit2$coefficients$est[3]
    } else {
      di0 <- fit2$coefficients$est[3]
      de0 <- fit2$coefficients$est[4]
    }

    model2b <- Model$new(
      as.formula(paste0(form," + (1|fexp(X,Y))")),
      data=dfanal,
      covariance = model2$covariance$parameters,
      mean = model2$mean$parameters,
      family = gaussian()
    )
    model2b$update_y(dfanal$sim_y)
    M <- model2b$information_matrix()
    se <- sqrt(diag(solve(M)))

    dfci$lower_ml[i] = fit2$coefficients$est[2] - qnorm(0.975)*se[2]
    dfci$upper_ml[i] = fit2$coefficients$est[2] + qnorm(0.975)*se[2]
    if(oneway){
      dfcid$lower_ml_e[i] = fit2$coefficients$est[3] -qnorm(0.975)*se[3]
      dfcid$upper_ml_e[i] = fit2$coefficients$est[3] +qnorm(0.975)*se[3]
    } else {
      dfcid$lower_ml_i[i] = fit2$coefficients$est[3] -qnorm(0.975)*se[3]
      dfcid$upper_ml_i[i] = fit2$coefficients$est[3] +qnorm(0.975)*se[3]
      dfcid$lower_ml_e[i] = fit2$coefficients$est[4] -qnorm(0.975)*se[4]
      dfcid$upper_ml_e[i] = fit2$coefficients$est[4] +qnorm(0.975)*se[4]
    }
    
    f1 <- model2$fitted()
    
    if(oneway){
      fitn <- tryCatch(nls(sim_y ~ fn2(distance, d2, int,b, del_e, 0, b_d1),
                           data = dfanal, 
                           start = list(int = 0, b = -0.2, del_e = 0.2, b_d1 = 0.0),
                           lower = c(-10,-3,0.0,-10), upper = c(10,3,1.0,1.0,10), algorithm = "port"), error = function(i)return(NA))
      
    } else {
      fitn <- tryCatch(nls(sim_y ~ fn2(distance, d2, int,b, del_e, del_i, b_d1),
                           data = dfanal, 
                           start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_d1 = 0.0),
                           lower = c(-10,-3,0.01,0.0,-10), upper = c(10,3,1.0,1.0,10), algorithm = "port"), error = function(i)return(NA))
    }
    
    if(is(fitn,"nls")){
      #f1 <- fitted(fitn)
      np <- fitn$m$getPars()
      dfci$bp[i] <- np[2]
      dfcid$dp_e[i] <- np[3]
      dfcid$dp_i[i] <- np[4]
      
      parallel::clusterExport(cl,c('dfanal','f1','oneway'))
      res <- pbapply::pbreplicate(200, genrep(dfanal,f1,as.matrix(L),oneway), cl = cl)
      if(is(res,"list"))res <- t(Reduce(rbind,res))
      
      dfci$lower[i] = fit2$coefficients$est[2] - qnorm(0.975)*sd(res[2,],na.rm=TRUE)
      dfci$upper[i] = fit2$coefficients$est[2] + qnorm(0.975)*sd(res[2,],na.rm=TRUE)
      dfci$lower2[i] = np[2] - qnorm(0.975)*sd(res[2,],na.rm=TRUE)
      dfci$upper2[i] = np[2] + qnorm(0.975)*sd(res[2,],na.rm=TRUE)
      if(!oneway){
        dfcid$lower_i[i] = fit2$coefficients$est[3] - qnorm(0.975)*sd(res[4,],na.rm=TRUE)
        dfcid$upper_i[i] = fit2$coefficients$est[3] + qnorm(0.975)*sd(res[4,],na.rm=TRUE)
        dfcid$lower_e[i] = fit2$coefficients$est[4] - qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        dfcid$upper_e[i] = fit2$coefficients$est[4] + qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        
        dfcid$lower_i2[i] = np[3] - qnorm(0.975)*sd(res[4,],na.rm=TRUE)
        dfcid$upper_i2[i] = np[3] + qnorm(0.975)*sd(res[4,],na.rm=TRUE)
        dfcid$lower_e2[i] = np[4] - qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        dfcid$upper_e2[i] = np[4] + qnorm(0.975)*sd(res[3,],na.rm=TRUE)
      } else {
        dfcid$lower_e[i] = fit2$coefficients$est[3] - qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        dfcid$upper_e[i] = fit2$coefficients$est[3] + qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        dfcid$lower_e2[i] = np[3] - qnorm(0.975)*sd(res[3,],na.rm=TRUE)
        dfcid$upper_e2[i] = np[3] + qnorm(0.975)*sd(res[3,],na.rm=TRUE)
      }
      
    }
  }
    
  dfci$b[i] <- b0
  dfcid$d_i[i] <- di0
  dfcid$d_e[i] <- de0
  
  if(oneway){
    dfci$pval[i] <- tryCatch(permute_p_value_b(pt,n_locs,0.3,200,de0), error = function(i)return(NA))
  } else {
    dfci$pval[i] <- tryCatch(permute_p_value_b(pt,n_locs,0.3,200,c(de0,di0)), error = function(i)return(NA))
  }
  
  #rm(pt)
  if(save_data & i %% 10 == 0){
    saveRDS(dfci,paste0("dfci102_b",gsub("\\.|-","",as.character(beta)),ifelse(adjust,"adj",""),"_dis.RDS"))
    saveRDS(dfcid,paste0("dfcid102_b",gsub("\\.|-","",as.character(beta)),ifelse(adjust,"adj",""),"_dis.RDS"))
  }
}

## SUMMARISE RESULTS
mean(dfcid$lower_ml_e < del_e & dfcid$upper_ml_e > del_e, na.rm=TRUE)
mean(dfcid$lower_ml_i < del_i & dfcid$upper_ml_i > del_i, na.rm=TRUE)
mean(dfci$lower_ml < beta & dfci$upper_ml > beta, na.rm=TRUE)
mean(dfci$lower < beta & dfci$upper > beta, na.rm=TRUE)

mean(dfcid$lower_e < del_e & dfcid$upper_e > del_e, na.rm=TRUE)
mean(dfcid$lower_i < del_e & dfcid$upper_i > del_e, na.rm=TRUE)

# bias
mean(dfcid$d_e - del_e, na.rm=TRUE)
mean(dfcid$dp_e - del_e, na.rm=TRUE)

mean(dfcid$d_i - del_i, na.rm=TRUE)
