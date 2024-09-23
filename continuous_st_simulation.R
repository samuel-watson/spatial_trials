# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
require(Matrix)
source(paste0(getwd(),"/solarized.R"))

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- TRUE # save newly created data

# load required functions
Rcpp::sourceCpp("src/perm_test.cpp")

#SIMULATION PARAMETERS

# Upper bound on Del_E 
max_de <- function(n_locs){
  sqrt(8/(n_locs*3*sqrt(3)))
}

# set simulation parameters
# if data already exists with these parameters and use_data == TRUE then it will instead be loaded
n_seed <- 10
n_child <- 100
nT <- 3 # number of time periods
cov_pars <- c(0.25,0.5) # G.P. variance, length scale (distance is Euclidean (x,y,t))

# GENERATE BASE DATA INCL. SAMPLE POINTS AND LATENT SURFACE

# if prior data exists
if(use_data){
  if(!file.exists(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_spt.RDS"))){
    stop("Data for this combination of parameters does not exist")
  } else {
    dfp <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_spt.RDS"))
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
  dfp$idx <- 1:nrow(dfp)
  
  dfp_t <- dfp
  dfp$t <- 1
  for(t in 2:nT){
    dfp_t$t <- t
    dfp <- rbind(dfp,dfp_t)
  }
  dfp$t_z <- (dfp$t - 1)/nT
  rm(dfp_t)
  for(t in 1:nT){
    dfp$tmp <- I(dfp$t == t)*1
    colnames(dfp)[ncol(dfp)] <- paste0("t",t)
  }
  
  mod <- Model$new(
    ~ (1|fexp(X,Y,t)),
    data = as.data.frame(dfp)[,c("X","Y","t")],
    covariance = cov_pars,
    mean = c(0),
    family = gaussian()
  )
  
  sim_data <- mod$sim_data(type="all")
  dfp$u <- sim_data$u
  rm(mod)
  
  
  if(save_data) saveRDS(dfp,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_spt.RDS"))
}

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp[dfp$t == 1, ], size = 0.5, alpha = 0.5)+
  theme_bw()+
  ggtitle("Simulated locations"); p_loc

# function to generate new

all_dists <- st_distance(dfp[dfp$t==1, ])

generate_intervention <- function(data, max_dist, beta, n_locs, rho_delta = 0.4,
                                  rho_beta = 0.8, plot = TRUE){
  
  # spatially regulated sampling scheme
  
  nT <- length(unique(data$t))
  data_t1 <- data[data$t == 1, ]
  
  int_idx <- sample(data_t1$idx,1)
  while(length(int_idx) < n_locs){
    int_idx_new <- sample(data_t1$idx,1)
    dists <- c(all_dists[int_idx,int_idx_new]) #c(st_distance(sampled_locs,dfp[int_idx_new,]))
    if(min(dists) > max_dist*1.5){
      int_idx <- c(int_idx, int_idx_new)
    } 
  }
  
  data$intervention <- 0
  data[data$idx %in% int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- rep(apply(all_dists[,int_idx],1,min), nT)
  data$fn <- NA
  data$y_true <- NA
  del_values <- c()
  
  for(j in 1:nT){
    del_values <- c(del_values, c(max_dist)*rho_delta^(1-j/nT))
  }
  
  data$fn <- fun(c(data$distance),50,4,8,del_values,nT,data$t)
  
  for(j in 1:nT){
    data[data$t == j, 'y_true'] <- data[data$t == j,]$fn * beta * rho_beta^(1-j/nT)
  }
  data$sim_p <- data$y_true + data$u
  data$sim_y <- rnorm(nrow(data),data$sim_p)
  
  
  if(plot){
    p_dist <- ggplot()+
      geom_sf(data=data[data$t == 1, ], aes(color = distance), size = 0.1)+
      geom_sf(data=data[data$intervention==1&data$t == 1,],color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "Distance")+
      theme_solar()+
      ggtitle("Distance")
    
    p_int <- ggplot()+
      geom_sf(data=data, aes(color = y_true), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      facet_wrap(~t)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Intervention effect")
    
    p_u <- ggplot()+
      geom_sf(data=data, aes(color = u), size = 0.1)+
      facet_wrap(~t)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Latent spatial effect")
    
    p_p <- ggplot()+
      geom_sf(data=data, aes(color = sim_p), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      facet_wrap(~t)+
      scico::scale_color_scico(palette = "roma", name = "Simulated\nprobability")+
      theme_solar()+
      ggtitle("Probability of outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

# test function & visualise
dfp <- generate_intervention(dfp, 0.35, -0.6, 8, plot = TRUE)


##########################################################
########### PERMUTATION TEST #############################
## TEST CPP
dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]

# simulation parameter values
beta <- 0
del_e <- 0.35
n_locs <- 8
max_del <- 0.5

model2b <- Model$new(
  ~ (1|fexp(X,Y,t)),
  data=dfanal,
  covariance = cov_pars,
  mean = 0,
  family = gaussian()
)
S <- model2b$Sigma()
E <- eigen(S)
B <- E$vectors%*%diag(1/sqrt(E$values))
rm(model2b)

# CREATE FORMULAE FOR MODELS

form <- "~ ("
for(t in 1:nT){
  form <- paste0(form,"t",t,"*b_",t)
  if(t < nT)form <- paste0(form,"+")
} 
form <- paste0(form,")*((1-(-0.02*log(exp(-50*(distance/(")
for(t in 1:nT){
  form <- paste0(form,"t",t,"*d_",t)
  if(t < nT)form <- paste0(form,"+")
} 
form <- paste0(form, ")))+exp(-50)))^(4))^(8))+(1|hsgp_fexp(X,Y,t))")

form2 <- "~ ("
for(t in 1:nT){
  form2 <- paste0(form2,"t",t,"*b_",t)
  if(t < nT)form2 <- paste0(form2,"+")
} 
form2 <- paste0(form2,")*fn+(1|hsgp_fexp(X,Y,t))")


pvals <- c()
pvals_ml <- c()
pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),t(B),dfanal$distance,all_dists,rep(0.01,nT),rep(max_del,nT),c(0),c(0),50,4,8,dfanal$t)
for(zz in 1:500){
  cat("\nITER: ",zz,"\n")
  dfp <- generate_intervention(dfp,0.3, 0.0, n_locs, plot=FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  update_y(pt,dfanal$sim_y-mean(dfanal$sim_y),dfanal$distance)
  
  # # maximum likelihood model
  
  # This is very slow - one iteration including the model fitting even with the approximation to the GP takes between 3 and 10 minutes.
  
  # model2 <- Model$new(
  #   as.formula(form),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01,rep(-0.01,nT), rep(0.1,nT)),
  #   family = gaussian()
  # )
  # 
  # # model2$set_trace(1)
  # model2$covariance$hsgp(m = c(15,15,3), L = c(1.2,1.2,1.2))
  # model2$update_parameters(cov.pars = cov_pars)
  # 
  # fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
  #                              lower.bound = c(-10,rep(-10,nT),rep(0.01,nT)),
  #                              upper.bound = c(10,rep(10,nT),rep(10,nT))),
  #                  error = function(i)return(list()))
  # 
  # model2b <- Model$new(
  #   ~ (1|hsgp_fexp(X,Y,t)),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01),
  #   family = gaussian()
  # )
  # 
  # # model2$set_trace(1)
  # model2b$covariance$hsgp(m = c(15,15,3), L = c(1.2,1.2,1.2))
  # model2b$update_parameters(cov.pars = cov_pars)
  # 
  # fit2b <- tryCatch(model2b$MCML(y = dfp$sim_y),
  #                  error = function(i)return(list()))
  # 
  # ll1 <- model2$log_likelihood()
  # ll0 <- model2b$log_likelihood()
  # lr <- -2*(ll0 - ll1)
  # pvals_ml <- c(pvals_ml, pchisq(lr,df = 4,lower.tail = FALSE))
  
  pval_new <- permute_p_value_b(pt,n_locs,0.45,200,c(0.1,0.2,0.3))
  pvals <- c(pvals, pval_new)
  
  if(save_data & i %% 10){
    saveRDS(pvals,paste0(getwd(),"/results/pvals_b0_st.RDS"))
    saveRDS(pvals_ml,paste0(getwd(),"/results/pvals_ml_b0_st.RDS"))
  }
}

mean(pvals < 0.05)
mean(pvals_ml < 0.05)

# VISUALISE THE CORRELATION SURFACE

# didf <- expand.grid(d1 = seq(0.01,0.5,length.out = 20), d2 = seq(0.01,0.5,length.out = 20), d3 = seq(0.01,0.5,length.out = 3), R = NA)
# for(i in 1:nrow(didf)){
#   didf$R[i] <- get_r_stat(pt,c(didf$d1[i],didf$d2[i],didf$d3[i]))
# }
# optim_r(pt,c(didf$d1[i],didf$d2[i],didf$d3[i]))
# 
# ggplot(data= didf, aes(x = d1, y = d2, fill = R))+
#   facet_wrap(~d3)+
#   geom_tile()+
#   scale_fill_viridis_c()

###############################################################
# CONFIDENCE INTERVALS
###############################################################

dfci_df <- data.frame(iter = 1:1000, lower = NA, upper = NA, lower_ml = NA, upper_ml = NA, b= NA, bp = NA, pval = NA)
dfcid_df <- data.frame(iter = 1:1000, lower = NA, upper = NA, lower_ml = NA, upper_ml = NA, d = NA, dp =NA)

dfci <- dfcid <- list()
for(i in 1:nT){
  dfci[[i]] <- dfci_df
  dfcid[[i]] <- dfcid_df
}

beta <- -0.6
spde <- make_mesh(dfanal, xy_cols = c("X","Y"), cutoff = 0.05)

for(i in 1:nrow(dfci[[1]])){
  cat("\nITER: ",i,"\n")
  
  dfp <- generate_intervention(dfp, 0.3, beta, n_locs, plot = FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),t(B),dfanal$distance,all_dists,rep(0.01,nT),rep(max_del,nT),c(0),c(0),50,4,8,dfanal$t)
  
  # Max likelihood
  # This is very slow 
  
  # model2 <- Model$new(
  #   as.formula(form),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01,rep(-0.01,nT), rep(0.1,nT)),
  #   family = gaussian()
  # )
  # 
  # # model2$set_trace(1)
  # model2$covariance$hsgp(m = c(15,15,3), L = c(1.2,1.2,1.2))
  # model2$update_parameters(cov.pars = cov_pars)
  # 
  # fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
  #                              lower.bound = c(-10,rep(-10,nT),rep(0.01,nT)),
  #                              upper.bound = c(10,rep(10,nT),rep(10,nT))),
  #                  error = function(i)return(list()))
  # 
  # if(is(fit2,"mcml")){
  #   b0 <- c()
  #   d0 <- c()
  #   
  #   for(t in 1:nT){
  #     b0 <- c(b0, fit2$coefficients$est[1+t])
  #     d0 <- c(d0, fit2$coefficients$est[nT+1+t])
  #     dfci[[t]]$lower_ml[i] = fit2$coefficients$lower[1+t]
  #     dfci[[t]]$upper_ml[i] = fit2$coefficients$upper[1+t]
  #     dfcid[[t]]$lower_ml[i] = fit2$coefficients$lower[nT+1+t]
  #     dfcid[[t]]$upper_ml[i] = fit2$coefficients$upper[nT+1+t]
  #     dfci[[t]]$b[i] <- fit2$coefficients$est[1+t]
  #     dfcid[[t]]$d[i] <- fit2$coefficients$est[nT+1+t]
  #   }
  # }
  
  best_del <- optim_r(pt,c(0.1,0.2,0.3))
  for(t in 1:nT) dfcid[[t]]$dp[i] <- best_del[t]
  dfanal$fn <- fun(dfanal$distance,50,4,8,best_del[1:nT],nT,dfanal$t)
  
  fitTMB <- sdmTMB(sim_y ~ 0, data = dfanal, mesh = spde, time = "t", spatiotemporal = "ar1", time_varying = ~ fn)
  
  # model2b <- Model$new(
  #   as.formula(form2),
  #   data=dfanal,
  #   covariance = cov_pars,
  #   mean = c(0.01,rep(-0.1,nT)),
  #   family = gaussian()
  # )
  # 
  # model2b$set_trace(1)
  # model2b$covariance$hsgp(m = c(15,15,3), L = c(1.2,1.2,1.2))
  # model2b$update_parameters(cov.pars = cov_pars)
  # 
  # fit2b <- tryCatch(model2b$MCML(y = dfp$sim_y, max.iter = 10),
  #                   error = function(i)return(list()))
  b0p <- c()
  for(t in 1:nT){
    b0p <- c(b0p, fitTMB$last.par.best[1+t])#fit2b$coefficients$est[1+t])
    dfci[[t]]$bp[i] <- fitTMB$last.par.best[1+t] # fit2b$coefficients$est[1+t]
  }
  
  ub <- confint_b(pt,n_locs,0.45,250,b0p+0.2,b0p,best_del[1:nT],TRUE)
  lb <- confint_b(pt,n_locs,0.45,250,b0p-0.2,b0p,best_del[1:nT],TRUE)
  ld <- confint_del(pt,n_locs,0.45,250,rep(0.05,nT),b0p,best_del[1:nT],TRUE)
  ud <- confint_del(pt,n_locs,0.45,250,best_del[1:nT]+0.1,b0p,best_del[1:nT],TRUE)
  
  for(t in 1:nT){
    dfci[[t]]$upper[i] <- ub[t]
    dfci[[t]]$lower[i] <- lb[t]
    dfcid[[t]]$lower[i] <- ld[t]
    dfcid[[t]]$upper[i] <- ud[t]
  }
  #dfci[[1]]$pval[i] <- permute_p_value_b(pt,n_locs,0.3,100,best_del[1:nT])
  
  if(save_data & i %% 10 == 0){
    cat("\n\n####################### ST ######################## \n\n")
    saveRDS(dfci,paste0(getwd(),"/results/dfci_b",gsub("\\.|-","",as.character(beta)),"_st.RDS"))
    saveRDS(dfcid,paste0(getwd(),"/results/dfcid_b",gsub("\\.|-","",as.character(beta)),"_st.RDS"))
  }
}


