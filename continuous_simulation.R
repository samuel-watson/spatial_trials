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
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
misspec <- TRUE

# GENERATE BASE DATA INCL. SAMPLE POINTS AND LATENT SURFACE

# if prior data exists
if(use_data){
  if(!file.exists(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_sp.RDS"))){
    stop("Data for this combination of parameters does not exist")
  } else {
    dfp <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_sp.RDS"))
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
  
  if(save_data) saveRDS(dfp,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_sp.RDS"))
}

mod <- Model$new(
  ~ (1|fexp(X,Y)),
  data = as.data.frame(dfp)[,c("X","Y")],
  covariance = cov_pars,
  mean = c(0),
  family = gaussian()
)

L <- t(chol(mod$Sigma())) # get cholesky decomposition of coviarance matrix
rm(mod)

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  theme_bw()+
  ggtitle("Simulated locations"); p_loc

# generate distance matrix 

all_dists <- st_distance(dfp)

# DRF function

fn <- function(x,int,b,del){
  int + b*((1-((-1/50)*log(exp(-50*x/del) + exp(-50)))^4)^8)
}

generate_intervention <- function(data, max_dist, beta, n_locs, plot = TRUE){
  
  # spatially regulated sampling scheme
  iter <- 1
  int_idx <- sample(1:nrow(data),1)
  while(length(int_idx) < n_locs){
    int_idx_new <- sample(1:nrow(data),1)
    dists <- c(all_dists[int_idx,int_idx_new]) #c(st_distance(sampled_locs,dfp[int_idx_new,]))
    if(min(dists) > max_dist*0.01){ # change this line to implement a randomisation scheme where locations are spatially regulated
      int_idx <- c(int_idx, int_idx_new)
    } 
    iter <- iter + 1
    if(iter > 500) stop("Iterations exceed max")
  }
  
  data$intervention <- 0
  data[int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- apply(all_dists[,which(data$intervention==1)],1,min)
  
  # generate intervention effect
  data$fn <- fun(data$distance,50,4,8,c(max_dist),1,data$t,misspec)
  data$y_true <- data$fn * beta
  # simulate outcome data
  data$sim_y <- data$y_true + L%*%rnorm(nrow(data))
  
  if(plot){
    p_dist <- ggplot()+
      geom_sf(data=data, aes(color = distance), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "Distance")+
      theme_solar()+
      ggtitle("Distance")
    
    p_int <- ggplot()+
      geom_sf(data=data, aes(color = y_true), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Intervention effect")
    
    p_u <- ggplot()+
      geom_sf(data=data, aes(color = u), size = 0.1)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Latent spatial effect")
    
    p_p <- ggplot()+
      geom_sf(data=data, aes(color = sim_y), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      scico::scale_color_scico(palette = "roma", name = "Value")+
      theme_solar()+
      ggtitle("Simulated outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

# test function & visualise
dfp <- generate_intervention(dfp, 0.3, -0.3, 15, TRUE)


##########################################################
########### PERMUTATION TEST #############################
# simulation parameter values
beta <- 0
del_e <- 0.25
n_locs <- 16
max_del <- max_de(n_locs)
Li <- solve(L)

pvals <- c()
pvals_ml <- c()
for(zz in 1:1000){
  cat("\nITER: ",zz,"\n")
  dfp <- generate_intervention(dfp, 0.3, beta, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),Li,dfanal$distance,all_dists,0.01,max_del,c(0),c(0),50,4,8,dfanal$t)
  
  # # # maximum likelihood model
  model2 <- Model$new(
    ~ twoway0(distance,8,4,50) + (1|hsgp_fexp(X,Y)),
    data=dfanal,
    covariance = cov_pars,
    mean = c(0.01,-0.01, 0.10),
    family = gaussian()
  )

  # model2$set_trace(1)
  model2$covariance$hsgp(m = c(10,10), L = c(1.05,1.05))
  model2$update_parameters(cov.pars = cov_pars)

  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
                               lower.bound = c(-10,-10,0.01),
                               upper.bound = c(10,10,10)),
                   error = function(i)return(list()))
  se <- tryCatch(sqrt(diag(solve(model2$information_matrix())))[2], error = function(e)return(NA))
  pvals_ml <- c(pvals_ml, 2*(1-pnorm(abs(fit2$coefficients$est[2]/se))))
  
  pval_new <- permute_p_value_b(pt,n_locs,0.01,200,0.3)
  pvals <- c(pvals, pval_new)
  
  if(save_data & zz %% 100){
    saveRDS(pvals,paste0(getwd(),"/pvals_b0_10",n_child,".RDS"))
    saveRDS(pvals_ml,paste0(getwd(),"/pvals_ml_b0_10",n_child,".RDS"))
  }
}

mean(pvals < 0.05)
mean(pvals_ml < 0.05)

##########################
# CONFIDENCE INTERVALS   #

dfci <- data.frame(iter = 1:1000, lower = NA, upper = NA,lower2 = NA, upper2 = NA,lower_ml = NA, upper_ml = NA, b= NA, bp = NA, pval = NA)
dfcid <- data.frame(iter = 1:1000, lower = NA, upper = NA,lower2 = NA, upper2 = NA,lower_ml = NA, upper_ml = NA, d = NA, dp =NA)

fn <- function(x,int,b,del){
  int + b*((1-((-1/50)*log(exp(-50*x/del) + exp(-50)))^4)^8)
}

genrep <- function(dfanal,f1,L,r){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn(distance, int, b, del),data = dfanal, 
                       start = list(int = 0, b = -0.2, del = 0.2),
                       lower = c(-10,-10,0.01), upper = c(10,10,1.0), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 3)
  }
  return(np)
}

cl <- parallel::makeCluster(6)
parallel::clusterExport(cl,c('fn','L','genrep'))

# set parameters of this specific simulation
beta <- -0.3
del_e <- 0.25
n_locs <- 16

for(i in 1:300){
  cat("\rITER: ",i)
  
  dfp <- generate_intervention(dfp, del_e, beta, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y - mean(dfanal$sim_y),Li,dfanal$distance,all_dists,0.01,0.44,0,0,50,4,8,dfanal$t) 
  
  # fit the model using HSGP approximation then calculate proper SEs below
  
  model2 <- Model$new(
    ~ twoway0(distance,8,4,10) + (1|hsgp_fexp(X,Y)),
    data=dfanal,
    covariance = cov_pars,
    mean = c(0.01,-0.01, del_e),
    family = gaussian()
  )

  # model2$set_trace(1)
  model2$covariance$hsgp(m = c(10,10), L = c(1.05,1.05))
  model2$update_parameters(cov.pars = cov_pars)

  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
                               lower.bound = c(-10,-10,0.01),
                               upper.bound = c(10,10,1.0)),
                   error = function(i)return(list()))
  
  if(is(fit2,"mcml")){
    
    b0 <- fit2$coefficients$est[2]
    d0 <- fit2$coefficients$est[3]
    dfci$b[i] <- b0
    dfcid$d[i] <- d0
    
    model2b <-  Model$new(
      ~ twoway0(distance,8,4,10) + (1|fexp(X,Y)),
      data=dfanal,
      covariance = model2$covariance$parameters, #model2$covariance$parameters
      mean = model2$mean$parameters,
      family = gaussian()
    )
    M <- model2b$information_matrix()
    se <- sqrt(diag(solve(M)))
    
    dfci$lower_ml[i] = fit2$coefficients$est[2] - qnorm(0.975)*se[2] #fit2$coefficients$est[2] - qt(0.975, desfac*nrow(dfanal))*se[2]
    dfci$upper_ml[i] = fit2$coefficients$est[2] + qnorm(0.975)*se[2] 
    dfcid$lower_ml[i] = fit2$coefficients$est[3] - qnorm(0.975)*se[3] 
    dfcid$upper_ml[i] = fit2$coefficients$est[3] + qnorm(0.975)*se[3] 
  }
  
  dfci$pval[i] <- tryCatch(permute_p_value_b(pt,n_locs,0.3,200,d0), error = function(i)return(NA))
  
  fitn <- tryCatch(nls(sim_y ~ fn(distance, int, b, del),data = dfanal, 
                       start = list(int = 0, b = -0.2, del = 0.3),
                       lower = c(-10,-10,0.01), upper = c(10,10,1.0), algorithm = "port"), error = function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
    dfci$bp[i] <- np[2]
    dfcid$dp[i] <- np[3]
    
    f1 <- model2b$fitted()
    parallel::clusterExport(cl,c('dfanal','f1'))
    res <- pbapply::pbreplicate(500, genrep(dfanal,f1,as.matrix(L)), cl = cl)
    
    dfci$lower[i] = np[2]-qnorm(0.975)*sd(res[2,],na.rm=TRUE)
    dfci$upper[i] = np[2]+qnorm(0.975)*sd(res[2,],na.rm=TRUE)
    dfcid$lower[i] = np[3]-qnorm(0.975)*sd(res[3,],na.rm=TRUE)
    dfcid$upper[i] = np[3]+qnorm(0.975)*sd(res[3,],na.rm=TRUE)
    
    dfci$lower2[i] = fit2$coefficients$est[2]-qnorm(0.975)*sd(res[2,],na.rm=TRUE)
    dfci$upper2[i] = fit2$coefficients$est[2]+qnorm(0.975)*sd(res[2,],na.rm=TRUE)
    dfcid$lower2[i] = fit2$coefficients$est[3]-qnorm(0.975)*sd(res[3,],na.rm=TRUE)
    dfcid$upper2[i] = fit2$coefficients$est[3]+qnorm(0.975)*sd(res[3,],na.rm=TRUE)
    
  }
  
  # rm(pt, model2, model2b)
  
  if(save_data & i %% 10 == 0){
    saveRDS(dfci,paste0(getwd(),"/dfci2_b",gsub("\\.|-","",as.character(beta)),"_cont_",n_child,"_",ifelse(misspec,"mis",""),".RDS"))
    saveRDS(dfcid,paste0(getwd(),"/dfcid2_b",gsub("\\.|-","",as.character(beta)),"_cont",n_child,"_",ifelse(misspec,"mis",""),".RDS"))
  }
}
