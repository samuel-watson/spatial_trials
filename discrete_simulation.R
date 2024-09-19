# PACKAGES AND PRE-REQUISITES


require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
source(paste0(getwd(),"/solarized.R"))

Rcpp::sourceCpp("src/perm_test.cpp")

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- FALSE # save newly created data

#SIMULATION PARAMETERS

# Upper bound on Del_E 
max_de <- function(n_locs){
  sqrt(8/(n_locs*3*sqrt(3)))
}

beta <- -0.3 # max absolute intervention effect
max_dist <- 0.3
radius <- 0.2
n_locs <- 8
n_seed <- 10
n_child <- 100
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
del_e_ul <- max_de(n_locs)

# GENERATE BASE DATA INCL. SAMPLE POINTS AND LATENT SURFACE


# if prior data exists
if(use_data){
  if(!file.exists(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))){
    stop("Data for this combination of parameters does not exist")
  } else {
    dfp <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))
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
  
  mod <- Model$new(
    ~ (1|fexp(X,Y)),
    data = as.data.frame(dfp)[,c("X","Y")],
    covariance = cov_pars,
    mean = c(0),
    family = gaussian()
  )
  
  sim_data <- mod$sim_data(type="all")
  dfp$u <- sim_data$u
  rm(mod)
  
  
  if(save_data) saveRDS(dfp,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"dis_sp.RDS"))
}

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  theme_bw()+
  ggtitle("Simulated locations"); p_loc

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

generate_intervention <- function(data, del_e, del_i, beta, n_locs, plot = TRUE){
  
  # spatially regulated sampling scheme
  
  idx <- sort(sample(1:nrow(dfi),n_locs))
  data$distance <- apply(dists_i[,idx],1,min)
  data$fn <- fun(data$distance,50,4,8,c(del_e,-del_i),1,data$t)
  
  # generate intervention effect
  data$y_true <- data$fn * beta
  # simulate outcome data
  data$sim_p <- data$y_true + data$u
  data$sim_y <- rnorm(nrow(data),data$sim_p)
  
  if(plot){
    p_dist <- ggplot()+
      geom_sf(data=data, aes(color = distance), size = 0.1)+
      geom_point(data=dfi[idx,],aes(x=x,y=y),color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "Distance")+
      theme_solar()+
      ggtitle("Distance")
    
    p_int <- ggplot()+
      geom_sf(data=data, aes(color = y_true), size = 0.1)+
      geom_point(data=dfi[idx,],aes(x=x,y=y),color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Intervention effect")
    
    p_u <- ggplot()+
      geom_sf(data=data, aes(color = u), size = 0.1)+
      scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
      theme_solar()+
      ggtitle("Latent spatial effect")
    
    p_p <- ggplot()+
      geom_sf(data=data, aes(color = sim_p), size = 0.1)+
      geom_point(data=dfi[idx,],aes(x=x,y=y),color="red",size=2)+
      scico::scale_color_scico(palette = "roma", name = "Simulated\nprobability")+
      theme_solar()+
      ggtitle("Probability of outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

##########################################################
########### PERMUTATION TEST #############################


# simulation parameter values
del_e <- 0.2
del_i <- 0.1
n_locs <- 8
max_del <- max_de(n_locs)

#test
dfp <- generate_intervention(dfp, del_e, del_i, -0.3, n_locs, TRUE)
dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]

model2b <- Model$new(
  ~ (1|fexp(X,Y)),
  data=dfanal,
  covariance = cov_pars,
  mean = 0,
  family = gaussian()
)
S <- model2b$Sigma()
E <- eigen(S)
B <- E$vectors%*%diag(1/sqrt(E$values))

beta <- 0
pvals <- c()
pvals_ml <- c()
for(zz in 1:1000){
  cat("\nITER: ",zz,"\n")
  dfp <- generate_intervention(dfp, del_e, del_i, beta, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),as.matrix(B),dfanal$distance,dists_i,c(0.01),c(0.4),c(-0.2),c(0.0),50,4,8,dfanal$t)
  
  # # maximum likelihood model
  model2 <- Model$new(
    ~ twoway2(distance,8,4,50) + (1|hsgp_fexp(X,Y)),
    data=dfanal,
    covariance = cov_pars,
    mean = c(0.01,-0.01, 0.10, 0.10),
    family = gaussian()
  )
  
  # model2$set_trace(1)
  model2$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
  model2$update_parameters(cov.pars = cov_pars)
  
  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y,
                               lower.bound = c(-10,-10,0.01,0.01),
                               upper.bound = c(10,10,10,1.0)),
                   error = function(i)return(list()))
  se <- tryCatch(sqrt(diag(solve(model2$information_matrix())))[2], error = function(e)return(NA))
  pvals_ml <- c(pvals_ml, 2*(1-pnorm(abs(fit2$coefficients$est[2]/se))))
  
  pval_new <- permute_p_value_b(pt,n_locs,0.3,200,c(0.2,-0.1))
  pvals <- c(pvals, pval_new)
  rm(pt)
  if(save_data & i %% 100){
    saveRDS(pvals,paste0(getwd(),"/results/pvals_b0_dis.RDS"))
    saveRDS(pvals_ml,paste0(getwd(),"/results/pvals_ml_b0_dis.RDS"))
  }
}

mean(pvals < 0.05)

# beta 0 non zero

dfci <- data.frame(iter = 1:1000, lower = NA, upper = NA, lower_ml = NA, upper_ml = NA, b= NA, bp = NA, pval = NA)
dfcid <- data.frame(iter = 1:1000, lower_e = NA, upper_e = NA, lower_ml_e = NA, upper_ml_e = NA, lower_i = NA, upper_i = NA, lower_ml_i = NA, upper_ml_i = NA, d_e = NA, dp_e =NA, d_i = NA, dp_i =NA)

model2b <- Model$new(
  ~ (1|fexp(X,Y)),
  data=dfanal,
  covariance = cov_pars,
  mean = 0,
  family = gaussian()
)
S <- model2b$Sigma()
E <- eigen(S)
B <- E$vectors%*%diag(1/sqrt(E$values))



for(i in 1:1000){
  cat("\nITER: ",i,"\n")
  dfp <- generate_intervention(dfp, 0.2, 0.1, -0.3, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),as.matrix(B),dfanal$distance,dists_i,c(0.01),c(0.4),c(-0.2),c(0.0),50,4,8,dfanal$t)
  
  model2 <- Model$new(
    ~ twoway2(distance,8,4,50) + (1|hsgp_fexp(X,Y)),
    data=dfanal,
    covariance = cov_pars,
    mean = c(0.01,-0.01, 0.10, 0.10),
    family = gaussian()
  )
  
  # model2$set_trace(1)
  model2$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
  model2$update_parameters(cov.pars = cov_pars)
  
  fit2 <- tryCatch(model2$MCML(y = dfp$sim_y, 
                               lower.bound = c(-10,-10,0.01,0.0), 
                               upper.bound = c(10,10,10,1)),
                   error = function(i)return(list()))
  
  if(is(fit2,"mcml")){
    b0 <- fit2$coefficients$est[2]
    di0 <- fit2$coefficients$est[3]
    de0 <- fit2$coefficients$est[4]
    
    dfci$lower_ml[i] = fit2$coefficients$lower[2]
    dfci$upper_ml[i] = fit2$coefficients$upper[2]
    dfcid$lower_ml_i[i] = fit2$coefficients$lower[3]
    dfcid$upper_ml_i[i] = fit2$coefficients$upper[3]
    dfcid$lower_ml_e[i] = fit2$coefficients$lower[4]
    dfcid$upper_ml_e[i] = fit2$coefficients$upper[4]
    dfci$b[i] <- b0
    dfcid$d_i[i] <- di0
    dfcid$d_e[i] <- de0
  }
  
  
  best_del <- optim_r(pt,c(0.2,0.1))
  dfcid$dp_e[i] <- best_del[1]
  dfcid$dp_i[i] <- best_del[3]

  dfanal$fn <- fun(dfanal$distance,50,4,8,best_del[c(1,3)],1,dfanal$t,TRUE)
  
  if(!any(!is.finite(dfanal$fn)|is.na(dfanal$fn))){
    
    model2b <- Model$new(
      ~ fn + (1|hsgp_fexp(X,Y)),
      data=dfanal,
      covariance = cov_pars,
      mean = c(0.01,-0.4),
      family = gaussian()
    )
    
    # model2$set_trace(1)
    model2b$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
    model2b$update_parameters(cov.pars = cov_pars)
    
    fit2b <- tryCatch(model2b$MCML(y = dfp$sim_y, 
                                   lower.bound = c(-10,-10), 
                                   upper.bound = c(10,10)),
                      error = function(i)return(list()))
    
    if(is(fit2b,"mcml")){
      dfci$bp[i] <- fit2b$coefficients$est[2]
      dfci$upper[i] <- confint_b(pt,n_locs,0.3,250,dfci$bp[i]+0.2,dfci$bp[i],best_del[c(1,3)])
      dfci$lower[i] <- confint_b(pt,n_locs,0.3,250,dfci$bp[i]-0.2,dfci$bp[i],best_del[c(1,3)])
      d_lower <- confint_del(pt,n_locs,0.3,250,c(0.05,-0.4),dfci$bp[i],best_del[c(1,3)])
      d_upper <- confint_del(pt,n_locs,0.3,250,c(best_del[1]+0.1,0.01),dfci$bp[i],best_del[c(1,3)])
      dfcid$lower_e[i] <- d_lower[1]
      dfcid$upper_e[i] <- d_upper[1]
      dfcid$lower_i[i] <- d_lower[2]
      dfcid$upper_i[i] <- d_upper[2]
    }
    
    dfci$pval[i] <- permute_p_value_b(pt,n_locs,0.3,200,best_del[c(1,3)])
  }
  
  rm(pt)
  if(save_data & i %% 10 == 0){
    saveRDS(dfci,paste0(getwd(),"/results/dfci_b",gsub("\\.|-","",as.character(beta)),"_dis.RDS"))
    saveRDS(dfcid,paste0(getwd(),"/results/dfcid_b",gsub("\\.|-","",as.character(beta)),"_dis.RDS"))
  }
}

