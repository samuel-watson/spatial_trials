# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
source(paste0(getwd(),"/solarized.R"))

# flags for data
use_data <- FALSE # use the simulated saved data
save_data <- TRUE # save newly created data

#SIMULATION PARAMETERS

# Upper bound on Del_E 
max_de <- function(n_locs){
  sqrt(8/(n_locs*3*sqrt(3)))
}

beta <- -0.3 # max absolute intervention effect
max_dist <- 0.3
n_locs <- 8
n_seed <- 10
n_child <- 100
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
del_e_ul <- max_de(n_locs)

# FUNCTIONS

# for simulating data
fn <- function(x,l,kappa,b,nu,del){
  b*((1-((-1/l)*log(exp(-l*x/del) + exp(-l)))^kappa)^nu)
}

fn2 <- function(d,theta,eff){
  x <- cos(d*pi/(2*theta))
  x[d > theta] <- 0
  x <- eff*x
  return(x)
}

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
  
  
  if(save_data) saveRDS(dfp,paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_sp.RDS"))
}

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  scale_color_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  theme_solar()+
  ggtitle("Simulated locations"); p_loc

# function to generate new

generate_intervention <- function(data, max_dist, beta, n_locs, plot = TRUE){
  
  # spatially regulated sampling scheme
  int_idx <- sample(1:nrow(data),1)
  sampled_locs <- data[int_idx,]
  while(nrow(sampled_locs) < n_locs){
    int_idx_new <- sample(1:nrow(data),1)
    dists <- c(st_distance(sampled_locs,dfp[int_idx_new,]))
    if(min(dists) > max_dist*1.5){
      sampled_locs <- rbind(sampled_locs, dfp[int_idx_new,])
      int_idx <- c(int_idx, int_idx_new)
    } 
  }
  
  data$intervention <- 0
  data[int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- NA
  for(i in 1:nrow(data)){
    data$distance[i] <- min(st_distance(data[i,],sampled_locs))
    cat("\rRow ",i," of ",nrow(data))
  }
  
  # generate intervention effect
  data$y_true <- fn(data$distance,50,4,beta,8,max_dist)
    #fn2(data$distance,max_dist,b1)
  # simulate outcome data
  data$sim_p <- data$y_true + data$u
  data$sim_y <- rnorm(nrow(data),data$sim_p)
  
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
      geom_sf(data=data, aes(color = sim_p), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      scico::scale_color_scico(palette = "roma", name = "Simulated\nprobability")+
      theme_solar()+
      ggtitle("Probability of outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

# test function & visualise
dfp <- generate_intervention(dfp, 0.3, -0.3, 8, TRUE)

## SIMULATION

result <- list(
  b_eff = c(),
  del_e = c(),
  cover_eff = c(),
  cover_del = c(),
  cover_eff_kr = c(),
  cover_del_kr = c(),
  ext = c(),
  lr = c()
)

for(i in 1:100){
  cat("\n\nITER: ",i," OF 100\n\n")
  
  fitted <- FALSE
  
  while(!fitted){
    dfp <- generate_intervention(dfp, 0.3, -0.3, 10, FALSE)
    dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
    
    # twoway0() is the function (12) in the article
    # twoway1() is function (11) with del_I = del_E
    # twoway2() is function (11) in the article
    
    model2 <- Model$new(
      ~ twoway0(distance,8,4,50) + (1|hsgp_fexp(X,Y)),
      data=dfanal,
      covariance = c(0.2,0.25),
      mean = c(0.01,-0.01, 0.10),
      family = gaussian()
    )
    
    model2$set_trace(1)
    model2$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
    model2$update_parameters(cov.pars = c(0.2,0.2))
    
    fit2 <- tryCatch(model2$MCML(y = dfp$sim_y, 
                        lower.bound = c(-10,-10,0.01), 
                        upper.bound = c(10,10,10)),
             error = function(i)return(list()))
    if(is(fit2,"mcml"))fitted <- TRUE
  }
  
  model2b <- Model$new(
    ~ twoway0(distance,8,4,50) + (1|fexp(X,Y)),
    data=dfanal,
    covariance = model2$covariance$parameters,
    mean = model2$mean$parameters,
    family = gaussian()
  )
  model2b$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
  model2b$update_parameters(cov.pars = model2$covariance$parameters)
  model2b$update_y(dfp$sim_y)
  kr <- model2b$small_sample_correction("KR")
  se <- sqrt(diag(kr$vcov_beta))
  
  yfit <- fn(dfp$distance,50,4,model2$mean$parameters[2],8,model2$mean$parameters[3])
  beta_est <- model2$mean$parameters
  yvar <- dfp$sim_y - yfit + model2$mean$X[,2:3] %*% beta_est[2:3] - model2$u()[,1]
  
  lfit1 <- lm(yvar ~ model2$mean$X - 1)
  lfit2 <- lm(yvar ~ 1)
  an1 <- anova(lfit1, lfit2)
  
  result$b_eff <- c(result$b_eff, fit2$coefficients$est[2])
  result$del_e <- c(result$del_e, fit2$coefficients$est[3])
  result$cover_eff <- c(result$cover_eff, fit2$coefficients$lower[2] < beta & fit2$coefficients$upper[2] > beta)
  result$cover_del <- c(result$cover_eff, fit2$coefficients$lower[3] < max_dist & fit2$coefficients$upper[3] > max_dist)
  result$cover_eff_kr <- c(result$cover_eff_kr, fit2$coefficients$est[2] - qt(0.975,df=kr$dof[2])*se[2] < beta & fit2$coefficients$est[2] + qt(0.975,df=kr$dof[2])*se[2] > beta)
  result$cover_del_kr <- c(result$cover_eff_kr, fit2$coefficients$est[3] - qt(0.975,df=kr$dof[3])*se[3] < max_dist & fit2$coefficients$est[3] + qt(0.975,df=kr$dof[3])*se[3] > max_dist)
  result$ext <- c(result$ext, 1-pnorm((model2$mean$parameters[3]-del_e_ul)/fit2$coefficients$SE[3]))
  result$lr <- c(result$lr, an1$`Pr(>F)`[2])
  
  if(i %% 100 == 0) saveRDS(result, paste0(getwd(),"/results/sim_",cov_pars[1],"_",cov_pars[2],"_",n_seed,n_child,"cont_sp.RDS"))
}



qplot(result$b_eff)
mean(result$b_eff)

qplot(result$del_e)
mean(result$del_e)

mean(result$cover_eff)
mean(result$cover_eff_kr)

mean(result$cover_del)
mean(result$cover_del_kr)

mean(result$lr <= 0.05)

