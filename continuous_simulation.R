# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
require(Matrix)
source(paste0(getwd(),"/solarized.R"))

# flags for data
use_data <- FALSE # use the simulated saved data
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
misspec <- FALSE

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
  theme_bw()+
  ggtitle("Simulated locations"); p_loc

# function to generate new

all_dists <- st_distance(dfp)

generate_intervention <- function(data, max_dist, beta, n_locs, plot = TRUE){
  
  # spatially regulated sampling scheme
  int_idx <- sample(1:nrow(data),1)
  while(length(int_idx) < n_locs){
    int_idx_new <- sample(1:nrow(data),1)
    dists <- c(all_dists[int_idx,int_idx_new]) #c(st_distance(sampled_locs,dfp[int_idx_new,]))
    if(min(dists) > max_dist*1.5){
      int_idx <- c(int_idx, int_idx_new)
    } 
  }
  
  data$intervention <- 0
  data[int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- apply(all_dists[,which(data$intervention==1)],1,min)
  
  # generate intervention effect
  data$fn <- fun(data$distance,50,4,8,c(max_dist),1,data$t,misspec)
  data$y_true <- data$fn * beta
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
dfp <- generate_intervention(dfp, 0.3, -0.3, 15, TRUE)


##########################################################
########### PERMUTATION TEST #############################
## TEST CPP
dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]

# simulation parameter values
beta <- 0
del_e <- 0.3
n_locs <- 8
max_del <- max_de(n_locs)

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

Si <- model2b$Sigma(TRUE)
L <- chol(Si)
rm(model2b)

pvals <- c()
pvals_ml <- c()
for(zz in 1:1000){
  cat("\nITER: ",zz,"\n")
  dfp <- generate_intervention(dfp, 0.3, beta, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y-mean(dfanal$sim_y),t(B),dfanal$distance,all_dists,0.01,max_del,c(0),c(0),50,4,8,dfanal$t)
  
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
  
  pval_new <- permute_p_value_b(pt,n_locs,0.3,200,0.3)
  pvals <- c(pvals, pval_new)
  
  if(save_data & zz %% 100){
    saveRDS(pvals,paste0(getwd(),"/results/pvals_b0.RDS"))
    saveRDS(pvals_ml,paste0(getwd(),"/results/pvals_ml_b0.RDS"))
  }
}

mean(pvals < 0.05)
mean(pvals_ml < 0.05)

#################
# CONFIDENCE INTERVALS

dfci <- data.frame(iter = 1:1000, lower = NA, upper = NA,lower2 = NA, upper2 = NA, lower_ml = NA, upper_ml = NA, b= NA, bp = NA, pval = NA)
dfcid <- data.frame(iter = 1:1000, lower = NA, upper = NA,lower2 = NA, upper2 = NA, lower_ml = NA, upper_ml = NA, d = NA, dp =NA)


# require(sdmTMB)
# spde <- make_mesh(dfanal, xy_cols = c("X","Y"), cutoff = 0.05)

beta <- -0.3
n_locs <- 8
#0.967

for(i in 77:nrow(dfci)){
  cat("\nITER: ",i,"\n")
  
  dfp <- generate_intervention(dfp, del_e, beta, n_locs, FALSE)
  dfanal <- as.data.frame(dfp)[,-which(colnames(dfp)=="dp")]
  pt <- new_r_stat(dfanal$sim_y - mean(dfanal$sim_y),t(B),dfanal$distance,all_dists,0.01,0.44,0,0,50,4,8,dfanal$t) # update_y(pt,dfanal$sim_y-mean(dfanal$sim_y),dfanal$distance)
  
  model2 <- Model$new(
    ~ twoway0(distance,8,4,50) + (1|hsgp_fexp(X,Y)),
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
                               upper.bound = c(10,10,0.5)),
                   error = function(i)return(list()))
  
  if(is(fit2,"mcml")){
    b0 <- fit2$coefficients$est[2]
    d0 <- fit2$coefficients$est[3]
    M <- model2$information_matrix(hessian.corr = "none")
    se <- sqrt(diag(solve(M)))
    
    dfci$lower_ml[i] = fit2$coefficients$lower[2] #fit2$coefficients$est[2] - qt(0.975, desfac*nrow(dfanal))*se[2]
    dfci$upper_ml[i] = fit2$coefficients$upper[2]
    dfcid$lower_ml[i] = fit2$coefficients$lower[3]
    dfcid$upper_ml[i] = fit2$coefficients$upper[3]
    dfci$b[i] <- b0
    dfcid$d[i] <- d0
    
    dfci$pval[i] <- tryCatch(permute_p_value_b(pt,n_locs,0.3,200,d0), error = function(i)return(NA))
    
    dfci$lower2[i] = fit2$coefficients$est[2] - qnorm(0.975)*se[2] #fit2$coefficients$est[2] - qt(0.975, desfac*nrow(dfanal))*se[2]
    dfci$upper2[i] = fit2$coefficients$est[2] + qnorm(0.975)*se[2] 
    dfcid$lower2[i] = fit2$coefficients$est[3] - qnorm(0.975)*se[3] 
    dfcid$upper2[i] = fit2$coefficients$est[3] + qnorm(0.975)*se[3] 
  }
  
  rm(pt)
  
  if(save_data & i %% 10 == 0){
    saveRDS(dfci,paste0(getwd(),"/results/dfci_b",gsub("\\.|-","",as.character(beta)),"_cont_",n_child,"_",ifelse(misspec,"mis",""),".RDS"))
    saveRDS(dfcid,paste0(getwd(),"/results/dfcid_b",gsub("\\.|-","",as.character(beta)),"_cont",n_child,"_",ifelse(misspec,"mis",""),".RDS"))
  }
}









require(rstan)
mod <- stan_model(paste0(getwd(),"/hsgp_gaussian.stan"))

dat <- list(
 D = 2,
 Q = 1,
 L = c(1.1,1.1),
 M = 10,
 M_nD = 100,
 Nsample = nrow(dfanal),
 nT = 1,
 y = dfanal$sim_y,
 x_grid = as.matrix(dfanal[,c("X","Y")]),
 indices = as.matrix(expand.grid(1:10,1:10)),
 X = matrix(1,nrow=nrow(dfanal),ncol=1),
 dists = dfanal$distance,
 prior_lscale = c(0,0.5),
 prior_var = c(0,0.5),
 prior_linpred_mean = array(0),
 prior_linpred_sd = array(1),
 mod = 1,
 known_cov = 0,
 sigma_data = c(1),
 phi_data = c(1),
 l = -50,
 kappa = 4,
 nu = 8
)

fit <- sampling(mod, data = dat, chains = 4, iter = 500, cores = 4)

print(fitl,c("gamma","b_fn","del_e"))

fitl <- optimizing(mod, data = dat, hessian =TRUE)


confint_both(pt,n_locs,0.4,5000,dfci$bp[i]+0.01,dfcid$dp[i]+0.01,dfcid$dp[i],dfci$bp[i])
confint_both(pt,n_locs,0.4,5000,dfci$bp[i]-0.01,dfcid$dp[i]-0.01,dfcid$dp[i],dfci$bp[i])

ggplot(data = dfanal, aes(x = distance, y = sim_y))+
  geom_point()+
  geom_smooth()

dfanal$pdist <- permute_distances(all_dists,8,0.45)
mean(dfanal$sim_y[dfanal$pdist < 0.1])
lfit <- loess(sim_y ~ distance, data = dfanal,span = 0.7)
x <- seq(0,0.6,length.out = 100)
qplot(x,predict(lfit, x))
cor(dfanal$distance,dfanal$sim_y)

kfit <- ksmooth(dfanal$distance, dfanal$sim_y)
qplot(x,predict(kfit, x))

gfit <- gam(sim_y ~ s(distance), data = dfanal)

require(mgcv)
help(gam)

confint_both(pt,n_locs,0.45,500,-0.2,0.5,best_del[1],dfci$bp[i])

mean(dfcid$lower2 < 0.3 & dfcid$upper2 > 0.3, na.rm=T)
mean(dfci$lower_ml < beta & dfci$upper_ml > beta, na.rm=T)
mean(tmp$lower_ml < beta & tmp$upper_ml > beta, na.rm=T)
mean(tmp$lower2 < 0.3 & tmp$upper2 > 0.3, na.rm=T)

# #permute_distances(all_dists,8,0.45)
fn2 <- fun(dfanal$distance,50,4,8,0.7,1,dfanal$t)
fn3 <- fun(dfanal$distance,50,4,8,dfcid$dp[i],1,dfanal$t)

resid <- (dfanal$sim_y - fit2$coefficients$est[1] - fit2$coefficients$est[2] * fn2)
residb <- (dfanal$sim_y - fit2$coefficients$est[1] - fit2$coefficients$est[2] * fn3)

tstat <- t(fn3)%*%Si%*%resid
tstat

fn2b <- fun(permute_distances(all_dists,8,0.45),50,4,8,best_del[1],1,dfanal$t)
tstat0 <- t(fn2b)%*%Si%*%resid; tstat0

Bx0 <- t(B)%*%(fn3)
Bx <- t(B)%*%(fitTMB$last.par.best[2] * fn2b)

(sum((L%*%resid)^2) - sum((L%*%residb)^2))/sum((L%*%residb)^2)

# Bx0 <- Bx0/norm(Bx0,"E")
# By <- By/norm(By,"E")
# Bx <- Bx/norm(Bx,"E")

t(Bx0)%*%By
t(Bx)%*%By

rs <- c()
for(z in 1:500){
  fn2b <- fun(permute_distances(all_dists,8,0.4),50,4,8,dfcid$dp[i],1,dfanal$t)
  #resid2 <- (dfanal$sim_y - fit2$coefficients$est[1] - fit2$coefficients$est[2] * fn2b)
  rs <- c(rs,t(fn2b)%*%Si%*%resid)
  #rs <- c(rs, (sum((L%*%resid)^2) - sum((L%*%resid2)^2))/sum((L%*%resid2)^2))
}
qplot(rs)
mean(abs(rs) > c(abs(tstat)))

resid2 <- (dfanal$sim_y - fit2$coefficients$est[1] - fit2$coefficients$est[2] * fn3)
oresid <- resid[order(ypred)]


sum(abs(resid))
sum(abs(resid2))

L%*%resid
L%*%resid2

ts <- Si %*% (dfanal$sim_y - 0.15 - fitTMB$last.par.best[2] * fn2)
sum(ts)

ts2 <- Si %*% (dfanal$sim_y - 0.15 - fitTMB$last.par.best[2] * fn2b)
sum(ts2)

fitl <- lm(dfanal$sim_y ~ fn2)
ypred <- fitted(fitl)


By <- L%*%(resid)
Bx0 <- L%*%(fn3)
Bx <- L%*%(fn2b)

# Bx0 <- Bx0/norm(Bx0,"E")
# By <- By/norm(By,"E")
# Bx <- Bx/norm(Bx,"E")

t(Bx0)%*%By
t(Bx)%*%By


didf <- data.frame(dist = seq(0.01,1,length.out = 100), R = NA)
for(i in 1:nrow(didf)){
  fn3 <- fun(dfanal$distance,50,4,8,didf$dist[i],1,dfanal$t)
  Bx0 <- t(B)%*%(fn3)
  Bx0 <- Bx0/norm(Bx0,"E")
  didf$R[i] <- -1.0*abs(t(Bx0)%*%By)
}

ggplot(data= didf, aes(x = dist, y = R))+
  geom_line()

# cor(By,Bx0)
# cor(By,Bx)
# qplot(Bx,By)
