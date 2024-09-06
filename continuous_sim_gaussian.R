# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
source(paste0(getwd(),"/solarized.R"))

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- FALSE # save newly created data


#SIMULATION PARAMETERS

baseline <- 0.3
beta <- -0.15 # max absolute intervention effect
max_dist <- 0.4
b0 <- log(baseline/(1-baseline))
b1 <- log((baseline+beta)/(1-baseline-beta)) - b0
n_locs <- 10
n_seed <- 10
n_child <- 300
cov_pars <- c(0.25,0.5) # G.P. variance, length scale

# FUNCTIONS

# for simulating data
fn2 <- function(d,theta,eff){
  x <- cos(d*pi/(2*theta))
  x[d > theta] <- 0
  x <- eff*x
  return(x)
}

# GENERATE BASE DATA INCL. SAMPLE POINTS AND LATENT SURFACE

# if prior data exists
if(use_data){
  if(!file.exists(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_cont_sp.RDS"))){
    stop("Data for this combination of parameters does not exist")
  } else {
    dfp <- readRDS(paste0(getwd(),"/data/data_",cov_pars[1],"_",cov_pars[2],"_cont_sp.RDS"))
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
    mean = c(b0-exp(cov_pars[1]/2)),
    family = binomial()
  )
  
  sim_data <- mod$sim_data(type="all")
  dfp$u <- sim_data$u
  rm(mod)
  
  
  if(save_data) saveRDS(dfp,paste0(getwd(),"/data/data_",gp_var,"_",lscale,"_cont_sp.RDS"))
  
}

# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  scale_color_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  theme_solar()+
  ggtitle("Simulated locations"); p_loc

# function to generate new

generate_intervention <- function(data, plot = TRUE){
  int_idx <- sample(1:nrow(data),n_locs)
  sampled_locs <- data[int_idx,]
  data$intervention <- 0
  data[int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- NA
  for(i in 1:nrow(data)){
    data$distance[i] <- min(st_distance(data[i,],sampled_locs))
    cat("\rRow ",i," of ",nrow(data))
  }
  
  # generate intervention effect
  data$y_true <- fn2(data$distance,max_dist,b1)
  # simulate outcome data
  data$sim_p <- data$y_true + data$u
  data$sim_p <- exp(data$sim_p)/(1+exp(data$sim_p))
  data$sim_y <- rbinom(nrow(data),1,data$sim_p)
  
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
dfp <- generate_intervention(dfp, TRUE)

## SIMULATION

dfanal <- as.data.frame(dfp)
dfanal <- dfanal[,-which(colnames(dfanal)=="dp")]

# twoway0() is the function (12) in the article
# twoway1() is function (11) with del_I = del_E
# twoway2() is function (11) in the article

model2 <- Model$new(
  ~ twoway0(distance,16,2,50) + (1|hsgp_fexp(X,Y)),
  data=dfanal,
  covariance = c(0.2,0.25),
  mean = c(b0,b1,max_dist),
  family = binomial()
)

model2$set_trace(1)
model2$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
model2$update_parameters(cov.pars = c(0.2,0.2))

# method = "saem", 
# algo = 2, 
# pr.average = FALSE,
# se.theta = FALSE,
# alpha = 0.6,
# convergence.prob = 0.99,
# conv.criterion = 2,

fit2 <- model2$MCML(y = dfp$sim_y, 
                    lower.bound = c(-10,-10,0.01), 
                    upper.bound = c(10,10,10))

saveRDS(fit2,paste0(file_path,"Dropbox/My PC (DESKTOP-NAUP832)/Documents/article_case_study_fit0.RDS"))
fit2 <- readRDS(paste0(file_path,"Dropbox/My PC (DESKTOP-NAUP832)/Documents/article_case_study_fit0.RDS"))


fn <- function(x,l,kappa,b,nu,del){
  b*((1-((-1/l)*log(exp(-l*x/del) + exp(-l)))^kappa)^nu)
}

df1 <- data.frame(distance = seq(0,0.75,length.out=100))
df1$y <- fn(df1$distance,
            50,
            4,
            model2$mean$parameters[2],
            16,
            model2$mean$parameters[3])

df1$cl <- sample(1:10,nrow(df1),replace=TRUE)

modeld <- Model$new(
  ~ twoway0(distance,16,2,50) + (1|gr(cl)),
  data=df1,
  covariance = c( 0.05),
  mean = model2$mean$parameters,
  family = poisson()
)

X0 <- modeld$mean$X
X0[,1] <- 0
M <- solve(model2$information_matrix())
rm(modeld)

df1$se <- sqrt(diag(X0%*%M%*%t(X0)))
df1$lci <- df1$y - qnorm(0.975)*df1$se
df1$uci <- df1$y + qnorm(0.975)*df1$se

df1$true <- fn2(df1$distance, max_dist, b1)


p4 <- ggplot()+
  geom_hline(yintercept = 0,lty=3)+
  geom_ribbon(data = df1, aes(x = distance, ymin=lci,ymax= uci), fill = unname(solar_color[14]), alpha = 0.2)+
  geom_line(data = df1, aes(x = distance, y = y))+
  geom_line(data=df1, aes(x=distance, y=true), lty =2)+
  theme_solar()+
  scale_x_continuous(expand = c(0.01,0))+
  labs(x="Distance (km)", y = "Log relative risk")+
  ggtitle("Estimated treatment effect")

p4

require(patchwork)

(p_int + p_u) / (p_p + p4)

# function
# 1 => 1
# 0 => 0
exp(0)
exp(-Inf)
fn2 <- function(x,theta){exp(-x/theta)}
x <- c(0,0.5,1,1.5,2)
fn2(x,0.2)
fn3 <- function(x,theta,b_del){
  exp(0.02*log(exp(-50*(x/b_del)) + exp(-50))/theta)
}
