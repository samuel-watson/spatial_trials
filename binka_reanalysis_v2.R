# binka re-analysis

require(glmmrBase)
require(ggplot2)
require(sf)
file_path <- "C:/Users/watsonsi/"
file_path <- "D:/Documents/"
source(paste0(file_path,"Dropbox/solarized.R"))

df <- read.csv(paste0(file_path,"Dropbox/binka_compounds.csv"))
df <- df[df$expected > 0 ,]
df$t <- 1
dfp <- rts2::create_points(df,pos_vars = c('x','y'), t_var = "t")
dfp$cl <- df$cluster
dfp$arm <- df$arm
dfp$deaths <- df$deaths
dfp$expected <- df$expected
dfp$nets <- df$nets
cl <- unique(df$cluster)
# create convex hull

for(i in cl){
  p1 <- concaveman::concaveman(dfp[dfp$cl==i,], concavity = 5) # otherwise
  p1$cl <- i
  p1$arm <- df[df$cluster == i, 'arm'][1]
  if(exists("dfpoly")){
    dfpoly <- rbind(dfpoly, p1)
  } else {
    dfpoly <- p1
  }
}

plot(dfpoly[dfpoly$cl > 0,"cl"])

# lets do some plotting

ggplot()+
  geom_sf(data=dfpoly[dfpoly$cl > 0,], aes(fill = factor(arm)),alpha = 0.2)+
  scale_fill_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  theme_solar()

# plot as spatial trial

dfpoly_s <- st_sf(st_sfc(st_convex_hull(st_union(dfp))))
dfpoly_t <- concaveman::concaveman(dfp, concavity = 20)

p1 <- ggplot()+
  geom_sf(data=dfpoly_t, alpha = 0.4, color = NA)+
  geom_sf(data = dfp, aes(color=factor(arm)), size=0.2, alpha=0.2)+
  geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
  scale_color_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  theme_solar()+
  theme(axis.text = element_blank())

p1b <- ggplot()+
  geom_sf(data=dfpoly_t, alpha = 0.4, color = NA)+
  geom_sf(data = dfp[dfp$deaths>0,], aes(color=factor(deaths)), size=0.3, alpha=0.2)+
  geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
  scale_color_manual(values = unname(solar_color[c(11:14)]),name = "N deaths")+
  theme_solar()+
  theme(axis.text = element_blank())

## now process into a "spatial trial"

## calculate distances
# if the point is inside an intervention area then it has distance 0 otherwise calculate minimum distance
dfp$distance <- NA
dfp$id <- 1:nrow(dfp)
int_area <- st_union(dfpoly[dfpoly$arm == "intervention",])
int_area_boundary <- st_boundary(int_area)
dfp_int <- st_filter(dfp,int_area)

for(i in 1:nrow(dfp)){
  if(i %in% dfp_int$id){
    dfp$distance[i] <- -1 * st_distance(dfp[i,],int_area_boundary)
  } else {
    dfp$distance[i] <- min(st_distance(dfp[i,],int_area), 1000) # set to lower value for plotting
  }
  cat("\rRow: ",i," of ",nrow(dfp))
}

dfp$distance_re <- dfp$distance + abs(min(dfp$distance))

## distance to any boundary
all_area <- st_union(dfpoly)
all_area_boundary <- st_boundary(all_area)
dfp$distance_all <- NA

for(i in 1:nrow(dfp)){
  dfp$distance_all[i] <- -1 * st_distance(dfp[i,],all_area_boundary) # set to lower value for plotting
  cat("\rRow: ",i," of ",nrow(dfp))
}

p2 <- ggplot()+
      geom_sf(data=dfpoly_t, alpha = 0.4, color = NA)+
      geom_sf(data = dfp, aes(color=distance), size=0.2, alpha=0.2)+
      geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
      scico::scale_color_scico(palette="roma")+
      theme_solar()+
      theme(axis.text = element_blank())

require(patchwork)
p1 + p2

ggplot(data=dfp,aes(x=distance))+
  geom_histogram()+
  theme_solar()+
  ggtitle("Distance from intervention area")

# lets try and fit the model

dfanal <- as.data.frame(dfp)
dfanal <- cbind(dfanal, df[,c('x','y')])
dfanal <- dfanal[,5:13]
# rescale x and y
xrange <- range(dfanal$x)
yrange <- range(dfanal$y)
scale_f <- max(diff(xrange),diff(yrange))
dfanal$x_re <- -1 + 2*(dfanal$x - min(dfanal$x))/scale_f  #-1 + (2 / diff(range(dfanal$x)))*(dfanal$x - min(dfanal$x))
dfanal$y_re <- -1 + 2*(dfanal$y - min(dfanal$y))/scale_f # -1 + (2 / diff(range(dfanal$y)))*(dfanal$x - min(dfanal$y))
dfanal <- dfanal[order(dfanal$y),]
# jitter the duplicated location
locs <- paste0(dfanal$x_re, dfanal$y_re)
dfanal[duplicated(locs),'y_re'] <- dfanal[duplicated(locs),'y_re'] + 1e-6
saveRDS(dfanal,"D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_analysis_data.RDS")
dfanal <- readRDS(paste0(file_path,"Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_analysis_data.RDS"))

model <- Model$new(
  ~ twoway1(distance,16,4,50) + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05,0.05),
  mean = c(0.07,-0.3,0.3),
  offset = log(dfanal$expected),
  family = poisson()
)

model3 <- Model$new(
  ~ twoway0(distance,8,50) + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05,0.05),
  mean = c(0.07,-0.3,0.3,4),
  offset = log(dfanal$expected),
  family = poisson()
)

# model <- Model$new(
#   ~ b_eff*((1 - (-0.02 *sign0(distance)* log(exp(-50*sign0(distance) * (((distance)+b_del_i)/(b_del_i + b_del_e))) + exp(-25*(sign0(distance)+1))))^(4))^(16)) + (1|hsgp_fexp(x_re,y_re)),
#   data=dfanal,
#   covariance = c(0.05,0.05),
#   mean = c(0.07,-0.01,0.4,0.8),
#   offset = log(dfanal$expected),
#   family = poisson()
# )

M <- model$information_matrix()
A <- model$information_matrix(hessian.corr = "return")

model$covariance$hsgp(m = c(25,25), L = c(1.2,1.2))
model$set_trace(1)
model$mcmc_options$warmup <- 100
model$mcmc_options$samps <- 30

fit0 <- model$MCML(y = dfanal$deaths, 
                  method = "saem", 
                  algo = 2, 
                  pr.average = FALSE,
                  se.theta = FALSE,
                  conv.criterion = 2,
                  lower.bound = c(-10,-10,0,0), 
                  upper.bound = c(10,10,2,2))

#0.3495881 -0.2829673 0.01 3.344217
# 0.0001125864 0.005914345
u <- model$u(scaled = FALSE)
saveRDS(u,"D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/u_binka.RDS")
saveRDS(fit0,"D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit_smooth.RDS")

# two-way model

model2 <- Model$new(
  ~ twoway2(distance,16,4,50) + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.08944, 0.09932),
  mean = c(0.37704, -0.29246 ,0.41069 ,0.76303),
  offset = log(dfanal$expected),
  family = poisson()
)

model2$covariance$hsgp(m = c(25,25), L = c(1.5,1.5))
model2$set_trace(1)
model2$mcmc_options$warmup <- 100
model2$mcmc_options$samps <- 30

fit2 <- model2$MCML(y = dfanal$deaths, 
                   method = "mcem", 
                   algo = 2, 
                   pr.average = FALSE,
                   se.theta = FALSE,
                   conv.criterion = 2,
                   lower.bound = c(-10,-10,0,0), 
                   upper.bound = c(10,10,2,2))

#  0.37704 -0.29246 0.41069 0.76303
# 0.08944 0.09932

model2$update_y(dfanal$deaths)
model2$update_parameters(mean.pars = c(0.37704 ,-0.29246, 0.41069, 0.76303),
                         cov.pars = c(0.08944, 0.09932))
model2$mcmc_sample()

M <- model2$information_matrix(hessian.corr = "none")
A <- model2$information_matrix(hessian.corr = "return")
M <- M+A
M <- readRDS("D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/M2.RDS")

saveRDS(M,"D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/M2.RDS")
saveRDS(fit2,"C:/Users/watsonsi/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2.RDS")
fit2 <- readRDS("D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2.RDS")


# need to add model that controls for potential distance - either: 
# i) distance indicators; ii) smooth function

for(i in 0:8){
  dfanal$tmp <- I(abs(dfanal$distance_all) > (i-1)*0.2 & abs(dfanal$distance_all) <= i*0.2)*1
  colnames(dfanal)[ncol(dfanal)] <- paste0("distance_all",i)
}
dfanal$distance_all_sq <- dfanal$distance_all^2

model2a <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all1 + distance_all2 + distance_all3 + distance_all4 + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.08016, 0.16235),
  mean = c(0.45605, -0.28641, 0.41107, 0.95507, -0.07335, -0.05357, -0.02834, -0.17584),
  offset = log(dfanal$expected),
  family = poisson()
)

model2a$covariance$hsgp(m = c(25,25), L = c(1.2,1.2))
model2a$set_trace(1)
model2a$mcmc_options$warmup <- 100
model2a$mcmc_options$samps <- 30
model2a$update_parameters(cov.pars = c(0.05161, 0.15722))
#model2a$update_y(dfanal$deaths)
#model2a$mcmc_sample()

fit2a <- model2a$MCML(y = dfanal$deaths, 
                      method = "saem", 
                      algo = 2, 
                      pr.average = FALSE,
                      se.theta = FALSE,
                      conv.criterion = 2,
                      lower.bound = c(-10,-10,0,0,-10,-10,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10,10,10))

#  0.45605 -0.28641 0.41107 0.95507 -0.07335 -0.05357 -0.02834 -0.17584
# 0.08016 0.16235

M2a <- model2a$information_matrix()
M2a0 <- model2a$information_matrix(hessian.corr = "none")
A <- model2a$information_matrix(hessian.corr = "return")
M2a <- M2a0 + A
saveRDS(M2a,"C:/Users/watsonsi/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/M2a.RDS")
saveRDS(fit2a,"C:/Users/watsonsi/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2a.RDS")
fit2a <- readRDS("D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2a.RDS")

# model2b <- Model$new(
#   ~ twoway2(distance,16,50) + distance_all + distance_all_sq + (1|hsgp_fexp(x_re,y_re)),
#   data=dfanal,
#   covariance = c(0.01,0.05),
#   mean = c(0.46069, -0.29128, 0.15855, 1.46512, 2.06206,0.05,-0.05),
#   offset = log(dfanal$expected),
#   family = poisson()
# )

model2b <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all + distance_all_sq + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.07511, 0.14205),
  mean = c(0.35819, -0.29053, 0.38594, 0.9419, 0.00151, 0.03722),
  offset = log(dfanal$expected),
  family = poisson()
)

model2b$update_y(dfanal$deaths)
model2b$covariance$hsgp(m = c(25,25), L = c(1.2,1.2))
model2b$set_trace(1)
model2b$mcmc_options$warmup <- 100
model2b$mcmc_options$samps <- 30
model2b$update_parameters(cov.pars = c(0.04702, 0.16173))
#model2b$mcmc_sample()

fit2b <- model2b$MCML(y = dfanal$deaths, 
                      method = "saem", 
                      algo = 2, 
                      pr.average = FALSE,
                      se.theta = FALSE,
                      conv.criterion = 2,
                      lower.bound = c(-10,-10,0,0,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10))

#  0.35819 -0.29053 0.38594 0.9419 0.00151 0.03722
# 0 0.07511 0.14205

M2b <- model2b$information_matrix(hessian.corr = "none")
Ab <- model2b$information_matrix(hessian.corr = "return")
M2b <- M2b + Ab
# M2b <- readRDS("D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/M2b.RDS")
saveRDS(M2b,"C:/Users/watsonsi/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/M2b.RDS")
saveRDS(fit2b,"C:/Users/watsonsi/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2b.RDS")
fit2b <- readRDS("D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit2b.RDS")

# model with distance markers
for(i in 0:10){
  if(i == 0){
    dfanal$tmp <- I(dfanal$distance == 0)*1
  } else {
    dfanal$tmp <- I(dfanal$distance > (i-1)*0.1 & dfanal$distance <= i*0.1)*1
  }
  colnames(dfanal)[ncol(dfanal)] <- paste0("distance",i)
}

# model with distance markers
for(i in 0:10){
  if(i == 0){
    dfanal$tmp <- I(dfanal$distance == 0)*1
  } else {
    dfanal$tmp <- I(dfanal$distance > (i-1)*0.1 & dfanal$distance <= i*0.1)*1
  }
  colnames(dfanal)[ncol(dfanal)] <- paste0("distance",i)
}

model <- Model$new(
  ~ distance0 + distance1 + distance2 + distance3 + distance4 + distance5 + distance6 + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.1,0.05),
  mean = c(0.38416 ,-0.30265 ,-0.12292, -0.33724, -0.18239, 0.03, -0.17747, 0.01),
  offset = log(dfanal$expected),
  family = poisson()
)

model$covariance$hsgp(m = c(25,25), L = c(1.2,1.2))
model$set_trace(1)
model$mcmc_options$warmup <- 100
model$mcmc_options$samps <- 30

fit6 <- model$MCML(y = dfanal$deaths, 
                   method = "saem", 
                   algo = 2, 
                   pr.average = FALSE,
                   se.theta = FALSE,
                   conv.criterion = 2)

saveRDS(fit6,"D:/Documents/Dropbox/My PC (DESKTOP-NAUP832)/Documents/spatial/binka_fit6.RDS")

# RESULTS
dfres <- data.frame(
  n = c(2,3,4,5,6),
  aic = c(10361.52,10299.62,10302.74,10406.23,10407.69)
)



# plot the individual res
u2 <- model$u()
dfanal$pred <- rowMeans(u)

ggplot(data = dfanal,aes(y=x_re,x=y_re,color=pred))+geom_point()

model$update_parameters(mean.pars = c(0.40319, -0.1, 0.8, 2),
                        cov.pars = c( 0.112, 0.3))


dfgrid <- expand.grid(x_re = seq(-1,1,length.out = 40), y_re = seq(-1,1,length.out = 40))
dfgrid$x <- dfgrid$x_re
dfgrid$y <- dfgrid$y_re
dfgrid$deaths <- 0
dfgrid$expected <- 0
dfgrid$nets <- 0
dfgrid$distance <- 3
dfgrid$id <- 0
dfgrid$death_ind <- 0
preds <- model$predict(newdata = dfgrid)

dfgrid <- expand.grid(x = seq(-1,1,length.out = 40), y = seq(-1,1,length.out = 40),t =1)
dfgrid$x <- (dfgrid$x + 1)/(2 / diff(range(dfanal$x))) + min(dfanal$x)
dfgrid$y <- (dfgrid$y + 1)/(2 / diff(range(dfanal$y))) + min(dfanal$y)
dfpgrid <- rts2::create_points(dfgrid,pos_vars = c('x','y'), t_var = "t")
dfpgrid$pred <- preds$re_parameters$vec
dfpgrid <- st_filter(dfpgrid,dfpoly_t)
dfgrid <- as.data.frame(dfpgrid)
dfgrid <- cbind(dfgrid,sf::st_coordinates(dfpgrid))

p3 <- ggplot(data = dfgrid)+
  geom_sf(data = dfpoly_t)+
  geom_tile(aes(x=X,y=Y,fill=pred))+
  geom_sf(data = dfp[dfp$deaths>0,], size=0.3)+
  geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
  scico::scale_fill_scico(palette = "bamako",name="Smoothed\nlatent\nrisk")+
  theme_solar()+
  theme(axis.text = element_blank())

### plot the function

fn <- function(x,l,kappa,nu,del_e, del_i, b){
  b*((1-((sign(x)/l)*log(exp(sign(x)*l*(x + del_i)/(del_e + del_i)) + exp(l*(sign(x)+1)/2)))^kappa)^nu) 
}

# fn <- function(x,l,kappa,b,nu,del){
#   b*((1-((-1/l)*log(exp(-l*x/del) + exp(-l)))^kappa)^nu)
# }

# -0.3119 0.16388 1.58156 2.06841
df1 <- data.frame(distance = rep(seq(-1,1.5,length.out=100),3),
                  model = rep(1:3,each=100))
df1 <- df1[df1$distance!=0,]
df1$expected <- 0
df1$distance_all1 <- 0
df1$distance_all2 <- 0
df1$distance_all3 <- 0
df1$distance_all4 <- 0
df1$distance_all_sq <- 0
df1$distance_all <- 0

#df1$y <- fn(df1$distance[1:100],-50,2.05801,16, 1.58349, 0.15759 ,-0.29911)
df1$y <- c(fn(df1$distance[1:100],-50,4,16,0.76303, 0.41069 ,-0.29246), 
  fn(df1$distance[1:100],-50,4,16, 0.95507 , 0.41107 ,-0.28641) , 
 fn(df1$distance[1:100],-50,4,16, 0.9419, 0.38594 ,-0.29053)) 

df1$cl <- sample(1:10,nrow(df1),replace=TRUE)
df1$se <- NA
df1$lci <- NA
df1$uci <- NA

modeld <- Model$new(
  ~ twoway2(distance,16,4,50) + (1|gr(cl)),
  data=df1[1:100,],
  covariance = c( 0.05),
  mean = c(0.37704, -0.29246, 0.41069, 0.76303),
  family = poisson()
)
X0 <- modeld$mean$X
X0[,1] <- 0
M <- solve(M)
df1$se[1:100] <- sqrt(diag(X0%*%M%*%t(X0)))
df1$lci[1:100] <- df1$y[1:100] - qnorm(0.975)*df1$se[1:100]
df1$uci[1:100] <- df1$y[1:100] + qnorm(0.975)*df1$se[1:100]
rm(modeld)

# fit 2 a
modeld <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all1 + distance_all2 + distance_all3 + distance_all4 + (1|gr(cl)),
  data=df1[101:200,],
  covariance = c( 0.05),
  mean = c(0.45605, -0.28641, 0.41107, 0.95507, -0.07335, -0.05357, -0.02834, -0.17584),
  family = poisson()
)
X0 <- modeld$mean$X
X0[,1] <- 0
M2a <- solve(M2a)
df1$se[101:200] <- sqrt(diag(X0%*%M2a%*%t(X0)))
df1$lci[101:200] <- df1$y[101:200] - qnorm(0.975)*df1$se[101:200]
df1$uci[101:200] <- df1$y[101:200] + qnorm(0.975)*df1$se[101:200]
rm(modeld)

# fit 2 b
modeld <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all + distance_all_sq + (1|gr(cl)),
  data=df1[201:300,],
  covariance = c( 0.05),
  mean = c(0.35819, -0.29053, 0.38594, 0.9419, 0.00151 ,0.03722),
  family = poisson()
)
X0 <- modeld$mean$X
X0[,1] <- 0
M2b <- solve(M2b)
df1$se[201:300] <- sqrt(diag(X0%*%M2b%*%t(X0)))
df1$lci[201:300] <- df1$y[201:300] - qnorm(0.975)*df1$se[201:300]
df1$uci[201:300] <- df1$y[201:300] + qnorm(0.975)*df1$se[201:300]
rm(modeld)

## discrete indicators
df2 <- fit3$coefficients[2:5,]
df2$distance <- c(0,0.05,0.15,0.25)

summary(glm(deaths ~ arm, data = df, family = poisson()))
saveRDS(df1,"D:/Documents/df1_pois.RDS")

# p4 <- ggplot()+
#   geom_hline(yintercept = 0,lty=2)+
#   geom_ribbon(data = df1[df1$distance < 0.6, ], aes(x = distance, ymin=lci,ymax= uci), fill = unname(solar_color[14]), alpha = 0.2)+
#   geom_line(data = df1[df1$distance < 0.6, ], aes(x = distance, y = y))+
#   geom_point(data = df2, aes(x=distance,y=est), color= unname(solar_color[11]))+
#   geom_linerange(data = df2, aes(x=distance, ymin = lower,ymax = upper), color= unname(solar_color[11]))+
#   theme_solar()+
#   scale_x_continuous(expand = c(0.01,0))+
#   labs(x="Distance (km)", y = "Log relative risk")

p4 <- ggplot()+
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  geom_ribbon(data = df1, aes(x = distance, ymin=lci,ymax= uci, fill= factor(model)), alpha = 0.2)+
  geom_line(data = df1, aes(x = distance, y = y, color = factor(model)))+
  theme_solar()+
  scale_x_continuous(expand = c(0.01,0))+
  labs(x="Distance (km)", y = "Log relative risk")+
  scale_fill_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial"), values = unname(solar_color[c(11,14,16)]))+
  scale_color_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial"), values = unname(solar_color[c(11,14,16)])); p4

p5 <- ggplot(data = dfres,aes(x=n,y=aic))+
  geom_point()+
  geom_line()+
  theme_solar()+
  labs(x="Number of 100m distance indicators", y = "AIC")



p5 + p4
