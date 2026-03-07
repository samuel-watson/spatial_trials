# binka re-analysis

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
source("src/solarized.R")

# load required functions
Rcpp::sourceCpp("src/perm_test.cpp")
source("src/reanalysis_fn.R")

# if USE_DATA is false it will regenerate/process the data (which can be slow), otherwise it'll load it
# note, to generate Figure S4 in the Supplementary Information, USE_DATA must be FALSE
# to generate the figure separately see the file "src/process_binka_data.R"
USE_DATA <- TRUE 

source("src/process_binka_data.R")

## standard cluster trial analysis

dfanal$treated <- (dfanal$arm == "intervention")*1
model_std <- Model$new(
  deaths ~ treated + (1|gr(cl)),
  data=dfanal,
  offset = dfanal$expected,
  family = gaussian()
)

model_std$fit()

# first model, no adjustment
# null model for permutation test

model_null <- Model$new(
  ~ (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05,0.05),
  mean = c(0.01),
  offset = dfanal$expected,
  family = gaussian()
)

model_null$covariance$hsgp(m = c(15,15), L = c(1.1,1.1))
model_null$update_parameters(cov.pars = c(0.05,0.05))
model_null$set_trace(1)

fit_null <- model_null$MCML(y = dfanal$deaths, reml = FALSE)

S <- model_null$Sigma()
Li <- solve(t(chol(S)))
rm(model_null,S)

n_locs <- length(unique(df[df$arm=="intervention","cluster"]))
pt <- new_r_stat(dfanal$deaths-dfanal$expected,Li,dfanal$distance,dists_i,0.01,2,c(-1),c(0),50,4,8,1)
permute_p_value_b(pt,n_locs,0,200,c(0.5,-0.02))

# fit full model, no adjustment

model <- Model$new(
  ~ twoway2(distance,8,4,50) + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05,0.05),
  mean = c(0.07,-0.3,0.3,0.2),
  offset = dfanal$expected,
  family = gaussian()
)

model$covariance$hsgp(m = c(15,15), L = c(1.1,1.1))
model$update_parameters(cov.pars = c(0.05,0.05))

fit0 <- model$MCML(y = dfanal$deaths, 
                   reml = FALSE,
                  lower.bound = c(-10,-10,0,0), 
                  upper.bound = c(10,10,2,2))

fit0

# bootstrapped confidence intervals
# get full covariance matrix
# this is slow due to the size of the matrices, requires ~8 GB of memory

model00 <- Model$new(
  ~ twoway2(distance,8,4,50) + (1|fexp(x_re,y_re)),
  data=dfanal,
  covariance = model$covariance$parameters,
  mean = model$mean$parameters,
  offset = dfanal$expected,
  family = gaussian()
)

L <- t(chol(model$Sigma()))
f0 <- model00$fitted()
rm(model00)

# check NLS

fitn <- nls(deaths ~ fn0(distance, int,b, del_e, del_i),data = dfanal, 
    start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1),
    lower = c(-10,-10,0.01,0.0), upper = c(10,10,1.0,0.5), algorithm = "port")

f01 <- fitted(fitn)

genrep0(dfanal, f01, L)

res <- pbapply::pbreplicate(1000, genrep0(dfanal,f01, L))
res <- Reduce(rbind, res)

fit0$coefficients$est[2] + sd(res[,2], na.rm=T)*qnorm(0.975)
fit0$coefficients$est[2] - sd(res[,2], na.rm=T)*qnorm(0.975)

fit0$coefficients$est[3] + sd(res[,4], na.rm=T)*qnorm(0.975)
fit0$coefficients$est[3] - sd(res[,4], na.rm=T)*qnorm(0.975)

fit0$coefficients$est[4] + sd(res[,3], na.rm=T)*qnorm(0.975)
fit0$coefficients$est[4] - sd(res[,3], na.rm=T)*qnorm(0.975)

#extract the information matrix for later plotting

M <- model$information_matrix()
rm(model)

# need to add model that controls for potential distance - either: 
# i) distance indicators; ii) smooth function

for(i in 0:8){
  dfanal$tmp <- I(abs(dfanal$distance_all) > (i-1)*0.2 & abs(dfanal$distance_all) <= i*0.2)*1
  colnames(dfanal)[ncol(dfanal)] <- paste0("distance_all",i)
}
dfanal$distance_all_sq <- dfanal$distance_all^2

# distance indicator model - set easy starting values to make it a bit quicker!

model2a <- Model$new(
  ~ twoway2(distance,8,4,50) + distance_all1 + distance_all2 + distance_all3 + distance_all4 + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.08016, 0.16235),
  mean = c(0.45605, -0.28641, 0.41107, 0.95507, -0.07335, -0.05357, -0.02834, -0.17584),
  offset = dfanal$expected,
  family = gaussian()
)

model2a$covariance$hsgp(m = c(15,15), L = c(1.2,1.2))
model2a$update_parameters(cov.pars = c(0.05161, 0.15722))

fit2a <- model2a$MCML(y = dfanal$deaths,
                      reml = FALSE,
                      lower.bound = c(-10,-10,0,0,-10,-10,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10,10,10))

fit2a

# check NLS

fitn <- nls(deaths ~ fn2a(distance,distance_all1,distance_all2,distance_all3,distance_all4, int,b, del_e, del_i,b_1, b_2,b_3,b_4),data = dfanal, 
            start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_1 = 0, b_2 = 0,b_3 = 0, b_4 = 0),
            lower = c(-10,-10,0.01,0.0,rep(-10,4)), upper = c(10,10,1.0,0.5,rep(10,4)), algorithm = "port")

f02a <- fitted(fitn)
f2a <- model2a$fitted()

genrep2a(dfanal, f2a, L)

res <- pbapply::pbreplicate(1000, genrep2a(dfanal,f2a, L))
res <- Reduce(rbind, res)

fit2a$coefficients$est[2] + sd(res[,2], na.rm=T)*qnorm(0.975)
fit2a$coefficients$est[2] - sd(res[,2], na.rm=T)*qnorm(0.975)

fit2a$coefficients$est[3] + sd(res[,4], na.rm=T)*qnorm(0.975)
fit2a$coefficients$est[3] - sd(res[,4], na.rm=T)*qnorm(0.975)

fit2a$coefficients$est[4] + sd(res[,3], na.rm=T)*qnorm(0.975)
fit2a$coefficients$est[4] - sd(res[,3], na.rm=T)*qnorm(0.975)

M2a <- model2a$information_matrix()
rm(model2a)

# degree-2 polynomial adjustment

model2b <- Model$new(
  ~ twoway2(distance,8,4,50) + distance_all + distance_all_sq + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05, 0.16),
  mean = c(0.3, -0.3, 0.2, 0.98, 0.01, 0.01),
  offset = dfanal$expected,
  family = gaussian()
)

model2b$covariance$hsgp(m = c(15,15), L = c(1.1,1.1))
model2b$update_parameters(cov.pars = c(0.05, 0.16))

fit2b <- model2b$MCML(y = dfanal$deaths,
                      reml = FALSE,
                      lower.bound = c(-10,-10,0,0,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10))

fit2b

fitn <- nls(deaths ~ fn2b(distance,distance_all,distance_all_sq,int,b, del_e, del_i,b_1, b_2),data = dfanal, 
            start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_1 = 0, b_2 = 0),
            lower = c(-10,-10,0.01,0.0,rep(-10,2)), upper = c(10,10,1.0,0.5,rep(10,2)), algorithm = "port")

f02b <- fitted(fitn)
f2b <- model2b$fitted()

genrep2b(dfanal, f2b, L)

res <- pbapply::pbreplicate(100, genrep2b(dfanal,f2b, L))
# res <- Reduce(rbind, res)

fit2b$coefficients$est[2] + sd(res[2,], na.rm=T)*qnorm(0.975)
fit2b$coefficients$est[2] - sd(res[2,], na.rm=T)*qnorm(0.975)

fit2b$coefficients$est[3] + sd(res[4,], na.rm=T)*qnorm(0.975)
fit2b$coefficients$est[3] - sd(res[4,], na.rm=T)*qnorm(0.975)

fit2b$coefficients$est[4] + sd(res[3,], na.rm=T)*qnorm(0.975)
fit2b$coefficients$est[4] - sd(res[3,], na.rm=T)*qnorm(0.975)

M2b <- model2b$information_matrix()
rm(model2b)

## final model - one-way only with distance polynomial


model2c <- Model$new(
  ~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*((distance)/(del_e))) + exp((-25)*(1+sign0(distance))))))^(4))^(8)) + distance_all + distance_all_sq + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.05, 0.16),
  mean = c(0.3, -0.3, 0.98, 0.01, 0.01),
  offset = dfanal$expected,
  family = gaussian()
)

model2c$covariance$hsgp(m = c(15,15), L = c(1.1,1.1))
model2c$update_parameters(cov.pars = c(0.05, 0.16))

fit2c <- model2c$MCML(y = dfanal$deaths,
                      reml = FALSE,
                      lower.bound = c(-10,-10,0,-10,-10), 
                      upper.bound = c(10,10,2,10,10))

fit2c

fitn <- nls(deaths ~ fn2b(distance,distance_all,distance_all_sq,int,b, del_e, 0,b_1, b_2),data = dfanal, 
            start = list(int = 0, b = -0.2, del_e = 0.2, b_1 = 0, b_2 = 0),
            lower = c(-10,-10,0.01,rep(-10,2)), upper = c(10,10,1.0,rep(10,2)), algorithm = "port")

f02c <- fitted(fitn)
f2c <- model2c$fitted()

genrep(dfanal, f2c, L)

res <- pbapply::pbreplicate(200, genrep(dfanal,f2c, L))
res <- Reduce(rbind, res)

fit2c$coefficients$est[2] + sd(res[,2], na.rm=T)*qnorm(0.975)
fit2c$coefficients$est[2] - sd(res[,2], na.rm=T)*qnorm(0.975)

fit2c$coefficients$est[3] + sd(res[,3], na.rm=T)*qnorm(0.975)
fit2c$coefficients$est[3] - sd(res[,3], na.rm=T)*qnorm(0.975)

M2c <- model2c$information_matrix()
rm(model2c)

### plot the function

df1 <- data.frame(distance = rep(seq(-1,1.5,length.out=100),4),
                  model = rep(1:4,each=100))

df1 <- df1[df1$distance!=0,]
df1$expected <- 0
df1$distance_all1 <- 0
df1$distance_all2 <- 0
df1$distance_all3 <- 0
df1$distance_all4 <- 0
df1$distance_all_sq <- 0
df1$distance_all <- 0

#df1$y <- fn(df1$distance[1:100],-50,2.05801,16, 1.58349, 0.15759 ,-0.29911)
df1$y <- c(fn(df1$distance[1:100],-50,4,8,fit0$coefficients$est[4], fit0$coefficients$est[3] ,fit0$coefficients$est[2]), 
  fn(df1$distance[1:100],-50,4,8, fit2a$coefficients$est[4] , fit2a$coefficients$est[3] ,fit2a$coefficients$est[2]) , 
 fn(df1$distance[1:100],-50,4,8, fit2b$coefficients$est[4], fit2b$coefficients$est[3] ,fit2b$coefficients$est[2]),
 fn(df1$distance[1:100],-50,4,8, fit2c$coefficients$est[3], 0 ,fit2c$coefficients$est[2])) 

df1$cl <- sample(1:10,nrow(df1),replace=TRUE)
df1$se <- NA
df1$lci <- NA
df1$uci <- NA

modeld <- Model$new(
  ~ twoway2(distance,8,4,50) + (1|gr(cl)),
  data=df1[1:100,],
  covariance = c( 0.05),
  mean = fit0$coefficients$est[1:4],
  family = poisson()
)
X0 <- modeld$mean$X
X0[,1] <- 0
Mi <- solve(M)
df1$se[1:100] <- sqrt(diag(X0%*%Mi%*%t(X0)))
df1$lci[1:100] <- df1$y[1:100] - qnorm(0.975)*df1$se[1:100]
df1$uci[1:100] <- df1$y[1:100] + qnorm(0.975)*df1$se[1:100]
rm(modeld)

# fit 2 a
modeld <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all1 + distance_all2 + distance_all3 + distance_all4 + (1|gr(cl)),
  data=df1[101:200,],
  covariance = c( 0.05),
  mean = fit2a$coefficients$est[1:8],
  family = poisson()
)
X0 <- modeld$mean$X
X0[,1] <- 0
M2ai <- solve(M2a)
df1$se[101:200] <- sqrt(diag(X0%*%M2ai%*%t(X0)))
df1$lci[101:200] <- df1$y[101:200] - qnorm(0.975)*df1$se[101:200]
df1$uci[101:200] <- df1$y[101:200] + qnorm(0.975)*df1$se[101:200]
rm(modeld)

# fit 2 b
modeld <- Model$new(
  ~ twoway2(distance,16,4,50) + distance_all + distance_all_sq + (1|gr(cl)),
  data=df1[201:300,],
  covariance = c( 0.05),
  mean = fit2b$coefficients$est[1:6],
  family = poisson()
)

X0 <- modeld$mean$X
X0[,1] <- 0
M2bi <- solve(M2b)
df1$se[201:300] <- sqrt(diag(X0%*%M2bi%*%t(X0)))
df1$lci[201:300] <- df1$y[201:300] - qnorm(0.975)*df1$se[201:300]
df1$uci[201:300] <- df1$y[201:300] + qnorm(0.975)*df1$se[201:300]
rm(modeld)

# fit 2 c
modeld <- Model$new(
  ~ b_eff * ((1 - (sign0(distance)*(-0.02)*(log(exp((-50)*sign0(distance)*((distance)/(del_e))) + exp((-25)*(1+sign0(distance))))))^(4))^(8)) + distance_all + distance_all_sq + (1|gr(cl)),
  data=df1[301:400,],
  covariance = c( 0.05),
  mean = fit2c$coefficients$est[1:5],
  family = poisson()
)

X0 <- modeld$mean$X
X0[,1] <- 0
M2ci <- solve(M2c)
df1$se[301:400] <- sqrt(diag(X0%*%M2ci%*%t(X0)))
df1$lci[301:400] <- df1$y[301:400] - qnorm(0.975)*df1$se[301:400]
df1$uci[301:400] <- df1$y[301:400] + qnorm(0.975)*df1$se[301:400]
rm(modeld)

p4 <- ggplot()+
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  geom_ribbon(data = df1[df1$distance<=1,], aes(x = distance, ymin=lci,ymax= uci, fill= factor(model)), alpha = 0.2)+
  geom_line(data = df1[df1$distance<=1,], aes(x = distance, y = y, color = factor(model)))+
  theme_solar()+
  scale_x_continuous(expand = c(0.01,0))+
  labs(x="Distance (km)", y = "Risk difference")+
  scale_fill_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial", "One-way"), values = unname(solar_color[c(9,11,14,16)]))+
  scale_color_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial", "One-way"), values = unname(solar_color[c(9,11,14,16)])); p4

p4 # this is Figure 3 in the article
