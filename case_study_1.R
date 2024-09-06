# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
source(paste0(getwd(),"/solarized.R"))

dfp <- readRDS(paste0(getwd(),"/data/data_case_study_1.RDS"))

# PLOTS
# plot the locations

p_loc <- ggplot()+
  geom_sf(data=dfp, size = 0.5, alpha = 0.5)+
  scale_color_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  theme_solar()+
  ggtitle("Simulated locations"); p_loc

p_dist <- ggplot()+
  geom_sf(data=dfp, aes(color = distance), size = 0.1)+
  geom_sf(data=dfp[dfp$intervention==1,],color="red",size=2)+
  scico::scale_color_scico(palette = "batlow", name = "Distance")+
  theme_solar()+
  ggtitle("Distance"); p_dist

p_int <- ggplot()+
  geom_sf(data=dfp, aes(color = y_true), size = 0.1)+
  geom_sf(data=dfp[dfp$intervention==1,],color="red",size=2)+
  scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
  theme_solar()+
  ggtitle("Intervention effect"); p_int

p_u <- ggplot()+
  geom_sf(data=dfp, aes(color = u), size = 0.1)+
  scico::scale_color_scico(palette = "batlow", name = "True\neffect")+
  theme_solar()+
  ggtitle("Latent spatial effect"); p_u

p_p <- ggplot()+
  geom_sf(data=dfp, aes(color = sim_p), size = 0.1)+
  geom_sf(data=dfp[dfp$intervention==1,],color="red",size=2)+
  scico::scale_color_scico(palette = "roma", name = "Simulated\nprobability")+
  theme_solar()+
  ggtitle("Probability of outcome"); p_p


## lets look at sampling

dfanal <- as.data.frame(dfp)
dfanal <- dfanal[,c(3:6,8:10)]

model2 <- Model$new(
  ~ b_eff*((1 - (-0.02 * log(exp(-50*(distance/b_del)) + exp(-50)))^(4))^(16)) + (1|hsgp_fexp(X,Y)),
  data=dfanal,
  covariance = c(0.25,0.25),
  mean = c(b0,b1,max_dist-0.01),
  family = binomial()
)

model2 <- Model$new(
  ~ twoway0(distance,16,100) + (1|hsgp_fexp(X,Y)),
  data=dfanal,
  covariance = c(0.2,0.25),
  mean = c(b0,b1,max_dist,2.01),
  family = binomial()
)

model2$set_trace(1)
model2$mcmc_options$warmup <- 100
model2$mcmc_options$samps <- 50
model2$covariance$hsgp(m = c(20,20), L = c(1.2,1.2))
model2$update_parameters(cov.pars = c(0.2,0.2))

fit2 <- model2$MCML(y = dfp$sim_y, 
                    method = "saem", 
                    algo = 2, 
                    pr.average = FALSE,
                    se.theta = FALSE,
                    alpha = 0.6,
                    convergence.prob = 0.99,
                    conv.criterion = 2,
                    lower.bound = c(-10,-10,0.01), 
                    upper.bound = c(10,10,2))

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
  ~ b_eff*((1 - (-0.02 * log(exp(-50*(distance/b_del)) + exp(-50)))^(4))^(16)) + (1|gr(cl)),
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
