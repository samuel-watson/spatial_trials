# binka re-analysis

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
source(paste0(getwd(),"/solarized.R"))

# load required functions
Rcpp::sourceCpp("src/perm_test.cpp")

df <- read.csv(paste0(getwd(),"/data/binka_compounds.csv"))
df <- df[df$expected > 0 ,]
df$t <- 1
dfp <- rts2::create_points(df,pos_vars = c('x','y'), t_var = "t")
dfp$cl <- df$cluster
dfp$arm <- df$arm
dfp$deaths <- df$deaths
dfp$expected <- df$expected
dfp$nets <- df$nets
cl <- unique(df$cluster)

# create convex hull of cluster shapes

for(i in cl){
  if(!require(concaveman))install.packages("concaveman")
  p1 <- concaveman::concaveman(dfp[dfp$cl==i,], concavity = 5) # change the concavity parameter for different levels of smoothing
  p1$cl <- i
  p1$arm <- df[df$cluster == i, 'arm'][1]
  if(exists("dfpoly")){
    dfpoly <- rbind(dfpoly, p1)
  } else {
    dfpoly <- p1
  }
}

# plot the cluster areas

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

p1
p1b

## now process into a "spatial trial"

## calculate distances

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

## distance to any boundary

all_area <- st_union(dfpoly)
all_area_boundary <- st_boundary(all_area)
dfp$distance_all <- NA

for(i in 1:nrow(dfp)){
  dfp$distance_all[i] <- -1 * st_distance(dfp[i,],all_area_boundary) 
  cat("\rRow: ",i," of ",nrow(dfp))
}

p2 <- ggplot()+
  geom_sf(data=dfpoly_t, alpha = 0.4, color = NA)+
  geom_sf(data = dfp, aes(color=distance), size=0.2, alpha=0.2)+
  geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
  scico::scale_color_scico(palette="roma")+
  theme_solar()+
  theme(axis.text = element_blank())


p1 + p2

ggplot(data=dfp,aes(x=distance))+
  geom_histogram()+
  theme_solar()+
  ggtitle("Distance from intervention area")

# prepare final fitting dataset

dfanal <- as.data.frame(dfp)
dfanal <- cbind(dfanal, df[,c('x','y')])
dfanal <- dfanal[,4:12]
# rescale x and y to [-1,1] for approx GP 
xrange <- range(dfanal$x)
yrange <- range(dfanal$y)
scale_f <- max(diff(xrange),diff(yrange))
dfanal$x_re <- -1 + 2*(dfanal$x - min(dfanal$x))/scale_f  #-1 + (2 / diff(range(dfanal$x)))*(dfanal$x - min(dfanal$x))
dfanal$y_re <- -1 + 2*(dfanal$y - min(dfanal$y))/scale_f # -1 + (2 / diff(range(dfanal$y)))*(dfanal$x - min(dfanal$y))
dfanal <- dfanal[order(dfanal$y),]
# jitter the duplicated location
locs <- paste0(dfanal$x_re, dfanal$y_re)
dfanal[duplicated(locs),'y_re'] <- dfanal[duplicated(locs),'y_re'] + 1e-6

# create distance matrix for observations to potential intervention areas for permutation test
dists_i <- matrix(NA,nrow=nrow(dfp),ncol=nrow(dfpoly))

for(i in 1:nrow(dfp)){
  for(j in 1:nrow(dfpoly)){
    if(dfp$cl[i] == (j-1)){
      dists_i[i,j] <- -1 * st_distance(dfp[i,],dfpoly[j,])
    } else {
      dists_i[i,j] <- st_distance(dfp[i,],dfpoly[j,])
    }
  }
  cat("\rRow: ",i, " of ",nrow(dfp))
}

saveRDS(dfanal,paste0(getwd(),"/data/binka_analysis_data.RDS"))
saveRDS(dists_i,paste0(getwd(),"/data/binka_dists.RDS"))



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

fit_null <- model_null$MCML(y = dfanal$deaths)

S <- model_null$Sigma()
E <- eigen(S)
B <- E$vectors%*%diag(1/sqrt(E$values))
rm(model_null,S,E)

n_locs <- sum(dfpoly$arm=="intervention")
pt <- new_r_stat(dfanal$sim_y-0.02165-dfanal$expected,t(B),dfanal$distance,dists_i,0.01,2,c(-1),c(0),50,4,8,1)
permute_p_value_b(pt,n_locs,0,200,c(0.5,-0.02))

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
model$set_trace(1)

fit0 <- model$MCML(y = dfanal$deaths, 
                  lower.bound = c(-10,-10,0,0), 
                  upper.bound = c(10,10,2,2))

fit0

#extract the information matrix for later plotting

M <- model$information_matrix()
rm(model)

# use these standard errors until the package is updated:
se <- sqrt(diag(solve(model$information_matrix())))

fit0$coefficients$est[1:4] + qnorm(0.975)*se[1:4]
fit0$coefficients$est[1:4] - qnorm(0.975)*se[1:4]
2*(1-pnorm(abs(fit0$coefficients$est[2]/se[2])))

# need to add model that controls for potential distance - either: 
# i) distance indicators; ii) smooth function

for(i in 0:8){
  dfanal$tmp <- I(abs(dfanal$distance_all) > (i-1)*0.2 & abs(dfanal$distance_all) <= i*0.2)*1
  colnames(dfanal)[ncol(dfanal)] <- paste0("distance_all",i)
}
dfanal$distance_all_sq <- dfanal$distance_all^2

# distance indicator model

model2a <- Model$new(
  ~ twoway2(distance,8,4,50) + distance_all1 + distance_all2 + distance_all3 + distance_all4 + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.08016, 0.16235),
  mean = c(0.45605, -0.28641, 0.41107, 0.95507, -0.07335, -0.05357, -0.02834, -0.17584),
  offset = dfanal$expected,
  family = gaussian()
)

model2a$covariance$hsgp(m = c(15,15), L = c(1.2,1.2))
model2a$set_trace(1)
model2a$update_parameters(cov.pars = c(0.05161, 0.15722))

fit2a <- model2a$MCML(y = dfanal$deaths,
                      lower.bound = c(-10,-10,0,0,-10,-10,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10,10,10))

fit2a

M2a <- model2a$information_matrix()
rm(model2a)
# degree-2 polynomial adjustment

model2b <- Model$new(
  ~ twoway2(distance,8,4,50) + distance_all + distance_all_sq + (1|hsgp_fexp(x_re,y_re)),
  data=dfanal,
  covariance = c(0.07511, 0.14205),
  mean = c(0.3, -0.3, 0.2, 0.98, 0.01, 0.01),
  offset = dfanal$expected,
  family = gaussian()
)

model2b$covariance$hsgp(m = c(15,15), L = c(1.1,1.1))
model2b$set_trace(1)
model2b$update_parameters(cov.pars = c(0.05, 0.16))

fit2b <- model2b$MCML(y = dfanal$deaths,
                      lower.bound = c(-10,-10,0,0,-10,-10), 
                      upper.bound = c(10,10,2,2,10,10))

fit2b

M2b <- model2b$information_matrix()
rm(model2b)

### plot the function

fn <- function(x,l,kappa,nu,del_e, del_i, b){
  b*((1-((sign(x)/l)*log(exp(sign(x)*l*(x + del_i)/(del_e + del_i)) + exp(l*(sign(x)+1)/2)))^kappa)^nu) 
}

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
df1$y <- c(fn(df1$distance[1:100],-50,4,8,fit0$coefficients$est[3], fit0$coefficients$est[4] ,fit0$coefficients$est[2]), 
  fn(df1$distance[1:100],-50,4,8, fit2a$coefficients$est[3] , fit2a$coefficients$est[4] ,fit2a$coefficients$est[2]) , 
 fn(df1$distance[1:100],-50,4,8, fit2b$coefficients$est[3], fit2b$coefficients$est[4] ,fit2b$coefficients$est[2])) 

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

## discrete indicators

p4 <- ggplot()+
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  geom_ribbon(data = df1[df1$distance<=1,], aes(x = distance, ymin=lci,ymax= uci, fill= factor(model)), alpha = 0.2)+
  geom_line(data = df1[df1$distance<=1,], aes(x = distance, y = y, color = factor(model)))+
  theme_solar()+
  scale_x_continuous(expand = c(0.01,0))+
  labs(x="Distance (km)", y = "Risk difference")+
  scale_fill_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial"), values = unname(solar_color[c(11,14,16)]))+
  scale_color_manual(name = "Model", labels = c("No adjustment", "Indicators", "Polynomial"), values = unname(solar_color[c(11,14,16)])); p4

p4