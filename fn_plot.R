# PLOTTING CONTINUOUS FUNCTIONS

require(ggplot2)
source(paste0(getwd(),"/solarized.R"))

fn1 <- function(x,l,kappa,nu,del_e, del_i, b){
  b*((1-((sign(x)/l)*log(exp(sign(x)*l*(x - del_i)/(del_e - del_i)) + exp(l*(sign(x)+1)/2)))^kappa)^nu) 
}

dfn <- expand.grid(d = seq(-1.5,1.5,length.out = 100), nu = c(1,2,4,8,16), l = -50,
                   kappa = c(0.5,1,2,4,8,16), beta = -1, del_e = 1, del_i = c(0,-0.5,-1),
                   fn = NA)

for(i in 1:nrow(dfn)){
  dfn$fn[i] <- fn1(dfn$d[i],dfn$l[i], dfn$kappa[i], dfn$nu[i], dfn$del_e[i], dfn$del_i[i], -1 )
}

ggplot(data=dfn[dfn$kappa <= dfn$nu,],aes(x=d,y=fn,color=factor(kappa)))+
  geom_hline(yintercept = c(-1,-0.5,0), lty= 3)+
  geom_vline(xintercept = c(0), lty = 2)+
  geom_line()+
  facet_grid(del_i~nu, scale= "free_y")+
  theme_solar()

# original version used the colour palette below
# +
#   scale_color_manual(name = expression(kappa), values = unname(solar_color[c(10,12,13,14,15,16)]))
