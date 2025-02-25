### FUNCTIONS - DISCRETE SIMULATION

generate_intervention <- function(data, del_e, del_i, beta, n_locs, plot = TRUE){
  
  # sample locations  
  idx <- sort(sample(1:nrow(dfi),n_locs))
  data$distance <- apply(dists_i[,idx],1,min)
  data$fn <- fun(data$distance,50,4,8,c(del_e,-del_i),1,data$t)
  
  # generate intervention effect
  data$y_true <- data$fn * beta
  # simulate outcome data
  data$u <- drop(L%*%rnorm(nrow(L)))
  data$sim_y <- data$y_true + data$u
  
  if(plot){
    p_dist <- ggplot()+
      geom_sf(data=data, aes(color = distance), size = 0.1)+
      geom_point(data=dfi[idx,],aes(x=x,y=y),color="red",size=2)+
      scico::scale_color_scico(palette = "batlow", name = "Distance")+
      theme_solar()+
      ggtitle("Distance")
    
    p_int <- ggplot()+
      geom_sf(data=data, aes(color = y_true), size = 0.1)+
      ggforce::geom_circle(data=dfi[idx,],aes(x0=x,y0=y,r=radius),fill="red",alpha = 0.1,color= NA)+
      ggforce::geom_circle(data=dfi[-idx,],aes(x0=x,y0=y,r=radius),fill="light blue",alpha = 0.3, color = NA)+
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
      geom_point(data=dfi[idx,],aes(x=x,y=y),color="red",size=2)+
      scico::scale_color_scico(palette = "roma", name = "y")+
      theme_solar()+
      ggtitle("Outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

# fn2 <- function(x,d1,d2,d3,int,b,del_e, del_i, b_d1, b_d2, b_d3){
#   int + d1*b_d1 + d2*b_d2+ d3*b_d3 + b*((1-((sign(x)/-50)*log(exp(sign(x)*-50*(x + del_i)/(del_e + del_i)) + exp(-50*(sign(x)+1)/2)))^4)^8) 
# }

fn2 <- function(x,d1,int,b,del_e, del_i, b_d1){
  int + d1*b_d1 + b*((1-((sign(x)/-50)*log(exp(sign(x)*-50*(x + del_i)/(del_e + del_i)) + exp(-50*(sign(x)+1)/2)))^4)^8) 
}

genrep <- function(dfanal,f1,L,oneway){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  if(oneway){
    fitn <- tryCatch(nls(ystar ~ fn2(distance, d2,int,b, del_e, 0, b_d1),data = dfanal, 
                         start = list(int = 0, b = -0.2, del_e = 0.2, b_d1 = 0.0),
                         lower = c(-10,-3,0.0,-10), upper = c(10,3,1.0,10), algorithm = "port"), error= function(i)return(NA))
  } else {
    fitn <- tryCatch(nls(ystar ~ fn2(distance, d2, int,b, del_e, del_i, b_d1),data = dfanal, 
                         start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_d1 = 0.0),
                         lower = c(-10,-3,0.01,0.0,-10), upper = c(10,3,1.0,1.0,10), algorithm = "port"), error= function(i)return(NA))
  }
  
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 3)
  }
  return(np)
}

# Upper bound on Del_E 
max_de <- function(n_locs){
  sqrt(8/(n_locs*3*sqrt(3)))
}