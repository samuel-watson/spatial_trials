## Functions for continuous randomisation simulations


fn <- function(x,int,b,del){
  int + b*((1-((-1/50)*log(exp(-50*x/del) + exp(-50)))^4)^8)
}

genrep <- function(dfanal,f1,L,r){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn(distance, int, b, del),data = dfanal, 
                       start = list(int = 0, b = -0.2, del = 0.2),
                       lower = c(-10,-10,0.01), upper = c(10,10,1.0), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 3)
  }
  return(np)
}

generate_intervention <- function(data, max_dist, beta, n_locs, misspec, plot = TRUE){
  
  # spatially regulated sampling scheme
  iter <- 1
  int_idx <- sample(1:nrow(data),1)
  while(length(int_idx) < n_locs){
    int_idx_new <- sample(1:nrow(data),1)
    dists <- c(all_dists[int_idx,int_idx_new]) #c(st_distance(sampled_locs,dfp[int_idx_new,]))
    if(min(dists) > max_dist*0.01){ # change this line to implement a randomisation scheme where locations are spatially regulated
      int_idx <- c(int_idx, int_idx_new)
    } 
    iter <- iter + 1
    if(iter > 500) stop("Iterations exceed max")
  }
  
  data$intervention <- 0
  data[int_idx,'intervention'] <- 1
  
  # generate distances from intervention effect
  data$distance <- apply(all_dists[,which(data$intervention==1)],1,min)
  
  # generate intervention effect
  data$fn <- fun(data$distance,50,4,8,c(max_dist),1,data$t,misspec)
  data$y_true <- data$fn * beta
  # simulate outcome data
  data$sim_y <- data$y_true + L%*%rnorm(nrow(data))
  
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
      geom_sf(data=data, aes(color = sim_y), size = 0.1)+
      geom_sf(data=data[data$intervention==1,],color="red",size=2)+
      scico::scale_color_scico(palette = "roma", name = "Value")+
      theme_solar()+
      ggtitle("Simulated outcome")
    
    print( (p_dist + p_int) / (p_u + p_p) )
  }
  return(data)
}

# Upper bound on Del_E 
max_de <- function(n_locs){
  sqrt(8/(n_locs*3*sqrt(3)))
}