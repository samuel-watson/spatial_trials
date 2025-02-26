if(!USE_DATA){
  
  df <- read.csv("data/binka_compounds.csv")
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
  
  # ggplot()+
  #   geom_sf(data=dfpoly[dfpoly$cl > 0,], aes(fill = factor(arm)),alpha = 0.2)+
  #   scale_fill_manual(values = unname(solar_color[c(11,14)]),name = "Arm")+
  #   theme_solar()
  
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
  
  # p1b <- ggplot()+
  #   geom_sf(data=dfpoly_t, alpha = 0.4, color = NA)+
  #   geom_sf(data = dfp[dfp$deaths>0,], aes(color=factor(deaths)), size=0.3, alpha=0.2)+
  #   geom_sf(data=dfpoly[dfpoly$arm == "intervention",],lty=1,alpha=0.2)+
  #   scale_color_manual(values = unname(solar_color[c(11:14)]),name = "N deaths")+
  #   theme_solar()+
  #   theme(axis.text = element_blank())
  # 
  # p1
  # p1b
  
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
  
  # this is figure S4 (Supplementary Information)
  p1 + p2
  
  ggplot(data=dfp,aes(x=distance))+
    geom_histogram()+
    theme_solar()+
    ggtitle("Distance from intervention area")
  
  # prepare final fitting dataset
  
  dfanal <- as.data.frame(dfp)
  dfanal <- cbind(dfanal, df[,c('x','y')])
  dfanal <- dfanal[,2:12]
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
  
  saveRDS(dfanal,"data/binka_analysis_data.RDS")
  saveRDS(dists_i,"data/binka_dists.RDS")
} else {
  df <- read.csv("data/binka_compounds.csv")
  dfanal <- readRDS("data/binka_analysis_data.RDS")
  dists_i <- readRDS("data/binka_dists.RDS")
}